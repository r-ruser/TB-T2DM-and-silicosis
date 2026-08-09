source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "DESeq2", "SummarizedExperiment", "WGCNA",
                  "AnnotationDbi", "org.Hs.eg.db", "ggplot2", "patchwork",
                  "svglite", "ragg"))
suppressPackageStartupMessages(library(WGCNA))

options(stringsAsFactors = FALSE)
WGCNA::allowWGCNAThreads(nThreads = 2)

acc <- "GSE114192"
model_file <- file.path(path_result, "02_bulk", acc, "models", "M01_GSE114192_primary_models.rds")
meta_file <- file.path(path_data, "02_GEO_bulk", acc, "processed", paste0(acc, "_metadata_curated.csv"))
if (!file.exists(model_file) || !file.exists(meta_file)) stop("Run R/24_gse114192_primary_deseq2.R first")

primary <- readRDS(model_file)
dds <- primary$dds
meta <- data.table::fread(meta_file, encoding = "UTF-8")[include_primary == TRUE]
meta[, group := stats::relevel(factor(group), ref = "TB_only")]
meta[, site := factor(site)]
rownames(meta) <- meta$sample_id

vsd <- DESeq2::vst(dds, blind = FALSE)
expr <- SummarizedExperiment::assay(vsd)
expr <- expr[, meta$sample_id, drop = FALSE]
gene_mad <- apply(expr, 1, stats::mad, na.rm = TRUE)
n_keep <- min(5000L, sum(is.finite(gene_mad) & gene_mad > 0))
selected <- names(sort(gene_mad, decreasing = TRUE))[seq_len(n_keep)]
dat_expr <- as.data.frame(t(expr[selected, , drop = FALSE]))

qc <- WGCNA::goodSamplesGenes(dat_expr, verbose = 0)
if (!qc$allOK) dat_expr <- dat_expr[qc$goodSamples, qc$goodGenes, drop = FALSE]
meta <- meta[match(rownames(dat_expr), sample_id)]
if (anyNA(meta$sample_id)) stop("WGCNA sample metadata reordering failed")
rownames(meta) <- meta$sample_id
if (nrow(dat_expr) < 20L || ncol(dat_expr) < 1000L) stop("Insufficient data after WGCNA QC")

powers <- c(1:10, seq(12, 20, 2))
sft <- WGCNA::pickSoftThreshold(dat_expr, powerVector = powers, networkType = "signed",
                                corFnc = "cor", verbose = 1)
fit <- data.table::as.data.table(sft$fitIndices)
data.table::setnames(fit, old = intersect(names(fit), c("Power", "SFT.R.sq", "slope", "mean.k.")),
                     new = intersect(c("power", "scale_free_r2", "slope", "mean_connectivity"),
                                     c("power", "scale_free_r2", "slope", "mean_connectivity")))
if (!all(c("power", "scale_free_r2", "slope") %in% names(fit))) {
  stop("Unexpected pickSoftThreshold output columns: ", paste(names(fit), collapse = ", "))
}
candidate <- fit[scale_free_r2 >= 0.80 & slope < 0]
if (nrow(candidate)) {
  soft_power <- min(candidate$power)
} else if (nrow(fit[slope < 0])) {
  soft_power <- fit[slope < 0][which.max(scale_free_r2), power]
} else {
  soft_power <- fit[which.max(scale_free_r2), power]
}

net <- WGCNA::blockwiseModules(
  dat_expr, power = soft_power, networkType = "signed", TOMType = "signed",
  corType = "pearson", minModuleSize = 30,
  mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE,
  maxBlockSize = 6000, saveTOMs = FALSE, verbose = 2
)
module_colors <- WGCNA::labels2colors(net$colors)
mes <- WGCNA::orderMEs(net$MEs)

assoc <- data.table::rbindlist(lapply(names(mes), function(me_name) {
  d <- data.frame(me = mes[[me_name]], group = meta$group, site = meta$site)
  fit_lm <- stats::lm(me ~ group + site, data = d)
  cf <- summary(fit_lm)$coefficients
  coef_name <- grep("^groupTB_DM$", rownames(cf), value = TRUE)
  if (length(coef_name) != 1L) stop("TB-DM coefficient missing for ", me_name)
  est <- cf[coef_name, "Estimate"]
  se <- cf[coef_name, "Std. Error"]
  data.table::data.table(
    module = sub("^ME", "", me_name), eigengene = me_name,
    estimate = est, std_error = se, ci_low = est - 1.96 * se,
    ci_high = est + 1.96 * se, p_value = cf[coef_name, "Pr(>|t|)"]
  )
}))
assoc <- assoc[module != "0"]
assoc[, fdr := stats::p.adjust(p_value, method = "BH")]
sizes <- data.table::data.table(module = as.character(net$colors))[, .(module_size = .N), by = module]
assoc <- merge(assoc, sizes, by = "module", all.x = TRUE)
data.table::setorder(assoc, fdr, p_value)

membership <- WGCNA::signedKME(dat_expr, mes, corFnc = "cor")
membership <- data.table::as.data.table(membership, keep.rownames = "gene_id")
membership_long <- data.table::melt(membership, id.vars = "gene_id",
                                    variable.name = "eigengene", value.name = "kME")
membership_long[, module := sub("^kME", "", eigengene)]
membership_long[, gene_symbol := unname(AnnotationDbi::mapIds(
  org.Hs.eg.db::org.Hs.eg.db, keys = sub("[.][0-9]+$", "", gene_id),
  keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first"))]

result_dir <- file.path(path_result, "03_cross_disease", "WGCNA", acc)
invisible(lapply(file.path(result_dir, c("tables", "source_data", "models", "figures")),
                 dir.create, recursive = TRUE, showWarnings = FALSE))
safe_write_csv(fit, file.path(result_dir, "source_data", "SD01_GSE114192_soft_threshold.csv"))
safe_write_csv(assoc, file.path(result_dir, "tables", "T01_GSE114192_module_TBDM_association.csv"))
safe_write_csv(membership_long, file.path(result_dir, "source_data", "SD02_GSE114192_module_membership.csv"))
saveRDS(list(network = net, module_colors = module_colors, eigengenes = mes,
             selected_power = soft_power, genes = colnames(dat_expr), sample_metadata = meta),
        file.path(result_dir, "models", "M01_GSE114192_WGCNA_network.rds"))

fit[, selected := power == soft_power]
p_a <- ggplot2::ggplot(fit, ggplot2::aes(power, scale_free_r2)) +
  ggplot2::geom_hline(yintercept = 0.80, linetype = 2, linewidth = 0.3, colour = "#767676") +
  ggplot2::geom_line(linewidth = 0.45, colour = "#416A9A") +
  ggplot2::geom_point(ggplot2::aes(fill = selected), shape = 21, size = 1.8,
                      colour = "#272727", stroke = 0.25) +
  ggplot2::scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = "#B84A4A"), guide = "none") +
  ggplot2::scale_x_continuous(breaks = powers) +
  ggplot2::labs(title = "Signed-network soft threshold", subtitle = paste("Selected power:", soft_power),
                x = "Power", y = expression("Scale-free topology " * R^2))

plot_assoc <- assoc[order(fdr, p_value)][seq_len(min(.N, 20L))]
plot_assoc[, module_label := sprintf("Module %s (n=%d)", module, module_size)]
plot_assoc[, module_label := factor(module_label, levels = rev(module_label))]
plot_assoc[, evidence := ifelse(fdr < 0.05, "BH-FDR < 0.05", "Not FDR-significant")]
p_b <- ggplot2::ggplot(plot_assoc, ggplot2::aes(estimate, module_label, colour = evidence)) +
  ggplot2::geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3, colour = "#767676") +
  ggplot2::geom_errorbar(ggplot2::aes(xmin = ci_low, xmax = ci_high),
                         orientation = "y", width = 0, linewidth = 0.45) +
  ggplot2::geom_point(size = 1.7) +
  ggplot2::scale_colour_manual(values = c("BH-FDR < 0.05" = "#B84A4A",
                                          "Not FDR-significant" = "#767676")) +
  ggplot2::labs(title = "Module association with TB-DM", subtitle = "Linear model adjusted for study site",
                x = "Adjusted eigengene difference", y = NULL, colour = NULL)

theme_wgcna <- ggplot2::theme_classic(base_size = 7, base_family = "Arial") +
  ggplot2::theme(axis.line = ggplot2::element_line(linewidth = 0.3),
    axis.ticks = ggplot2::element_line(linewidth = 0.3), legend.position = "top",
    plot.title = ggplot2::element_text(size = 8, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 6, colour = "#767676"),
    plot.tag = ggplot2::element_text(size = 8, face = "bold"))
fig <- (p_a | p_b) + patchwork::plot_layout(widths = c(0.8, 1.2)) +
  patchwork::plot_annotation(tag_levels = "a", title = "GSE114192 co-expression network",
                             subtitle = "Independent signed WGCNA; 5,000 most variable genes") & theme_wgcna

out <- file.path(result_dir, "figures", "F05_GSE114192_WGCNA")
w <- 183 / 25.4; h <- 92 / 25.4
svglite::svglite(paste0(out, ".svg"), width = w, height = h); print(fig); grDevices::dev.off()
assert_english_only_svg(paste0(out, ".svg"))
grDevices::cairo_pdf(paste0(out, ".pdf"), width = w, height = h, family = "Arial"); print(fig); grDevices::dev.off()
ragg::agg_tiff(paste0(out, ".tiff"), width = w, height = h, units = "in", res = 600,
               compression = "lzw"); print(fig); grDevices::dev.off()
ragg::agg_png(paste0(out, ".png"), width = w, height = h, units = "in", res = 300);
print(fig); grDevices::dev.off()

write_log("GSE114192 WGCNA completed: power=", soft_power, "; modules=", nrow(assoc),
          "; FDR-significant module associations=", sum(assoc$fdr < 0.05),
          "; English-only SVG QA passed")
