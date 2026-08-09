source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "DESeq2", "SummarizedExperiment", "WGCNA",
                  "ggplot2", "svglite", "ragg"))
suppressPackageStartupMessages(library(WGCNA))
WGCNA::allowWGCNAThreads(nThreads = 2)

tb_model_file <- file.path(path_result, "03_cross_disease", "WGCNA", "GSE114192",
                           "models", "M01_GSE114192_WGCNA_network.rds")
sil_model_file <- file.path(path_result, "03_cross_disease", "WGCNA", "GSE165489",
                            "models", "M01_GSE165489_WGCNA_network.rds")
if (!all(file.exists(c(tb_model_file, sil_model_file)))) stop("Independent WGCNA networks missing")
tb_net <- readRDS(tb_model_file)
sil_net <- readRDS(sil_model_file)

# Reconstruct the exact GSE114192 selected-gene expression matrix.
tb_primary <- readRDS(file.path(path_result, "02_bulk", "GSE114192", "models",
                                "M01_GSE114192_primary_models.rds"))
tb_vsd <- DESeq2::vst(tb_primary$dds, blind = FALSE)
tb_expr_all <- SummarizedExperiment::assay(tb_vsd)
tb_mad <- apply(tb_expr_all, 1, stats::mad, na.rm = TRUE)
tb_selected <- names(sort(tb_mad[is.finite(tb_mad) & tb_mad > 0], decreasing = TRUE))[1:5000]
if (!identical(tb_selected, tb_net$genes)) stop("GSE114192 WGCNA gene reconstruction mismatch")
tb_dat <- as.data.frame(t(tb_expr_all[tb_selected, , drop = FALSE]))

tb_map <- unique(data.table::fread(file.path(path_result, "03_cross_disease", "WGCNA",
  "GSE114192", "source_data", "SD02_GSE114192_module_membership.csv"))[, .(gene_id, gene_symbol)])
sil_map <- unique(data.table::fread(file.path(path_data, "02_GEO_bulk", "GSE165489", "processed",
  "GSE165489_gene_map.csv"))[, .(entrez_id, gene_symbol)])

standardize_symbols <- function(dat_expr, gene_ids, colors, map_ids, map_symbols) {
  if (length(colors) != ncol(dat_expr)) stop("Module color length mismatch")
  info <- data.table::data.table(index = seq_along(gene_ids), gene_id = gene_ids,
    gene_symbol = map_symbols[match(gene_ids, map_ids)],
    variability = apply(dat_expr, 2, stats::mad, na.rm = TRUE), color = as.character(colors))
  info <- info[!is.na(gene_symbol) & gene_symbol != ""]
  data.table::setorder(info, gene_symbol, -variability)
  info <- info[, .SD[1], by = gene_symbol]
  z <- dat_expr[, info$index, drop = FALSE]
  colnames(z) <- info$gene_symbol
  cvec <- info$color; names(cvec) <- info$gene_symbol
  list(data = z, colors = cvec, info = info)
}
tb_std <- standardize_symbols(tb_dat, tb_net$genes, tb_net$module_colors,
                              tb_map$gene_id, tb_map$gene_symbol)
sil_colors <- WGCNA::labels2colors(sil_net$network$colors)
sil_std <- standardize_symbols(sil_net$dat_expr, sil_net$genes, sil_colors,
                               sil_map$entrez_id, sil_map$gene_symbol)
common <- intersect(colnames(tb_std$data), colnames(sil_std$data))
if (length(common) < 2000L) stop("Too few common WGCNA genes for preservation: ", length(common))
multi_expr <- list(TB_DM = list(data = tb_std$data[, common, drop = FALSE]),
                   Silicosis = list(data = sil_std$data[, common, drop = FALSE]))
multi_color <- list(TB_DM = tb_std$colors[common], Silicosis = sil_std$colors[common])

result_dir <- file.path(path_result, "03_cross_disease", "WGCNA", "module_preservation")
invisible(lapply(file.path(result_dir, c("tables", "source_data", "models", "figures")),
                 dir.create, recursive = TRUE, showWarnings = FALSE))
mp_file <- file.path(result_dir, "models", "M01_bidirectional_module_preservation.rds")
if (file.exists(mp_file)) {
  mp <- readRDS(mp_file)
} else {
  set.seed(20260808)
  mp <- WGCNA::modulePreservation(multi_expr, multi_color, referenceNetworks = 1:2,
    nPermutations = 100, randomSeed = 20260808, quickCor = 0,
    networkType = "signed", verbose = 2)
  saveRDS(mp, mp_file, compress = FALSE)
}

extract_preservation <- function(component, statistic_type) {
  out <- list()
  for (ref_name in names(component)) {
    for (test_name in names(component[[ref_name]])) {
      z <- data.table::as.data.table(component[[ref_name]][[test_name]], keep.rownames = "module")
      z[, `:=`(reference_network = ref_name, test_network = test_name,
               statistic_type = statistic_type)]
      out[[paste(ref_name, test_name, sep = "__")]] <- z
    }
  }
  data.table::rbindlist(out, fill = TRUE)
}
z_tab <- extract_preservation(mp$preservation$Z, "Z")
obs_tab <- extract_preservation(mp$preservation$observed, "observed")
z_col <- grep("^Zsummary", names(z_tab), value = TRUE)[1]
rank_col <- grep("^medianRank", names(obs_tab), value = TRUE)[1]
if (is.na(z_col) || is.na(rank_col)) stop("Expected module-preservation statistics not found")
summary_z <- z_tab[, .(module, reference_network, test_network,
                       module_size = get(grep("^moduleSize", names(z_tab), value = TRUE)[1]),
                       z_summary = get(z_col))]
summary_rank <- obs_tab[, .(module, reference_network, test_network,
                           median_rank = get(rank_col))]
pres <- merge(summary_z, summary_rank,
              by = c("module", "reference_network", "test_network"), all = TRUE)
pres <- pres[!is.na(module) & module != "" & is.finite(z_summary) &
               !grepl("gold|grey", module, ignore.case = TRUE)]
pres[, preservation_class := data.table::fcase(
  z_summary > 10, "Strong (Z > 10)", z_summary >= 2, "Moderate (2 <= Z <= 10)",
  default = "Weak (Z < 2)"
)]

safe_write_csv(pres, file.path(result_dir, "tables", "T01_bidirectional_module_preservation.csv"))
safe_write_csv(data.table::data.table(gene_symbol = common),
               file.path(result_dir, "source_data", "SD01_common_WGCNA_genes.csv"))

pres[, direction := ifelse(grepl("TB_DM", reference_network),
  "TB-DM modules tested in silicosis", "Silicosis modules tested in TB-DM")]
pres[, module_label := factor(paste0(module, " (n=", module_size, ")"),
  levels = rev(unique(paste0(module, " (n=", module_size, ")"))))]
p <- ggplot2::ggplot(pres, ggplot2::aes(z_summary, module_label, colour = preservation_class)) +
  ggplot2::geom_vline(xintercept = c(2, 10), linetype = 2, linewidth = 0.3, colour = "#767676") +
  ggplot2::geom_point(size = 2) +
  ggplot2::facet_wrap(~ direction, scales = "free_y", ncol = 2) +
  ggplot2::scale_colour_manual(values = c("Weak (Z < 2)" = "#AFAFAF",
    "Moderate (2 <= Z <= 10)" = "#416A9A", "Strong (Z > 10)" = "#B84A4A")) +
  ggplot2::labs(title = "Bidirectional co-expression module preservation",
                subtitle = paste("100 permutations;", format(length(common), big.mark = ","),
                                 "common genes; lower median rank is better"),
                x = "Z-summary preservation statistic", y = NULL, colour = NULL) +
  ggplot2::theme_classic(base_size = 7, base_family = "Arial") +
  ggplot2::theme(axis.line = ggplot2::element_line(linewidth = 0.3),
    axis.ticks = ggplot2::element_line(linewidth = 0.3), legend.position = "top",
    strip.text = ggplot2::element_text(size = 7, face = "bold"),
    plot.title = ggplot2::element_text(size = 8, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 6, colour = "#767676"))
out <- file.path(result_dir, "figures", "F11_bidirectional_WGCNA_module_preservation")
w <- 183 / 25.4; h <- 90 / 25.4
svglite::svglite(paste0(out, ".svg"), width = w, height = h); print(p); grDevices::dev.off()
assert_english_only_svg(paste0(out, ".svg"))
grDevices::cairo_pdf(paste0(out, ".pdf"), width = w, height = h, family = "Arial"); print(p); grDevices::dev.off()
ragg::agg_tiff(paste0(out, ".tiff"), width = w, height = h, units = "in", res = 600,
               compression = "lzw"); print(p); grDevices::dev.off()
ragg::agg_png(paste0(out, ".png"), width = w, height = h, units = "in", res = 300);
print(p); grDevices::dev.off()
write_log("Bidirectional WGCNA preservation completed: common genes=", length(common),
          "; modules tested=", nrow(pres), "; English-only SVG QA passed")
