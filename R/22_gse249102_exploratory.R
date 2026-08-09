source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "limma", "ggplot2", "patchwork", "svglite", "ragg"))

acc <- "GSE249102"
matrix_file <- file.path(path_data, "02_GEO_bulk", acc, "metadata",
                         paste0(acc, "_series_matrix.txt.gz"))
metadata_file <- file.path(path_result, "02_bulk", acc,
                           paste0("A_", acc, "_sample_metadata_raw.csv"))
if (!file.exists(matrix_file) || !file.exists(metadata_file)) stop("GSE249102 inputs missing")

con <- gzfile(matrix_file, "rt")
lines <- readLines(con, warn = FALSE, encoding = "UTF-8")
close(con)
a <- match("!series_matrix_table_begin", lines)
b <- match("!series_matrix_table_end", lines)
if (is.na(a) || is.na(b) || b <= a + 2L) stop("No expression table in GSE249102 series matrix")
tab <- data.table::fread(text = paste(lines[(a + 1L):(b - 1L)], collapse = "\n"),
                         sep = "\t", header = TRUE, check.names = FALSE)
probe_id <- as.character(tab[[1]])
expr <- as.matrix(tab[, -1, with = FALSE])
storage.mode(expr) <- "double"
rownames(expr) <- probe_id

meta <- data.table::fread(metadata_file, encoding = "UTF-8")
if (!all(meta$sample_id %in% colnames(expr))) stop("Expression/metadata sample mismatch")
expr <- expr[, meta$sample_id, drop = FALSE]
meta[, group := factor(group, levels = c("CTRL", "TB", "TBDM", "DM2", "PDM2"))]
if (anyNA(meta$group)) stop("Unmapped GSE249102 group")

design <- stats::model.matrix(~ 0 + group, data = meta)
colnames(design) <- sub("^group", "", colnames(design))
fit <- limma::lmFit(expr, design)
contrasts <- limma::makeContrasts(
  TBDM_vs_TB = TBDM - TB,
  PDM2_vs_DM2 = PDM2 - DM2,
  levels = design
)
fit2 <- limma::eBayes(limma::contrasts.fit(fit, contrasts), trend = TRUE, robust = TRUE)

result_dir <- file.path(path_result, "02_bulk", acc)
dir.create(file.path(result_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(result_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(result_dir, "source_data"), recursive = TRUE, showWarnings = FALSE)

de_list <- list()
for (contrast in colnames(contrasts)) {
  tt <- data.table::as.data.table(limma::topTable(fit2, coef = contrast, number = Inf,
                                                  sort.by = "P", adjust.method = "BH"), keep.rownames = "probe_id")
  tt[, contrast := contrast]
  tt[, evidence_class := data.table::fcase(
    adj.P.Val < 0.05 & abs(logFC) >= 1, "FDR<0.05 and |logFC|>=1",
    adj.P.Val < 0.05, "FDR<0.05",
    default = "Not FDR-significant")]
  safe_write_csv(tt, file.path(result_dir, "tables", paste0("T_", contrast, "_full.csv")))
  de_list[[contrast]] <- tt
}

pca <- stats::prcomp(t(expr), center = TRUE, scale. = FALSE)
var_exp <- 100 * pca$sdev^2 / sum(pca$sdev^2)
pca_df <- data.frame(sample_id = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2])
pca_df <- merge(pca_df, meta[, .(sample_id, group)], by = "sample_id", sort = FALSE)

palette_groups <- c(CTRL = "#767676", TB = "#416A9A", TBDM = "#B84A4A",
                    DM2 = "#7A9E6D", PDM2 = "#C49A55")
theme_nature <- ggplot2::theme_classic(base_size = 7, base_family = "Arial") +
  ggplot2::theme(axis.line = ggplot2::element_line(linewidth = 0.3),
                 axis.ticks = ggplot2::element_line(linewidth = 0.3),
                 legend.position = "top", legend.title = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(size = 8, face = "bold"),
                 plot.subtitle = ggplot2::element_text(size = 6, colour = "#767676"),
                 plot.tag = ggplot2::element_text(size = 8, face = "bold"))

p_a <- ggplot2::ggplot(pca_df, ggplot2::aes(PC1, PC2, colour = group)) +
  ggplot2::geom_hline(yintercept = 0, colour = "#D8D8D8", linewidth = 0.25) +
  ggplot2::geom_vline(xintercept = 0, colour = "#D8D8D8", linewidth = 0.25) +
  ggplot2::geom_point(size = 2, alpha = 0.9) +
  ggplot2::scale_colour_manual(values = palette_groups, drop = FALSE) +
  ggplot2::labs(title = "Expression-space overview",
                subtitle = "GSE249102; n = 4 per group",
                x = sprintf("PC1 (%.1f%%)", var_exp[1]), y = sprintf("PC2 (%.1f%%)", var_exp[2])) + theme_nature

volcano_plot <- function(tt, title) {
  tt <- data.table::copy(tt)
  tt[, neglog10_fdr := -log10(pmax(adj.P.Val, .Machine$double.xmin))]
  ggplot2::ggplot(tt, ggplot2::aes(logFC, neglog10_fdr, colour = evidence_class)) +
    ggplot2::geom_point(size = 0.45, alpha = 0.55) +
    ggplot2::geom_vline(xintercept = c(-1, 1), linetype = 2, linewidth = 0.3, colour = "#767676") +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.3, colour = "#767676") +
    ggplot2::scale_colour_manual(values = c("FDR<0.05 and |logFC|>=1" = "#B84A4A",
      "FDR<0.05" = "#416A9A", "Not FDR-significant" = "#C9C9C9"), drop = FALSE) +
    ggplot2::labs(title = title,
                  subtitle = "Red: FDR<0.05 & |logFC|>=1; blue: FDR<0.05",
                  x = "log2 fold change", y = expression(-log[10](FDR))) +
    theme_nature + ggplot2::theme(legend.position = "none")
}
p_b <- volcano_plot(de_list$TBDM_vs_TB, "TBDM versus TB")
p_c <- volcano_plot(de_list$PDM2_vs_DM2, "Poorly controlled versus DM")

fig <- (p_a | p_b | p_c) + patchwork::plot_annotation(tag_levels = "a",
  title = "Exploratory glycaemic-gradient validation",
  subtitle = "Small cohort (n = 20); effect direction and pathway/module validation only")

source_data <- data.table::rbindlist(de_list, fill = TRUE)
safe_write_csv(source_data, file.path(result_dir, "source_data", "SD_GSE249102_limma_results.csv"))
safe_write_csv(pca_df, file.path(result_dir, "source_data", "SD_GSE249102_PCA.csv"))

out <- file.path(result_dir, "figures", "F02_GSE249102_exploratory")
w <- 183 / 25.4; h <- 65 / 25.4
svglite::svglite(paste0(out, ".svg"), width = w, height = h); print(fig); grDevices::dev.off()
assert_english_only_svg(paste0(out, ".svg"))
grDevices::cairo_pdf(paste0(out, ".pdf"), width = w, height = h, family = "Arial"); print(fig); grDevices::dev.off()
ragg::agg_tiff(paste0(out, ".tiff"), width = w, height = h, units = "in", res = 600, compression = "lzw"); print(fig); grDevices::dev.off()
ragg::agg_png(paste0(out, ".png"), width = w, height = h, units = "in", res = 300); print(fig); grDevices::dev.off()

saveRDS(list(fit = fit2, design = design, contrasts = contrasts),
        file.path(result_dir, "GSE249102_limma_model.rds"))
write_log("GSE249102 exploratory analysis completed; English-only SVG QA passed")
