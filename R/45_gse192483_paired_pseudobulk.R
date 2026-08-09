source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "SingleCellExperiment", "SummarizedExperiment", "Matrix", "edgeR",
  "ggplot2", "patchwork", "svglite", "ragg"))

acc <- "GSE192483"
result_dir <- file.path(path_result, "04_scRNA", acc)
sce_file <- file.path(path_data, "03_scRNA", acc, "processed", paste0(acc, "_sce_qc.rds"))
ann_file <- file.path(result_dir, "source_data", "SD03_GSE192483_cell_annotation_UMAP.csv")
sce <- readRDS(sce_file)
ann <- data.table::fread(ann_file)
idx <- match(colnames(sce), ann$cell_id)
if (anyNA(idx)) stop("Annotation and SCE cell identifiers do not match")
ann <- ann[idx]
if (!identical(colnames(sce), ann$cell_id)) stop("Cell ordering failed")

paired_patients <- ann[, data.table::uniqueN(region), by = patient][V1 == 2L, patient]
if (length(paired_patients) != 5L) stop("Expected five paired patients")
ann[, inferential_sample := patient %in% paired_patients &
  !cell_type %in% c("Mixed epithelial-myeloid", "Unresolved")]
sample_meta <- unique(ann[patient %in% paired_patients,
  .(sample_id, sample_code, patient, region)])
data.table::setorder(sample_meta, patient, region)
sample_meta[, patient := factor(patient)]
sample_meta[, region := factor(region, levels = c("Less-involved lung", "TB lesion"))]
if (nrow(sample_meta) != 10L) stop("Paired sample metadata is incomplete")

cell_n <- ann[inferential_sample == TRUE, .N, by = .(cell_type, sample_id)]
eligible <- cell_n[, .(n_samples = .N, minimum_cells = min(N)), by = cell_type][
  n_samples == nrow(sample_meta) & minimum_cells >= 20L, cell_type]
if (!length(eligible)) stop("No cell type meets paired sample and cell-count thresholds")
write_log("GSE192483 eligible paired pseudobulk cell types: ", paste(eligible, collapse = ", "))

counts <- SummarizedExperiment::assay(sce, "counts")
gene_symbol <- as.character(SummarizedExperiment::rowData(sce)$gene_symbol)
all_res <- list(); audit <- list()
for (ct in eligible) {
  keep_cells <- ann$inferential_sample & ann$cell_type == ct
  sample_factor <- factor(ann$sample_id[keep_cells], levels = sample_meta$sample_id)
  mm <- Matrix::sparse.model.matrix(~ 0 + sample_factor)
  colnames(mm) <- sample_meta$sample_id
  pb <- round(as.matrix(counts[, keep_cells, drop = FALSE] %*% mm))
  design <- stats::model.matrix(~ patient + region, data = sample_meta)
  y <- edgeR::DGEList(counts = pb)
  keep_gene <- edgeR::filterByExpr(y, design = design)
  y <- edgeR::calcNormFactors(y[keep_gene, , keep.lib.sizes = FALSE])
  y <- edgeR::estimateDisp(y, design, robust = TRUE)
  fit <- edgeR::glmQLFit(y, design, robust = TRUE)
  qlf <- edgeR::glmQLFTest(fit, coef = "regionTB lesion")
  tab <- data.table::as.data.table(edgeR::topTags(qlf, n = Inf, sort.by = "none")$table,
    keep.rownames = "gene_id")
  data.table::setnames(tab, c("logFC", "logCPM", "F", "PValue", "FDR"),
    c("log2FC", "logCPM", "F_statistic", "p_value", "FDR"), skip_absent = TRUE)
  tab[, gene_symbol := gene_symbol[match(gene_id, rownames(sce))]]
  tab[, cell_type := ct]
  tab[, direction := data.table::fcase(FDR < 0.05 & log2FC > 0, "Higher in TB lesion",
    FDR < 0.05 & log2FC < 0, "Lower in TB lesion", default = "Not FDR-significant")]
  all_res[[ct]] <- tab
  audit[[ct]] <- data.table::data.table(cell_type = ct, paired_patients = length(paired_patients),
    tissue_samples = ncol(pb), genes_tested = nrow(tab), FDR_genes = sum(tab$FDR < 0.05),
    nominal_p05_genes = sum(tab$p_value < 0.05), minimum_cells_per_sample = min(cell_n[cell_type == ct]$N))
}
res <- data.table::rbindlist(all_res)
audit_dt <- data.table::rbindlist(audit)
res[, global_FDR := stats::p.adjust(p_value, method = "BH")]
safe_write_csv(res, file.path(result_dir, "source_data", "SD05_GSE192483_paired_pseudobulk_edgeR.csv"))
safe_write_csv(audit_dt, file.path(result_dir, "tables", "T08_GSE192483_paired_pseudobulk_summary.csv"))
safe_write_csv(cell_n, file.path(result_dir, "tables", "T09_GSE192483_cells_per_type_sample.csv"))

res[, neglog10P := -log10(pmax(p_value, 1e-300))]
res[, evidence := data.table::fcase(FDR < 0.05 & log2FC > 0, "FDR: higher in lesion",
  FDR < 0.05 & log2FC < 0, "FDR: lower in lesion", p_value < 0.05, "Nominal P < 0.05",
  default = "P >= 0.05")]
res[, facet_label := paste0(cell_type, "\n", audit_dt$FDR_genes[match(cell_type, audit_dt$cell_type)], " FDR")]
cols <- c("FDR: higher in lesion" = "#B84A4A", "FDR: lower in lesion" = "#416A9A",
  "Nominal P < 0.05" = "#D08A3E", "P >= 0.05" = "#B8B8B8")
p_a <- ggplot2::ggplot(res, ggplot2::aes(log2FC, neglog10P, colour = evidence)) +
  ggplot2::geom_point(size = 0.34, alpha = 0.58) +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.3) +
  ggplot2::scale_colour_manual(values = cols) +
  ggplot2::facet_wrap(~facet_label, scales = "free_y", ncol = 3) +
  ggplot2::labs(title = "Paired lesion-versus-lung pseudobulk",
    subtitle = "edgeR quasi-likelihood with patient blocking; five paired patients",
    x = "log2 fold change", y = "-log10 P", colour = NULL) +
  ggplot2::theme(legend.position = "bottom")
if (sum(audit_dt$FDR_genes) > 0L) {
  audit_dt[, display_count := FDR_genes]
  bar_title <- "FDR-significant genes"; bar_subtitle <- "BH FDR within each cell type"
  bar_fill <- "#3E6D8E"
} else {
  audit_dt[, display_count := nominal_p05_genes]
  bar_title <- "Nominal signals only"; bar_subtitle <- "Unadjusted P < 0.05; no FDR hits"
  bar_fill <- "#D08A3E"
}
p_b <- ggplot2::ggplot(audit_dt, ggplot2::aes(reorder(cell_type, display_count), display_count)) +
  ggplot2::geom_col(width = 0.65, fill = bar_fill) + ggplot2::coord_flip() +
  ggplot2::geom_text(ggplot2::aes(label = display_count), hjust = -0.12, size = 2.1) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.16))) +
  ggplot2::labs(title = bar_title, subtitle = bar_subtitle, x = NULL, y = "Genes")
theme_sc <- ggplot2::theme_classic(base_size = 7, base_family = "Arial") +
  ggplot2::theme(axis.line = ggplot2::element_line(linewidth = 0.3),
    axis.ticks = ggplot2::element_line(linewidth = 0.3),
    legend.key.width = grid::unit(4, "mm"), legend.text = ggplot2::element_text(size = 5.3),
    strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_text(size = 6.2, face = "bold"),
    plot.title = ggplot2::element_text(size = 8, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 6, colour = "#767676"),
    plot.tag = ggplot2::element_text(size = 8, face = "bold"))
fig <- (p_a | p_b) + patchwork::plot_layout(widths = c(2.5, 1)) +
  patchwork::plot_annotation(tag_levels = "a", title = "GSE192483 paired TB lung analysis",
    subtitle = "Patient is the biological replicate; SP020 excluded from paired inference") & theme_sc
out <- file.path(result_dir, "figures", "F16_GSE192483_paired_pseudobulk_DE")
w <- 183 / 25.4; h <- 142 / 25.4
svglite::svglite(paste0(out, ".svg"), width = w, height = h); print(fig); grDevices::dev.off()
assert_english_only_svg(paste0(out, ".svg"))
grDevices::cairo_pdf(paste0(out, ".pdf"), width = w, height = h, family = "Arial"); print(fig); grDevices::dev.off()
ragg::agg_tiff(paste0(out, ".tiff"), width = w, height = h, units = "in", res = 600,
  compression = "lzw"); print(fig); grDevices::dev.off()
ragg::agg_png(paste0(out, ".png"), width = w, height = h, units = "in", res = 300)
print(fig); grDevices::dev.off()
write_log("GSE192483 paired pseudobulk completed: eligible cell types=", paste(eligible, collapse = ", "),
  "; English-only SVG QA passed")
