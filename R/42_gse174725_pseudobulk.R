source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "SingleCellExperiment", "SummarizedExperiment", "Matrix",
  "edgeR", "ggplot2", "patchwork", "svglite", "ragg"))

acc <- "GSE174725"
result_dir <- file.path(path_result, "04_scRNA", acc)
sce_file <- file.path(path_data, "03_scRNA", acc, "processed", paste0(acc, "_sce_qc.rds"))
ann_file <- file.path(result_dir, "source_data", "SD03_GSE174725_cell_annotation_UMAP.csv")
sce <- readRDS(sce_file)
ann <- data.table::fread(ann_file)
idx <- match(colnames(sce), ann$cell_id)
if (anyNA(idx)) stop("Annotation and SCE cell identifiers do not match")
ann <- ann[idx]
if (!identical(colnames(sce), ann$cell_id)) stop("Cell ordering failed")

min_cells <- 20L
donors <- unique(ann$sample_id)
cell_n <- ann[, .N, by = .(cell_type, sample_id)]
eligible <- cell_n[, .(n_donors = .N, minimum_cells = min(N)), by = cell_type][
  n_donors == length(donors) & minimum_cells >= 20L, cell_type]
if (!length(eligible)) stop("No cell type meets the donor and cell-count threshold")

counts <- SummarizedExperiment::assay(sce, "counts")
group_by_sample <- unique(ann[, .(sample_id, group)])
group_by_sample[, group := factor(group, levels = c("Silica-exposed control", "Silicosis"))]
all_res <- list(); audit <- list()
for (ct in eligible) {
  keep_cells <- ann$cell_type == ct
  sample_factor <- factor(ann$sample_id[keep_cells], levels = donors)
  mm <- Matrix::sparse.model.matrix(~ 0 + sample_factor)
  colnames(mm) <- donors
  pb <- counts[, keep_cells, drop = FALSE] %*% mm
  pb <- round(as.matrix(pb))
  meta <- group_by_sample[match(colnames(pb), sample_id)]
  design <- stats::model.matrix(~ group, data = meta)
  y <- edgeR::DGEList(counts = pb, samples = as.data.frame(meta))
  keep_gene <- edgeR::filterByExpr(y, design = design)
  y <- edgeR::calcNormFactors(y[keep_gene, , keep.lib.sizes = FALSE])
  y <- edgeR::estimateDisp(y, design, robust = TRUE)
  fit <- edgeR::glmQLFit(y, design, robust = TRUE)
  qlf <- edgeR::glmQLFTest(fit, coef = "groupSilicosis")
  tab <- data.table::as.data.table(edgeR::topTags(qlf, n = Inf, sort.by = "none")$table, keep.rownames = "gene_symbol")
  data.table::setnames(tab, c("logFC", "logCPM", "F", "PValue", "FDR"),
    c("log2FC", "logCPM", "F_statistic", "p_value", "FDR"), skip_absent = TRUE)
  tab[, cell_type := ct]
  tab[, direction := data.table::fcase(FDR < 0.05 & log2FC > 0, "Higher in silicosis",
    FDR < 0.05 & log2FC < 0, "Lower in silicosis", default = "Not FDR-significant")]
  all_res[[ct]] <- tab
  audit[[ct]] <- data.table::data.table(cell_type = ct, donors = ncol(pb),
    control_donors = sum(meta$group == "Silica-exposed control"),
    silicosis_donors = sum(meta$group == "Silicosis"), genes_tested = nrow(tab),
    FDR_genes = sum(tab$FDR < 0.05), nominal_p05_genes = sum(tab$p_value < 0.05),
    minimum_cells_per_donor = min(cell_n[cell_type == ct]$N))
}
res <- data.table::rbindlist(all_res)
audit_dt <- data.table::rbindlist(audit)
res[, global_FDR := stats::p.adjust(p_value, method = "BH")]
safe_write_csv(res, file.path(result_dir, "source_data", "SD05_GSE174725_pseudobulk_edgeR.csv"))
safe_write_csv(audit_dt, file.path(result_dir, "tables", "T08_GSE174725_pseudobulk_summary.csv"))
safe_write_csv(cell_n, file.path(result_dir, "tables", "T09_GSE174725_cells_per_type_donor.csv"))

plot_types <- audit_dt[order(cell_type)]$cell_type
res[, neglog10P := -log10(pmax(p_value, 1e-300))]
res[, label := FALSE]
res[, label := data.table::frank(p_value, ties.method = "first") <= 3, by = cell_type]
res[, evidence := ifelse(p_value < 0.05, "Nominal P < 0.05", "P >= 0.05")]
res[, facet_label := paste0(cell_type, "\n0 FDR; ",
  audit_dt$nominal_p05_genes[match(cell_type, audit_dt$cell_type)], " nominal")]
cols <- c("Nominal P < 0.05" = "#D08A3E", "P >= 0.05" = "#B8B8B8")
p_a <- ggplot2::ggplot(res, ggplot2::aes(log2FC, neglog10P, colour = evidence)) +
  ggplot2::geom_point(size = 0.45, alpha = 0.60) +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.3) +
  ggplot2::scale_colour_manual(values = cols) +
  ggplot2::facet_wrap(~facet_label, scales = "free_y", nrow = 2) +
  ggplot2::labs(title = "Donor-level pseudobulk differential expression",
    subtitle = "edgeR quasi-likelihood; two exposed controls versus three silicosis donors",
    x = "log2 fold change", y = "-log10 P", colour = NULL) +
  ggplot2::theme(legend.position = "bottom")
p_b <- ggplot2::ggplot(audit_dt, ggplot2::aes(reorder(cell_type, nominal_p05_genes), nominal_p05_genes)) +
  ggplot2::geom_col(width = 0.65, fill = "#D08A3E") +
  ggplot2::coord_flip() +
  ggplot2::geom_text(ggplot2::aes(label = nominal_p05_genes), hjust = -0.15, size = 2.2) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
  ggplot2::labs(title = "Nominal signals only", subtitle = "Unadjusted P < 0.05; no FDR hits",
    x = NULL, y = "Genes")
theme_sc <- ggplot2::theme_classic(base_size = 7, base_family = "Arial") +
  ggplot2::theme(axis.line = ggplot2::element_line(linewidth = 0.3),
    axis.ticks = ggplot2::element_line(linewidth = 0.3),
    legend.key.width = grid::unit(4, "mm"), legend.text = ggplot2::element_text(size = 5.5),
    strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_text(size = 6.5, face = "bold"),
    plot.title = ggplot2::element_text(size = 8, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 6, colour = "#767676"),
    plot.tag = ggplot2::element_text(size = 8, face = "bold"))
fig <- (p_a | p_b) + patchwork::plot_layout(widths = c(2.3, 1)) +
  patchwork::plot_annotation(tag_levels = "a", title = "GSE174725 donor-level cell-type analysis",
    subtitle = "Exploratory inference; donor is the biological replicate") & theme_sc
out <- file.path(result_dir, "figures", "F13_GSE174725_pseudobulk_DE")
w <- 183 / 25.4; h <- 112 / 25.4
svglite::svglite(paste0(out, ".svg"), width = w, height = h); print(fig); grDevices::dev.off()
assert_english_only_svg(paste0(out, ".svg"))
grDevices::cairo_pdf(paste0(out, ".pdf"), width = w, height = h, family = "Arial"); print(fig); grDevices::dev.off()
ragg::agg_tiff(paste0(out, ".tiff"), width = w, height = h, units = "in", res = 600,
  compression = "lzw"); print(fig); grDevices::dev.off()
ragg::agg_png(paste0(out, ".png"), width = w, height = h, units = "in", res = 300)
print(fig); grDevices::dev.off()
write_log("GSE174725 pseudobulk completed: eligible cell types=", paste(eligible, collapse = ", "),
  "; English-only SVG QA passed")
