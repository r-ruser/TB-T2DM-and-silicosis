source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "SingleCellExperiment", "SummarizedExperiment", "Matrix",
  "scran", "batchelor", "BiocParallel", "uwot", "igraph", "ggplot2", "patchwork",
  "svglite", "ragg"))

acc <- "GSE192483"
sce_file <- file.path(path_data, "03_scRNA", acc, "processed", paste0(acc, "_sce_qc.rds"))
result_dir <- file.path(path_result, "04_scRNA", acc)
if (!file.exists(sce_file)) stop("Run R/43_gse192483_prepare_qc.R first")
sce <- readRDS(sce_file)
if (ncol(sce) < 50000L || data.table::uniqueN(sce$sample_id) != 11L) stop("Unexpected QC object")

cache_file <- file.path(result_dir, "models", "GSE192483_mnn_umap_cluster.rds")
if (file.exists(cache_file)) {
  write_log("GSE192483 annotation: loading cached MNN/UMAP/clusters")
  cache <- readRDS(cache_file)
  hvg <- cache$hvg; cell_order <- cache$cell_order; umap <- cache$umap; cluster <- cache$cluster
} else {
  write_log("GSE192483 annotation: sample-blocked HVG modeling")
  dec <- scran::modelGeneVar(sce, block = factor(sce$sample_id), BPPARAM = BiocParallel::SerialParam())
  hvg <- scran::getTopHVGs(dec, n = 2500)
  split_sce <- lapply(unique(sce$sample_id), function(sid) sce[, sce$sample_id == sid])
  names(split_sce) <- unique(sce$sample_id)
  set.seed(20260808)
  mnn_args <- c(unname(split_sce), list(subset.row = hvg, d = 30, k = 20,
    cos.norm = TRUE, auto.merge = TRUE, BPPARAM = BiocParallel::SerialParam()))
  mnn <- do.call(batchelor::fastMNN, mnn_args)
  corrected <- SingleCellExperiment::reducedDim(mnn, "corrected")
  cell_order <- match(colnames(mnn), colnames(sce))
  if (anyNA(cell_order)) stop("fastMNN cell ordering failed")
  set.seed(20260808)
  umap <- uwot::umap(corrected, n_neighbors = 30, min_dist = 0.30, metric = "cosine",
    n_threads = 4, verbose = TRUE, ret_model = FALSE)
  rownames(umap) <- colnames(mnn)
  g <- scran::buildSNNGraph(mnn, use.dimred = "corrected", k = 20,
    BPPARAM = BiocParallel::SerialParam())
  cluster <- as.character(igraph::membership(igraph::cluster_louvain(g)))
  names(cluster) <- colnames(mnn)
  saveRDS(list(hvg = hvg, cell_order = cell_order, umap = umap, cluster = cluster),
    cache_file, compress = FALSE)
}
cell_ids <- names(cluster)

cluster_factor <- factor(cluster, levels = sort(unique(cluster)))
design_cluster <- Matrix::sparse.model.matrix(~ 0 + cluster_factor)
colnames(design_cluster) <- levels(cluster_factor)
log_expr <- SummarizedExperiment::assay(sce, "logcounts")[, cell_order, drop = FALSE]
cluster_sum <- log_expr %*% design_cluster
cluster_mean <- sweep(as.matrix(cluster_sum), 2, as.numeric(table(cluster_factor)), "/")
symbol <- as.character(SummarizedExperiment::rowData(sce)$gene_symbol)
symbol_to_row <- split(seq_along(symbol), symbol)

marker_sets <- list(
  Macrophage = c("C1QA", "C1QB", "C1QC", "APOC1", "MARCO", "FABP4", "PPARG", "MRC1", "MSR1", "FCER1G"),
  Monocyte = c("S100A8", "S100A9", "S100A10", "S100A12", "FCN1", "VCAN", "CCR2", "CTSS", "LILRB1"),
  Neutrophil = c("FCGR3B", "CSF3R", "CXCR2", "FPR1", "S100A8", "S100A9", "MNDA"),
  Dendritic = c("FCER1A", "CD1C", "CLEC10A", "CLEC9A", "CST3"),
  `CD4 T` = c("CD3D", "CD3E", "TRAC", "IL7R", "LTB", "MAL", "CCR7"),
  `CD8 T` = c("CD3D", "CD3E", "TRAC", "CD8A", "CD8B", "CCL5", "GZMK"),
  NK = c("NKG7", "GNLY", "KLRD1", "PRF1", "CTSW", "GZMB"),
  B = c("MS4A1", "CD79A", "CD37", "CD74", "HLA-DRA", "CD22"),
  Plasma = c("MZB1", "JCHAIN", "SDC1", "DERL3", "IGKC"),
  Epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT5", "KRT17"),
  Endothelial = c("PECAM1", "VWF", "EMCN", "KDR", "CLDN5", "RAMP2"),
  Fibroblast = c("COL1A1", "COL1A2", "DCN", "COL3A1", "LUM", "PDGFRA"),
  Mast = c("TPSAB1", "TPSB2", "KIT", "CPA3", "MS4A2"),
  Platelet = c("PPBP", "PF4", "GP9", "GP6", "CLEC1B", "PEAR1")
)
marker_ids <- lapply(marker_sets, function(gs) {
  rows <- unique(unlist(symbol_to_row[intersect(gs, names(symbol_to_row))]))
  rownames(sce)[rows]
})
if (any(lengths(marker_ids) < 3L)) stop("Too few markers for: ",
  paste(names(marker_ids)[lengths(marker_ids) < 3L], collapse = ", "))
marker_union <- unique(unlist(marker_ids))
marker_z <- t(scale(t(cluster_mean[marker_union, , drop = FALSE])))
marker_z[!is.finite(marker_z)] <- 0
score <- vapply(marker_ids, function(ids) colMeans(marker_z[ids, , drop = FALSE]),
  numeric(ncol(marker_z)))
score <- t(score); rownames(score) <- names(marker_ids); colnames(score) <- colnames(cluster_mean)
cluster_type <- apply(score, 2, function(z) names(which.max(z)))
score_sorted <- apply(score, 2, sort, decreasing = TRUE)
score_margin <- score_sorted[1, ] - score_sorted[2, ]

gene_z <- t(scale(t(cluster_mean)))
gene_z[!is.finite(gene_z)] <- -Inf
top_markers <- data.table::rbindlist(lapply(colnames(gene_z), function(cl) {
  candidate <- which(rowMeans(cluster_mean) > 0.05)
  ord <- candidate[order(gene_z[candidate, cl], cluster_mean[candidate, cl], decreasing = TRUE)][1:20]
  data.table::data.table(cluster = cl, rank = seq_along(ord), gene_id = rownames(gene_z)[ord],
    gene_symbol = symbol[match(rownames(gene_z)[ord], rownames(sce))], specificity_z = gene_z[ord, cl])
}))
score_dt <- data.table::as.data.table(as.table(score))
data.table::setnames(score_dt, c("cell_type", "cluster", "marker_score"))
score_dt[, assigned_cell_type := cluster_type[as.character(cluster)]]
score_dt[, assignment_margin := score_margin[as.character(cluster)]]

cell_meta <- data.table::data.table(cell_id = cell_ids, UMAP1 = umap[, 1], UMAP2 = umap[, 2],
  cluster = cluster, cell_type = unname(cluster_type[cluster]),
  annotation_method = "Canonical marker score with top-marker audit")
base_meta <- data.table::as.data.table(as.data.frame(SummarizedExperiment::colData(sce)))
cell_meta <- merge(cell_meta, base_meta[, .(cell_id, sample_id, sample_code, patient, region, paired_patient)],
  by = "cell_id", all.x = TRUE, sort = FALSE)
if (anyNA(cell_meta$sample_id)) stop("Annotation metadata merge failed")
abundance <- cell_meta[, .(cells = .N), by = .(sample_id, sample_code, patient, region, cell_type)]
abundance[, percent := 100 * cells / sum(cells), by = sample_id]

safe_write_csv(cell_meta, file.path(result_dir, "source_data", "SD03_GSE192483_cell_annotation_UMAP.csv"))
safe_write_csv(data.table::data.table(gene_id = hvg,
  gene_symbol = symbol[match(hvg, rownames(sce))]), file.path(result_dir, "source_data", "SD04_GSE192483_HVGs.csv"))
safe_write_csv(score_dt, file.path(result_dir, "tables", "T04_GSE192483_cluster_marker_scores.csv"))
safe_write_csv(top_markers, file.path(result_dir, "tables", "T05_GSE192483_cluster_top_markers.csv"))
safe_write_csv(abundance, file.path(result_dir, "tables", "T06_GSE192483_cell_type_abundance.csv"))

cell_type_order <- c("Macrophage", "Monocyte", "Neutrophil", "Dendritic", "CD4 T", "CD8 T", "NK", "B",
  "Plasma", "Epithelial", "Endothelial", "Fibroblast", "Mast", "Platelet")
pal <- c("Macrophage" = "#3E6D8E", "Monocyte" = "#6A9FB5", "Neutrophil" = "#B84A4A",
  "Dendritic" = "#7B6FA6", "CD4 T" = "#58A27C", "CD8 T" = "#82B36A", "NK" = "#D08A3E",
  "B" = "#4B8C86", "Plasma" = "#A8648A", "Epithelial" = "#C7A35A",
  "Endothelial" = "#4D7C8A", "Fibroblast" = "#9B7E62", "Mast" = "#8B6E55", "Platelet" = "#9C755F")
cell_meta[, cell_type := factor(cell_type, levels = cell_type_order)]
abundance[, cell_type := factor(cell_type, levels = cell_type_order)]
sample_levels <- unique(cell_meta$sample_code)[order(match(unique(cell_meta$sample_code),
  c("SP019H", "SP019L", "SP020L", "SP021H", "SP021L", "SP023H", "SP023L", "SP024H", "SP024L", "SP025H", "SP025L")))]
sample_pal <- stats::setNames(c("#B84A4A", "#6A9FB5", "#9BC1D0", "#D1736E", "#477C9B", "#A94B4B",
  "#7AA6B8", "#C86460", "#557F9F", "#9E5A58", "#86B2C3"), sample_levels)

p_a <- ggplot2::ggplot(cell_meta, ggplot2::aes(UMAP1, UMAP2, colour = cell_type)) +
  ggplot2::geom_point(size = 0.09, alpha = 0.55) + ggplot2::scale_colour_manual(values = pal, drop = TRUE) +
  ggplot2::guides(colour = ggplot2::guide_legend(nrow = 3, byrow = TRUE,
    override.aes = list(size = 1.5, alpha = 1))) +
  ggplot2::labs(title = "Lung cell identities", subtitle = "Marker-score annotation",
    x = "UMAP 1", y = "UMAP 2", colour = NULL)
p_b <- ggplot2::ggplot(cell_meta, ggplot2::aes(UMAP1, UMAP2, colour = sample_code)) +
  ggplot2::geom_point(size = 0.09, alpha = 0.48) + ggplot2::scale_colour_manual(values = sample_pal) +
  ggplot2::guides(colour = ggplot2::guide_legend(nrow = 3, byrow = TRUE,
    override.aes = list(size = 1.5, alpha = 1))) +
  ggplot2::labs(title = "Sample mixing audit", subtitle = "Six patients and 11 tissue samples",
    x = "UMAP 1", y = "UMAP 2", colour = NULL)
p_c <- ggplot2::ggplot(abundance, ggplot2::aes(sample_code, percent, fill = cell_type)) +
  ggplot2::geom_col(width = 0.76) + ggplot2::scale_fill_manual(values = pal, drop = TRUE) +
  ggplot2::scale_x_discrete(limits = sample_levels) +
  ggplot2::labs(title = "Cell-type composition by tissue sample", subtitle = "H: TB lesion; L: less-involved lung",
    x = NULL, y = "Cells (%)", fill = NULL)
theme_sc <- ggplot2::theme_classic(base_size = 7, base_family = "Arial") +
  ggplot2::theme(axis.line = ggplot2::element_line(linewidth = 0.3),
    axis.ticks = ggplot2::element_line(linewidth = 0.3),
    legend.key.height = grid::unit(2.5, "mm"), legend.key.width = grid::unit(3.5, "mm"),
    legend.text = ggplot2::element_text(size = 5.4),
    plot.title = ggplot2::element_text(size = 8, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 6, colour = "#767676"),
    plot.tag = ggplot2::element_text(size = 8, face = "bold"))
p_a <- p_a + theme_sc + ggplot2::theme(legend.position = "bottom", aspect.ratio = 1)
p_b <- p_b + theme_sc + ggplot2::theme(legend.position = "bottom", aspect.ratio = 1)
p_c <- p_c + theme_sc + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 38, hjust = 1),
  legend.position = "none")
fig <- ((p_a | p_b) / p_c) + patchwork::plot_layout(heights = c(2.15, 1)) +
  patchwork::plot_annotation(tag_levels = "a", title = "GSE192483 TB lung cell atlas",
    subtitle = "Patient-resolved lesion and less-involved tissue profiling")
out <- file.path(result_dir, "figures", "F15_GSE192483_cell_atlas")
w <- 183 / 25.4; h <- 164 / 25.4
svglite::svglite(paste0(out, ".svg"), width = w, height = h); print(fig); grDevices::dev.off()
assert_english_only_svg(paste0(out, ".svg"))
grDevices::cairo_pdf(paste0(out, ".pdf"), width = w, height = h, family = "Arial"); print(fig); grDevices::dev.off()
ragg::agg_tiff(paste0(out, ".tiff"), width = w, height = h, units = "in", res = 600,
  compression = "lzw"); print(fig); grDevices::dev.off()
ragg::agg_png(paste0(out, ".png"), width = w, height = h, units = "in", res = 300)
print(fig); grDevices::dev.off()
write_log("GSE192483 annotation completed: cells=", nrow(cell_meta), "; clusters=",
  data.table::uniqueN(cell_meta$cluster), "; cell types=", data.table::uniqueN(cell_meta$cell_type),
  "; English-only SVG QA passed")
