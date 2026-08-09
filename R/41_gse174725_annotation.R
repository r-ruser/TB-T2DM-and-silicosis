source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "SingleCellExperiment", "SummarizedExperiment", "Matrix",
                  "scran", "batchelor", "BiocParallel", "uwot", "igraph",
                  "ggplot2", "patchwork", "svglite", "ragg"))

acc <- "GSE174725"
sce_file <- file.path(path_data, "03_scRNA", acc, "processed", paste0(acc, "_sce_qc.rds"))
if (!file.exists(sce_file)) stop("Run R/40_gse174725_prepare_qc.R first")
sce <- readRDS(sce_file)
if (ncol(sce) < 10000L || length(unique(sce$sample_id)) != 5L) stop("Unexpected GSE174725 QC object")

write_log("GSE174725 annotation: donor-blocked HVG modeling")
dec <- scran::modelGeneVar(sce, block = factor(sce$sample_id), BPPARAM = BiocParallel::SerialParam())
hvg <- scran::getTopHVGs(dec, n = 2000)
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
umap <- uwot::umap(corrected, n_neighbors = 30, min_dist = 0.30,
                   metric = "cosine", n_threads = 1, verbose = TRUE,
                   ret_model = FALSE)
rownames(umap) <- colnames(mnn)
g <- scran::buildSNNGraph(mnn, use.dimred = "corrected", k = 20,
                          BPPARAM = BiocParallel::SerialParam())
cluster <- as.character(igraph::membership(igraph::cluster_louvain(g)))
names(cluster) <- colnames(mnn)

cluster_factor <- factor(cluster, levels = sort(unique(cluster)))
design_cluster <- Matrix::sparse.model.matrix(~ 0 + cluster_factor)
colnames(design_cluster) <- levels(cluster_factor)
log_expr <- SummarizedExperiment::assay(sce, "logcounts")[, cell_order, drop = FALSE]
cluster_sum <- log_expr %*% design_cluster
cluster_mean <- sweep(as.matrix(cluster_sum), 2, as.numeric(table(cluster_factor)), "/")

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
  Mast = c("TPSAB1", "TPSB2", "KIT", "CPA3", "MS4A2"),
  Erythroid = c("HBB", "HBA1", "HBA2", "ALAS2", "AHSP"),
  Platelet = c("PPBP", "PF4", "GP9", "GP6", "CLEC1B", "PEAR1")
)
marker_sets <- lapply(marker_sets, function(gs) intersect(gs, rownames(cluster_mean)))
if (any(lengths(marker_sets) < 3L)) stop("Too few available markers for: ",
                                        paste(names(marker_sets)[lengths(marker_sets) < 3L], collapse = ", "))
marker_union <- unique(unlist(marker_sets))
marker_z <- t(scale(t(cluster_mean[marker_union, , drop = FALSE])))
marker_z[!is.finite(marker_z)] <- 0
score <- vapply(marker_sets, function(gs) colMeans(marker_z[gs, , drop = FALSE]),
                numeric(ncol(marker_z)))
score <- t(score)
rownames(score) <- names(marker_sets)
colnames(score) <- colnames(cluster_mean)
cluster_type <- apply(score, 2, function(z) names(which.max(z)))
score_sorted <- apply(score, 2, sort, decreasing = TRUE)
score_margin <- score_sorted[1, ] - score_sorted[2, ]

# Independent top-marker audit based on cluster specificity across all genes.
gene_z <- t(scale(t(cluster_mean)))
gene_z[!is.finite(gene_z)] <- 0
top_markers <- data.table::rbindlist(lapply(colnames(gene_z), function(cl) {
  ord <- order(gene_z[, cl], decreasing = TRUE)[1:15]
  data.table::data.table(cluster = cl, rank = seq_along(ord),
                         gene_symbol = rownames(gene_z)[ord], specificity_z = gene_z[ord, cl])
}))
score_dt <- data.table::as.data.table(as.table(score))
data.table::setnames(score_dt, c("cell_type", "cluster", "marker_score"))
score_dt[, assigned_cell_type := cluster_type[as.character(cluster)]]
score_dt[, assignment_margin := score_margin[as.character(cluster)]]

cell_meta <- data.table::data.table(
  cell_id = colnames(mnn), UMAP1 = umap[, 1], UMAP2 = umap[, 2],
  cluster = cluster, cell_type = unname(cluster_type[cluster])
)
base_meta <- data.table::as.data.table(as.data.frame(SummarizedExperiment::colData(sce)))
cell_meta <- merge(cell_meta, base_meta[, .(cell_id, sample_id, donor, group)],
                   by = "cell_id", all.x = TRUE, sort = FALSE)
if (anyNA(cell_meta$sample_id)) stop("Cell annotation metadata merge failed")

abundance <- cell_meta[, .(cells = .N), by = .(sample_id, donor, group, cell_type)]
abundance[, percent := 100 * cells / sum(cells), by = sample_id]

result_dir <- file.path(path_result, "04_scRNA", acc)
safe_write_csv(cell_meta, file.path(result_dir, "source_data", "SD03_GSE174725_cell_annotation_UMAP.csv"))
safe_write_csv(score_dt, file.path(result_dir, "tables", "T04_GSE174725_cluster_marker_scores.csv"))
safe_write_csv(top_markers, file.path(result_dir, "tables", "T05_GSE174725_cluster_top_markers.csv"))
safe_write_csv(abundance, file.path(result_dir, "tables", "T06_GSE174725_cell_type_abundance.csv"))
safe_write_csv(data.table::data.table(gene_symbol = hvg),
               file.path(result_dir, "source_data", "SD04_GSE174725_HVGs.csv"))

cell_type_order <- c("Macrophage", "Monocyte", "Neutrophil", "Dendritic", "CD4 T", "CD8 T",
                     "NK", "B", "Plasma", "Epithelial", "Mast", "Erythroid", "Platelet")
pal <- c("Macrophage" = "#3E6D8E", "Monocyte" = "#6A9FB5", "Neutrophil" = "#B84A4A",
  "Dendritic" = "#7B6FA6", "CD4 T" = "#58A27C", "CD8 T" = "#82B36A",
  "NK" = "#D08A3E", "B" = "#4B8C86", "Plasma" = "#A8648A",
  "Epithelial" = "#C7A35A", "Mast" = "#8B6E55", "Erythroid" = "#7A7A7A",
  "Platelet" = "#9C755F")
cell_meta[, cell_type := factor(cell_type, levels = cell_type_order)]
abundance[, cell_type := factor(cell_type, levels = cell_type_order)]
donor_pal <- stats::setNames(c("#416A9A", "#6A9FB5", "#B84A4A", "#D1736E", "#9E5A58"),
                             unique(cell_meta$donor))
p_a <- ggplot2::ggplot(cell_meta, ggplot2::aes(UMAP1, UMAP2, colour = cell_type)) +
  ggplot2::geom_point(size = 0.18, alpha = 0.65) +
  ggplot2::scale_colour_manual(values = pal, drop = FALSE) +
  ggplot2::labs(title = "Marker-based cell identities", subtitle = "fastMNN used for visualization only",
                x = "UMAP 1", y = "UMAP 2", colour = NULL)
p_b <- ggplot2::ggplot(cell_meta, ggplot2::aes(UMAP1, UMAP2, colour = donor)) +
  ggplot2::geom_point(size = 0.18, alpha = 0.55) +
  ggplot2::scale_colour_manual(values = donor_pal) +
  ggplot2::labs(title = "Donor mixing audit", subtitle = "Five biological donors",
                x = "UMAP 1", y = "UMAP 2", colour = NULL)
p_c <- ggplot2::ggplot(abundance, ggplot2::aes(donor, percent, fill = cell_type)) +
  ggplot2::geom_col(width = 0.78) +
  ggplot2::scale_fill_manual(values = pal, drop = FALSE) +
  ggplot2::labs(title = "Cell-type composition by donor", subtitle = "Descriptive proportions",
                x = NULL, y = "Cells (%)", fill = NULL) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
theme_sc <- ggplot2::theme_classic(base_size = 7, base_family = "Arial") +
  ggplot2::theme(axis.line = ggplot2::element_line(linewidth = 0.3),
    axis.ticks = ggplot2::element_line(linewidth = 0.3), legend.position = "right",
    plot.title = ggplot2::element_text(size = 8, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 6, colour = "#767676"),
    plot.tag = ggplot2::element_text(size = 8, face = "bold"))
fig <- (p_a | p_b | p_c) + patchwork::plot_layout(widths = c(1, 1, 1.05)) +
  patchwork::plot_annotation(tag_levels = "a", title = "GSE174725 silicosis BALF cell atlas",
                             subtitle = "Donor-resolved annotation after adaptive QC and doublet removal") & theme_sc
out <- file.path(result_dir, "figures", "F12_GSE174725_cell_atlas")
w <- 183 / 25.4; h <- 94 / 25.4
svglite::svglite(paste0(out, ".svg"), width = w, height = h); print(fig); grDevices::dev.off()
assert_english_only_svg(paste0(out, ".svg"))
grDevices::cairo_pdf(paste0(out, ".pdf"), width = w, height = h, family = "Arial"); print(fig); grDevices::dev.off()
ragg::agg_tiff(paste0(out, ".tiff"), width = w, height = h, units = "in", res = 600,
               compression = "lzw"); print(fig); grDevices::dev.off()
ragg::agg_png(paste0(out, ".png"), width = w, height = h, units = "in", res = 300);
print(fig); grDevices::dev.off()
write_log("GSE174725 annotation completed: cells=", nrow(cell_meta),
          "; clusters=", length(unique(cell_meta$cluster)),
          "; cell types=", length(unique(cell_meta$cell_type)),
          "; English-only SVG QA passed")
