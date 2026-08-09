source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "ggplot2", "patchwork", "svglite", "ragg"))

acc <- "GSE192483"
result_dir <- file.path(path_result, "04_scRNA", acc)
ann_file <- file.path(result_dir, "source_data", "SD03_GSE192483_cell_annotation_UMAP.csv")
score_file <- file.path(result_dir, "tables", "T04_GSE192483_cluster_marker_scores.csv")
cell_meta <- data.table::fread(ann_file)
score_dt <- data.table::fread(score_file)

review_map <- data.table::data.table(cluster = as.character(sort(unique(as.integer(cell_meta$cluster)))))
review_map[, automatic_cell_type := cell_meta$cell_type[match(cluster, cell_meta$cluster)]]
review_map[, reviewed_cell_type := automatic_cell_type]
review_map[cluster == "15", reviewed_cell_type := "Unresolved"]
review_map[cluster == "16", reviewed_cell_type := "Mixed epithelial-myeloid"]
review_map[cluster == "17", reviewed_cell_type := "Epithelial"]
review_map[, review_reason := data.table::fcase(
  cluster == "15", "Low marker-score margin and no stable lineage-specific top markers",
  cluster == "16", "Concurrent SCGB1A1/SFTPC and FABP4/MARCO expression",
  cluster == "17", "EPCAM/KRT8/MUC1/SFTA2 top-marker support",
  default = "Automatic label retained after top-marker review")]
cell_meta[, automatic_cell_type := cell_type]
cell_meta[, cell_type := review_map$reviewed_cell_type[match(as.character(cluster), review_map$cluster)]]
cell_meta[, annotation_method := "Canonical marker score plus top-marker review"]
score_dt[, reviewed_cell_type := review_map$reviewed_cell_type[match(as.character(cluster), review_map$cluster)]]
abundance <- cell_meta[, .(cells = .N), by = .(sample_id, sample_code, patient, region, cell_type)]
abundance[, percent := 100 * cells / sum(cells), by = sample_id]
safe_write_csv(cell_meta, ann_file)
safe_write_csv(score_dt, score_file)
safe_write_csv(abundance, file.path(result_dir, "tables", "T06_GSE192483_cell_type_abundance.csv"))
safe_write_csv(review_map, file.path(result_dir, "tables", "T07_GSE192483_reviewed_cluster_annotation.csv"))

cell_type_order <- c("Macrophage", "Monocyte", "Neutrophil", "Dendritic", "CD4 T", "CD8 T", "NK", "B",
  "Plasma", "Epithelial", "Endothelial", "Fibroblast", "Mast", "Platelet",
  "Mixed epithelial-myeloid", "Unresolved")
pal <- c("Macrophage" = "#3E6D8E", "Monocyte" = "#6A9FB5", "Neutrophil" = "#B84A4A",
  "Dendritic" = "#7B6FA6", "CD4 T" = "#58A27C", "CD8 T" = "#82B36A", "NK" = "#D08A3E",
  "B" = "#4B8C86", "Plasma" = "#A8648A", "Epithelial" = "#C7A35A",
  "Endothelial" = "#4D7C8A", "Fibroblast" = "#9B7E62", "Mast" = "#8B6E55", "Platelet" = "#9C755F",
  "Mixed epithelial-myeloid" = "#B39B8A", "Unresolved" = "#A8A8A8")
cell_meta[, cell_type := factor(cell_type, levels = cell_type_order)]
abundance[, cell_type := factor(cell_type, levels = cell_type_order)]
sample_levels <- c("SP019H", "SP019L", "SP020L", "SP021H", "SP021L", "SP023H", "SP023L",
  "SP024H", "SP024L", "SP025H", "SP025L")
sample_pal <- stats::setNames(c("#B84A4A", "#6A9FB5", "#9BC1D0", "#D1736E", "#477C9B", "#A94B4B",
  "#7AA6B8", "#C86460", "#557F9F", "#9E5A58", "#86B2C3"), sample_levels)
p_a <- ggplot2::ggplot(cell_meta, ggplot2::aes(UMAP1, UMAP2, colour = cell_type)) +
  ggplot2::geom_point(size = 0.09, alpha = 0.55) + ggplot2::scale_colour_manual(values = pal, drop = TRUE) +
  ggplot2::guides(colour = ggplot2::guide_legend(nrow = 3, byrow = TRUE,
    override.aes = list(size = 1.5, alpha = 1))) +
  ggplot2::labs(title = "Reviewed lung cell identities", subtitle = "Marker scores and top-marker audit",
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
    legend.text = ggplot2::element_text(size = 5.2),
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
write_log("GSE192483 reviewed annotation figure completed; English-only SVG QA passed")
