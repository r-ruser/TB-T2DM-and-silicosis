source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "ggplot2", "patchwork", "svglite", "ragg"))

acc <- "GSE174725"
result_dir <- file.path(path_result, "04_scRNA", acc)
ann_file <- file.path(result_dir, "source_data", "SD03_GSE174725_cell_annotation_UMAP.csv")
score_file <- file.path(result_dir, "tables", "T04_GSE174725_cluster_marker_scores.csv")
if (!file.exists(ann_file) || !file.exists(score_file)) stop("Run R/41_gse174725_annotation.R first")

cell_meta <- data.table::fread(ann_file)
score_dt <- data.table::fread(score_file)

# Canonical top-marker review: cluster 13 expresses GP6, PEAR1 and CLEC1B and is platelet-like.
cell_meta[cluster == 13, cell_type := "Platelet"]
cell_meta[, annotation_method := "Canonical marker score plus top-marker review"]
score_dt[cluster == 13, `:=`(assigned_cell_type = "Platelet", assignment_margin = NA_real_)]

donor_labels <- c(
  exposure_patient_1 = "Exposure 1", exposure_patient_2 = "Exposure 2",
  silicosis_patient_1 = "Silicosis 1", silicosis_patient_2 = "Silicosis 2",
  silicosis_patient_3 = "Silicosis 3"
)
cell_meta[, donor_label := unname(donor_labels[donor])]
if (anyNA(cell_meta$donor_label)) stop("Missing donor display label")
abundance <- cell_meta[, .(cells = .N), by = .(sample_id, donor, donor_label, group, cell_type)]
abundance[, percent := 100 * cells / sum(cells), by = sample_id]

safe_write_csv(cell_meta, ann_file)
safe_write_csv(score_dt, score_file)
safe_write_csv(abundance, file.path(result_dir, "tables", "T06_GSE174725_cell_type_abundance.csv"))
safe_write_csv(data.table::data.table(cluster = sort(unique(cell_meta$cluster)),
  reviewed_cell_type = cell_meta[, unique(cell_type), by = cluster][order(cluster)]$V1,
  method = "Canonical marker score plus top-marker review"),
  file.path(result_dir, "tables", "T07_GSE174725_reviewed_cluster_annotation.csv"))

cell_type_order <- c("Macrophage", "Monocyte", "Neutrophil", "Dendritic", "CD4 T", "CD8 T",
  "NK", "B", "Plasma", "Epithelial", "Mast", "Erythroid", "Platelet")
pal <- c("Macrophage" = "#3E6D8E", "Monocyte" = "#6A9FB5", "Neutrophil" = "#B84A4A",
  "Dendritic" = "#7B6FA6", "CD4 T" = "#58A27C", "CD8 T" = "#82B36A",
  "NK" = "#D08A3E", "B" = "#4B8C86", "Plasma" = "#A8648A",
  "Epithelial" = "#C7A35A", "Mast" = "#8B6E55", "Erythroid" = "#7A7A7A",
  "Platelet" = "#9C755F")
cell_meta[, cell_type := factor(cell_type, levels = cell_type_order)]
abundance[, cell_type := factor(cell_type, levels = cell_type_order)]
donor_pal <- c("Exposure 1" = "#416A9A", "Exposure 2" = "#6A9FB5",
  "Silicosis 1" = "#B84A4A", "Silicosis 2" = "#D1736E", "Silicosis 3" = "#9E5A58")

p_a <- ggplot2::ggplot(cell_meta, ggplot2::aes(UMAP1, UMAP2, colour = cell_type)) +
  ggplot2::geom_point(size = 0.16, alpha = 0.62) +
  ggplot2::scale_colour_manual(values = pal, drop = TRUE) +
  ggplot2::guides(colour = ggplot2::guide_legend(nrow = 3, byrow = TRUE,
    override.aes = list(size = 1.5, alpha = 1))) +
  ggplot2::labs(title = "Reviewed cell identities", subtitle = "Marker scores and top-marker audit",
    x = "UMAP 1", y = "UMAP 2", colour = NULL) +
  ggplot2::theme(legend.position = "bottom", aspect.ratio = 1)
p_b <- ggplot2::ggplot(cell_meta, ggplot2::aes(UMAP1, UMAP2, colour = donor_label)) +
  ggplot2::geom_point(size = 0.16, alpha = 0.52) +
  ggplot2::scale_colour_manual(values = donor_pal) +
  ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2, byrow = TRUE,
    override.aes = list(size = 1.5, alpha = 1))) +
  ggplot2::labs(title = "Donor mixing audit", subtitle = "Five biological donors",
    x = "UMAP 1", y = "UMAP 2", colour = NULL) +
  ggplot2::theme(legend.position = "bottom", aspect.ratio = 1)
p_c <- ggplot2::ggplot(abundance, ggplot2::aes(donor_label, percent, fill = cell_type)) +
  ggplot2::geom_col(width = 0.76) +
  ggplot2::scale_fill_manual(values = pal, drop = TRUE) +
  ggplot2::scale_x_discrete(limits = names(donor_pal)) +
  ggplot2::labs(title = "Cell-type composition by donor", subtitle = "Descriptive proportions",
    x = NULL, y = "Cells (%)", fill = NULL) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 32, hjust = 1), legend.position = "none")

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
p_c <- p_c + theme_sc + ggplot2::theme(
  axis.text.x = ggplot2::element_text(angle = 32, hjust = 1), legend.position = "none")
fig <- ((p_a | p_b) / p_c) + patchwork::plot_layout(heights = c(2.15, 1)) +
  patchwork::plot_annotation(tag_levels = "a", title = "GSE174725 silicosis BALF cell atlas",
    subtitle = "Donor-resolved annotation after adaptive QC and doublet removal")
out <- file.path(result_dir, "figures", "F12_GSE174725_cell_atlas")
w <- 183 / 25.4; h <- 164 / 25.4
svglite::svglite(paste0(out, ".svg"), width = w, height = h); print(fig); grDevices::dev.off()
assert_english_only_svg(paste0(out, ".svg"))
grDevices::cairo_pdf(paste0(out, ".pdf"), width = w, height = h, family = "Arial"); print(fig); grDevices::dev.off()
ragg::agg_tiff(paste0(out, ".tiff"), width = w, height = h, units = "in", res = 600,
  compression = "lzw"); print(fig); grDevices::dev.off()
ragg::agg_png(paste0(out, ".png"), width = w, height = h, units = "in", res = 300)
print(fig); grDevices::dev.off()
write_log("GSE174725 reviewed annotation figure completed; English-only SVG QA passed")
