source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "ggplot2", "patchwork", "svglite", "ragg"))

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))

cat("=== Nature-Figure Style UMAP (Fixed Font Size) ===\n")

# ============================================================
# Nature Figure Theme & Functions (from skill)
# ============================================================
theme_nature <- function(base_size = 6.5, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.5, colour = "black"),
      legend.title = element_text(size = base_size - 0.3, face = "bold"),
      legend.text = element_text(size = base_size - 0.7, colour = "black"),
      strip.text = element_text(size = base_size - 0.3, face = "bold"),
      plot.title = element_text(size = base_size + 1, face = "bold"),
      plot.subtitle = element_text(size = base_size - 0.5, colour = "#666666"),
      panel.grid = element_blank(),
      plot.margin = margin(5, 5, 5, 5)
    )
}

save_pub_r <- function(plot, filename, width_mm = 183, height_mm = 120, dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4

  svglite::svglite(paste0(filename, ".svg"), width = w, height = h)
  print(plot)
  dev.off()

  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = "Arial")
  print(plot)
  dev.off()

  ragg::agg_tiff(paste0(filename, ".tiff"), width = w, height = h, units = "in", res = dpi, compression = "lzw")
  print(plot)
  dev.off()

  ragg::agg_png(paste0(filename, ".png"), width = w, height = h, units = "in", res = 300)
  print(plot)
  dev.off()

  cat("Saved:", filename, "\n")
}

# ============================================================
# Color palettes (Nature style)
# ============================================================
myeloid_colors <- c(
  "Classical monocyte" = "#E41A1C",
  "Inflammatory monocyte" = "#FF7F00",
  "Non-classical monocyte" = "#377EB8",
  "IFN-stimulated monocyte" = "#4DAF4A",
  "Alveolar macrophage" = "#984EA3",
  "C1QC macrophage" = "#A65628",
  "SPP1 macrophage" = "#F781BF",
  "Inflammatory macrophage" = "#66C2A5",
  "IFN-stimulated macrophage" = "#FC8D62"
)

# ============================================================
# Load secondary annotation data
# ============================================================
sec_dir <- file.path(path_result, "04_scRNA", "secondary_annotation", "source_data")
fig_dir <- file.path(path_result, "06_final", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

cohorts <- list(
  list(
    name = "GSE174725",
    file = file.path(sec_dir, "SD_GSE174725_secondary_cell_annotation.csv"),
    label = "Silicosis BALF",
    short = "Silicosis"
  ),
  list(
    name = "GSE192483",
    file = file.path(sec_dir, "SD_GSE192483_secondary_cell_annotation.csv"),
    label = "TB lung",
    short = "TB"
  ),
  list(
    name = "GSE268210",
    file = file.path(sec_dir, "SD_GSE268210_secondary_cell_annotation.csv"),
    umap_file = file.path(path_result, "04_scRNA", "GSE268210", "source_data",
      "SD_GSE268210_UMAP_coordinates.csv"),
    label = "T2DM PBMC",
    short = "T2DM"
  )
)

# ============================================================
# Create UMAP for each cohort
# ============================================================
plots <- list()

for (cohort in cohorts) {
  cat("\n--- Processing:", cohort$label, "---\n")

  ann <- fread(cohort$file)

  # Load UMAP if not in main file
  if (!"UMAP1" %in% names(ann) && !is.null(cohort$umap_file) && file.exists(cohort$umap_file)) {
    umap_coords <- fread(cohort$umap_file)
    idx <- match(ann$cell_id, umap_coords$cell_id)
    ann[, UMAP1 := umap_coords$UMAP1[idx]]
    ann[, UMAP2 := umap_coords$UMAP2[idx]]
    ann <- ann[!is.na(UMAP1) & !is.na(UMAP2)]
  }

  # Filter out unresolved
  ann_resolved <- ann[!grepl("^Unresolved", secondary_cell_type)]
  cat("Resolved cells:", nrow(ann_resolved), "\n")

  # Get myeloid cells for detailed view
  myeloid_types <- names(myeloid_colors)
  ann_myeloid <- ann_resolved[secondary_cell_type %in% myeloid_types]
  cat("Myeloid cells:", nrow(ann_myeloid), "\n")

  # Create UMAP plot - all cells
  p_all <- ggplot(ann_resolved, aes(UMAP1, UMAP2, colour = secondary_cell_type)) +
    geom_point(size = 0.3, alpha = 0.6) +
    labs(title = cohort$label,
      subtitle = paste("All cell types (n =", nrow(ann_resolved), ")"),
      colour = NULL) +
    theme_nature() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_text(size = 6),
      legend.position = "right",
      legend.key.size = unit(2, "mm"),
      legend.text = element_text(size = 5)
    ) +
    guides(colour = guide_legend(override.aes = list(size = 1.5, alpha = 1), ncol = 1))

  # Create UMAP plot - myeloid only
  if (nrow(ann_myeloid) > 0) {
    p_myeloid <- ggplot(ann_myeloid, aes(UMAP1, UMAP2, colour = secondary_cell_type)) +
      geom_point(size = 0.4, alpha = 0.7) +
      scale_colour_manual(values = myeloid_colors, drop = TRUE) +
      labs(title = paste(cohort$label, "- Myeloid"),
        subtitle = paste("Myeloid subtypes (n =", nrow(ann_myeloid), ")"),
        colour = NULL) +
      theme_nature() +
      theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 6),
        legend.position = "right",
        legend.key.size = unit(2, "mm"),
        legend.text = element_text(size = 5)
      ) +
      guides(colour = guide_legend(override.aes = list(size = 1.5, alpha = 1), ncol = 1))
  } else {
    p_myeloid <- NULL
  }

  # Save individual plots
  out_base <- file.path(fig_dir, paste0("Figure_UMAP_nature_", cohort$name))
  save_pub_r(p_all, out_base, width_mm = 183, height_mm = 120)

  if (!is.null(p_myeloid)) {
    out_myeloid <- file.path(fig_dir, paste0("Figure_UMAP_myeloid_nature_", cohort$name))
    save_pub_r(p_myeloid, out_myeloid, width_mm = 183, height_mm = 120)
  }

  plots[[cohort$name]] <- list(all = p_all, myeloid = p_myeloid)
}

# ============================================================
# Create combined 3-disease UMAP
# ============================================================
cat("\n--- Creating combined 3-disease UMAP ---\n")

# Prepare data for combined plot
combined_list <- list()
for (cohort in cohorts) {
  ann <- fread(cohort$file)
  if (!"UMAP1" %in% names(ann) && !is.null(cohort$umap_file) && file.exists(cohort$umap_file)) {
    umap_coords <- fread(cohort$umap_file)
    idx <- match(ann$cell_id, umap_coords$cell_id)
    ann[, UMAP1 := umap_coords$UMAP1[idx]]
    ann[, UMAP2 := umap_coords$UMAP2[idx]]
  }
  ann <- ann[!grepl("^Unresolved", secondary_cell_type)]
  cohort_short <- cohort$short
  ann[, cohort := cohort_short]
  combined_list[[cohort$name]] <- ann
}

combined <- rbindlist(combined_list, fill = TRUE)
combined <- combined[!is.na(UMAP1) & !is.na(UMAP2)]

# Create myeloid-focused combined plot
myeloid_types <- names(myeloid_colors)
combined_myeloid <- combined[secondary_cell_type %in% myeloid_types]

p_combined <- ggplot(combined_myeloid, aes(UMAP1, UMAP2, colour = secondary_cell_type)) +
  geom_point(size = 0.2, alpha = 0.5) +
  scale_colour_manual(values = myeloid_colors, drop = TRUE) +
  facet_wrap(~ cohort, ncol = 3, scales = "free") +
  labs(title = "Myeloid Subtypes Across Three Diseases",
    subtitle = "Secondary cell type annotation",
    colour = NULL) +
  theme_nature() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 6),
    strip.text = element_text(size = 7, face = "bold"),
    legend.position = "right",
    legend.key.size = unit(2, "mm"),
    legend.text = element_text(size = 5)
  ) +
  guides(colour = guide_legend(override.aes = list(size = 1.5, alpha = 1), ncol = 1))

out_combined <- file.path(fig_dir, "Figure_UMAP_myeloid_combined_three_disease")
save_pub_r(p_combined, out_combined, width_mm = 183, height_mm = 90)

cat("\n=== Nature-figure UMAP completed ===\n")
