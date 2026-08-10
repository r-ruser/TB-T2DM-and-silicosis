source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "ggplot2", "patchwork", "svglite", "ragg"))

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))

cat("=== Circular UMAP Visualization ===\n")

# ============================================================
# Configuration
# ============================================================
circlize_dir <- file.path(path_result, "06_final", "figures")
dir.create(circlize_dir, recursive = TRUE, showWarnings = FALSE)

# Nature-figure theme (using sans-serif for Windows compatibility)
theme_nature <- function(base_size = 7) {
  theme_classic(base_size = base_size) +
    theme(
      axis.line = element_line(linewidth = 0.3, color = "black"),
      axis.ticks = element_line(linewidth = 0.3, color = "black"),
      axis.text = element_text(size = 6, color = "black"),
      axis.title = element_text(size = 7, color = "black"),
      legend.text = element_text(size = 5.5),
      legend.title = element_text(size = 6, face = "bold"),
      plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 6, color = "#666666", hjust = 0.5)
    )
}

# Color palette for cell types
cell_type_colors <- c(
  "Monocyte" = "#E41A1C",
  "Macrophage" = "#377EB8",
  "Myeloid" = "#E41A1C",
  "T" = "#4DAF4A",
  "CD4 T" = "#4DAF4A",
  "CD8 T" = "#984EA3",
  "NK" = "#A65628",
  "B" = "#FF7F0E",
  "Plasma" = "#F781BF",
  "Dendritic" = "#999999",
  "Mast" = "#66C2A5",
  "Epithelial" = "#FC8D62",
  "Neutrophil" = "#D53E4F",
  "Platelet" = "#8DA0CB",
  "Other" = "#BBBBBB"
)

# ============================================================
# Helper function: Create circular UMAP using coord_polar
# ============================================================
create_circular_umap <- function(umap_data, cohort_name, color_by = NULL) {
  cat("\n--- Creating circular UMAP for", cohort_name, "---\n")

  umap_dt <- as.data.table(umap_data)

  # Auto-detect cell type column
  if (is.null(color_by)) {
    if ("cell_type" %in% names(umap_dt)) {
      color_by <- "cell_type"
    } else if ("broad_cell_type" %in% names(umap_dt)) {
      color_by <- "broad_cell_type"
    } else {
      stop("No cell type column found in data")
    }
  }
  cat("Using color column:", color_by, "\n")

  # Scale coordinates to [-1, 1]
  umap_dt[, UMAP1_scaled := (UMAP1 - mean(UMAP1, na.rm = TRUE)) / sd(UMAP1, na.rm = TRUE)]
  umap_dt[, UMAP2_scaled := (UMAP2 - mean(UMAP2, na.rm = TRUE)) / sd(UMAP2, na.rm = TRUE)]

  # Clip to [-2, 2] to avoid extreme outliers
  umap_dt[, UMAP1_scaled := pmax(pmin(UMAP1_scaled, 2), -2)]
  umap_dt[, UMAP2_scaled := pmax(pmin(UMAP2_scaled, 2), -2)]

  # Scale to [0, 1] for polar coordinates
  umap_dt[, UMAP1_polar := (UMAP1_scaled + 2) / 4]
  umap_dt[, UMAP2_polar := (UMAP2_scaled + 2) / 4]

  # Get available colors
  available_types <- unique(umap_dt[[color_by]])
  colors_to_use <- cell_type_colors[names(cell_type_colors) %in% available_types]

  # Create circular plot
  p <- ggplot(umap_dt, aes(x = UMAP1_polar, y = UMAP2_polar, color = get(color_by))) +
    geom_point(size = 0.15, alpha = 0.5) +
    coord_polar() +
    scale_color_manual(values = colors_to_use, name = "Cell Type") +
    labs(title = paste(cohort_name, "- Circular UMAP"),
      subtitle = "Polar coordinate visualization") +
    theme_nature() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      legend.position = "right"
    )

  # Save
  out_base <- file.path(circlize_dir, paste0("Figure_circlize_UMAP_", cohort_name))

  # Subsample for saving if too large
  if (nrow(umap_dt) > 50000) {
    set.seed(20260810)
    save_idx <- sample(nrow(umap_dt), 50000)
    p_save <- ggplot(umap_dt[save_idx], aes(x = UMAP1_polar, y = UMAP2_polar, color = get(color_by))) +
      geom_point(size = 0.1, alpha = 0.4) +
      coord_polar() +
      scale_color_manual(values = colors_to_use, name = "Cell Type") +
      labs(title = paste(cohort_name, "- Circular UMAP"),
        subtitle = paste("Polar coordinate visualization (n =", 50000, "cells)")) +
      theme_nature() +
      theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "right"
      )
  } else {
    p_save <- p
  }

  ggsave(paste0(out_base, ".pdf"), p_save, width = 8, height = 8)
  ggsave(paste0(out_base, ".png"), p_save, width = 8, height = 8, dpi = 300)

  cat("Circular UMAP saved:", out_base, "\n")
  return(p)
}

# ============================================================
# Load UMAP data for each cohort
# ============================================================
cat("\n--- Loading UMAP data ---\n")

# GSE174725 Silicosis BALF
ann174 <- fread(file.path(path_result, "04_scRNA", "GSE174725", "source_data",
  "SD03_GSE174725_cell_annotation_UMAP.csv"))
cat("GSE174725:", nrow(ann174), "cells\n")

# GSE192483 TB lung
ann192 <- fread(file.path(path_result, "04_scRNA", "GSE192483", "source_data",
  "SD03_GSE192483_cell_annotation_UMAP.csv"))
cat("GSE192483:", nrow(ann192), "cells\n")

# GSE268210 T2DM PBMC
ann268 <- fread(file.path(path_result, "04_scRNA", "GSE268210", "source_data",
  "SD_GSE268210_UMAP_coordinates.csv"))
cat("GSE268210:", nrow(ann268), "cells\n")
cat("GSE268210 columns:", paste(names(ann268), collapse = ", "), "\n")

# ============================================================
# Create circular UMAPs
# ============================================================
p1 <- create_circular_umap(ann174, "Silicosis_BALF")
p2 <- create_circular_umap(ann192, "TB_lung")
p3 <- create_circular_umap(ann268, "T2DM_PBMC")

# ============================================================
# Create combined circular UMAP with facets
# ============================================================
cat("\n--- Creating combined circular UMAP ---\n")

# Prepare data for combined plot
ann174_combined <- copy(ann174)[, cohort := "Silicosis BALF"]
ann192_combined <- copy(ann192)[, cohort := "TB lung"]
ann268_combined <- copy(ann268)[, cohort := "T2DM PBMC"]

# Standardize column names
setnames(ann174_combined, "cell_type", "cell_type_main", skip_absent = TRUE)
setnames(ann192_combined, "cell_type", "cell_type_main", skip_absent = TRUE)
setnames(ann268_combined, "broad_cell_type", "cell_type_main", skip_absent = TRUE)

# Combine
combined <- rbindlist(list(
  ann174_combined[, .(UMAP1, UMAP2, cell_type = cell_type_main, cohort)],
  ann192_combined[, .(UMAP1, UMAP2, cell_type = cell_type_main, cohort)],
  ann268_combined[, .(UMAP1, UMAP2, cell_type = cell_type_main, cohort)]
), fill = TRUE)

# Scale coordinates
combined[, UMAP1_scaled := (UMAP1 - mean(UMAP1, na.rm = TRUE)) / sd(UMAP1, na.rm = TRUE)]
combined[, UMAP2_scaled := (UMAP2 - mean(UMAP2, na.rm = TRUE)) / sd(UMAP2, na.rm = TRUE)]
combined[, UMAP1_scaled := pmax(pmin(UMAP1_scaled, 2), -2)]
combined[, UMAP2_scaled := pmax(pmin(UMAP2_scaled, 2), -2)]
combined[, UMAP1_polar := (UMAP1_scaled + 2) / 4]
combined[, UMAP2_polar := (UMAP2_scaled + 2) / 4]

# Get available colors
available_types <- unique(combined$cell_type)
colors_to_use <- cell_type_colors[names(cell_type_colors) %in% available_types]

# Create combined plot
p_combined <- ggplot(combined, aes(x = UMAP1_polar, y = UMAP2_polar, color = cell_type)) +
  geom_point(size = 0.1, alpha = 0.4) +
  coord_polar() +
  facet_wrap(~ cohort, ncol = 2) +
  scale_color_manual(values = colors_to_use, name = "Cell Type") +
  labs(title = "Combined Circular UMAP",
    subtitle = "Three-disease single-cell atlas") +
  theme_nature() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    strip.text = element_text(size = 7, face = "bold"),
    legend.position = "right"
  )

# Save combined plot
out_file <- file.path(circlize_dir, "Figure_circlize_UMAP_combined")
ggsave(paste0(out_file, ".pdf"), p_combined, width = 12, height = 10)
ggsave(paste0(out_file, ".png"), p_combined, width = 12, height = 10, dpi = 300)
svglite::svglite(paste0(out_file, ".svg"), width = 12, height = 10)
print(p_combined)
dev.off()

cat("Combined circular UMAP saved\n")

# ============================================================
# Summary
# ============================================================
cat("\n--- Circular UMAP summary ---\n")
cat("Silicosis BALF:", nrow(ann174), "cells\n")
cat("TB lung:", nrow(ann192), "cells\n")
cat("T2DM PBMC:", nrow(ann268), "cells\n")
cat("Total:", nrow(combined), "cells\n")

cat("\n=== Circlize circular UMAP completed ===\n")
write_log("Circlize circular UMAP generated")
