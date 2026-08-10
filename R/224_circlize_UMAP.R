source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "ggplot2", "patchwork", "svglite", "ragg", "circlize"))

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(circlize))

cat("=== Circlize Circular UMAP Visualization ===\n")

# ============================================================
# Configuration
# ============================================================
circlize_dir <- file.path(path_result, "06_final", "figures")
dir.create(circlize_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Helper function: Create circular UMAP
# ============================================================
create_circular_umap <- function(umap_data, cohort_name, color_by = "cell_type") {
  cat("\n--- Creating circular UMAP for", cohort_name, "---\n")

  # Get UMAP coordinates
  umap_dt <- as.data.table(umap_data)

  # Scale coordinates to [0, 1]
  umap_dt[, UMAP1_scaled := (UMAP1 - min(UMAP1)) / (max(UMAP1) - min(UMAP1))]
  umap_dt[, UMAP2_scaled := (UMAP2 - min(UMAP2)) / (max(UMAP2) - min(UMAP2))]

  # Convert to polar coordinates
  umap_dt[, angle := atan2(UMAP2_scaled - 0.5, UMAP1_scaled - 0.5)]
  umap_dt[, radius := sqrt((UMAP1_scaled - 0.5)^2 + (UMAP2_scaled - 0.5)^2)]

  # Color palette
  cell_types <- unique(umap_dt[[color_by]])
  n_types <- length(cell_types)

  # Use nature-figure color palette
  if (n_types <= 10) {
    colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F0E",
      "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62")
  } else {
    colors <- colorRampPalette(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"))(n_types)
  }
  color_map <- setNames(colors[1:n_types], cell_types)

  # Create circular plot
  out_file <- file.path(circlize_dir, paste0("Figure_circlize_UMAP_", cohort_name))

  # PDF
  pdf(paste0(out_file, ".pdf"), width = 8, height = 8)

  # Initialize circular layout
  circos.clear()
  circos.par(
    start.degree = 90,
    gap.degree = 2,
    track.margin = c(0.01, 0.01),
    cell.padding = c(0, 0, 0, 0)
  )

  # Create empty track
  circos.initialize(factors = "a", xlim = c(0, 1))

  # Add background track
  circos.track(
    factors = "a",
    ylim = c(0, 1),
    panel.fun = function(x, y) {
      # Draw grid lines
      for (r in seq(0.1, 0.9, by = 0.2)) {
        circos.lines(
          h = r * cos(seq(0, 2 * pi, length.out = 100)),
          v = r * sin(seq(0, 2 * pi, length.out = 100)),
          col = "#F0F0F0",
          lwd = 0.5
        )
      }
    },
    track.height = 0.8,
    bg.border = NA
  )

  # Plot cells as points
  for (ct in cell_types) {
    ct_data <- umap_dt[get(color_by) == ct]
    if (nrow(ct_data) == 0) next

    # Subsample if too many points
    if (nrow(ct_data) > 1000) {
      set.seed(20260810)
      ct_data <- ct_data[sample(.N, 1000)]
    }

    # Convert to Cartesian coordinates on circle
    x <- ct_data$radius * cos(ct_data$angle)
    y <- ct_data$radius * sin(ct_data$angle)

    # Add points
    circos.points(
      x = x,
      y = y,
      sector.index = "a",
      track.index = 1,
      col = alpha(color_map[ct], 0.6),
      pch = 16,
      cex = 0.3
    )
  }

  # Add legend
  legend(
    "bottomright",
    legend = names(color_map),
    col = color_map,
    pch = 16,
    cex = 0.7,
    bty = "n",
    title = "Cell Type"
  )

  # Add title
  title(main = paste(cohort_name, "- Circular UMAP"),
    sub = "Circlize circular layout")

  dev.off()

  # PNG
  png(paste0(out_file, ".png"), width = 2400, height = 2400, res = 300)

  circos.clear()
  circos.par(
    start.degree = 90,
    gap.degree = 2,
    track.margin = c(0.01, 0.01),
    cell.padding = c(0, 0, 0, 0)
  )

  circos.initialize(factors = "a", xlim = c(0, 1))

  circos.track(
    factors = "a",
    ylim = c(0, 1),
    panel.fun = function(x, y) {
      for (r in seq(0.1, 0.9, by = 0.2)) {
        circos.lines(
          h = r * cos(seq(0, 2 * pi, length.out = 100)),
          v = r * sin(seq(0, 2 * pi, length.out = 100)),
          col = "#F0F0F0",
          lwd = 0.5
        )
      }
    },
    track.height = 0.8,
    bg.border = NA
  )

  for (ct in cell_types) {
    ct_data <- umap_dt[get(color_by) == ct]
    if (nrow(ct_data) == 0) next

    if (nrow(ct_data) > 1000) {
      set.seed(20260810)
      ct_data <- ct_data[sample(.N, 1000)]
    }

    x <- ct_data$radius * cos(ct_data$angle)
    y <- ct_data$radius * sin(ct_data$angle)

    circos.points(
      x = x,
      y = y,
      sector.index = "a",
      track.index = 1,
      col = alpha(color_map[ct], 0.6),
      pch = 16,
      cex = 0.3
    )
  }

  legend(
    "bottomright",
    legend = names(color_map),
    col = color_map,
    pch = 16,
    cex = 0.7,
    bty = "n",
    title = "Cell Type"
  )

  title(main = paste(cohort_name, "- Circular UMAP"),
    sub = "Circlize circular layout")

  dev.off()

  # SVG
  svglite::svglite(paste0(out_file, ".svg"), width = 8, height = 8)

  circos.clear()
  circos.par(
    start.degree = 90,
    gap.degree = 2,
    track.margin = c(0.01, 0.01),
    cell.padding = c(0, 0, 0, 0)
  )

  circos.initialize(factors = "a", xlim = c(0, 1))

  circos.track(
    factors = "a",
    ylim = c(0, 1),
    panel.fun = function(x, y) {
      for (r in seq(0.1, 0.9, by = 0.2)) {
        circos.lines(
          h = r * cos(seq(0, 2 * pi, length.out = 100)),
          v = r * sin(seq(0, 2 * pi, length.out = 100)),
          col = "#F0F0F0",
          lwd = 0.5
        )
      }
    },
    track.height = 0.8,
    bg.border = NA
  )

  for (ct in cell_types) {
    ct_data <- umap_dt[get(color_by) == ct]
    if (nrow(ct_data) == 0) next

    if (nrow(ct_data) > 1000) {
      set.seed(20260810)
      ct_data <- ct_data[sample(.N, 1000)]
    }

    x <- ct_data$radius * cos(ct_data$angle)
    y <- ct_data$radius * sin(ct_data$angle)

    circos.points(
      x = x,
      y = y,
      sector.index = "a",
      track.index = 1,
      col = alpha(color_map[ct], 0.6),
      pch = 16,
      cex = 0.3
    )
  }

  legend(
    "bottomright",
    legend = names(color_map),
    col = color_map,
    pch = 16,
    cex = 0.7,
    bty = "n",
    title = "Cell Type"
  )

  title(main = paste(cohort_name, "- Circular UMAP"),
    sub = "Circlize circular layout")

  dev.off()

  circos.clear()

  cat("Circular UMAP saved:", out_file, "\n")
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

# ============================================================
# Create circular UMAPs
# ============================================================
create_circular_umap(ann174, "Silicosis_BALF")
create_circular_umap(ann192, "TB_lung")
create_circular_umap(ann268, "T2DM_PBMC")

# ============================================================
# Create combined circular UMAP
# ============================================================
cat("\n--- Creating combined circular UMAP ---\n")

# Combine all cohorts
ann174_combined <- copy(ann174)[, cohort := "Silicosis BALF"]
ann192_combined <- copy(ann192)[, cohort := "TB lung"]
ann268_combined <- copy(ann268)[, cohort := "T2DM PBMC"]

# Ensure column names match
setnames(ann174_combined, "cell_type", "cell_type_main", skip_absent = TRUE)
setnames(ann192_combined, "cell_type", "cell_type_main", skip_absent = TRUE)
setnames(ann268_combined, "broad_cell_type", "cell_type_main", skip_absent = TRUE)

# Combine
combined <- rbindlist(list(
  ann174_combined[, .(UMAP1, UMAP2, cell_type = cell_type_main, cohort)],
  ann192_combined[, .(UMAP1, UMAP2, cell_type = cell_type_main, cohort)],
  ann268_combined[, .(UMAP1, UMAP2, cell_type = cell_type_main, cohort)]
), fill = TRUE)

# Create facetted circular UMAP
out_file <- file.path(circlize_dir, "Figure_circlize_UMAP_combined")

# For combined plot, use ggplot2 with coord_polar for better faceting
p_combined <- ggplot(combined, aes(x = UMAP1, y = UMAP2, color = cell_type)) +
  geom_point(size = 0.2, alpha = 0.5) +
  facet_wrap(~ cohort, ncol = 2) +
  coord_polar() +
  scale_color_manual(values = c(
    "Monocyte" = "#E41A1C",
    "Macrophage" = "#377EB8",
    "T" = "#4DAF4A",
    "CD4 T" = "#4DAF4A",
    "CD8 T" = "#984EA3",
    "B" = "#FF7F0E",
    "NK" = "#A65628",
    "Other" = "#999999"
  ), na.value = "#CCCCCC") +
  labs(title = "Combined Circular UMAP",
    subtitle = "Three-disease single-cell atlas",
    color = "Cell Type") +
  theme_classic(base_size = 7) +
  theme(
    plot.title = element_text(size = 8, face = "bold"),
    strip.text = element_text(size = 7, face = "bold"),
    legend.position = "right",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )

# Save
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
