source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "ggplot2", "patchwork", "svglite", "ragg"))

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))

cat("=== Circular UMAP with Secondary Annotation ===\n")

# ============================================================
# Configuration
# ============================================================
circlize_dir <- file.path(path_result, "06_final", "figures")
dir.create(circlize_dir, recursive = TRUE, showWarnings = FALSE)

# Nature-figure theme
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

# Secondary cell type color palette (myeloid-focused)
secondary_colors <- c(
  # Monocyte subtypes
  "Classical monocyte" = "#E41A1C",
  "Inflammatory monocyte" = "#FF7F00",
  "Non-classical monocyte" = "#377EB8",
  "IFN-stimulated monocyte" = "#4DAF4A",
  # Macrophage subtypes
  "Alveolar macrophage" = "#984EA3",
  "C1QC macrophage" = "#A65628",
  "SPP1 macrophage" = "#F781BF",
  "Inflammatory macrophage" = "#66C2A5",
  "IFN-stimulated macrophage" = "#FC8D62",
  # T cell subtypes
  "CD4 naive" = "#8DA0CB",
  "CD4 memory" = "#E78AC3",
  "CD8 naive" = "#D53E4F",
  "CD8 cytotoxic" = "#333333",
  "CD8 effector-memory" = "#1F78B4",
  "Treg" = "#B2DF8A",
  "Cycling T" = "#FB9A99",
  # B cell subtypes
  "Naive B" = "#CCEBC5",
  "Memory B" = "#BC80BD",
  "Activated B" = "#FED9A6",
  "Plasma cell" = "#B3B3B3",
  # NK subtypes
  "NK FCGR3A cytotoxic" = "#E5C494",
  "NK XCL1 chemokine" = "#8DD3C7",
  # Other
  "Mast" = "#FFFFB3",
  "Resting mast" = "#BEBADA",
  "Activated mast" = "#FB8072",
  "Neutrophil" = "#80B1D3",
  "Activated neutrophil" = "#FDB462",
  "Mature neutrophil" = "#B3DE69",
  "Epithelial" = "#FCCDE5",
  "AT1" = "#D9D9D9",
  "AT2" = "#BC80BD",
  "Club" = "#CCEBC5",
  "Ciliated" = "#FED9A6",
  "Basal" = "#FFED6F",
  "Platelet" = "#E41A1C"
)

# ============================================================
# Load secondary annotation data
# ============================================================
cat("\n--- Loading secondary annotation data ---\n")

sec_dir <- file.path(path_result, "04_scRNA", "secondary_annotation", "source_data")

# Load each cohort
cohorts <- list(
  list(
    name = "Silicosis_BALF",
    file = file.path(sec_dir, "SD_GSE174725_secondary_cell_annotation.csv"),
    umap_file = NULL,  # UMAP coords in same file
    label = "Silicosis BALF"
  ),
  list(
    name = "TB_lung",
    file = file.path(sec_dir, "SD_GSE192483_secondary_cell_annotation.csv"),
    umap_file = NULL,  # UMAP coords in same file
    label = "TB lung"
  ),
  list(
    name = "T2DM_PBMC",
    file = file.path(sec_dir, "SD_GSE268210_secondary_cell_annotation.csv"),
    umap_file = file.path(path_result, "04_scRNA", "GSE268210", "source_data",
      "SD_GSE268210_UMAP_coordinates.csv"),
    label = "T2DM PBMC"
  )
)

# ============================================================
# Helper function: Create circular UMAP with secondary annotation
# ============================================================
create_circular_umap_secondary <- function(ann_file, umap_file, cohort_name, label) {
  cat("\n--- Creating circular UMAP for", label, "---\n")

  if (!file.exists(ann_file)) {
    cat("File not found:", ann_file, "\n")
    return(NULL)
  }

  ann <- fread(ann_file)
  cat("Cells:", nrow(ann), "\n")

  # Load UMAP coordinates if not in main file
  if (!"UMAP1" %in% names(ann) && !is.null(umap_file) && file.exists(umap_file)) {
    cat("Loading UMAP coordinates from separate file\n")
    umap_coords <- fread(umap_file)
    # Match by cell_id
    idx <- match(ann$cell_id, umap_coords$cell_id)
    ann[, UMAP1 := umap_coords$UMAP1[idx]]
    ann[, UMAP2 := umap_coords$UMAP2[idx]]
    # Remove cells without UMAP coords
    ann <- ann[!is.na(UMAP1) & !is.na(UMAP2)]
    cat("Cells with UMAP:", nrow(ann), "\n")
  }

  # Use secondary_cell_type
  ann[, cell_type := secondary_cell_type]

  # Filter out "Unresolved" types for cleaner visualization
  ann_resolved <- ann[!grepl("^Unresolved", secondary_cell_type)]
  cat("Resolved cells:", nrow(ann_resolved), "\n")

  # Scale coordinates to [-1, 1]
  ann_resolved[, UMAP1_scaled := (UMAP1 - mean(UMAP1, na.rm = TRUE)) / sd(UMAP1, na.rm = TRUE)]
  ann_resolved[, UMAP2_scaled := (UMAP2 - mean(UMAP2, na.rm = TRUE)) / sd(UMAP2, na.rm = TRUE)]

  # Clip to [-2, 2]
  ann_resolved[, UMAP1_scaled := pmax(pmin(UMAP1_scaled, 2), -2)]
  ann_resolved[, UMAP2_scaled := pmax(pmin(UMAP2_scaled, 2), -2)]

  # Scale to [0, 1] for polar coordinates
  ann_resolved[, UMAP1_polar := (UMAP1_scaled + 2) / 4]
  ann_resolved[, UMAP2_polar := (UMAP2_scaled + 2) / 4]

  # Get available colors
  available_types <- unique(ann_resolved$cell_type)
  colors_to_use <- secondary_colors[names(secondary_colors) %in% available_types]

  # Add any missing colors
  missing_types <- setdiff(available_types, names(colors_to_use))
  if (length(missing_types) > 0) {
    extra_colors <- colorRampPalette(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"))(length(missing_types))
    names(extra_colors) <- missing_types
    colors_to_use <- c(colors_to_use, extra_colors)
  }

  # Create circular plot
  p <- ggplot(ann_resolved, aes(x = UMAP1_polar, y = UMAP2_polar, color = cell_type)) +
    geom_point(size = 0.2, alpha = 0.6) +
    coord_polar() +
    scale_color_manual(values = colors_to_use, name = "Cell Type") +
    labs(title = paste(label, "- Secondary Annotation"),
      subtitle = paste("Polar coordinate (n =", nrow(ann_resolved), "resolved cells)")) +
    theme_nature() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      legend.position = "right",
      legend.key.size = unit(2, "mm")
    )

  # Subsample for saving if too large
  if (nrow(ann_resolved) > 50000) {
    set.seed(20260810)
    save_idx <- sample(nrow(ann_resolved), 50000)
    p_save <- ggplot(ann_resolved[save_idx], aes(x = UMAP1_polar, y = UMAP2_polar, color = cell_type)) +
      geom_point(size = 0.15, alpha = 0.5) +
      coord_polar() +
      scale_color_manual(values = colors_to_use, name = "Cell Type") +
      labs(title = paste(label, "- Secondary Annotation"),
        subtitle = paste("Polar coordinate (n =", 50000, "sampled cells)")) +
      theme_nature() +
      theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "right",
        legend.key.size = unit(2, "mm")
      )
  } else {
    p_save <- p
  }

  # Save
  out_base <- file.path(circlize_dir, paste0("Figure_circlize_UMAP_secondary_", cohort_name))

  ggsave(paste0(out_base, ".pdf"), p_save, width = 12, height = 10)
  ggsave(paste0(out_base, ".png"), p_save, width = 12, height = 10, dpi = 300)

  cat("Saved:", out_base, "\n")
  return(p_save)
}

# ============================================================
# Create circular UMAPs for each cohort
# ============================================================
plots <- list()
for (cohort in cohorts) {
  p <- create_circular_umap_secondary(cohort$file, cohort$umap_file, cohort$name, cohort$label)
  if (!is.null(p)) {
    plots[[cohort$name]] <- p
  }
}

# ============================================================
# Create combined plot
# ============================================================
if (length(plots) >= 2) {
  cat("\n--- Creating combined circular UMAP ---\n")

  combined <- wrap_plots(plots, ncol = 2) +
    plot_annotation(
      title = "Secondary Cell Type Annotation - Circular UMAP",
      subtitle = "Three-disease single-cell atlas with myeloid subtypes",
      theme = theme(
        plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8)
      )
    )

  out_file <- file.path(circlize_dir, "Figure_circlize_UMAP_secondary_combined")
  ggsave(paste0(out_file, ".pdf"), combined, width = 20, height = 15)
  ggsave(paste0(out_file, ".png"), combined, width = 20, height = 15, dpi = 300)

  cat("Combined saved:", out_file, "\n")
}

# ============================================================
# Summary
# ============================================================
cat("\n--- Summary ---\n")
for (cohort in cohorts) {
  if (file.exists(cohort$file)) {
    ann <- fread(cohort$file)
    n_resolved <- sum(!grepl("^Unresolved", ann$secondary_cell_type))
    n_types <- length(unique(ann[!grepl("^Unresolved", secondary_cell_type)]$secondary_cell_type))
    cat(cohort$label, ":", n_resolved, "resolved cells,", n_types, "subtypes\n")
  }
}

cat("\n=== Circular UMAP with secondary annotation completed ===\n")
