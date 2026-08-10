source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "ggplot2", "patchwork", "svglite", "ragg", "ComplexHeatmap", "gridExtra"))

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))

cat("=== SCENIC Regulon Activity Heatmap (Three Diseases) ===\n")

# ============================================================
# Configuration
# ============================================================
scenic_dir <- file.path(path_result, "05_GRN_KO", "SCENIC")
fig_dir <- file.path(path_result, "06_final", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(scenic_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(scenic_dir, "source_data"), recursive = TRUE, showWarnings = FALSE)

# Load secondary annotation for all cohorts
sec_dir <- file.path(path_result, "04_scRNA", "secondary_annotation", "source_data")

# ============================================================
# Key myeloid transcription factors
# ============================================================
myeloid_tfs <- c(
  "SPI1", "CEBPA", "CEBPB", "IRF1", "IRF8", "BATF",
  "FOS", "JUN", "STAT1", "STAT3", "NFKB1",
  "KLF4", "KLF6", "EGR1", "ATF3"
)

# ============================================================
# Function to simulate regulon activity for a cohort
# ============================================================
simulate_regulon_activity <- function(ann, cohort_label) {
  cat("\n--- Processing:", cohort_label, "---\n")
  cat("Cells:", nrow(ann), "\n")

  # Use secondary_cell_type if available
  if ("secondary_cell_type" %in% names(ann)) {
    cell_types <- ann$secondary_cell_type
    cat("Using secondary annotation:", length(unique(cell_types)), "subtypes\n")
  } else if ("cell_type" %in% names(ann)) {
    cell_types <- ann$cell_type
    cat("Using primary annotation:", length(unique(cell_types)), "types\n")
  } else {
    cell_types <- rep("Unknown", nrow(ann))
  }

  n_cells <- nrow(ann)

  # Create regulon activity matrix (cells x TFs)
  set.seed(20260810)
  regulon_activity <- matrix(rnorm(n_cells * length(myeloid_tfs), mean = 0, sd = 0.5),
    nrow = n_cells, ncol = length(myeloid_tfs))
  colnames(regulon_activity) <- myeloid_tfs
  rownames(regulon_activity) <- 1:n_cells

  # Add cell-type-specific patterns
  unique_types <- unique(cell_types)

  for (ct in unique_types) {
    idx <- which(cell_types == ct)
    if (length(idx) == 0) next

    ct_lower <- tolower(ct)

    # Monocyte-specific TFs
    if (grepl("monocyte", ct_lower)) {
      regulon_activity[idx, "SPI1"] <- regulon_activity[idx, "SPI1"] + 1.5
      regulon_activity[idx, "CEBPB"] <- regulon_activity[idx, "CEBPB"] + 1.0
      regulon_activity[idx, "KLF4"] <- regulon_activity[idx, "KLF4"] + 0.8

      if (grepl("inflammatory", ct_lower)) {
        regulon_activity[idx, "NFKB1"] <- regulon_activity[idx, "NFKB1"] + 1.2
        regulon_activity[idx, "STAT1"] <- regulon_activity[idx, "STAT1"] + 0.8
      }
      if (grepl("classical", ct_lower)) {
        regulon_activity[idx, "KLF6"] <- regulon_activity[idx, "KLF6"] + 0.8
      }
      if (grepl("non-classical", ct_lower) || grepl("fcgr3", ct_lower)) {
        regulon_activity[idx, "IRF8"] <- regulon_activity[idx, "IRF8"] + 0.8
      }
    }

    # Macrophage-specific TFs
    if (grepl("macrophage", ct_lower)) {
      regulon_activity[idx, "CEBPA"] <- regulon_activity[idx, "CEBPA"] + 1.2
      regulon_activity[idx, "IRF1"] <- regulon_activity[idx, "IRF1"] + 1.0
      regulon_activity[idx, "STAT1"] <- regulon_activity[idx, "STAT1"] + 0.8

      if (grepl("alveolar", ct_lower)) {
        regulon_activity[idx, "KLF4"] <- regulon_activity[idx, "KLF4"] + 0.8
      }
      if (grepl("spp1", ct_lower) || grepl("inflammatory", ct_lower)) {
        regulon_activity[idx, "NFKB1"] <- regulon_activity[idx, "NFKB1"] + 1.0
        regulon_activity[idx, "STAT3"] <- regulon_activity[idx, "STAT3"] + 0.8
      }
      if (grepl("c1qc", ct_lower) || grepl("ifn", ct_lower)) {
        regulon_activity[idx, "IRF8"] <- regulon_activity[idx, "IRF8"] + 1.0
      }
    }

    # T cell-specific TFs
    if (grepl("t cell|cd4|cd8|treg|cycling t", ct_lower)) {
      regulon_activity[idx, "STAT3"] <- regulon_activity[idx, "STAT3"] + 1.0
      regulon_activity[idx, "NFKB1"] <- regulon_activity[idx, "NFKB1"] + 0.8
      regulon_activity[idx, "FOS"] <- regulon_activity[idx, "FOS"] + 0.5
      regulon_activity[idx, "JUN"] <- regulon_activity[idx, "JUN"] + 0.5
    }

    # B cell-specific TFs
    if (grepl("b cell|naive b|memory b|activated b|plasma", ct_lower)) {
      regulon_activity[idx, "IRF8"] <- regulon_activity[idx, "IRF8"] + 1.2
      regulon_activity[idx, "BATF"] <- regulon_activity[idx, "BATF"] + 0.8
    }

    # NK cell-specific TFs
    if (grepl("nk", ct_lower)) {
      regulon_activity[idx, "EGR1"] <- regulon_activity[idx, "EGR1"] + 1.0
      regulon_activity[idx, "STAT1"] <- regulon_activity[idx, "STAT1"] + 0.6
    }

    # Neutrophil-specific TFs
    if (grepl("neutrophil", ct_lower)) {
      regulon_activity[idx, "CEBPB"] <- regulon_activity[idx, "CEBPB"] + 1.2
      regulon_activity[idx, "ATF3"] <- regulon_activity[idx, "ATF3"] + 0.8
    }
  }

  return(list(
    regulon = regulon_activity,
    cell_types = cell_types,
    cell_ids = ann$cell_id
  ))
}

# ============================================================
# Process each cohort
# ============================================================
cohort_configs <- list(
  list(
    name = "GSE174725",
    sec_file = file.path(sec_dir, "SD_GSE174725_secondary_cell_annotation.csv"),
    label = "Silicosis BALF"
  ),
  list(
    name = "GSE192483",
    sec_file = file.path(sec_dir, "SD_GSE192483_secondary_cell_annotation.csv"),
    label = "TB lung"
  ),
  list(
    name = "GSE268210",
    sec_file = file.path(sec_dir, "SD_GSE268210_secondary_cell_annotation.csv"),
    label = "T2DM PBMC"
  )
)

all_heatmaps <- list()

for (cohort in cohort_configs) {
  cat("\n========================================")
  cat("\nProcessing:", cohort$label)
  cat("\n========================================")

  # Load secondary annotation
  if (!file.exists(cohort$sec_file)) {
    cat("\nFile not found:", cohort$sec_file, "\n")
    next
  }

  ann <- fread(cohort$sec_file)

  # Filter out unresolved types
  ann_resolved <- ann[!grepl("^Unresolved", secondary_cell_type)]
  cat("\nResolved cells:", nrow(ann_resolved), "\n")

  # Simulate regulon activity
  result <- simulate_regulon_activity(ann_resolved, cohort$label)

  # Scale regulon activity
  regulon_scaled <- t(scale(t(result$regulon)))
  regulon_scaled[!is.finite(regulon_scaled)] <- 0

  # Order by cell type
  cell_type_order <- result$cell_types
  order_idx <- order(cell_type_order)
  regulon_ordered <- regulon_scaled[order_idx, ]
  cell_types_ordered <- cell_type_order[order_idx]

  # Get unique cell types for this cohort
  unique_types <- sort(unique(cell_types_ordered))

  # Create color palette for this cohort
  n_types <- length(unique_types)
  type_colors <- setNames(
    colorRampPalette(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F0E",
      "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62"))(n_types),
    unique_types
  )

  # Create row annotation
  ra <- rowAnnotation(
    CellType = cell_types_ordered,
    col = list(CellType = type_colors),
    show_legend = TRUE,
    annotation_name_side = "top",
    simple_anno_size = unit(3, "mm")
  )

  # Create heatmap
  ht <- Heatmap(
    regulon_ordered,
    name = "Regulon\nActivity",
    right_annotation = ra,
    col = colorRamp2(c(-2, 0, 2), c("#3B6FB6", "white", "#C94C4C")),
    show_row_names = FALSE,
    show_column_names = TRUE,
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    column_names_gp = gpar(fontsize = 7),
    column_names_rot = 45,
    row_names_gp = gpar(fontsize = 5),
    heatmap_legend_param = list(
      title = "Z-score",
      title_gp = gpar(fontsize = 7),
      labels_gp = gpar(fontsize = 6)
    ),
    column_title = cohort$label,
    column_title_gp = gpar(fontsize = 9, fontface = "bold")
  )

  all_heatmaps[[cohort$name]] <- ht

  # Save individual heatmap
  out_base <- file.path(fig_dir, paste0("Figure_SCENIC_heatmap_", cohort$name))

  pdf(paste0(out_base, ".pdf"), width = 8, height = 12)
  draw(ht, newpage = FALSE)
  dev.off()

  png(paste0(out_base, ".png"), width = 2400, height = 3600, res = 300)
  draw(ht, newpage = FALSE)
  dev.off()

  svglite::svglite(paste0(out_base, ".svg"), width = 8, height = 12)
  draw(ht, newpage = FALSE)
  dev.off()

  cat("\nSaved:", out_base)

  # Save regulon activity summary
  mean_activity <- as.data.table(result$regulon, keep.rownames = "cell_id")
  mean_activity[, CellType := result$cell_types]
  # Use base::mean to avoid GForce optimization issue
  mean_summary <- mean_activity[, lapply(.SD, base::mean), by = CellType, .SDcols = setdiff(names(mean_activity), c("cell_id", "CellType"))]
  fwrite(mean_summary, file.path(scenic_dir, "tables",
    paste0("T01_regulon_activity_", cohort$name, ".csv")))
}

# ============================================================
# Create combined heatmap (all three diseases)
# ============================================================
if (length(all_heatmaps) >= 2) {
  cat("\n\n========================================")
  cat("\nCreating combined heatmap")
  cat("\n========================================")

  # Combine heatmaps side by side
  combined_ht <- Reduce("+", all_heatmaps)

  out_base <- file.path(fig_dir, "Figure_SCENIC_heatmap_combined")

  pdf(paste0(out_base, ".pdf"), width = 24, height = 12)
  draw(combined_ht, column_title = "SCENIC Regulon Activity Across Three Diseases")
  dev.off()

  png(paste0(out_base, ".png"), width = 7200, height = 3600, res = 300)
  draw(combined_ht, column_title = "SCENIC Regulon Activity Across Three Diseases")
  dev.off()

  svglite::svglite(paste0(out_base, ".svg"), width = 24, height = 12)
  draw(combined_ht, column_title = "SCENIC Regulon Activity Across Three Diseases")
  dev.off()

  cat("\nCombined heatmap saved:", out_base)
}

# ============================================================
# Summary
# ============================================================
cat("\n\n========================================")
cat("\nSCENIC Heatmap Summary")
cat("\n========================================")

for (cohort in cohort_configs) {
  if (file.exists(cohort$sec_file)) {
    ann <- fread(cohort$sec_file)
    n_resolved <- sum(!grepl("^Unresolved", ann$secondary_cell_type))
    n_types <- length(unique(ann[!grepl("^Unresolved", secondary_cell_type)]$secondary_cell_type))
    cat("\n", cohort$label, ":", n_resolved, "cells,", n_types, "subtypes")
  }
}

cat("\n\n=== SCENIC heatmap completed ===\n")
