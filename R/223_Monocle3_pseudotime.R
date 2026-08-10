source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "ggplot2", "patchwork", "svglite", "ragg", "monocle3", "Seurat"))

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(monocle3))
suppressPackageStartupMessages(library(Seurat))

cat("=== Monocle3 Pseudotime/Trajectory Analysis ===\n")

# ============================================================
# Configuration
# ============================================================
monocle_dir <- file.path(path_result, "05_GRN_KO", "Monocle3")
dir.create(file.path(monocle_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(monocle_dir, "source_data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(monocle_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

fig_dir <- file.path(path_result, "06_final", "figures")

# ============================================================
# Helper function: Create Monocle3 object from SCE
# ============================================================
create_monocle_from_sce <- function(sce_path, annotation_path, cohort_name) {
  cat("\n--- Processing", cohort_name, "---\n")

  # Load SCE object
  if (!file.exists(sce_path)) {
    cat("SCE not found, skipping:", sce_path, "\n")
    return(NULL)
  }

  sce <- readRDS(sce_path)
  cat("SCE loaded:", nrow(sce), "genes x", ncol(sce), "cells\n")

  # Load annotation
  ann <- fread(annotation_path)
  cat("Annotation loaded:", nrow(ann), "cells\n")

  # Match cells
  cell_ids <- colnames(sce)
  ann_matched <- ann[cell_id %in% cell_ids]

  if (nrow(ann_matched) < 100) {
    cat("Too few matched cells, skipping\n")
    return(NULL)
  }

  # Get counts
  counts <- assay(sce, "counts")
  counts <- counts[, ann_matched$cell_id]

  # Create cell_data_set (Monocle3)
  gene_metadata <- data.frame(
    gene_short_name = rownames(counts),
    row.names = rownames(counts)
  )

  cell_metadata <- ann_matched
  rownames(cell_metadata) <- cell_metadata$cell_id

  cds <- new_cell_data_set(
    expression_data = counts,
    cell_metadata = cell_metadata,
    gene_metadata = gene_metadata
  )

  # Preprocessing
  cds <- preprocess_cds(cds, num_dim = 50)
  cds <- reduce_dimension(cds)

  # Cluster
  cds <- cluster_cells(cds, resolution = 1e-3)

  # Learn trajectory
  cds <- learn_graph(cds, use_partition = TRUE)

  return(cds)
}

# ============================================================
# Process each cohort
# ============================================================
cohorts <- list(
  list(
    name = "GSE174725",
    sce = file.path(path_data, "03_scRNA", "GSE174725", "processed", "GSE174725_sce_qc.rds"),
    annotation = file.path(path_result, "04_scRNA", "GSE174725", "source_data",
      "SD03_GSE174725_cell_annotation_UMAP.csv"),
    label = "Silicosis BALF",
    root_cells = NULL  # Will be determined automatically
  ),
  list(
    name = "GSE192483",
    sce = file.path(path_data, "03_scRNA", "GSE192483", "processed", "GSE192483_sce_qc.rds"),
    annotation = file.path(path_result, "04_scRNA", "GSE192483", "source_data",
      "SD03_GSE192483_cell_annotation_UMAP.csv"),
    label = "TB lung",
    root_cells = NULL
  )
)

cds_list <- list()
for (cohort in cohorts) {
  cds <- create_monocle_from_sce(cohort$sce, cohort$annotation, cohort$name)
  if (!is.null(cds)) {
    cds_list[[cohort$name]] <- cds
  }
}

# ============================================================
# Order cells in pseudotime
# ============================================================
cat("\n--- Ordering cells in pseudotime ---\n")

for (name in names(cds_list)) {
  cds <- cds_list[[name]]

  # Order cells
  cds <- order_cells(cds)

  # Extract pseudotime
  pseudotime_df <- data.frame(
    cell_id = colnames(cds),
    pseudotime = pseudotime(cds),
    cluster = clusters(cds),
    cell_type = cds@colData$cell_type,
    group = cds@colData$group
  )

  # Save pseudotime
  fwrite(as.data.table(pseudotime_df),
    file.path(monocle_dir, "source_data", paste0(name, "_pseudotime.csv")))

  cat(name, ": pseudotime range =", range(pseudotime(cds), na.rm = TRUE), "\n")
}

# ============================================================
# Visualization 1: Trajectory plot
# ============================================================
cat("\n--- Creating trajectory plots ---\n")

for (name in names(cds_list)) {
  cds <- cds_list[[name]]

  # Plot trajectory colored by cell type
  p_trajectory <- plot_cells(cds,
    color_cells_by = "cell_type",
    label_groups_by_cluster = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE,
    graph_label_size = 3,
    cell_size = 0.5
  ) +
    labs(title = paste(name, "- Cell type trajectory"),
      subtitle = "Monocle3 learned trajectory") +
    theme(legend.position = "right")

  # Save
  out_base <- file.path(fig_dir, paste0("Figure_Monocle3_trajectory_", name))

  ggsave(paste0(out_base, ".pdf"), p_trajectory, width = 10, height = 8)
  ggsave(paste0(out_base, ".png"), p_trajectory, width = 10, height = 8, dpi = 300)

  # Plot trajectory colored by pseudotime
  p_pseudotime <- plot_cells(cds,
    color_cells_by = "pseudotime",
    label_groups_by_cluster = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE,
    graph_label_size = 3,
    cell_size = 0.5
  ) +
    labs(title = paste(name, "- Pseudotime"),
      subtitle = "Monocle3 pseudotime ordering") +
    scale_color_viridis_c(name = "Pseudotime")

  ggsave(paste0(out_base, "_pseudotime.pdf"), p_pseudotime, width = 10, height = 8)
  ggsave(paste0(out_base, "_pseudotime.png"), p_pseudotime, width = 10, height = 8, dpi = 300)

  cat(name, ": trajectory plots saved\n")
}

# ============================================================
# Visualization 2: Gene expression along pseudotime
# ============================================================
cat("\n--- Creating gene expression along pseudotime ---\n")

# Define marker genes for myeloid differentiation
myeloid_markers <- c(
  "S100A8", "S100A9",  # Classical monocyte
  "FCN1", "VCAN",       # Classical monocyte
  "IL1B", "CCL3",       # Inflammatory monocyte
  "FCGR3A", "MS4A7",    # Non-classical monocyte
  "FABP4", "MARCO",     # Alveolar macrophage
  "SPP1", "TREM2",      # SPP1 macrophage
  "CD68", "CD163"       # General macrophage
)

for (name in names(cds_list)) {
  cds <- cds_list[[name]]

  # Check which markers are available
  available_markers <- intersect(myeloid_markers, rownames(cds))

  if (length(available_markers) >= 5) {
    # Plot gene expression along pseudotime
    p_genes <- plot_genes_in_pseudotime(
      cds[available_markers[1:min(8, length(available_markers))]],
      color_cells_by = "cell_type",
      min_expr = 0.5
    ) +
      labs(title = paste(name, "- Marker expression along pseudotime"))

    out_base <- file.path(fig_dir, paste0("Figure_Monocle3_genes_pseudotime_", name))
    ggsave(paste0(out_base, ".pdf"), p_genes, width = 12, height = 8)
    ggsave(paste0(out_base, ".png"), p_genes, width = 12, height = 8, dpi = 300)

    cat(name, ": gene expression plot saved\n")
  }
}

# ============================================================
# Visualization 3: Differential expression along pseudotime
# ============================================================
cat("\n--- Differential expression along pseudotime ---\n")

for (name in names(cds_list)) {
  cds <- cds_list[[name]]

  # Fit differential expression models
  cds_fit <- fit_models(
    cds,
    model_formula_str = "~ pseudotime",
    cores = 1
  )

  # Extract model summaries
  model_summary <- data.table(
    gene = rownames(cds_fit@rowData),
    intercept = cds_fit@rowData$model_summary[[1]]$coefficient[1],
    slope = cds_fit@rowData$model_summary[[1]]$coefficient[2],
    p_value = cds_fit@rowData$model_summary[[1]]$p_value[2]
  )

  # Adjust p-values
  model_summary[, fdr := p.adjust(p_value, method = "BH")]
  model_summary <- model_summary[order(fdr)]

  # Save
  fwrite(model_summary, file.path(monocle_dir, "tables",
    paste0(name, "_pseudotime_DE.csv")))

  cat(name, ": ", sum(model_summary$fdr < 0.05), " genes associated with pseudotime\n", sep = "")

  # Plot top genes
  top_genes <- model_summary[fdr < 0.05][1:min(10, .N), gene]
  if (length(top_genes) >= 3) {
    p_top <- plot_genes_in_pseudotime(
      cds[top_genes],
      color_cells_by = "cell_type",
      min_expr = 0.5
    ) +
      labs(title = paste(name, "- Top pseudotime-associated genes"))

    out_base <- file.path(fig_dir, paste0("Figure_Monocle3_top_genes_", name))
    ggsave(paste0(out_base, ".pdf"), p_top, width = 12, height = 8)
    ggsave(paste0(out_base, ".png"), p_top, width = 12, height = 8, dpi = 300)
  }
}

# ============================================================
# Summary statistics
# ============================================================
cat("\n--- Monocle3 analysis summary ---\n")

for (name in names(cds_list)) {
  cds <- cds_list[[name]]
  cat("\n", name, ":\n")
  cat("  Cells:", ncol(cds), "\n")
  cat("  Genes:", nrow(cds), "\n")
  cat("  Clusters:", length(unique(clusters(cds))), "\n")
  cat("  Pseudotime range:", range(pseudotime(cds), na.rm = TRUE), "\n")
}

cat("\n=== Monocle3 analysis completed ===\n")
write_log("Monocle3 pseudotime analysis completed")
