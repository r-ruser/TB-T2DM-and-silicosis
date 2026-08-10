source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "ggplot2", "patchwork", "svglite", "ragg", "Seurat", "CellChat"))

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(CellChat))

cat("=== CellChat Cell-Cell Communication Analysis ===\n")

# ============================================================
# Configuration
# ============================================================
cellchat_dir <- file.path(path_result, "05_GRN_KO", "CellChat")
dir.create(file.path(cellchat_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cellchat_dir, "source_data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cellchat_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

fig_dir <- file.path(path_result, "06_final", "figures")

# ============================================================
# Helper function: Create CellChat object from SCE
# ============================================================
create_cellchat_from_sce <- function(sce_path, annotation_path, cohort_name) {
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

  # Create Seurat object
  seu <- CreateSeuratObject(counts = counts, meta.data = ann_matched)
  seu <- NormalizeData(seu)

  # Create CellChat object
  cellchat <- createCellChat(object = seu, group.by = "cell_type")

  # Set database
  CellChatDB <- CellChatDB.human
  cellchat@DB <- CellChatDB

  # Preprocessing
  tryCatch({
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat, do.fast = TRUE)
    cellchat <- identifyOverExpressedInteractions(cellchat)
  }, error = function(e) {
    cat("  Warning in preprocessing:", e$message, "\n")
    cat("  Attempting without do.fast...\n")
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
    cellchat <- identifyOverExpressedInteractions(cellchat)
  })

  # Run CellChat
  tryCatch({
    cellchat <- computeCommunProb(cellchat, type = "triMean")
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
  }, error = function(e) {
    cat("  Error in CellChat pipeline:", e$message, "\n")
    return(NULL)
  })

  # Save
  tryCatch({
    saveRDS(cellchat, file.path(cellchat_dir, "source_data",
      paste0(cohort_name, "_cellchat.rds")))
  }, error = function(e) {
    cat("  Warning: Could not save CellChat object:", e$message, "\n")
  })

  return(cellchat)
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
    label = "Silicosis BALF"
  ),
  list(
    name = "GSE192483",
    sce = file.path(path_data, "03_scRNA", "GSE192483", "processed", "GSE192483_sce_qc.rds"),
    annotation = file.path(path_result, "04_scRNA", "GSE192483", "source_data",
      "SD03_GSE192483_cell_annotation_UMAP.csv"),
    label = "TB lung"
  )
)

cellchat_list <- list()
for (cohort in cohorts) {
  cc <- create_cellchat_from_sce(cohort$sce, cohort$annotation, cohort$name)
  if (!is.null(cc)) {
    cellchat_list[[cohort$name]] <- cc
  }
}

# ============================================================
# Merge CellChat objects for comparison
# ============================================================
if (length(cellchat_list) >= 2) {
  cat("\n--- Merging CellChat objects ---\n")

  # Merge
  cellchat_merged <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))

  # Save merged object
  saveRDS(cellchat_merged, file.path(cellchat_dir, "source_data", "merged_cellchat.rds"))

  # ============================================================
  # Visualization 1: Interaction count and weight
  # ============================================================
  cat("Creating interaction count heatmap...\n")

  pdf(file.path(fig_dir, "Figure_CellChat_interaction_count.pdf"), width = 10, height = 8)
  print(compareInteractions(cellchat_merged, show.legend = TRUE, group = c(1, 2),
    measure = "count"))
  dev.off()

  png(file.path(fig_dir, "Figure_CellChat_interaction_count.png"), width = 3000, height = 2400, res = 300)
  print(compareInteractions(cellchat_merged, show.legend = TRUE, group = c(1, 2),
    measure = "count"))
  dev.off()

  # ============================================================
  # Visualization 2: Differential communication
  # ============================================================
  cat("Creating differential communication heatmap...\n")

  pdf(file.path(fig_dir, "Figure_CellChat_differential_heatmap.pdf"), width = 10, height = 8)
  print(netVisual_heatmap(cellchat_merged, comparison = c(1, 2), color.use = NULL,
    title.name = "Differential communication"))
  dev.off()

  png(file.path(fig_dir, "Figure_CellChat_differential_heatmap.png"), width = 3000, height = 2400, res = 300)
  print(netVisual_heatmap(cellchat_merged, comparison = c(1, 2), color.use = NULL,
    title.name = "Differential communication"))
  dev.off()

  # ============================================================
  # Visualization 3: Bubble plot
  # ============================================================
  cat("Creating bubble plot...\n")

  pdf(file.path(fig_dir, "Figure_CellChat_bubble_plot.pdf"), width = 12, height = 10)
  print(netVisual_bubble(cellchat_merged, comparison = c(1, 2), angle.x = 45,
    remove.isolate = TRUE))
  dev.off()

  png(file.path(fig_dir, "Figure_CellChat_bubble_plot.png"), width = 3600, height = 3000, res = 300)
  print(netVisual_bubble(cellchat_merged, comparison = c(1, 2), angle.x = 45,
    remove.isolate = TRUE))
  dev.off()
}

# ============================================================
# Per-cohort network visualization
# ============================================================
for (name in names(cellchat_list)) {
  cat("\n--- Network visualization for", name, "---\n")

  cc <- cellchat_list[[name]]

  # Circle plot
  pdf(file.path(fig_dir, paste0("Figure_CellChat_network_", name, ".pdf")), width = 10, height = 10)
  netVisual_circle(cc@net$count, vertex.weight = table(cc@idents),
    weight.scale = TRUE, label.edge = FALSE, title.name = paste(name, "- Interaction count"))
  dev.off()

  png(file.path(fig_dir, paste0("Figure_CellChat_network_", name, ".png")), width = 3000, height = 3000, res = 300)
  netVisual_circle(cc@net$count, vertex.weight = table(cc@idents),
    weight.scale = TRUE, label.edge = FALSE, title.name = paste(name, "- Interaction count"))
  dev.off()
}

# ============================================================
# Summary statistics
# ============================================================
cat("\n--- CellChat summary ---\n")

for (name in names(cellchat_list)) {
  cc <- cellchat_list[[name]]
  cat("\n", name, ":\n")
  cat("  Cell types:", paste(unique(cc@idents), collapse = ", "), "\n")
  cat("  Interactions:", sum(cc@net$count), "\n")
  cat("  Pathways:", length(unique(cc@netP$pathways)), "\n")
}

cat("\n=== CellChat analysis completed ===\n")
write_log("CellChat analysis completed")
