source(file.path("R", "00_config.R"), encoding = "UTF-8")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(hdf5r))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))

cat("=== GSE326212 TB BAL Analysis (Simplified) ===\n\n")

# ============================================================
# 1. Load and process samples individually
# ============================================================
cat("1. Loading samples...\n")

h5_dir <- file.path(path_data, "03_scRNA", "GSE326212", "raw", "selected_h5")
h5_files <- list.files(h5_dir, pattern = "\\.h5$", full.names = TRUE)

# Load metadata
metadata <- fread(file.path(path_result, "04_scRNA", "GSE326212", "tables", "T05_GSE326212_primary_sample_selection.csv"))

# Process first 5 samples
samples_to_load <- head(h5_files, 5)
obj_list <- list()

for (h5 in samples_to_load) {
  sample_name <- gsub("GSM\\d+_(.+)\\.h5", "\\1", basename(h5))
  cat("  Loading", sample_name, "...")

  tryCatch({
    h5file <- H5File$new(h5, mode = "r")
    counts <- h5file[["matrix/data"]]$read()
    barcodes <- h5file[["matrix/barcodes"]]$read()
    features <- h5file[["matrix/features/name"]]$read()
    indices <- h5file[["matrix/indices"]]$read()
    indptr <- h5file[["matrix/indptr"]]$read()
    shape <- h5file[["matrix/shape"]]$read()
    h5file$close_all()

    sparse_mat <- Matrix::sparseMatrix(
      i = indices + 1,
      j = rep(1:length(indptr[-1]), diff(indptr)),
      x = counts,
      dims = shape
    )
    colnames(sparse_mat) <- barcodes
    rownames(sparse_mat) <- features

    obj <- CreateSeuratObject(counts = sparse_mat, project = sample_name, min.cells = 3, min.features = 200)
    obj$sample_id <- sample_name
    obj_list[[sample_name]] <- obj
    cat(" ", ncol(obj), "cells\n")
  }, error = function(e) {
    cat(" ERROR:", e$message, "\n")
  })
}

cat("\n  Loaded", length(obj_list), "samples\n")

# ============================================================
# 2. Process each sample individually
# ============================================================
cat("\n2. Processing each sample...\n")

processed_list <- list()
for (name in names(obj_list)) {
  cat("  Processing", name, "...")
  obj <- obj_list[[name]]
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 10, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:10, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:10, verbose = FALSE)
  processed_list[[name]] <- obj
  cat(" done\n")
}

# ============================================================
# 3. Cell type annotation
# ============================================================
cat("\n3. Cell type annotation...\n")

# Merge all objects
obj <- merge(processed_list[[1]], y = processed_list[-1], add.cell.ids = names(processed_list))
cat("  Total cells:", ncol(obj), "\n")

# Add clinical metadata
obj$group <- NA_character_
for (i in 1:nrow(metadata)) {
  sample <- metadata$sample_id[i]
  idx <- which(obj$sample_id == sample)
  if (length(idx) > 0) {
    obj$group[idx] <- metadata$analysis_group[i]
  }
}

# Simple cell type annotation based on marker expression
cat("  Annotating cell types...\n")

# Function to get mean expression of a gene across cells
get_mean_expr <- function(obj, gene_pattern) {
  all_genes <- rownames(obj)
  idx <- grep(gene_pattern, all_genes, ignore.case = TRUE)
  if (length(idx) > 0) {
    # Use LayerData for Seurat v5
    expr_mat <- LayerData(obj, layer = "data")
    return(colMeans(expr_mat[idx, , drop = FALSE]))
  }
  return(rep(0, ncol(obj)))
}

# Calculate scores for each cell type
obj$monocyte_score <- get_mean_expr(obj, "^CD14[-_]") + get_mean_expr(obj, "^LYZ[-_]") + get_mean_expr(obj, "^S100A8[-_]")
obj$macrophage_score <- get_mean_expr(obj, "^CD68[-_]") + get_mean_expr(obj, "^MARCO[-_]")
obj$neutrophil_score <- get_mean_expr(obj, "^FCGR3B[-_]") + get_mean_expr(obj, "^CSF3R[-_]")
obj$epithelial_score <- get_mean_expr(obj, "^EPCAM[-_]") + get_mean_expr(obj, "^KRT19[-_]")

# Assign cell types based on highest score
obj$cell_type <- "Other"
obj$cell_type[obj$monocyte_score > obj$macrophage_score & obj$monocyte_score > obj$neutrophil_score & obj$monocyte_score > obj$epithelial_score] <- "Monocytes"
obj$cell_type[obj$macrophage_score > obj$monocyte_score & obj$macrophage_score > obj$neutrophil_score & obj$macrophage_score > obj$epithelial_score] <- "Macrophages"
obj$cell_type[obj$neutrophil_score > obj$monocyte_score & obj$neutrophil_score > obj$macrophage_score & obj$neutrophil_score > obj$epithelial_score] <- "Neutrophils"
obj$cell_type[obj$epithelial_score > obj$monocyte_score & obj$epithelial_score > obj$macrophage_score & obj$epithelial_score > obj$neutrophil_score] <- "Epithelial"

cat("  Cell type distribution:\n")
print(table(obj$cell_type))

# ============================================================
# 4. Run UMAP on merged object
# ============================================================
cat("\n4. Running UMAP...\n")

# Join layers for Seurat v5
obj <- JoinLayers(obj)

obj <- RunPCA(obj, npcs = 20, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:15, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:15, verbose = FALSE)

# ============================================================
# 5. UMAP visualization
# ============================================================
cat("\n5. Creating UMAP...\n")

fig_dir <- file.path(path_result, "04_scRNA", "GSE326212", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

p_umap_ct <- DimPlot(obj, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE) +
  ggtitle("Cell Type") + theme(plot.title = element_text(size = 10, face = "bold"))

p_umap_group <- DimPlot(obj, reduction = "umap", group.by = "group", label = TRUE) +
  ggtitle("Clinical Group") + theme(plot.title = element_text(size = 10, face = "bold"))

fig_umap <- p_umap_ct | p_umap_group

out_base <- file.path(fig_dir, "F19_GSE326212_UMAP")
w <- 183/25.4; h <- 90/25.4

svglite::svglite(paste0(out_base, ".svg"), width = w, height = h)
print(fig_umap)
dev.off()

grDevices::cairo_pdf(paste0(out_base, ".pdf"), width = w, height = h, family = "Arial")
print(fig_umap)
dev.off()

ragg::agg_tiff(paste0(out_base, ".tiff"), width = w, height = h, units = "in", res = 600, compression = "lzw")
print(fig_umap)
dev.off()

ragg::agg_png(paste0(out_base, ".png"), width = w, height = h, units = "in", res = 300)
print(fig_umap)
dev.off()

cat("  UMAP saved\n")

# ============================================================
# 5. Save results
# ============================================================
cat("\n5. Saving results...\n")

# Save Seurat object
model_dir <- file.path(path_result, "04_scRNA", "GSE326212", "models")
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(obj, file.path(model_dir, "GSE326212_seurat.rds"))

# Save cell type summary
ct_summary <- data.table(
  cell_type = names(table(obj$cell_type)),
  count = as.numeric(table(obj$cell_type)),
  percent = round(100 * as.numeric(table(obj$cell_type)) / ncol(obj), 1)
)
fwrite(ct_summary, file.path(path_result, "06_final", "tables", "scRNA_GSE326212_cell_type_summary.csv"), bom = TRUE)

cat("  Saved\n")

write_log("GSE326212 analysis completed: ", ncol(obj), " cells, ", length(unique(obj$cell_type)), " cell types")
