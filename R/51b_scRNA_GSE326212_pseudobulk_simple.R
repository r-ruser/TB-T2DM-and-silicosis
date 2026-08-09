source(file.path("R", "00_config.R"), encoding = "UTF-8")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(edgeR))

cat("=== GSE326212 Pseudobulk DE (Simplified) ===\n\n")

# ============================================================
# 1. Load all h5 files and create Seurat object
# ============================================================
cat("1. Loading all h5 files...\n")

h5_dir <- file.path(path_data, "03_scRNA", "GSE326212", "raw", "selected_h5")
h5_files <- list.files(h5_dir, pattern = "\\.h5$", full.names = TRUE)
cat("  H5 files:", length(h5_files), "\n")

# Load all samples
obj_list <- list()
for (h5 in h5_files) {
  sample_name <- gsub("GSM\\d+_(.+)\\.h5", "\\1", basename(h5))
  tryCatch({
    h5file <- hdf5r::H5File$new(h5, mode = "r")
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
  }, error = function(e) {
    cat("    Error loading", sample_name, "\n")
  })
}

cat("  Loaded", length(obj_list), "samples\n")

# Merge all objects
cat("  Merging objects...\n")
obj <- merge(obj_list[[1]], y = obj_list[-1], add.cell.ids = names(obj_list))
cat("  Total cells:", ncol(obj), "\n")

# Join layers for Seurat v5
obj <- JoinLayers(obj)

# ============================================================
# 2. Load metadata
# ============================================================
cat("\n2. Loading metadata...\n")
metadata_file <- file.path(path_result, "04_scRNA", "GSE326212", "tables", "T05_GSE326212_primary_sample_selection.csv")
metadata <- fread(metadata_file)

# Create mapping from h5 sample_id to analysis_group
# h5 files have names like "T017_B_filtered", metadata has "T017_B"
sample_map <- unique(data.table(
  sample_id = metadata$sample_id,
  group = metadata$analysis_group
))

cat("  Sample groups:\n")
print(table(sample_map$group))

# ============================================================
# 3. Create pseudobulk by sample
# ============================================================
cat("\n3. Creating pseudobulk counts...\n")

# Get counts from Seurat object
counts <- GetAssayData(obj, layer = "counts")

# Create sample identifier (remove _filtered suffix)
obj$sample_clean <- gsub("_filtered$", "", obj$sample_id)

# Aggregate by sample
cat("  Aggregating by sample...\n")
samples <- unique(obj$sample_clean)
counts_list <- list()
for (s in samples) {
  idx <- which(obj$sample_clean == s)
  if (length(idx) > 0) {
    counts_list[[s]] <- rowSums(counts[, idx, drop = FALSE])
  }
}
counts_pb <- do.call(cbind, counts_list)

cat("  Pseudobulk matrix:", nrow(counts_pb), "genes x", ncol(counts_pb), "samples\n")

# Create metadata
pb_meta <- data.table(
  sample_id = colnames(counts_pb)
)
pb_meta <- merge(pb_meta, sample_map, by = "sample_id", all.x = TRUE)

cat("  Matched samples:", sum(!is.na(pb_meta$group)), "/", nrow(pb_meta), "\n")
cat("  Groups:", paste(unique(pb_meta$group), collapse=", "), "\n")

# ============================================================
# 4. DE analysis: Active TB vs Stable controller
# ============================================================
cat("\n4. Running DE analysis...\n")

# Filter to only Active TB and Stable controller
pb_meta_compare <- pb_meta[group %in% c("Active TB", "Stable controller")]
cat("  Samples for comparison:", nrow(pb_meta_compare), "\n")

if (nrow(pb_meta_compare) < 4) {
  cat("  Too few samples for comparison\n")
} else {
  counts_compare <- counts_pb[, pb_meta_compare$sample_id, drop = FALSE]

  # Filter low expressed genes
  dge <- DGEList(counts = counts_compare)
  keep <- filterByExpr(dge, group = pb_meta_compare$group, min.count = 10, min.total.count = 50)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  cat("  Genes after filtering:", nrow(dge), "\n")

  # Design matrix
  group_factor <- factor(pb_meta_compare$group, levels = c("Stable controller", "Active TB"))
  design <- model.matrix(~ group_factor)

  # DE analysis
  dds <- DESeqDataSetFromMatrix(countData = dge$counts, colData = pb_meta_compare, design = design)
  dds <- DESeq(dds)

  # Get results
  res <- results(dds, contrast = c("group", "Active TB", "Stable controller"), alpha = 0.05)

  res_dt <- data.table(as.data.frame(res), keep.rownames = "gene_symbol")
  res_dt[, fdr := padj]

  # Save results
  result_dir <- file.path(path_result, "06_final", "tables")
  fwrite(res_dt, file.path(result_dir, "scRNA_GSE326212_pseudobulk_Active_TB_vs_Controller.csv"), bom = TRUE)

  n_deg <- sum(res_dt$fdr < 0.05, na.rm = TRUE)
  n_up <- sum(res_dt$log2FoldChange > 0 & res_dt$fdr < 0.05, na.rm = TRUE)
  n_down <- sum(res_dt$log2FoldChange < 0 & res_dt$fdr < 0.05, na.rm = TRUE)

  cat("\n  Results:\n")
  cat("  Genes tested:", nrow(res_dt), "\n")
  cat("  DEG (FDR<0.05):", n_deg, "\n")
  cat("  Up:", n_up, "\n")
  cat("  Down:", n_down, "\n")

  # Top 20 genes
  cat("\n  Top 20 genes:\n")
  print(res_dt[fdr < 0.05][order(fdr)][1:min(20, n_deg), .(gene_symbol, log2FoldChange, fdr)])
}

# ============================================================
# 5. Summary
# ============================================================
cat("\n5. Summary...\n")

summary_dt <- data.table(
  metric = c("Total cells", "Samples (Active TB vs Controller)", "Genes tested", "DEG (FDR<0.05)"),
  value = c(ncol(obj), nrow(pb_meta_compare), ifelse(exists("res_dt"), nrow(res_dt), 0),
            ifelse(exists("res_dt"), sum(res_dt$fdr < 0.05, na.rm = TRUE), 0))
)

fwrite(summary_dt, file.path(result_dir, "scRNA_GSE326212_pseudobulk_summary.csv"), bom = TRUE)

cat("\n  Summary:\n")
print(summary_dt)

write_log("GSE326212 pseudobulk DE completed")
