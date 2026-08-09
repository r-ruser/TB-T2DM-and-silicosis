source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "Seurat", "edgeR", "ggplot2", "patchwork", "svglite", "ragg"))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))

cat("=== GSE326212 TB BAL Complete Analysis ===\n\n")

# ============================================================
# 1. Load metadata and select samples
# ============================================================
cat("1. Loading metadata...\n")

metadata_file <- file.path(path_result, "04_scRNA", "GSE326212", "tables", "T05_GSE326212_primary_sample_selection.csv")
metadata <- fread(metadata_file)
cat("  Selected samples:", nrow(metadata), "\n")

# Get h5 files
h5_dir <- file.path(path_data, "03_scRNA", "GSE326212", "raw", "selected_h5")
h5_files <- list.files(h5_dir, pattern = "\\.h5$", full.names = TRUE)
cat("  H5 files:", length(h5_files), "\n")

# ============================================================
# 2. Load and process data
# ============================================================
cat("\n2. Loading 10x data...\n")

# Load fewer samples to fit in 16GB memory
samples_to_load <- head(h5_files, 5)
cat("  Loading", length(samples_to_load), "samples (memory constraint)...\n")

obj_list <- list()
for (h5 in samples_to_load) {
  sample_name <- gsub("GSM\\d+_(.+)\\.h5", "\\1", basename(h5))
  tryCatch({
    # Read from matrix group (standard 10x format)
    h5file <- hdf5r::H5File$new(h5, mode = "r")
    counts <- h5file[["matrix/data"]]$read()
    barcodes <- h5file[["matrix/barcodes"]]$read()
    features <- h5file[["matrix/features/name"]]$read()
    indices <- h5file[["matrix/indices"]]$read()
    indptr <- h5file[["matrix/indptr"]]$read()
    shape <- h5file[["matrix/shape"]]$read()
    h5file$close_all()

    # Create sparse matrix
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
    cat("    ", sample_name, ":", ncol(obj), "cells\n")
  }, error = function(e) {
    cat("    Error loading", sample_name, ":", e$message, "\n")
  })
}

# Merge objects later after QC
cat("\n  Loaded", length(obj_list), "samples\n")

# ============================================================
# 3. QC
# ============================================================
cat("\n3. Quality control...\n")

# Calculate QC metrics
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")

# QC plots
p_qc <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "sample_id", ncol = 3)

# Filter
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)
cat("  Cells after QC:", ncol(obj), "\n")

# ============================================================
# 4. Normalization and HVG
# ============================================================
cat("\n4. Normalization...\n")

# Merge objects first, then normalize
obj <- merge(obj_list[[1]], y = obj_list[-1], add.cell.ids = names(obj_list))

# Normalize
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj)

# ============================================================
# 5. PCA and Clustering
# ============================================================
cat("\n5. PCA and clustering...\n")

obj <- RunPCA(obj, npcs = 20)
obj <- FindNeighbors(obj, dims = 1:15)
obj <- FindClusters(obj, resolution = 0.5)
obj <- RunUMAP(obj, dims = 1:15)

cat("  Clusters:", length(unique(obj$seurat_clusters)), "\n")

# ============================================================
# 6. Cell type annotation
# ============================================================
cat("\n6. Cell type annotation...\n")

# Canonical markers
markers <- list(
  "T cells" = c("CD3D", "CD3E", "CD3G"),
  "B cells" = c("MS4A1", "CD79A"),
  "NK cells" = c("GNLY", "NKG7"),
  "Monocytes" = c("CD14", "LYZ", "S100A8", "S100A9"),
  "Macrophages" = c("CD68", "MARCO", "MSR1"),
  "Dendritic cells" = c("FCER1A", "CD1C"),
  "Neutrophils" = c("FCGR3B", "CSF3R"),
  "Mast cells" = c("TPSAB1", "TPSB2"),
  "Epithelial" = c("EPCAM", "KRT19")
)

# Score each cell type
for (ct in names(markers)) {
  genes <- intersect(markers[[ct]], rownames(obj))
  if (length(genes) > 0) {
    obj <- AddModuleScore(obj, features = list(genes), name = paste0("score_", ct))
  }
}

# Assign cell types based on highest score
score_cols <- grep("^score_", colnames(obj@meta.data), value = TRUE)
obj$cell_type <- apply(obj@meta.data[, score_cols], 1, function(x) {
  names(which.max(x))
})
obj$cell_type <- gsub("^score_", "", obj$cell_type)

cat("  Cell type distribution:\n")
print(table(obj$cell_type))

# ============================================================
# 7. Add clinical metadata
# ============================================================
cat("\n7. Adding clinical metadata...\n")

# Match sample IDs to metadata
obj$group <- NA_character_
for (i in 1:nrow(metadata)) {
  sample <- metadata$sample_id[i]
  idx <- which(obj$sample_id == sample)
  if (length(idx) > 0) {
    obj$group[idx] <- metadata$analysis_group[i]
  }
}

cat("  Group distribution:\n")
print(table(obj$group, useNA = "ifany"))

# ============================================================
# 8. UMAP visualization
# ============================================================
cat("\n8. Creating UMAP plots...\n")

fig_dir <- file.path(path_result, "04_scRNA", "GSE326212", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# UMAP by cell type
p_umap_ct <- DimPlot(obj, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE) +
  ggtitle("Cell Type") +
  theme(plot.title = element_text(size = 10, face = "bold"))

# UMAP by group
p_umap_group <- DimPlot(obj, reduction = "umap", group.by = "group", label = TRUE) +
  ggtitle("Clinical Group") +
  theme(plot.title = element_text(size = 10, face = "bold"))

# Combined UMAP
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
# 9. Pseudobulk DE analysis
# ============================================================
cat("\n9. Pseudobulk DE analysis...\n")

# Aggregate by sample and cell type
obj$sample_celltype <- paste0(obj$sample_id, "_", obj$cell_type)

# Create pseudobulk counts
cat("  Creating pseudobulk matrix...\n")
DefaultAssay(obj) <- "RNA"
counts_pb <- AggregateExpression(obj, group.by = "sample_celltype", assays = "RNA")[[1]]

# Get metadata
pb_meta <- data.table(
  sample_celltype = colnames(counts_pb),
  sample_id = gsub("_(.+)$", "", colnames(counts_pb)),
  cell_type = gsub("^(.+)_", "", colnames(counts_pb))
)

# Merge with clinical metadata
pb_meta <- merge(pb_meta, metadata[, .(sample_id, analysis_group)], by = "sample_id", all.x = TRUE)

# Filter to monocytes and macrophages
myeloid_types <- c("Monocytes", "Macrophages")
pb_myeloid <- counts_pb[, pb_meta$cell_type %in% myeloid_types]
meta_myeloid <- pb_meta[cell_type %in% myeloid_types]

cat("  Myeloid pseudobulk samples:", ncol(pb_myeloid), "\n")

# DE analysis for each cell type
for (ct in myeloid_types) {
  cat("\n  ---", ct, "---\n")

  idx <- which(meta_myeloid$cell_type == ct)
  if (length(idx) < 4) {
    cat("  Too few samples for", ct, "\n")
    next
  }

  counts_ct <- pb_myeloid[, idx]
  meta_ct <- meta_myeloid[idx]

  # Filter low expressed genes
  dge <- DGEList(counts = counts_ct)
  keep <- filterByExpr(dge, group = meta_ct$analysis_group, min.count = 10)
  dge <- dge[keep, , keep.lib.sizes = FALSE]

  # Design
  group_factor <- factor(meta_ct$analysis_group)
  if (length(levels(group_factor)) < 2) {
    cat("  Only one group for", ct, "\n")
    next
  }

  design <- model.matrix(~ group_factor)

  # DE analysis
  dds <- DESeqDataSetFromMatrix(countData = dge$counts, colData = meta_ct, design = design)
  dds <- DESeq(dds)

  # Get results for Active TB vs Controller
  if ("Active TB" %in% levels(group_factor) && "Stable controller" %in% levels(group_factor)) {
    res <- results(dds, contrast = c("group_factor", "Active TB", "Stable controller"), alpha = 0.05)

    res_dt <- data.table(as.data.frame(res), keep.rownames = "gene_symbol")
    res_dt[, fdr := padj]
    res_dt[, cell_type := ct]

    # Save
    ct_name <- gsub(" ", "_", ct)
    fwrite(res_dt, file.path(path_result, "06_final", "tables",
                             paste0("scRNA_GSE326212_", ct_name, "_Active_TB_vs_Controller.csv")), bom = TRUE)

    cat("  DEGs (FDR<0.05):", sum(res_dt$fdr < 0.05, na.rm = TRUE), "\n")
    cat("  Up:", sum(res_dt$log2FoldChange > 0 & res_dt$fdr < 0.05, na.rm = TRUE), "\n")
    cat("  Down:", sum(res_dt$log2FoldChange < 0 & res_dt$fdr < 0.05, na.rm = TRUE), "\n")
  }
}

# ============================================================
# 10. Save Seurat object
# ============================================================
cat("\n10. Saving Seurat object...\n")

model_dir <- file.path(path_result, "04_scRNA", "GSE326212", "models")
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(obj, file.path(model_dir, "GSE326212_seurat.rds"))
cat("  Saved\n")

write_log("GSE326212 complete analysis: ", ncol(obj), " cells, ", length(unique(obj$cell_type)), " cell types")
