source(file.path("R", "00_config.R"), encoding = "UTF-8")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))

cat("=== GSE326212 Pseudobulk DE Analysis ===\n\n")

# ============================================================
# 1. Load Seurat object
# ============================================================
cat("1. Loading Seurat object...\n")

obj <- readRDS(file.path(path_result, "04_scRNA", "GSE326212", "models", "GSE326212_seurat.rds"))
cat("  Loaded:", ncol(obj), "cells\n")
cat("  Cell types:", paste(unique(obj$cell_type), collapse=", "), "\n")
cat("  Groups:", paste(unique(obj$group), collapse=", "), "\n")

# ============================================================
# 2. Create pseudobulk counts
# ============================================================
cat("\n2. Creating pseudobulk counts...\n")

# Create sample-celltype identifier
obj$sample_celltype <- paste0(obj$sample_id, "_", obj$cell_type)

# Aggregate counts
cat("  Aggregating counts by sample and cell type...\n")
DefaultAssay(obj) <- "RNA"
counts_pb <- Seurat::AggregateExpression(obj, group.by = "sample_celltype", assays = "RNA")[[1]]

cat("  Pseudobulk matrix:", nrow(counts_pb), "genes x", ncol(counts_pb), "samples\n")

# Create metadata
# Seurat replaces underscores with dashes, so format is like "T017-B-filtered-Monocytes"
# Need to parse differently
pb_meta <- data.table(
  sample_celltype = colnames(counts_pb),
  cell_type = sapply(strsplit(colnames(counts_pb), "-"), function(x) x[length(x)])
)
# Extract sample_id by removing cell_type suffix
pb_meta[, sample_id := sub(paste0("-", cell_type, "$"), "", sample_celltype)]

# Add clinical metadata from the metadata file
metadata_file <- file.path(path_result, "04_scRNA", "GSE326212", "tables", "T05_GSE326212_primary_sample_selection.csv")
metadata <- fread(metadata_file)

# Create sample info with proper matching
sample_info <- unique(data.table(
  sample_id = metadata$sample_id,
  group = metadata$analysis_group
))

# Match sample_id in pb_meta
# pb_meta has sample_ids like "T017-B-filtered" (dashes from Seurat)
# metadata has sample_ids like "T017_B" (underscores)
# Need to: 1) replace dashes with underscores, 2) remove "_filtered" suffix
pb_meta[, sample_id_clean := gsub("-", "_", sample_id)]
pb_meta[, sample_id_clean := gsub("_filtered$", "", sample_id_clean)]
pb_meta <- merge(pb_meta, sample_info, by.x = "sample_id_clean", by.y = "sample_id", all.x = TRUE)

cat("  Matched samples:", sum(!is.na(pb_meta$group)), "/", nrow(pb_meta), "\n")
cat("  Sample IDs (first 5):", head(pb_meta$sample_id_clean, 5), "\n")
cat("  Groups:", paste(unique(pb_meta$group), collapse=", "), "\n")

cat("\n  Pseudobulk sample distribution:\n")
print(table(pb_meta$cell_type, pb_meta$group))

# ============================================================
# 3. DE analysis for each cell type
# ============================================================
cat("\n3. Running DE analysis...\n")

result_dir <- file.path(path_result, "06_final", "tables")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

# Focus on key cell types
key_types <- c("Monocytes", "Macrophages", "Neutrophils")
all_results <- list()

for (ct in key_types) {
  cat("\n  --- ", ct, " ---\n")

  # Get samples for this cell type
  idx <- which(pb_meta$cell_type == ct)
  if (length(idx) < 4) {
    cat("  Too few samples for", ct, "(", length(idx), ")\n")
    next
  }

  meta_ct <- pb_meta[idx]
  counts_ct <- counts_pb[, idx, drop = FALSE]

  # Check groups
  groups <- unique(meta_ct$group)
  cat("  Groups:", paste(groups, collapse=", "), "\n")

  if (length(groups) < 2) {
    cat("  Only one group, skipping\n")
    next
  }

  # Filter low expressed genes
  dge <- DGEList(counts = counts_ct)
  keep <- filterByExpr(dge, group = meta_ct$group, min.count = 10, min.total.count = 50)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  cat("  Genes after filtering:", nrow(dge), "\n")

  # Design matrix
  group_factor <- factor(meta_ct$group)

  # Ensure we have the right contrast
  if ("Active TB" %in% levels(group_factor) && "Stable controller" %in% levels(group_factor)) {
    group_factor <- relevel(group_factor, ref = "Stable controller")
  }

  design <- model.matrix(~ group_factor)

  # DE analysis
  dds <- DESeqDataSetFromMatrix(countData = dge$counts, colData = meta_ct, design = design)
  dds <- DESeq(dds)

  # Get results for each comparison
  coefs <- resultsNames(dds)
  cat("  Coefficients:", paste(coefs, collapse=", "), "\n")

  for (coef in coefs) {
    if (coef == "Intercept") next

    res <- results(dds, name = coef, alpha = 0.05)

    res_dt <- data.table(as.data.frame(res), keep.rownames = "gene_symbol")
    res_dt[, fdr := padj]
    res_dt[, cell_type := ct]
    res_dt[, comparison := coef]

    # Save
    coef_name <- gsub("group_factor", "", coef)
    coef_name <- gsub(" ", "_", coef_name)
    filename <- paste0("scRNA_GSE326212_", gsub(" ", "_", ct), "_", coef_name, ".csv")
    fwrite(res_dt, file.path(result_dir, filename), bom = TRUE)

    n_deg <- sum(res_dt$fdr < 0.05, na.rm = TRUE)
    n_up <- sum(res_dt$log2FoldChange > 0 & res_dt$fdr < 0.05, na.rm = TRUE)
    n_down <- sum(res_dt$log2FoldChange < 0 & res_dt$fdr < 0.05, na.rm = TRUE)

    cat("  ", coef, ": DEG =", n_deg, "(up:", n_up, ", down:", n_down, ")\n")

    all_results[[paste(ct, coef)]] <- res_dt
  }
}

# ============================================================
# 4. Summary
# ============================================================
cat("\n4. Summary...\n")

summary_dt <- rbindlist(lapply(names(all_results), function(x) {
  dt <- all_results[[x]]
  parts <- strsplit(x, " ")[[1]]
  data.table(
    cell_type = parts[1],
    comparison = paste(parts[-1], collapse=" "),
    n_tested = nrow(dt),
    n_deg = sum(dt$fdr < 0.05, na.rm = TRUE),
    n_up = sum(dt$log2FoldChange > 0 & dt$fdr < 0.05, na.rm = TRUE),
    n_down = sum(dt$log2FoldChange < 0 & dt$fdr < 0.05, na.rm = TRUE)
  )
}))

fwrite(summary_dt, file.path(result_dir, "scRNA_GSE326212_pseudobulk_summary.csv"), bom = TRUE)

cat("\n  Pseudobulk DE Summary:\n")
print(summary_dt)

write_log("GSE326212 pseudobulk DE completed")
