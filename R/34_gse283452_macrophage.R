source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "DESeq2", "edgeR", "readxl", "GEOquery",
                  "ggplot2", "patchwork", "svglite", "ragg"))
suppressPackageStartupMessages(library(GEOquery))

acc <- "GSE283452"

# ============================================================
# 1. Load count data from xlsx
# ============================================================
count_file <- file.path(path_data, "02_GEO_bulk", acc, "raw", paste0(acc, "_RawCounts.xlsx"))
counts_raw <- readxl::read_excel(count_file, sheet = 1)
cat("Loaded counts:", nrow(counts_raw), "genes x", ncol(counts_raw) - 2, "samples\n")

counts <- as.matrix(counts_raw[, -(1:2)])
rownames(counts) <- counts_raw$Name
colnames(counts) <- sub(" \\(GE\\) - Total counts", "", colnames(counts))

# ============================================================
# 2. Load metadata
# ============================================================
geo_file <- file.path(path_data, "02_GEO_bulk", acc, "metadata", paste0(acc, "_series_matrix.txt.gz"))
eset <- getGEO(filename = geo_file, GSEMatrix = TRUE, getGPL = FALSE)
if (is.list(eset)) eset <- eset[[1]]
meta <- data.table::as.data.table(Biobase::pData(eset))

# Standardize labels
meta[, celltype := `cell type:ch1`]
meta[, t2dm := NA_character_]
meta[grepl("T2D", `diabetes status:ch1`, ignore.case = TRUE), t2dm := "T2DM"]
meta[grepl("noT2D", `diabetes status:ch1`, ignore.case = TRUE), t2dm := "Control"]

meta[, mtb := NA_character_]
meta[grepl("Mtb", `treatment:ch1`, ignore.case = TRUE), mtb := "Mtb"]
meta[grepl("Uninfected", `treatment:ch1`, ignore.case = TRUE), mtb := "Media"]

meta[, timepoint := `time:ch1`]
meta[, donor := `donor number:ch1`]

# Match by order (both have 89 samples)
cat("Count cols:", ncol(counts), "; Meta rows:", nrow(meta), "\n")
if (ncol(counts) == nrow(meta)) {
  cat("Matching by order\n")
  meta$sample_id <- colnames(counts)
} else {
  stop("Sample count mismatch: ", ncol(counts), " vs ", nrow(meta))
}

cat("\nT2DM x Mtb distribution:\n")
print(table(meta$t2dm, meta$mtb, useNA = "ifany"))
cat("\nCell types:\n")
print(table(meta$celltype, useNA = "ifany"))
cat("\nTimepoints:\n")
print(table(meta$timepoint, useNA = "ifany"))

# ============================================================
# 3. DESeq2 analysis per cell type
# ============================================================
result_dir <- file.path(path_result, "02_bulk", acc)
invisible(lapply(file.path(result_dir, c("tables", "source_data", "figures")),
                 dir.create, recursive = TRUE, showWarnings = FALSE))

for (ct in unique(meta$celltype)) {
  cat("\n=== Cell type:", ct, "===\n")
  meta_ct <- meta[celltype == ct]
  cat("Samples:", nrow(meta_ct), "\n")
  print(table(meta_ct$t2dm, meta_ct$mtb))

  if (nrow(meta_ct) < 10) {
    cat("Too few samples, skipping\n")
    next
  }

  counts_ct <- counts[, meta_ct$sample_id, drop = FALSE]

  # Low-expression filter
  dge <- edgeR::DGEList(counts = counts_ct)
  keep <- edgeR::filterByExpr(dge, group = paste(meta_ct$t2dm, meta_ct$mtb), min.count = 5)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  cat("Genes after filtering:", nrow(dge), "\n")

  # Design: T2DM x Mtb interaction + timepoint
  t2dm_factor <- factor(meta_ct$t2dm, levels = c("Control", "T2DM"))
  mtb_factor <- factor(meta_ct$mtb, levels = c("Media", "Mtb"))
  time_factor <- factor(meta_ct$timepoint)
  donor_factor <- factor(meta_ct$donor)

  coldata <- data.frame(
    t2dm = t2dm_factor,
    mtb = mtb_factor,
    time = time_factor,
    donor = donor_factor,
    row.names = colnames(dge)
  )

  # Simplified design: T2DM x Mtb interaction
  design <- model.matrix(~ t2dm * mtb, data = coldata)

  ct_name <- gsub(" ", "_", ct)

  tryCatch({
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = dge$counts, colData = coldata, design = design)
    dds <- DESeq2::DESeq(dds)

    cat("Available results names:", paste(DESeq2::resultsNames(dds), collapse = ", "), "\n")

    # Extract T2DM x Mtb interaction
    int_coef <- grep("t2dmT2DM.mtbMtb|t2dmT2DM:mtbMtb", DESeq2::resultsNames(dds), value = TRUE)
    if (length(int_coef) > 0) {
      res_int <- DESeq2::results(dds, name = int_coef[1], alpha = 0.05)
      res_dt <- data.table::as.data.table(as.data.frame(res_int), keep.rownames = "gene_symbol")
      if ("rn" %in% colnames(res_dt)) data.table::setnames(res_dt, "rn", "gene_symbol")
      res_dt[, fdr := stats::p.adjust(pvalue, method = "BH")]

      safe_write_csv(res_dt, file.path(result_dir, "tables",
                      paste0("T01_", acc, "_", ct_name, "_T2DM_x_Mtb.csv")))

      cat("T2DM x Mtb interaction DEG (FDR<0.05):", sum(res_dt$fdr < 0.05, na.rm = TRUE), "\n")
    }

    # T2DM main effect
    t2dm_coef <- grep("^t2dmT2DM$", DESeq2::resultsNames(dds), value = TRUE)
    if (length(t2dm_coef) > 0) {
      res_t2dm <- DESeq2::results(dds, name = t2dm_coef[1], alpha = 0.05)
      res_t2dm_dt <- data.table::as.data.table(as.data.frame(res_t2dm), keep.rownames = "gene_symbol")
      if ("rn" %in% colnames(res_t2dm_dt)) data.table::setnames(res_t2dm_dt, "rn", "gene_symbol")
      res_t2dm_dt[, fdr := stats::p.adjust(pvalue, method = "BH")]

      safe_write_csv(res_t2dm_dt, file.path(result_dir, "tables",
                      paste0("T02_", acc, "_", ct_name, "_T2DM_main.csv")))

      cat("T2DM main effect DEG (FDR<0.05):", sum(res_t2dm_dt$fdr < 0.05, na.rm = TRUE), "\n")
    }

    # Mtb main effect
    mtb_coef <- grep("^mtbMtb$", DESeq2::resultsNames(dds), value = TRUE)
    if (length(mtb_coef) > 0) {
      res_mtb <- DESeq2::results(dds, name = mtb_coef[1], alpha = 0.05)
      res_mtb_dt <- data.table::as.data.table(as.data.frame(res_mtb), keep.rownames = "gene_symbol")
      if ("rn" %in% colnames(res_mtb_dt)) data.table::setnames(res_mtb_dt, "rn", "gene_symbol")
      res_mtb_dt[, fdr := stats::p.adjust(pvalue, method = "BH")]

      safe_write_csv(res_mtb_dt, file.path(result_dir, "tables",
                      paste0("T03_", acc, "_", ct_name, "_Mtb_main.csv")))

      cat("Mtb main effect DEG (FDR<0.05):", sum(res_mtb_dt$fdr < 0.05, na.rm = TRUE), "\n")
    }
  }, error = function(e) {
    cat("DESeq2 error:", e$message, "\n")
  })
}

# ============================================================
# 4. Module score analysis
# ============================================================
cat("\n=== Module score analysis ===\n")
discovery_file <- file.path(path_result, "02_bulk", "GSE114192", "tables", "T02_GSE114192_DESeq2_results.csv")
if (file.exists(discovery_file)) {
  discovery <- data.table::fread(discovery_file)
  shared_up <- discovery[gene_symbol != "" & fdr < 0.05 & log2FC > 0, gene_symbol]
  shared_down <- discovery[gene_symbol != "" & fdr < 0.05 & log2FC < 0, gene_symbol]

  shared_up_match <- intersect(shared_up, rownames(counts))
  shared_down_match <- intersect(shared_down, rownames(counts))

  cat("Shared up genes matched:", length(shared_up_match), "\n")
  cat("Shared down genes matched:", length(shared_down_match), "\n")

  if (length(shared_up_match) > 5) {
    # Calculate module scores
    dge_all <- edgeR::DGEList(counts = counts)
    dge_all <- edgeR::calcNormFactors(dge_all)
    logcpm <- edgeR::cpm(dge_all, log = TRUE)

    meta$module_up <- colMeans(logcpm[shared_up_match, , drop = FALSE], na.rm = TRUE)
    meta$module_down <- colMeans(logcpm[shared_down_match, , drop = FALSE], na.rm = TRUE)

    # Summary by group
    for (ct in unique(meta$celltype)) {
      meta_ct <- meta[celltype == ct]
      cat("\n", ct, "module scores:\n")
      for (t2dm_val in c("Control", "T2DM")) {
        for (mtb_val in c("Media", "Mtb")) {
          sub <- meta_ct[t2dm == t2dm_val & mtb == mtb_val]
          if (nrow(sub) > 0) {
            cat("  ", t2dm_val, mtb_val, ": up =", round(mean(sub$module_up), 2),
                "; down =", round(mean(sub$module_down), 2), "; n =", nrow(sub), "\n")
          }
        }
      }
    }

    safe_write_csv(meta, file.path(result_dir, "source_data", paste0("SD01_", acc, "_module_scores.csv")))
  }
}

write_log(acc, " macrophage Mtb challenge analysis completed")
