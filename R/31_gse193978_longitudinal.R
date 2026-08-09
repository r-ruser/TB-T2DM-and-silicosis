source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "DESeq2", "edgeR", "GEOquery",
                  "ggplot2", "patchwork", "svglite", "ragg"))
suppressPackageStartupMessages(library(GEOquery))

acc <- "GSE193978"

# ============================================================
# 1. Load count data
# ============================================================
count_file <- file.path(path_data, "02_GEO_bulk", acc, "raw", paste0(acc, "_TANDEM_longitudinal_paper1_rawdata.txt.gz"))
counts_raw <- data.table::fread(count_file)
cat("Loaded counts:", nrow(counts_raw), "genes x", ncol(counts_raw) - 1, "samples\n")

# First column is gene ID (V1)
counts <- as.matrix(counts_raw[, -1])
rownames(counts) <- counts_raw$V1
colnames(counts) <- paste0("RSEQ", sprintf("%03d", seq_len(ncol(counts))))

# ============================================================
# 2. Load metadata
# ============================================================
geo_file <- file.path(path_data, "02_GEO_bulk", acc, "metadata", paste0(acc, "_series_matrix.txt.gz"))
geo <- getGEO(filename = geo_file, GSEMatrix = TRUE, getGPL = FALSE)
if (is.list(geo)) geo <- geo[[1]]
meta <- data.table::as.data.table(Biobase::pData(geo))

# Standardize labels
meta[, disease_category := `disease_category:ch1`]
meta[, timepoint := `timepoint:ch1`]
meta[, patient_id := `patient id:ch1`]
meta[, site := `field_site:ch1`]

# Normalize timepoint labels
meta[, timepoint_norm := NA_character_]
meta[grepl("Diagnosis", timepoint, ignore.case = TRUE), timepoint_norm := "Diagnosis"]
meta[grepl("Week.?2", timepoint, ignore.case = TRUE), timepoint_norm := "Week_2"]
meta[grepl("Month.?2$", timepoint, ignore.case = TRUE), timepoint_norm := "Month_2"]
meta[grepl("Month.?6", timepoint, ignore.case = TRUE), timepoint_norm := "Month_6"]

cat("\nDisease categories:\n")
print(table(meta$disease_category, useNA = "ifany"))
cat("\nTimepoints:\n")
print(table(meta$timepoint_norm, useNA = "ifany"))

# Match by order: first 267 metadata rows = 267 RSEQ columns
n_counts <- ncol(counts)
cat("\nMatching", n_counts, "count columns to first", n_counts, "metadata rows\n")
meta_matched <- meta[seq_len(n_counts)]
meta_matched$rseq_id <- colnames(counts)

# Verify some matches
cat("First 5 metadata titles:", paste(head(meta_matched$title, 5), collapse = "\n  "), "\n")

# ============================================================
# 3. PRIMARY ANALYSIS: TB-DM vs TB-only at Diagnosis
# ============================================================
cat("\n=== PRIMARY: TB-DM vs TB-only at Diagnosis ===\n")
meta_base <- meta_matched[disease_category %in% c("TB-DM", "TB-only") & timepoint_norm == "Diagnosis"]
cat("Baseline samples:", nrow(meta_base), "\n")
print(table(meta_base$disease_category, meta_base$site))

if (nrow(meta_base) >= 10) {
  counts_base <- counts[, meta_base$rseq_id, drop = FALSE]

  # Low-expression filter
  dge <- edgeR::DGEList(counts = counts_base)
  keep <- edgeR::filterByExpr(dge, group = meta_base$disease_category, min.count = 10, min.total.count = 50)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  cat("Genes after filtering:", nrow(dge), "\n")

  # DESeq2
  group_factor <- factor(meta_base$disease_category, levels = c("TB-only", "TB-DM"))
  site_factor <- factor(meta_base$site)
  coldata <- data.frame(group = group_factor, site = site_factor, row.names = colnames(dge))

  # Check if site can be included
  if (length(levels(site_factor)) > 1) {
    design <- ~ site + group
  } else {
    design <- ~ group
  }

  dds <- DESeq2::DESeqDataSetFromMatrix(countData = dge$counts, colData = coldata, design = design)
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds, contrast = c("group", "TB-DM", "TB-only"), alpha = 0.05)

  res_dt <- data.table::as.data.table(as.data.frame(res), keep.rownames = "gene_symbol")
  if ("rn" %in% colnames(res_dt)) data.table::setnames(res_dt, "rn", "gene_symbol")
  res_dt[, fdr := stats::p.adjust(pvalue, method = "BH")]
  res_dt[, is_deg := fdr < 0.05]

  cat("\nDEG (FDR<0.05):", sum(res_dt$is_deg, na.rm = TRUE), "\n")

  # Save results
  result_dir <- file.path(path_result, "02_bulk", acc)
  dir.create(file.path(result_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
  safe_write_csv(res_dt, file.path(result_dir, "tables", paste0("T01_", acc, "_TB_DM_vs_TB_diagnosis.csv")))
}

# ============================================================
# 4. Longitudinal: TB-DM vs TB-only across timepoints
# ============================================================
cat("\n=== Longitudinal: TB-DM vs TB-only across timepoints ===\n")
meta_long <- meta_matched[disease_category %in% c("TB-DM", "TB-only")]
cat("Longitudinal samples:", nrow(meta_long), "\n")

for (tp in c("Diagnosis", "Week_2", "Month_2", "Month_6")) {
  meta_tp <- meta_long[timepoint_norm == tp]
  if (nrow(meta_tp) < 10) {
    cat(tp, ": too few samples (", nrow(meta_tp), "), skipping\n")
    next
  }

  counts_tp <- counts[, meta_tp$rseq_id, drop = FALSE]
  dge <- edgeR::DGEList(counts = counts_tp)
  keep <- edgeR::filterByExpr(dge, group = meta_tp$disease_category, min.count = 5)
  dge <- dge[keep, , keep.lib.sizes = FALSE]

  group_factor <- factor(meta_tp$disease_category, levels = c("TB-only", "TB-DM"))
  coldata <- data.frame(group = group_factor, row.names = colnames(dge))

  dds <- DESeq2::DESeqDataSetFromMatrix(countData = dge$counts, colData = coldata, design = ~ group)
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds, contrast = c("group", "TB-DM", "TB-only"), alpha = 0.05)

  n_deg <- sum(res$padj < 0.05, na.rm = TRUE)
  cat(tp, ": samples =", nrow(meta_tp), "; DEG(FDR<0.05) =", n_deg, "\n")
}

# ============================================================
# 5. Module score analysis
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
    # Calculate log-CPM for module scores
    dge_all <- edgeR::DGEList(counts = counts)
    dge_all <- edgeR::calcNormFactors(dge_all)
    logcpm <- edgeR::cpm(dge_all, log = TRUE)

    # Module scores per sample
    meta_matched$module_up <- colMeans(logcpm[shared_up_match, , drop = FALSE], na.rm = TRUE)
    meta_matched$module_down <- colMeans(logcpm[shared_down_match, , drop = FALSE], na.rm = TRUE)

    # Summary by disease category and timepoint
    cat("\nModule scores by disease category and timepoint:\n")
    for (dc in c("TB-DM", "TB-only")) {
      for (tp in c("Diagnosis", "Week_2", "Month_2", "Month_6")) {
        sub <- meta_matched[disease_category == dc & timepoint_norm == tp]
        if (nrow(sub) > 0) {
          cat(dc, tp, ": up =", round(mean(sub$module_up), 2),
              "; down =", round(mean(sub$module_down), 2), "; n =", nrow(sub), "\n")
        }
      }
    }

    safe_write_csv(meta_matched, file.path(result_dir, "source_data", paste0("SD01_", acc, "_module_scores.csv")))
  }
}

# ============================================================
# 6. Summary
# ============================================================
result_dir <- file.path(path_result, "02_bulk", acc)
summary_dt <- data.table::data.table(
  metric = c("Total samples (matched)", "TB-DM Diagnosis", "TB-only Diagnosis",
             "Unique patients", "Sites"),
  value = c(nrow(meta_matched),
            sum(meta_matched$disease_category == "TB-DM" & meta_matched$timepoint_norm == "Diagnosis", na.rm = TRUE),
            sum(meta_matched$disease_category == "TB-only" & meta_matched$timepoint_norm == "Diagnosis", na.rm = TRUE),
            length(unique(meta_matched$patient_id)),
            paste(unique(meta_matched$site), collapse = ", "))
)
safe_write_csv(summary_dt, file.path(result_dir, "tables", paste0("T00_", acc, "_summary.csv")))

write_log(acc, " longitudinal analysis completed: matched=", nrow(meta_matched),
          "; unique_patients=", length(unique(meta_matched$patient_id)))
