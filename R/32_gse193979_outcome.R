source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "DESeq2", "edgeR", "GEOquery",
                  "ggplot2", "patchwork", "svglite", "ragg"))
suppressPackageStartupMessages(library(GEOquery))

acc <- "GSE193979"

# ============================================================
# 1. Load count data
# ============================================================
count_file <- file.path(path_data, "02_GEO_bulk", acc, "raw", paste0(acc, "_TANDEM_longitudinal_paper2_rawdata.txt.gz"))
counts_raw <- data.table::fread(count_file)
cat("Loaded counts:", nrow(counts_raw), "genes x", ncol(counts_raw) - 1, "samples\n")

counts <- as.matrix(counts_raw[, -1])
rownames(counts) <- counts_raw$V1

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

cat("\nDisease categories:\n")
print(table(meta$disease_category, useNA = "ifany"))
cat("\nTimepoints:\n")
print(table(meta$timepoint, useNA = "ifany"))

# Match by order
n_counts <- ncol(counts)
cat("\nMatching", n_counts, "count columns to first", n_counts, "metadata rows\n")
meta_matched <- meta[seq_len(n_counts)]
meta_matched$rseq_id <- colnames(counts)

# ============================================================
# 3. TB-DM vs TB-only at Diagnosis
# ============================================================
cat("\n=== TB-DM vs TB-only at Diagnosis ===\n")
meta_base <- meta_matched[disease_category %in% c("TB-DM", "TB-only") & timepoint == "Diagnosis"]
cat("Baseline samples:", nrow(meta_base), "\n")

if (nrow(meta_base) >= 10) {
  counts_base <- counts[, meta_base$rseq_id, drop = FALSE]
  dge <- edgeR::DGEList(counts = counts_base)
  keep <- edgeR::filterByExpr(dge, group = meta_base$disease_category, min.count = 10, min.total.count = 50)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  cat("Genes after filtering:", nrow(dge), "\n")

  group_factor <- factor(meta_base$disease_category, levels = c("TB-only", "TB-DM"))
  site_factor <- factor(meta_base$site)
  coldata <- data.frame(group = group_factor, site = site_factor, row.names = colnames(dge))

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

  cat("DEG (FDR<0.05):", sum(res_dt$fdr < 0.05, na.rm = TRUE), "\n")

  result_dir <- file.path(path_result, "02_bulk", acc)
  dir.create(file.path(result_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
  safe_write_csv(res_dt, file.path(result_dir, "tables", paste0("T01_", acc, "_TB_DM_vs_TB_diagnosis.csv")))
}

# ============================================================
# 4. Summary
# ============================================================
result_dir <- file.path(path_result, "02_bulk", acc)
summary_dt <- data.table::data.table(
  metric = c("Total samples (matched)", "TB-DM Diagnosis", "TB-only Diagnosis",
             "Unique patients"),
  value = c(nrow(meta_matched),
            sum(meta_matched$disease_category == "TB-DM" & meta_matched$timepoint == "Diagnosis", na.rm = TRUE),
            sum(meta_matched$disease_category == "TB-only" & meta_matched$timepoint == "Diagnosis", na.rm = TRUE),
            length(unique(meta_matched$patient_id)))
)
safe_write_csv(summary_dt, file.path(result_dir, "tables", paste0("T00_", acc, "_summary.csv")))

write_log(acc, " outcome analysis completed: matched=", nrow(meta_matched),
          "; unique_patients=", length(unique(meta_matched$patient_id)))
