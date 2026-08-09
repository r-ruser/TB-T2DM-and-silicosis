source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table"))

cat("=== Final Integration ===\n\n")

# ============================================================
# 1. Collect all validation results
# ============================================================
results_list <- list()

# GSE114192 (discovery)
f <- file.path(path_result, "02_bulk", "GSE114192", "tables", "T02_GSE114192_DESeq2_results.csv")
if (file.exists(f)) {
  dt <- data.table::fread(f)
  results_list[["GSE114192"]] <- data.table::data.table(
    dataset = "GSE114192",
    role = "Discovery",
    contrast = "TB-DM vs TB-only",
    n_total = NA_integer_,
    n_tested = nrow(dt),
    n_deg_fdr005 = sum(dt$fdr < 0.05, na.rm = TRUE),
    n_deg_strong = sum(dt$fdr < 0.05 & abs(dt$log2FC) > 0.5, na.rm = TRUE)
  )
}

# GSE181143 (validation)
f <- file.path(path_result, "02_bulk", "GSE181143", "tables", "T00_GSE181143_summary.csv")
if (file.exists(f)) {
  dt <- data.table::fread(f)
  n_deg <- dt[metric == "DEG (FDR<0.05)", value]
  results_list[["GSE181143"]] <- data.table::data.table(
    dataset = "GSE181143",
    role = "External validation",
    contrast = "TB-DM vs TB-only",
    n_total = as.integer(dt[metric == "Total samples", value]),
    n_tested = as.integer(dt[metric == "Genes tested", value]),
    n_deg_fdr005 = as.integer(n_deg),
    n_deg_strong = as.integer(dt[metric == "Strong DEG (FDR<0.05 & |logFC|>0.5)", value])
  )
}

# GSE193978 (longitudinal)
f <- file.path(path_result, "02_bulk", "GSE193978", "tables", "T00_GSE193978_summary.csv")
if (file.exists(f)) {
  dt <- data.table::fread(f)
  results_list[["GSE193978"]] <- data.table::data.table(
    dataset = "GSE193978",
    role = "Longitudinal validation",
    contrast = "TB-DM vs TB-only (Diagnosis)",
    n_total = as.integer(dt[metric == "Total samples (matched)", value]),
    n_tested = NA_integer_,
    n_deg_fdr005 = 0L,
    n_deg_strong = 0L
  )
}

# GSE193979 (outcome)
f <- file.path(path_result, "02_bulk", "GSE193979", "tables", "T00_GSE193979_summary.csv")
if (file.exists(f)) {
  dt <- data.table::fread(f)
  results_list[["GSE193979"]] <- data.table::data.table(
    dataset = "GSE193979",
    role = "Outcome validation",
    contrast = "TB-DM vs TB-only (Diagnosis)",
    n_total = as.integer(dt[metric == "Total samples (matched)", value]),
    n_tested = NA_integer_,
    n_deg_fdr005 = 2L,
    n_deg_strong = NA_integer_
  )
}

# GSE249102 (glycaemic gradient)
f <- file.path(path_result, "02_bulk", "GSE249102", "tables", "T00_GSE249102_summary.csv")
if (file.exists(f)) {
  dt <- data.table::fread(f)
  results_list[["GSE249102"]] <- data.table::data.table(
    dataset = "GSE249102",
    role = "Glycaemic gradient",
    contrast = "TB-DM vs TB",
    n_total = as.integer(dt[metric == "Total samples", value]),
    n_tested = NA_integer_,
    n_deg_fdr005 = as.integer(dt[metric == "TB-DM vs TB DEG (FDR<0.05)", value]),
    n_deg_strong = NA_integer_
  )
}

# GSE283452 (macrophage Mtb)
f <- file.path(path_result, "02_bulk", "GSE283452", "tables", "T01_GSE283452_HAM_T2DM_x_Mtb.csv")
if (file.exists(f)) {
  dt <- data.table::fread(f)
  results_list[["GSE283452_HAM"]] <- data.table::data.table(
    dataset = "GSE283452",
    role = "Mechanistic bridge",
    contrast = "T2DM x Mtb interaction (HAM)",
    n_total = 66L,
    n_tested = nrow(dt),
    n_deg_fdr005 = sum(dt$fdr < 0.05, na.rm = TRUE),
    n_deg_strong = sum(dt$fdr < 0.05 & abs(dt$log2FoldChange) > 0.5, na.rm = TRUE)
  )
}

# GSE165489 (silicosis discovery)
f <- file.path(path_result, "02_bulk", "GSE165489", "tables", "T01_GSE165489_DESeq2_results.csv")
if (file.exists(f)) {
  dt <- data.table::fread(f)
  results_list[["GSE165489"]] <- data.table::data.table(
    dataset = "GSE165489",
    role = "Silicosis discovery",
    contrast = "Silicosis vs exposed-no-silicosis",
    n_total = NA_integer_,
    n_tested = nrow(dt),
    n_deg_fdr005 = sum(dt$fdr < 0.05, na.rm = TRUE),
    n_deg_strong = sum(dt$fdr < 0.05 & abs(dt$log2FC) > 0.5, na.rm = TRUE)
  )
}

# Combine all results
if (length(results_list) > 0) {
  final_table <- data.table::rbindlist(results_list, fill = TRUE)

  result_dir <- file.path(path_result, "06_final", "tables")
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(final_table, file.path(result_dir, "T01_all_dataset_summary.csv"), bom = TRUE)

  cat("Final integration table:\n")
  print(final_table)
}

# ============================================================
# 2. Cross-disease convergence summary
# ============================================================
cat("\n=== Cross-disease convergence summary ===\n")

# Check RRHO results
rrho_file <- file.path(path_result, "03_cross_disease", "GSE114192_GSE165489", "tables", "T01_RRHO_results.csv")
if (file.exists(rrho_file)) {
  rrho <- data.table::fread(rrho_file)
  cat("RRHO overlap regions:", nrow(rrho), "\n")
}

# Check GSEA concordance
gsea_file <- file.path(path_result, "03_cross_disease", "tables", "T02_GSEA_concordance.csv")
if (file.exists(gsea_file)) {
  gsea <- data.table::fread(gsea_file)
  cat("GSEA concordant pathways:", sum(gsea$concordant, na.rm = TRUE), "/", nrow(gsea), "\n")
}

# Check WGCNA preservation
wgcna_file <- file.path(path_result, "03_cross_disease", "WGCNA", "tables", "T01_WGCNA_module_preservation.csv")
if (file.exists(wgcna_file)) {
  wgcna <- data.table::fread(wgcna_file)
  cat("WGCNA preserved modules:", sum(wgcna$zsummary > 2, na.rm = TRUE), "/", nrow(wgcna), "\n")
}

# ============================================================
# 3. Update analysis decisions
# ============================================================
cat("\n=== Updating analysis decisions ===\n")

decisions_file <- file.path(path_result, "00_audit", "analysis_decisions.md")
decisions <- readLines(decisions_file)

# Add new decisions
new_decisions <- c(
  "",
  "## 2026-08-08 — validation completion",
  "",
  "1. GSE181143 (TB-DM validation, n=60 baseline): 1 DEG (FDR<0.05). All baseline samples from India.",
  "2. GSE193978 (TANDEM longitudinal, n=267 matched): 0 DEGs at Diagnosis, 1 at Month_2, 7 at Month_6.",
  "3. GSE193979 (TANDEM outcome, n=148 matched): 2 DEGs at Diagnosis.",
  "4. GSE249102 (glycaemic gradient, n=20): 0 DEGs. Small sample size, underpowered.",
  "5. GSE283452 (macrophage Mtb challenge, HAM n=66): 906 T2DM x Mtb interaction DEGs.",
  "6. MDM (n=23) had only T2DM samples; interaction model not identifiable.",
  "7. TANDEM datasets (GSE193978/193979) RSEQ-to-GSM mapping inferred by order.",
  "8. SCENIC on GSE174725 myeloid cells: GRNBoost2 run pending.",
  "9. All validation datasets show limited TB-DM DEGs, consistent with small effect sizes",
  "   and cross-population heterogeneity. Module-level and pathway-level evidence",
  "   is prioritized over single-gene DEG counts."
)

writeLines(c(decisions, new_decisions), decisions_file)
cat("Analysis decisions updated\n")

# ============================================================
# 4. Log completion
# ============================================================
write_log("Final integration completed: ",
          length(results_list), " datasets summarized; ",
          "analysis decisions updated")

cat("\n=== All analyses completed ===\n")
