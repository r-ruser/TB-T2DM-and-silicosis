source(file.path("R", "00_config.R"), encoding = "UTF-8")
suppressPackageStartupMessages(library(data.table))

cat("=== Final Cleanup and Compilation ===\n\n")

final_tables <- file.path(path_result, "06_final", "tables")
final_figures <- file.path(path_result, "06_final", "figures")

# ============================================================
# 1. Clean up old/incorrect figures
# ============================================================
cat("1. Cleaning up old figures...\n")

# Remove old proxy-based GBD figures (keep new prevalence ones)
old_figs <- list.files(final_figures, pattern = "GBD.*proxy|proxy.*GBD", full.names = TRUE)
if (length(old_figs) > 0) {
  file.remove(old_figs)
  cat("  Removed", length(old_figs), "old proxy GBD figures\n")
}

# Remove old shared-gene based figures
old_shared <- list.files(final_figures, pattern = "shared.*gene|overlap.*gene", full.names = TRUE, ignore.case = TRUE)
if (length(old_shared) > 0) {
  file.remove(old_shared)
  cat("  Removed", length(old_shared), "old shared gene figures\n")
}

# ============================================================
# 2. Compile final tables
# ============================================================
cat("\n2. Compiling final tables...\n")

# Create comprehensive summary
summary_list <- list()

# GBD
summary_list[["GBD"]] <- data.table(
  Analysis = "GBD 2023",
  Component = "Triple-high burden",
  Result = "15 locations (75th percentile)",
  Files = "T01-T04 GBD"
)

# TB-DM Meta
summary_list[["TBDM_meta"]] <- data.table(
  Analysis = "TB-DM meta-analysis",
  Component = "Cross-cohort validation",
  Result = "380 FDR<0.05 genes",
  Files = "T12 meta-analysis"
)

# Silicosis
summary_list[["Silicosis"]] <- data.table(
  Analysis = "Silicosis discovery",
  Component = "GSE165489",
  Result = "108 FDR<0.05 genes",
  Files = "GSE165489 tables"
)

# Cross-disease
summary_list[["Cross_disease"]] <- data.table(
  Analysis = "Cross-disease comparison",
  Component = "Gene-level overlap",
  Result = "1 shared gene (limited)",
  Files = "T14 cross-disease"
)

# Pathway convergence
summary_list[["Pathway"]] <- data.table(
  Analysis = "Pathway convergence",
  Component = "Hallmark + Reactome GSEA",
  Result = "36 concordant pathways",
  Files = "T15 pathway convergence"
)

# Single-cell
summary_list[["scRNA"]] <- data.table(
  Analysis = "Single-cell localization",
  Component = "Monocytes + Macrophages",
  Result = "94% + 83% positive for shared module",
  Files = "scRNA tables"
)

# T2DM x Mtb
summary_list[["T2DM_Mtb"]] <- data.table(
  Analysis = "Mechanistic bridge",
  Component = "GSE283452",
  Result = "906 interaction DEGs",
  Files = "GSE283452 tables"
)

summary_dt <- rbindlist(summary_list)
fwrite(summary_dt, file.path(final_tables, "SUMMARY_three_disease_integration.csv"), bom = TRUE)

cat("  Summary saved\n")

# ============================================================
# 3. Final file counts
# ============================================================
cat("\n3. Final file counts:\n")

tables <- list.files(final_tables, pattern = "\\.csv$")
figures <- list.files(final_figures, pattern = "\\.(svg|pdf|tiff|png)$")

cat("  Tables:", length(tables), "\n")
cat("  Figures:", length(figures), "\n")

# ============================================================
# 4. Key results summary
# ============================================================
cat("\n")
cat(paste(rep("=", 60), collapse=""), "\n")
cat("FINAL THREE-DISEASE INTEGRATION RESULTS\n")
cat(paste(rep("=", 60), collapse=""), "\n\n")

cat("1. GBD 2023:\n")
cat("   - 15 triple-high burden locations\n")
cat("   - TB + T2DM + Silicosis co-burden\n\n")

cat("2. TB-DM Meta-Analysis:\n")
cat("   - 2 cohorts (GSE114192, GSE181143)\n")
cat("   - 380 FDR-significant genes\n")
cat("   - 251 up, 129 down\n\n")

cat("3. Silicosis:\n")
cat("   - GSE165489 (n=86)\n")
cat("   - 108 FDR-significant genes\n")
cat("   - 99 up, 9 down\n\n")

cat("4. Cross-Disease:\n")
cat("   - Gene-level: Only 1 shared gene\n")
cat("   - Pathway-level: 36 concordant pathways\n")
cat("   - Key pathways: Neutrophil degranulation,\n")
cat("     Complement, ROS production, Heme metabolism\n\n")

cat("5. Single-Cell:\n")
cat("   - Monocytes: 94% positive for shared module\n")
cat("   - Macrophages: 83% positive for shared module\n\n")

cat("6. Mechanistic Bridge:\n")
cat("   - GSE283452: T2DM x Mtb in macrophages\n")
cat("   - 906 interaction DEGs\n\n")

cat("KEY FINDING:\n")
cat("TB-DM and Silicosis show LIMITED gene-level overlap\n")
cat("but STRONG pathway-level convergence in myeloid\n")
cat("immune programs.\n")

write_log("Final compilation completed")
