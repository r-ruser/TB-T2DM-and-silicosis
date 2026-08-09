source(file.path("R", "00_config.R"), encoding = "UTF-8")
suppressPackageStartupMessages(library(data.table))

cat("=== Figure Reorganization ===\n\n")

fig_dir <- file.path(path_result, "06_final", "figures")
supp_dir <- file.path(path_result, "06_final", "figures", "supplementary")
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Main Figures (Fig 1-5)
# ============================================================
cat("1. Main Figures:\n")
cat("   Fig 1: GBD Global Co-burden\n")
cat("   Fig 2: Cross-disease Transcriptomic Convergence\n")
cat("   Fig 3: Cross-cohort TB-DM Meta-analysis\n")
cat("   Fig 4: Single-cell Localization\n")
cat("   Fig 5: T2DM x Mtb Macrophage Interaction\n\n")

# Rename main figures
main_figures <- list(
  # Fig 1: GBD
  list(src = "GBD_F01_GBD2023_world_maps", dst = "Fig1_GBD_global_coburden"),
  list(src = "GBD_F02_GBD2023_additional_analysis", dst = "Fig1_GBD_coburden_index"),

  # Fig 2: Cross-disease convergence
  list(src = "Cross_RRHO_F09_cross_disease_RRHO_pathway_concordance", dst = "Fig2_RRHO_pathway_concordance"),
  list(src = "F10_pathway_convergence", dst = "Fig2_pathway_convergence"),

  # Fig 3: Meta-analysis
  list(src = "F11_meta_analysis", dst = "Fig3_TBDM_meta_analysis"),

  # Fig 4: Single-cell localization
  list(src = "scRNA_cross_disease_F17_RRHO_local_signature_single_cell", dst = "Fig4_single_cell_localization")
)

for (fig in main_figures) {
  for (ext in c("svg", "pdf", "tiff", "png")) {
    src_file <- file.path(fig_dir, paste0(fig$src, ".", ext))
    dst_file <- file.path(fig_dir, paste0(fig$dst, ".", ext))
    if (file.exists(src_file)) {
      file.rename(src_file, dst_file)
    }
  }
  cat("   Renamed:", fig$src, "->", fig$dst, "\n")
}

# ============================================================
# Supplementary Figures
# ============================================================
cat("\n2. Supplementary Figures:\n")

# Move QC and atlas figures to supplementary
supp_figures <- c(
  "scRNA_GSE174725_F06_GSE174725_QC",
  "scRNA_GSE174725_F12_GSE174725_cell_atlas",
  "scRNA_GSE174725_F13_GSE174725_pseudobulk_DE",
  "scRNA_GSE192483_F14_GSE192483_QC",
  "scRNA_GSE192483_F15_GSE192483_cell_atlas",
  "scRNA_GSE192483_F16_GSE192483_paired_pseudobulk_DE",
  "scRNA_GSE326212_F18_GSE326212_QC_summary"
)

# Rename supplementary figures
supp_mapping <- list(
  list(src = "scRNA_GSE174725_F06_GSE174725_QC", dst = "S1_GSE174725_QC"),
  list(src = "scRNA_GSE174725_F12_GSE174725_cell_atlas", dst = "S2_GSE174725_cell_atlas"),
  list(src = "scRNA_GSE174725_F13_GSE174725_pseudobulk_DE", dst = "S3_GSE174725_pseudobulk"),
  list(src = "scRNA_GSE192483_F14_GSE192483_QC", dst = "S4_GSE192483_QC"),
  list(src = "scRNA_GSE192483_F15_GSE192483_cell_atlas", dst = "S5_GSE192483_cell_atlas"),
  list(src = "scRNA_GSE192483_F16_GSE192483_paired_pseudobulk_DE", dst = "S6_GSE192483_pseudobulk"),
  list(src = "scRNA_GSE326212_F18_GSE326212_QC_summary", dst = "S7_GSE326212_QC")
)

for (fig in supp_mapping) {
  for (ext in c("svg", "pdf", "tiff", "png")) {
    src_file <- file.path(fig_dir, paste0(fig$src, ".", ext))
    dst_file <- file.path(supp_dir, paste0(fig$dst, ".", ext))
    if (file.exists(src_file)) {
      file.rename(src_file, dst_file)
    }
  }
  cat("   Moved:", fig$src, "-> supplementary/", fig$dst, "\n")
}

# ============================================================
# Summary
# ============================================================
cat("\n3. Final Structure:\n")
cat("   Main Figures (Fig1-5):\n")
main_files <- list.files(fig_dir, pattern = "^Fig[1-5]")
print(main_files)

cat("\n   Supplementary Figures:\n")
supp_files <- list.files(supp_dir, pattern = "^S[0-9]")
print(supp_files)

cat("\nFigure reorganization completed!\n")
