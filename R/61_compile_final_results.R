source(file.path("R", "00_config.R"), encoding = "UTF-8")

cat("=== Compiling Final Results ===\n\n")

final_tables <- file.path(path_result, "06_final", "tables")
final_figures <- file.path(path_result, "06_final", "figures")
dir.create(final_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(final_figures, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. Copy GBD results
# ============================================================
cat("1. Copying GBD results...\n")

# Tables
gbd_tables <- c(
  "T01_GBD2023_prevalence_CBI.csv",
  "T02_top30_prevalence_CBI.csv",
  "T03_threshold_sensitivity_prevalence.csv",
  "T04_prevalence_summary.csv"
)

for (f in gbd_tables) {
  src <- file.path(path_result, "01_GBD2023", "tables", f)
  if (file.exists(src)) {
    file.copy(src, file.path(final_tables, paste0("GBD_", f)), overwrite = TRUE)
    cat("  Copied:", f, "\n")
  }
}

# Figures
gbd_figs <- c(
  "F01_GBD2023_world_maps",
  "F02_GBD2023_additional_analysis"
)

for (base in gbd_figs) {
  for (ext in c(".svg", ".pdf", ".tiff", ".png")) {
    src <- file.path(path_result, "01_GBD2023", "figures", paste0(base, ext))
    if (file.exists(src)) {
      file.copy(src, file.path(final_figures, paste0("GBD_", base, ext)), overwrite = TRUE)
    }
  }
  cat("  Copied:", base, "\n")
}

# ============================================================
# 2. Copy Bulk transcriptomics results
# ============================================================
cat("\n2. Copying Bulk results...\n")

bulk_datasets <- c("GSE114192", "GSE165489", "GSE181143", "GSE193978", "GSE193979", "GSE249102", "GSE283452")

for (ds in bulk_datasets) {
  ds_dir <- file.path(path_result, "02_bulk", ds, "tables")
  if (dir.exists(ds_dir)) {
    files <- list.files(ds_dir, pattern = "^T\\d+.*\\.csv$")
    for (f in files) {
      src <- file.path(ds_dir, f)
      dst <- file.path(final_tables, paste0(ds, "_", f))
      file.copy(src, dst, overwrite = TRUE)
    }
    cat("  Copied", length(files), "files from", ds, "\n")
  }
}

# ============================================================
# 3. Copy scRNA results
# ============================================================
cat("\n3. Copying scRNA results...\n")

sc_datasets <- c("GSE174725", "GSE192483", "GSE326212", "cross_disease")

for (ds in sc_datasets) {
  ds_dir <- file.path(path_result, "04_scRNA", ds, "tables")
  if (dir.exists(ds_dir)) {
    files <- list.files(ds_dir, pattern = "^T\\d+.*\\.csv$")
    for (f in files) {
      src <- file.path(ds_dir, f)
      dst <- file.path(final_tables, paste0("scRNA_", ds, "_", f))
      file.copy(src, dst, overwrite = TRUE)
    }
    cat("  Copied", length(files), "files from", ds, "\n")
  }
}

# scRNA figures
for (ds in c("GSE174725", "GSE192483", "GSE326212", "cross_disease")) {
  fig_dir <- file.path(path_result, "04_scRNA", ds, "figures")
  if (dir.exists(fig_dir)) {
    files <- list.files(fig_dir, pattern = "\\.(svg|pdf|tiff|png)$")
    for (f in files) {
      src <- file.path(fig_dir, f)
      dst <- file.path(final_figures, paste0("scRNA_", ds, "_", f))
      file.copy(src, dst, overwrite = TRUE)
    }
    cat("  Copied", length(files), "figures from", ds, "\n")
  }
}

# ============================================================
# 4. Copy Cross-disease results
# ============================================================
cat("\n4. Copying Cross-disease results...\n")

cross_dirs <- c("RRHO", "GSEA", "WGCNA")
for (cd in cross_dirs) {
  cd_path <- file.path(path_result, "03_cross_disease", cd)
  if (dir.exists(cd_path)) {
    # Tables
    tbl_dir <- file.path(cd_path, "tables")
    if (dir.exists(tbl_dir)) {
      files <- list.files(tbl_dir, pattern = "^T\\d+.*\\.csv$")
      for (f in files) {
        src <- file.path(tbl_dir, f)
        dst <- file.path(final_tables, paste0("Cross_", cd, "_", f))
        file.copy(src, dst, overwrite = TRUE)
      }
      cat("  Copied", length(files), "tables from", cd, "\n")
    }
    # Figures
    fig_dir <- file.path(cd_path, "figures")
    if (dir.exists(fig_dir)) {
      files <- list.files(fig_dir, pattern = "\\.(svg|pdf|tiff|png)$")
      for (f in files) {
        src <- file.path(fig_dir, f)
        dst <- file.path(final_figures, paste0("Cross_", cd, "_", f))
        file.copy(src, dst, overwrite = TRUE)
      }
      cat("  Copied", length(files), "figures from", cd, "\n")
    }
  }
}

# ============================================================
# 5. Copy Virtual KO results
# ============================================================
cat("\n5. Copying Virtual KO results...\n")

ko_dir <- file.path(path_result, "05_GRN_KO", "tables")
if (dir.exists(ko_dir)) {
  files <- list.files(ko_dir, pattern = "^T\\d+.*\\.csv$")
  for (f in files) {
    src <- file.path(ko_dir, f)
    dst <- file.path(final_tables, paste0("KO_", f))
    file.copy(src, dst, overwrite = TRUE)
  }
  cat("  Copied", length(files), "files\n")
}

# ============================================================
# 6. Generate summary table
# ============================================================
cat("\n6. Generating summary table...\n")

summary_list <- list()

# GBD summary
summary_list[["GBD"]] <- data.table::data.table(
  Analysis = "GBD 2023",
  Dataset = "GBD 2023",
  Metric = "Triple-high burden countries",
  Value = "15",
  Notes = "75th percentile threshold, 204 locations"
)

# Bulk DEG counts
for (ds in bulk_datasets) {
  ds_dir <- file.path(path_result, "02_bulk", ds, "tables")
  if (dir.exists(ds_dir)) {
    # Try to find the main DEG file
    deg_files <- list.files(ds_dir, pattern = "DESeq2|_DE\\.csv|_primary", full.names = TRUE)
    if (length(deg_files) > 0) {
      tryCatch({
        dt <- data.table::fread(deg_files[1])
        n_deg <- sum(dt$padj < 0.05, na.rm = TRUE)
        summary_list[[ds]] <- data.table::data.table(
          Analysis = "Bulk transcriptomics",
          Dataset = ds,
          Metric = "DEGs (FDR<0.05)",
          Value = as.character(n_deg),
          Notes = paste(nrow(dt), "genes tested")
        )
      }, error = function(e) NULL)
    }
  }
}

# scRNA summary
summary_list[["scRNA_GSE174725"]] <- data.table::data.table(
  Analysis = "Single-cell",
  Dataset = "GSE174725",
  Metric = "Cells",
  Value = "20,646",
  Notes = "Silicosis BALF, 10 cell types"
)

summary_list[["scRNA_GSE326212"]] <- data.table::data.table(
  Analysis = "Single-cell",
  Dataset = "GSE326212",
  Metric = "QC completed",
  Value = "Yes",
  Notes = "TB BAL, in progress"
)

# Combine summary
if (length(summary_list) > 0) {
  summary_dt <- data.table::rbindlist(summary_list, fill = TRUE)
  data.table::fwrite(summary_dt, file.path(final_tables, "SUMMARY_all_analyses.csv"), bom = TRUE)
  cat("Summary table saved\n")
  print(summary_dt)
}

# ============================================================
# 7. List all files
# ============================================================
cat("\n=== Final files ===\n")
cat("\nTables:\n")
print(list.files(final_tables, pattern = "\\.csv$"))
cat("\nFigures:\n")
print(list.files(final_figures, pattern = "\\.(svg|pdf|tiff|png)$"))

write_log("Final results compiled to 06_final/")
