source(file.path("R", "00_config.R"), encoding = "UTF-8")
suppressPackageStartupMessages(library(data.table))

cat("=== GBD 2023 Figure 1: Full Audit & Rebuild ===\n\n")

# ============================================================
# STEP 1: Load and audit ALL raw GBD data
# ============================================================
cat("STEP 1: Loading raw GBD data for audit...\n")

# Check multiple possible locations for GBD data
raw_files <- c(
  list.files(file.path(path_data, "01_GBD2023", "raw"), pattern = "\\.csv$", full.names = TRUE),
  list.files(file.path(path_data, "gbd_raw"), pattern = "\\.csv$", full.names = TRUE),
  list.files(file.path(path_data, "gbd_prevalence_1"), pattern = "\\.csv$", full.names = TRUE),
  list.files(file.path(path_data, "gbd_prevalence_2"), pattern = "\\.csv$", full.names = TRUE),
  list.files(file.path(path_data), pattern = "IHME.*\\.csv$", full.names = TRUE)
)
raw_files <- unique(raw_files)
cat("Raw files found:", length(raw_files), "\n")

# Load all raw files
all_raw <- list()
for (f in raw_files) {
  dt <- fread(f, encoding = "UTF-8")
  dt[, source_file := basename(f)]
  all_raw[[basename(f)]] <- dt
  cat("  ", basename(f), ":", nrow(dt), "rows,", ncol(dt), "cols\n")
}

raw_all <- rbindlist(all_raw, fill = TRUE)

# ============================================================
# STEP 2: Complete indicator audit
# ============================================================
cat("\nSTEP 2: Complete indicator audit...\n")

# Create comprehensive audit table
audit_dt <- raw_all[, .(
  cause = cause_name,
  cause_id = cause_id,
  measure = measure_name,
  measure_id = measure_id,
  metric = metric_name,
  metric_id = metric_id,
  age = age_name,
  age_id = age_id,
  sex = sex_name,
  sex_id = sex_id,
  year = year,
  location = location_name,
  location_id = location_id,
  value = val,
  lower = lower,
  upper = upper,
  source_file = source_file
)]

# Save full audit
result_dir <- file.path(path_result, "01_GBD2023", "tables")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(audit_dt, file.path(result_dir, "T01_GBD2023_indicator_audit.csv"), bom = TRUE)

cat("  Total audit records:", nrow(audit_dt), "\n")
cat("  Unique causes:", paste(unique(audit_dt$cause), collapse = ", "), "\n")
cat("  Unique measures:", paste(unique(audit_dt$measure), collapse = ", "), "\n")
cat("  Unique metrics:", paste(unique(audit_dt$metric), collapse = ", "), "\n")
cat("  Years:", paste(unique(audit_dt$year), collapse = ", "), "\n")
cat("  Sex:", paste(unique(audit_dt$sex), collapse = ", "), "\n")
cat("  Age:", paste(unique(audit_dt$age), collapse = ", "), "\n")

# ============================================================
# STEP 3: TB Prevalence vs Incidence audit
# ============================================================
cat("\nSTEP 3: TB Prevalence vs Incidence audit...\n")

# Filter to TB data (Chinese: 结核病)
tb_data <- audit_dt[grepl("结核", cause, ignore.case = TRUE)]
cat("  TB records:", nrow(tb_data), "\n")
cat("  TB measures:", paste(unique(tb_data$measure), collapse = ", "), "\n")
cat("  TB metrics:", paste(unique(tb_data$metric), collapse = ", "), "\n")

# Check what's available
for (m in unique(tb_data$measure)) {
  for (met in unique(tb_data$metric)) {
    sub <- tb_data[measure == m & metric == met & sex == "Both" & age == "Age-standardized"]
    if (nrow(sub) > 0) {
      cat("\n  Measure:", m, ", Metric:", met, "\n")
      cat("    N locations:", nrow(sub), "\n")
      cat("    Year:", unique(sub$year), "\n")
      cat("    Value range:", range(sub$value, na.rm = TRUE), "\n")
      cat("    Median:", median(sub$value, na.rm = TRUE), "\n")
      cat("    IQR:", quantile(sub$value, 0.25, na.rm = TRUE), "-", quantile(sub$value, 0.75, na.rm = TRUE), "\n")
      cat("    P75:", quantile(sub$value, 0.75, na.rm = TRUE), "\n")
      cat("    P90:", quantile(sub$value, 0.90, na.rm = TRUE), "\n")
      cat("    Max:", max(sub$value, na.rm = TRUE), "\n")

      # Top 20 locations
      top20 <- sub[order(-value)][1:min(20, .N)]
      cat("    Top 20:\n")
      print(top20[, .(location, value, lower, upper)])
    }
  }
}

# Save TB prevalence audit (Chinese: 患病率)
tb_prev <- tb_data[measure == "患病率" & metric == "率" & sex == "合计" & age == "年龄标准化" & year == 2023]
if (nrow(tb_prev) > 0) {
  fwrite(tb_prev, file.path(result_dir, "T02_TB_prevalence_audit.csv"), bom = TRUE)
  cat("\n  TB Prevalence audit saved:", nrow(tb_prev), "locations\n")
  cat("    Value range:", range(tb_prev$value, na.rm = TRUE), "\n")
  cat("    Median:", median(tb_prev$value, na.rm = TRUE), "\n")
  cat("    P75:", quantile(tb_prev$value, 0.75, na.rm = TRUE), "\n")
  cat("    P90:", quantile(tb_prev$value, 0.90, na.rm = TRUE), "\n")
  cat("    Max:", max(tb_prev$value, na.rm = TRUE), "\n")
} else {
  cat("\n  TB Prevalence (Rate, Age-standardized, Both, 2023) NOT FOUND\n")
}

# Save TB incidence audit (Chinese: 发病率)
tb_inc <- tb_data[measure == "发病率" & metric == "率" & sex == "合计" & age == "年龄标准化" & year == 2023]
if (nrow(tb_inc) > 0) {
  fwrite(tb_inc, file.path(result_dir, "T03_TB_incidence_audit.csv"), bom = TRUE)
  cat("  TB Incidence audit saved:", nrow(tb_inc), "locations\n")
  cat("    Value range:", range(tb_inc$value, na.rm = TRUE), "\n")
  cat("    Median:", median(tb_inc$value, na.rm = TRUE), "\n")
  cat("    P75:", quantile(tb_inc$value, 0.75, na.rm = TRUE), "\n")
  cat("    P90:", quantile(tb_inc$value, 0.90, na.rm = TRUE), "\n")
  cat("    Max:", max(tb_inc$value, na.rm = TRUE), "\n")
}

# ============================================================
# STEP 4: T2DM and Silicosis audit
# ============================================================
cat("\nSTEP 4: T2DM and Silicosis audit...\n")

# T2DM (Chinese: 糖尿病, 2型糖尿病)
t2dm_data <- audit_dt[grepl("糖尿病", cause, ignore.case = TRUE)]
cat("  T2DM records:", nrow(t2dm_data), "\n")
cat("  T2DM measures:", paste(unique(t2dm_data$measure), collapse = ", "), "\n")

t2dm_prev <- t2dm_data[measure == "患病率" & metric == "率" & sex == "合计" & age == "年龄标准化" & year == 2023]
if (nrow(t2dm_prev) > 0) {
  fwrite(t2dm_prev, file.path(result_dir, "T04_T2DM_audit.csv"), bom = TRUE)
  cat("  T2DM Prevalence audit saved:", nrow(t2dm_prev), "locations\n")
  cat("    Value range:", range(t2dm_prev$value, na.rm = TRUE), "\n")
  cat("    Median:", median(t2dm_prev$value, na.rm = TRUE), "\n")
}

# Silicosis (Chinese: 尘肺病)
sil_data <- audit_dt[grepl("尘肺", cause, ignore.case = TRUE)]
cat("  Silicosis records:", nrow(sil_data), "\n")
cat("  Silicosis measures:", paste(unique(sil_data$measure), collapse = ", "), "\n")

sil_prev <- sil_data[measure == "患病率" & metric == "率" & sex == "合计" & age == "年龄标准化" & year == 2023]
if (nrow(sil_prev) > 0) {
  fwrite(sil_prev, file.path(result_dir, "T05_silicosis_audit.csv"), bom = TRUE)
  cat("  Silicosis Prevalence audit saved:", nrow(sil_prev), "locations\n")
  cat("    Value range:", range(sil_prev$value, na.rm = TRUE), "\n")
  cat("    Median:", median(sil_prev$value, na.rm = TRUE), "\n")
} else {
  cat("  Silicosis Prevalence (Rate, Age-standardized, Both, 2023) NOT FOUND\n")
  # Check what's available
  cat("  Available silicosis data:\n")
  print(unique(sil_data[, .(measure, metric, age, sex, year)]))
}

# ============================================================
# STEP 5: Decision on TB indicator
# ============================================================
cat("\nSTEP 5: TB indicator decision...\n")

if (nrow(tb_prev) > 0) {
  max_prev <- max(tb_prev$value, na.rm = TRUE)
  cat("  TB Prevalence max:", max_prev, "per 100,000\n")

  if (max_prev > 5000) {
    cat("  WARNING: TB prevalence max > 5,000/100k is suspicious for active TB\n")
    cat("  This may indicate LTBI data or unit error\n")
    cat("  Checking TB incidence instead...\n")

    if (nrow(tb_inc) > 0) {
      max_inc <- max(tb_inc$value, na.rm = TRUE)
      cat("  TB Incidence max:", max_inc, "per 100,000\n")
      cat("  Decision: USE TB INCIDENCE (more reliable for active TB)\n")
      tb_indicator <- "Incidence"
      tb_indicator_name <- "Age-standardized TB incidence rate per 100,000"
    } else {
      cat("  TB Incidence NOT available\n")
      cat("  Decision: USE TB PREVALENCE with caveat\n")
      tb_indicator <- "Prevalence"
      tb_indicator_name <- "Age-standardized TB prevalence rate per 100,000"
    }
  } else {
    cat("  TB prevalence values appear reasonable\n")
    cat("  Decision: USE TB PREVALENCE\n")
    tb_indicator <- "Prevalence"
    tb_indicator_name <- "Age-standardized TB prevalence rate per 100,000"
  }
} else {
  cat("  TB Prevalence not available, using Incidence\n")
  tb_indicator <- "Incidence"
  tb_indicator_name <- "Age-standardized TB incidence rate per 100,000"
}

cat("\n  FINAL TB INDICATOR:", tb_indicator, "\n")
cat("  Figure 1A title:", tb_indicator_name, "\n")

write_log("GBD audit completed: TB indicator =", tb_indicator)
