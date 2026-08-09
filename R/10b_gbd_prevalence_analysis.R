source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "ggplot2", "patchwork", "sf", "spdep",
                  "svglite", "ragg"))

cat("=== GBD 2023 Prevalence Analysis ===\n\n")

# ============================================================
# 1. Load and clean data
# ============================================================
cat("1. Loading prevalence data...\n")

# Load combined prevalence data
prev_file <- file.path(path_data, "01_GBD2023", "processed", "GBD2023_prevalence_combined.csv")
prev <- data.table::fread(prev_file, encoding = "UTF-8")

# Filter to 2023 and both sexes combined
prev_2023 <- prev[year == 2023 & sex_name == "合计"]
cat("  2023 data:", nrow(prev_2023), "rows\n")

# Map cause names (the CSV has encoding issues, use cause_id)
# 297 = TB, 976 = T2DM, 510 = Pneumoconiosis/Silicosis
prev_2023[, disease := NA_character_]
prev_2023[cause_id == 297, disease := "TB"]
prev_2023[cause_id == 976, disease := "T2DM"]
prev_2023[cause_id == 510, disease := "Silicosis"]

# Filter to diseases of interest
prev_diseases <- prev_2023[disease %in% c("TB", "T2DM", "Silicosis")]
cat("  Diseases:", paste(unique(prev_diseases$disease), collapse = ", "), "\n")
cat("  Locations:", length(unique(prev_diseases$location_name)), "\n")

# Reshape to wide format
prev_wide <- data.table::dcast(prev_diseases, location_id + location_name ~ disease,
                               value.var = c("val", "upper", "lower"))
cat("  Wide format:", nrow(prev_wide), "locations\n")

# ============================================================
# 2. Percentile normalization and combined burden index
# ============================================================
cat("\n2. Calculating combined burden index...\n")

# Use (rank - 0.5) / N to avoid 0 and 1
calc_percentile <- function(x) {
  (rank(x, na.last = "keep") - 0.5) / sum(!is.na(x))
}

prev_wide[, p_tb := calc_percentile(val_TB)]
prev_wide[, p_dm := calc_percentile(val_T2DM)]
prev_wide[, p_silicosis := calc_percentile(val_Silicosis)]

# Combined burden index (geometric mean)
prev_wide[, cbi := (p_tb * p_dm * p_silicosis)^(1/3)]

# Triple-high burden (>=75th percentile)
prev_wide[, triple_high := p_tb >= 0.75 & p_dm >= 0.75 & p_silicosis >= 0.75]

# Eight-type classification
prev_wide[, tb_high := p_tb >= 0.75]
prev_wide[, dm_high := p_dm >= 0.75]
prev_wide[, sil_high := p_silicosis >= 0.75]
prev_wide[, burden_class := "Low-Low-Low"]
prev_wide[tb_high & !dm_high & !sil_high, burden_class := "TB high only"]
prev_wide[!tb_high & dm_high & !sil_high, burden_class := "DM high only"]
prev_wide[!tb_high & !dm_high & sil_high, burden_class := "Silicosis high only"]
prev_wide[tb_high & dm_high & !sil_high, burden_class := "TB + DM high"]
prev_wide[tb_high & !dm_high & sil_high, burden_class := "TB + Silicosis high"]
prev_wide[!tb_high & dm_high & sil_high, burden_class := "DM + Silicosis high"]
prev_wide[tb_high & dm_high & sil_high, burden_class := "TB + DM + Silicosis high"]

cat("  Triple-high countries:", sum(prev_wide$triple_high, na.rm = TRUE), "\n")
cat("  Burden class distribution:\n")
print(table(prev_wide$burden_class))

# ============================================================
# 3. Save results
# ============================================================
cat("\n3. Saving results...\n")

result_dir <- file.path(path_result, "01_GBD2023", "tables")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

# Save combined burden index
data.table::fwrite(prev_wide, file.path(result_dir, "T01_GBD2023_prevalence_CBI.csv"), bom = TRUE)

# Save top 30 countries
top30 <- prev_wide[order(-cbi)][1:min(30, nrow(prev_wide))]
data.table::fwrite(top30, file.path(result_dir, "T02_top30_prevalence_CBI.csv"), bom = TRUE)

cat("  Top 10 countries by CBI:\n")
print(top30[1:min(10, nrow(top30)), .(location_name, val_TB, val_T2DM, val_Silicosis, cbi, burden_class)])

# ============================================================
# 4. Spatial analysis
# ============================================================
cat("\n4. Spatial analysis...\n")

# Load world map
world <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
# Use a simpler world map
tryCatch({
  world <- sf::st_read(file.path(path_data, "05_reference", "world_simple.shp"), quiet = TRUE)
}, error = function(e) {
  cat("  Using built-in world map\n")
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
})

# Merge with CBI data
# Need to match location names - this may need manual mapping
# For now, save the data for mapping

# ============================================================
# 5. Sensitivity analysis
# ============================================================
cat("\n5. Sensitivity analysis...\n")

# Threshold sensitivity
thresholds <- c(0.67, 0.75, 0.80, 0.90)
sens_results <- lapply(thresholds, function(thr) {
  n_triple <- sum(prev_wide$p_tb >= thr & prev_wide$p_dm >= thr & prev_wide$p_silicosis >= thr, na.rm = TRUE)
  data.table::data.table(threshold = thr, n_countries = n_triple)
})
sens_dt <- data.table::rbindlist(sens_results)
data.table::fwrite(sens_dt, file.path(result_dir, "T03_threshold_sensitivity_prevalence.csv"), bom = TRUE)
cat("  Threshold sensitivity:\n")
print(sens_dt)

# ============================================================
# 6. Summary
# ============================================================
cat("\n=== Summary ===\n")
cat("Total locations:", nrow(prev_wide), "\n")
cat("TB prevalence range:", range(prev_wide$val_TB, na.rm = TRUE), "\n")
cat("T2DM prevalence range:", range(prev_wide$val_T2DM, na.rm = TRUE), "\n")
cat("Silicosis prevalence range:", range(prev_wide$val_Silicosis, na.rm = TRUE), "\n")
cat("Triple-high countries (75th):", sum(prev_wide$triple_high, na.rm = TRUE), "\n")

write_log("GBD 2023 prevalence analysis completed: locations=", nrow(prev_wide),
          "; triple_high=", sum(prev_wide$triple_high, na.rm = TRUE))
