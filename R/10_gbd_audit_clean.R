source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "dplyr", "tidyr"))

csv_candidates <- list.files(file.path(path_data, "gbd_raw"), pattern = "\\.csv$",
                             full.names = TRUE, ignore.case = TRUE)
zip_candidates <- list.files(path_data, pattern = "^IHME-GBD_2023_DATA-.*\\.zip$",
                             full.names = TRUE, ignore.case = TRUE)
source_tables <- list()
for (f in csv_candidates) {
  x <- data.table::fread(f, encoding = "UTF-8")
  x[, source_archive := basename(f)]
  source_tables[[length(source_tables) + 1L]] <- x
}
for (z in zip_candidates) {
  members <- utils::unzip(z, list = TRUE)$Name
  members <- members[grepl("\\.csv$", members, ignore.case = TRUE)]
  for (member in members) {
    extract_dir <- tempfile("gbd_zip_")
    dir.create(extract_dir)
    extracted <- utils::unzip(z, files = member, exdir = extract_dir)
    x <- data.table::fread(extracted, encoding = "UTF-8")
    x[, source_archive := paste0(basename(z), "::", member)]
    source_tables[[length(source_tables) + 1L]] <- x
    unlink(extract_dir, recursive = TRUE)
  }
}
if (!length(source_tables)) stop("No GBD CSV found in data/gbd_raw or IHME ZIP in data/")
d_all <- data.table::rbindlist(source_tables, fill = TRUE, use.names = TRUE)
core_cols <- setdiff(names(d_all), "source_archive")
d <- unique(d_all, by = core_cols)

required <- c("measure_name", "location_id", "location_name", "sex_name", "age_name",
              "cause_name", "metric_name", "year", "val", "upper", "lower")
missing_cols <- setdiff(required, names(d))
if (length(missing_cols)) stop("GBD file missing columns: ", paste(missing_cols, collapse = ", "))

inventory <- data.table::rbindlist(lapply(required, function(v) {
  x <- d[[v]]
  data.table::data.table(variable = v, n_unique = data.table::uniqueN(x),
                         n_missing = sum(is.na(x)), examples = paste(head(unique(x), 8), collapse = " | "))
}))
safe_write_csv(inventory, file.path(path_result, "00_audit", "A02_GBD_variable_inventory.csv"))

expected_primary <- data.frame(
  disease = c("Tuberculosis", "Type 2 diabetes mellitus", "Silicosis"),
  required_measure = c("Incidence", "Prevalence", "Prevalence"),
  required_sex = "Both",
  required_age = "Age-standardized",
  required_year = 2023,
  available_cause = c("结核病", "2型糖尿病", "尘肺病"),
  available_measure = c("发病率", "发病率", "发病率"),
  status = c("primary-compatible metric", "disease-specific incidence proxy", "cause-total incidence proxy")
)
safe_write_csv(expected_primary, file.path(path_result, "00_audit", "A03_GBD_protocol_compatibility.csv"))

subset <- d[year == 2023 & sex_name == "合计" & age_name == "年龄标准化" &
              metric_name == "率" & measure_name == "发病率" &
              cause_name %chin% c("结核病", "2型糖尿病", "尘肺病")]
if (nrow(subset) == 0) stop("No eligible both-sex age-standardized incidence rows")

key_counts <- subset[, .N, by = .(location_id, location_name, cause_name)]
if (any(key_counts$N != 1L)) stop("Non-unique location-cause rows detected")

wide <- data.table::dcast(subset, location_id + location_name ~ cause_name,
                          value.var = c("val", "lower", "upper"))
data.table::setnames(wide,
  old = c("val_结核病", "val_2型糖尿病", "val_尘肺病",
          "lower_结核病", "lower_2型糖尿病", "lower_尘肺病",
          "upper_结核病", "upper_2型糖尿病", "upper_尘肺病"),
  new = c("tb_asir", "t2dm_asir_proxy", "pneumoconiosis_asir_proxy",
          "tb_lower", "t2dm_lower", "pneumoconiosis_lower",
          "tb_upper", "t2dm_upper", "pneumoconiosis_upper"), skip_absent = TRUE)

needed_rates <- c("tb_asir", "t2dm_asir_proxy", "pneumoconiosis_asir_proxy")
if (anyNA(wide[, ..needed_rates])) stop("Incomplete three-cause incidence matrix; no imputation performed")
if (any(unlist(wide[, ..needed_rates]) <= 0)) stop("Non-positive rate encountered; log model undefined")

safe_write_csv(wide, file.path(path_data, "01_GBD2023", "processed", "GBD2023_proxy_rates.csv"))
safe_write_csv(wide, file.path(path_result, "01_GBD2023", "source_data", "SD01_GBD2023_proxy_rates.csv"))

summary_audit <- data.frame(
  item = c("source_files", "rows_after_exact_deduplication", "locations", "year", "sex", "age", "metric",
           "causes", "primary_metric_compatible", "missing_rate_cells", "exact_duplicate_rows"),
  value = c(paste(unique(d_all$source_archive), collapse = "; "), nrow(d), data.table::uniqueN(d$location_id), "2023", "Both",
            "Age-standardized", "Rate", paste(sort(unique(d$cause_name)), collapse = "; "),
            "No: T2DM and pneumoconiosis are incidence proxies; pneumoconiosis is not silicosis", sum(is.na(wide[, ..needed_rates])),
            sum(duplicated(d)))
)
safe_write_csv(summary_audit, file.path(path_result, "00_audit", "A04_GBD_data_audit_summary.csv"))
write_log("GBD audit/clean completed: ", nrow(wide), " locations; proxy flag retained")
