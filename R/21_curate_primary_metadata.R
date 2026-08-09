source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table"))

read_meta <- function(accession) {
  f <- file.path(path_result, "02_bulk", accession, paste0("A_", accession, "_sample_metadata_raw.csv"))
  if (!file.exists(f)) stop("Missing raw metadata audit for ", accession, "; run R/20_geo_metadata_audit.R")
  data.table::fread(f, encoding = "UTF-8")
}

curate_gse114192 <- function() {
  x <- read_meta("GSE114192")
  group_col <- grep("disease.state", names(x), value = TRUE)[1]
  site_col <- grep("field_site", names(x), value = TRUE)[1]
  subject_col <- grep("patient.id", names(x), value = TRUE)[1]
  if (anyNA(c(group_col, site_col, subject_col))) stop("GSE114192 metadata fields not found")
  out <- x[, .(sample_id, subject_id = as.character(get(subject_col)),
               group_raw = get(group_col), site = get(site_col), title, source_name)]
  out[, group := data.table::fcase(
    group_raw == "TB_DM", "TB_DM",
    group_raw == "TB_only", "TB_only",
    group_raw == "TB_IH", "TB_IH",
    group_raw == "DM_only", "DM_only",
    group_raw %in% c("Healthy", "HC", "Control", "Healthy_Control"), "Healthy",
    group_raw == "IH", "IH",
    default = NA_character_)]
  out[, include_primary := group %in% c("TB_DM", "TB_only")]
  out[, exclusion_reason := data.table::fifelse(is.na(group), "unmapped group",
    data.table::fifelse(include_primary, "", "not in primary contrast"))]
  out
}

curate_gse165489 <- function() {
  x <- read_meta("GSE165489")
  rad_col <- grep("^radiography$", names(x), value = TRUE)[1]
  if (is.na(rad_col)) stop("GSE165489 radiography field not found")
  out <- x[, .(sample_id, subject_id = sample_id, group_raw = get(rad_col),
               exposure_years = suppressWarnings(as.numeric(get(grep("silica.exposure", names(x), value = TRUE)[1]))),
               title, source_name)]
  raw_low <- tolower(paste(out$group_raw, out$source_name))
  out[, group := data.table::fcase(
    grepl("silicosis", raw_low) & !grepl("without|non[- ]?silicosis", raw_low), "silicosis",
    grepl("without|non[- ]?silicosis|exposed.*no", raw_low), "exposed_no_silicosis",
    grepl("healthy|unexposed|control", raw_low), "unexposed",
    default = NA_character_)]
  out[, include_primary := group %in% c("silicosis", "exposed_no_silicosis")]
  out[, exclusion_reason := data.table::fifelse(is.na(group), "unmapped group",
    data.table::fifelse(include_primary, "", "not in primary contrast"))]
  out
}

datasets <- list(GSE114192 = curate_gse114192(), GSE165489 = curate_gse165489())
audit <- list()
for (acc in names(datasets)) {
  out <- datasets[[acc]]
  target_data <- file.path(path_data, "02_GEO_bulk", acc, "processed")
  target_result <- file.path(path_result, "02_bulk", acc)
  dir.create(target_data, recursive = TRUE, showWarnings = FALSE)
  dir.create(target_result, recursive = TRUE, showWarnings = FALSE)
  safe_write_csv(out, file.path(target_data, paste0(acc, "_metadata_curated.csv")))
  safe_write_csv(out, file.path(target_result, paste0("T_", acc, "_metadata_curated.csv")))
  out[, accession := acc]
  audit[[acc]] <- out[, .(n_samples = .N, n_subjects = data.table::uniqueN(subject_id),
                          n_in_primary = sum(include_primary), n_missing_group = sum(is.na(group))),
                      by = .(accession, group)]
}
audit_table <- data.table::rbindlist(audit, fill = TRUE)
safe_write_csv(audit_table, file.path(path_result, "00_audit", "A08_primary_metadata_group_counts.csv"))
write_log("Primary metadata curation completed for GSE114192 and GSE165489")
