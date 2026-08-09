source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "rhdf5"))

acc <- "GSE326212"
raw_dir <- file.path(path_data, "03_scRNA", acc, "raw")
result_dir <- file.path(path_result, "04_scRNA", acc)
extract_dir <- file.path(raw_dir, "selected_h5")
dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)
meta <- data.table::fread(file.path(raw_dir, paste0(acc, "_scRNAseq_tb_bal_metadata.csv.gz")))
meta[, time_code := sub("^.*_", "", sample_id)]
meta[, analysis_group := NA_character_]
meta[grepl("^TB_neut_", status), analysis_group := "Active TB"]

controller <- meta[status == "IGRApos_PETneg"]
controller[, priority := data.table::fcase(time_code == "B", 1L, time_code == "F", 2L, default = 9L)]
controller <- controller[order(subject_id, priority)][, .SD[1], by = subject_id]
controller[, analysis_group := "Stable controller"]

progressor <- meta[grepl("^Progressor_", status) & time_code %in% c("B", "F")]
progressor[, priority := data.table::fcase(time_code == "F", 1L, time_code == "B", 2L, default = 9L)]
progressor <- progressor[order(subject_id, priority)][, .SD[1], by = subject_id]
progressor[, analysis_group := "Pre-progression"]

selected <- data.table::rbindlist(list(meta[analysis_group == "Active TB"], controller, progressor),
  fill = TRUE, use.names = TRUE)
selected[, priority := NULL]
data.table::setorder(selected, analysis_group, subject_id)
if (nrow(selected) != 33L || data.table::uniqueN(selected$subject_id) != 33L) {
  stop("Selection must contain 33 unique subjects; found samples=", nrow(selected),
    ", subjects=", data.table::uniqueN(selected$subject_id))
}
expected_groups <- c("Active TB" = 21L, "Stable controller" = 6L, "Pre-progression" = 6L)
if (!all(table(selected$analysis_group)[names(expected_groups)] == expected_groups)) stop("Unexpected group counts")

tar_file <- file.path(raw_dir, paste0(acc, "_RAW.tar"))
members <- utils::untar(tar_file, list = TRUE)
selected[, archive_member := members[match(paste0(sample_id, "_filtered.h5"),
  sub("^GSM[0-9]+_", "", members))]]
if (anyNA(selected$archive_member)) stop("Selected H5 member matching failed")
missing <- selected$archive_member[!file.exists(file.path(extract_dir, selected$archive_member))]
if (length(missing)) {
  write_log("Extracting ", length(missing), " selected GSE326212 H5 files")
  utils::untar(tar_file, files = missing, exdir = extract_dir)
}
selected[, local_h5 := normalizePath(file.path(extract_dir, archive_member), winslash = "/", mustWork = TRUE)]
selected[, selection_rule := data.table::fcase(
  analysis_group == "Active TB", "All active TB baseline BAL samples",
  analysis_group == "Stable controller", "One sample per subject; baseline preferred over follow-up",
  analysis_group == "Pre-progression", "One pre-diagnosis sample per subject; follow-up preferred over baseline")]
safe_write_csv(selected, file.path(result_dir, "tables", "T05_GSE326212_primary_sample_selection.csv"))

h5_audit <- rhdf5::h5ls(selected$local_h5[1], recursive = TRUE)
safe_write_csv(data.table::as.data.table(h5_audit),
  file.path(result_dir, "tables", "T06_GSE326212_10x_h5_structure.csv"))
write_log("GSE326212 primary selection completed: ", paste(names(expected_groups), expected_groups,
  sep = "=", collapse = "; "), "; unique subjects=", data.table::uniqueN(selected$subject_id))
