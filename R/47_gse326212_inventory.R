source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table"))

acc <- "GSE326212"
raw_dir <- file.path(path_data, "03_scRNA", acc, "raw")
result_dir <- file.path(path_result, "04_scRNA", acc)
invisible(lapply(file.path(result_dir, c("tables", "source_data", "models", "figures")), dir.create,
  recursive = TRUE, showWarnings = FALSE))

tar_file <- file.path(raw_dir, paste0(acc, "_RAW.tar"))
meta_file <- file.path(raw_dir, paste0(acc, "_scRNAseq_tb_bal_metadata.csv.gz"))
if (!file.exists(tar_file) || !file.exists(meta_file)) stop("GSE326212 local inputs are incomplete")
members <- utils::untar(tar_file, list = TRUE)
h5 <- members[grepl("[.]h5$", members, ignore.case = TRUE)]
manifest <- data.table::data.table(
  archive_member = h5,
  gsm = sub("_.*$", "", h5),
  sample_code = sub("^GSM[0-9]+_", "", sub("_filtered[.]h5$", "", h5)),
  suffix = sub("^.*_([A-Za-z])$", "\\1", sub("_filtered[.]h5$", "", h5))
)
safe_write_csv(manifest, file.path(result_dir, "source_data", "SD01_GSE326212_h5_manifest.csv"))

write_log("Reading GSE326212 author cell metadata")
meta <- data.table::fread(meta_file, showProgress = TRUE)
name_audit <- data.table::data.table(column_index = seq_along(names(meta)), column_name = names(meta),
  storage_type = vapply(meta, typeof, character(1)), missing = vapply(meta, function(x) sum(is.na(x)), integer(1)),
  unique_values = vapply(meta, data.table::uniqueN, integer(1)))
safe_write_csv(name_audit, file.path(result_dir, "tables", "T01_GSE326212_metadata_column_audit.csv"))

candidate_pattern <- "sample|patient|donor|subject|group|status|disease|outcome|progress|time|visit|cluster|cell|type|annot"
candidate_cols <- names(meta)[grepl(candidate_pattern, names(meta), ignore.case = TRUE)]
candidate_summary <- data.table::rbindlist(lapply(candidate_cols, function(nm) {
  x <- as.character(meta[[nm]])
  tab <- sort(table(x, useNA = "ifany"), decreasing = TRUE)
  data.table::data.table(column = nm, value = names(tab)[seq_len(min(50L, length(tab)))],
    cells = as.integer(tab[seq_len(min(50L, length(tab)))]))
}), fill = TRUE)
safe_write_csv(candidate_summary, file.path(result_dir, "tables", "T02_GSE326212_metadata_candidate_values.csv"))

summary_dt <- data.table::data.table(
  archive_bytes = file.info(tar_file)$size,
  h5_files = length(h5), gsm_records = data.table::uniqueN(manifest$gsm),
  metadata_rows = nrow(meta), metadata_columns = ncol(meta),
  candidate_metadata_columns = length(candidate_cols)
)
safe_write_csv(summary_dt, file.path(result_dir, "tables", "T03_GSE326212_inventory_summary.csv"))
safe_write_csv(meta[seq_len(min(100L, nrow(meta))), ..candidate_cols],
  file.path(result_dir, "source_data", "SD02_GSE326212_metadata_candidate_preview.csv"))
write_log("GSE326212 inventory completed: h5 files=", length(h5), "; GSM records=",
  data.table::uniqueN(manifest$gsm), "; metadata rows=", nrow(meta), "; columns=", ncol(meta),
  "; candidate fields=", paste(candidate_cols, collapse = ", "))
