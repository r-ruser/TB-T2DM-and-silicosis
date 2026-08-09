source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table"))

manifest_file <- file.path(path_data, "00_manifest", "M06_GEO_selected_downloads.csv")
if (!file.exists(manifest_file)) stop("Selected-download manifest missing")
m <- data.table::fread(manifest_file, encoding = "UTF-8")
audit <- list()
for (i in seq_len(nrow(m))) {
  final <- m$final[i]
  candidates <- c(paste0(final, ".part"), paste0(final, ".download"),
                  paste0(final, ".aria_test"), paste0(final, ".aria"))
  candidates <- candidates[file.exists(candidates)]
  action <- "none"
  source <- NA_character_
  if (!file.exists(final) && length(candidates)) {
    sizes <- file.info(candidates)$size
    valid <- which(sizes == m$bytes[i])
    if (length(valid)) {
      source <- candidates[valid[1]]
      if (file.rename(source, final)) action <- "finalized" else action <- "rename_failed"
    }
  }
  final_bytes <- if (file.exists(final)) file.info(final)$size else NA_real_
  audit[[i]] <- data.table::data.table(accession = m$accession[i], file = m$file[i],
    source = source, final = final, expected_bytes = m$bytes[i], final_bytes = final_bytes,
    verified = !is.na(final_bytes) && final_bytes == m$bytes[i], action = action)
}
audit <- data.table::rbindlist(audit)
safe_write_csv(audit, file.path(path_data, "00_manifest", "M08_GEO_finalize_audit.csv"))
write_log("Verified staging finalization: ", sum(audit$verified), "/", nrow(audit), " selected files complete")
