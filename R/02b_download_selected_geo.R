source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "curl"))

manifest_file <- file.path(path_data, "00_manifest", "M04_GEO_supplementary_files.csv")
if (!file.exists(manifest_file)) stop("Run R/02_download_geo.R --mode=supp-manifest first")
m <- data.table::fread(manifest_file, encoding = "UTF-8")
m <- m[!grepl("^filelist[.]txt$", file, ignore.case = TRUE)]
m[, selected_reason := "required supplementary source"]

prefer_only <- function(accession_value, pattern, reason) {
  idx <- which(m$accession == accession_value)
  m[idx, `:=`(selected = grepl(pattern, file, ignore.case = TRUE),
              selected_reason = reason)]
}
m[, selected := TRUE]
prefer_only("GSE181143", "count_data", "analysis-ready count matrix; RAW archive is redundant")
prefer_only("GSE193978", "paper1_rawdata", "analysis-ready longitudinal matrix")
prefer_only("GSE193979", "paper2_rawdata", "analysis-ready outcome matrix")
prefer_only("GSE249102", "Normalized_signal", "analysis-ready expression matrix; series matrix also retained")
prefer_only("GSE174725", "RAW[.]tar|barcodes|features|matrix[.]mtx",
            "10x components plus donor-resolved RAW archive required for pseudobulk")
prefer_only("GSE326212", "RAW[.]tar|metadata[.]csv|Table_1", "R-readable RAW plus metadata; 12-GiB h5ad excluded as redundant")

selected <- m[selected == TRUE]
selected[, layer := geo_catalog$layer[match(accession, geo_catalog$accession)]]
selected[, final := file.path(path_data, ifelse(layer == "bulk", "02_GEO_bulk", "03_scRNA"),
                              accession, "raw", file)]
selected[, staging := paste0(final, ".download")]
selected[, existing_complete := file.exists(final) & file.info(final)$size == bytes]
safe_write_csv(selected, file.path(path_data, "00_manifest", "M06_GEO_selected_downloads.csv"))

jobs <- selected[existing_complete == FALSE]
if (!nrow(jobs)) {
  write_log("All selected GEO supplementary files already complete")
  quit(save = "no", status = 0)
}

worker <- function(i, urls, staging, expected) {
  dest <- staging[i]
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(dest) && file.info(dest)$size == expected[i]) {
    return(data.frame(i = i, status = 0L, downloaded = 0, message = "staging already complete"))
  }
  ans <- tryCatch({
    curl::curl_download(urls[i], dest, quiet = TRUE, mode = "wb")
    actual <- if (file.exists(dest)) file.info(dest)$size else NA_real_
    if (is.na(actual) || actual != expected[i]) {
      data.frame(i = i, status = 2L, downloaded = actual,
                 message = paste0("size mismatch: expected ", expected[i], ", got ", actual))
    } else data.frame(i = i, status = 0L, downloaded = actual, message = "ok")
  }, error = function(e) data.frame(i = i, status = 1L, downloaded = NA_real_,
                                    message = conditionMessage(e)))
  ans
}

n_workers <- min(4L, nrow(jobs))
cl <- parallel::makePSOCKcluster(n_workers)
on.exit(parallel::stopCluster(cl), add = TRUE)
res_list <- parallel::parLapply(cl, seq_len(nrow(jobs)), worker,
                               urls = jobs$url, staging = jobs$staging, expected = jobs$bytes)
parallel::stopCluster(cl)
on.exit(NULL, add = FALSE)
res <- data.table::rbindlist(res_list)
res <- cbind(jobs[, .(accession, file, url, final, staging, expected_bytes = bytes)], res)

for (i in seq_len(nrow(res))) {
  if (res$status[i] == 0L && file.exists(res$staging[i]) && !file.exists(res$final[i])) {
    if (!file.rename(res$staging[i], res$final[i])) {
      res$status[i] <- 3L
      res$message[i] <- "download valid but final rename failed"
    }
  }
  res$final_bytes[i] <- if (file.exists(res$final[i])) file.info(res$final[i])$size else NA_real_
  res$verified[i] <- !is.na(res$final_bytes[i]) && res$final_bytes[i] == res$expected_bytes[i]
}
safe_write_csv(res, file.path(path_data, "00_manifest", "M07_GEO_selected_download_status.csv"))
write_log("Selected GEO download completed: verified ", sum(res$verified), "/", nrow(res),
          "; selected total GiB=", round(sum(selected$bytes) / 1024^3, 2))
if (!all(res$verified)) quit(save = "no", status = 2)
