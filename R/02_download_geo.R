source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "curl"))

args <- commandArgs(trailingOnly = TRUE)
mode_arg <- grep("^--mode=", args, value = TRUE)
mode <- if (length(mode_arg)) sub("^--mode=", "", mode_arg[[1]]) else "manifest"
if (!mode %in% c("manifest", "matrix", "supp-manifest", "supp", "supp-parallel")) {
  stop("mode must be manifest, matrix, supp-manifest, supp, or supp-parallel")
}
parallel_jobs <- list()

geo_base_dir <- function(layer, accession) {
  file.path(path_data, if (layer == "bulk") "02_GEO_bulk" else "03_scRNA", accession)
}

series_prefix <- function(acc) paste0(substr(acc, 1, nchar(acc) - 3L), "nnn")
matrix_url <- function(acc) sprintf(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/%s_series_matrix.txt.gz",
  series_prefix(acc), acc, acc
)
supp_index_url <- function(acc) sprintf(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/",
  series_prefix(acc), acc
)

download_one <- function(url, dest) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(dest) && file.info(dest)$size > 0) {
    write_log("Skip existing: ", dest)
    return(TRUE)
  }
  tmp <- paste0(dest, ".part")
  if (file.exists(tmp)) unlink(tmp)
  ok <- tryCatch({
    curl::curl_download(url, tmp, quiet = FALSE, mode = "wb")
    if (!file.exists(tmp) || file.info(tmp)$size <= 0) stop("empty download")
    if (!file.rename(tmp, dest)) stop("cannot rename temporary file")
    TRUE
  }, error = function(e) {
    write_log("Download failed: ", url, " :: ", conditionMessage(e))
    FALSE
  })
  ok
}

manifest <- geo_catalog
manifest$matrix_url <- vapply(manifest$accession, matrix_url, character(1))
manifest$supp_index_url <- vapply(manifest$accession, supp_index_url, character(1))
manifest$download_mode <- mode
manifest$checked_at <- timestamp()
safe_write_csv(manifest, file.path(path_data, "00_manifest", "M02_GEO_series_catalog.csv"))

for (i in seq_len(nrow(geo_catalog))) {
  row <- geo_catalog[i, ]
  target <- geo_base_dir(row$layer, row$accession)
  dir.create(file.path(target, "metadata"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(target, "raw"), recursive = TRUE, showWarnings = FALSE)
  if (mode == "manifest") next
  if (mode == "matrix") {
    download_one(matrix_url(row$accession),
                 file.path(target, "metadata", paste0(row$accession, "_series_matrix.txt.gz")))
  }
  if (mode %in% c("supp-manifest", "supp", "supp-parallel")) {
    idx <- tryCatch(readLines(supp_index_url(row$accession), warn = FALSE, encoding = "UTF-8"),
                    error = function(e) character())
    if (!length(idx)) {
      write_log("No supplementary index retrieved: ", row$accession)
      next
    }
    links <- unique(sub('.*href="([^"]+)".*', '\\1', grep('href="[^"]+"', idx, value = TRUE)))
    links <- links[!grepl("^(\\?|/|#|ftp:|http:|https:)", links) & !links %in% c("../")]
    file_urls <- paste0(supp_index_url(row$accession), links)
    get_remote_size <- function(url) {
      tryCatch({
        h <- curl::new_handle(nobody = TRUE, followlocation = TRUE)
        res <- curl::curl_fetch_memory(url, handle = h)
        headers <- rawToChar(res$headers)
        hit <- regmatches(headers, regexpr("(?i)content-length:[[:space:]]*[0-9]+", headers, perl = TRUE))
        if (!length(hit) || !nzchar(hit)) return(NA_real_)
        as.numeric(gsub("[^0-9]", "", hit))
      }, error = function(e) NA_real_)
    }
    sizes <- vapply(file_urls, get_remote_size, numeric(1))
    supp_manifest <- data.frame(accession = row$accession, file = links, url = file_urls,
                                bytes = sizes, gib = sizes / 1024^3, checked_at = timestamp())
    safe_write_csv(supp_manifest,
                   file.path(target, "metadata", paste0(row$accession, "_supplementary_manifest.csv")))
    if (mode == "supp") {
      for (j in seq_along(links)) download_one(file_urls[j], file.path(target, "raw", links[j]))
    }
    if (mode == "supp-parallel") {
      parallel_jobs[[length(parallel_jobs) + 1L]] <- data.table::data.table(
        accession = row$accession, url = file_urls,
        final = file.path(target, "raw", links), expected_bytes = sizes)
    }
  }
}

if (mode == "supp-parallel" && length(parallel_jobs)) {
  jobs <- unique(data.table::rbindlist(parallel_jobs), by = "url")
  jobs[, complete := file.exists(final) & file.info(final)$size > 0 &
         (is.na(expected_bytes) | file.info(final)$size == expected_bytes)]
  jobs <- jobs[complete == FALSE]
  if (nrow(jobs)) {
    jobs[, partial := paste0(final, ".part")]
    for (i in seq_len(nrow(jobs))) {
      legacy <- paste0(jobs$partial[i], ".curltmp")
      if (!file.exists(jobs$partial[i]) && file.exists(legacy)) file.rename(legacy, jobs$partial[i])
    }
    result <- curl::multi_download(jobs$url, jobs$partial, resume = TRUE, progress = TRUE,
                                   multi_timeout = Inf, multiplex = TRUE)
    safe_write_csv(result, file.path(path_data, "00_manifest", "M05_GEO_parallel_download_status.csv"))
    for (i in seq_len(nrow(jobs))) {
      good <- file.exists(jobs$partial[i]) && file.info(jobs$partial[i])$size > 0 &&
        (is.na(jobs$expected_bytes[i]) || file.info(jobs$partial[i])$size == jobs$expected_bytes[i])
      if (good && !file.exists(jobs$final[i])) {
        if (!file.rename(jobs$partial[i], jobs$final[i])) write_log("Rename failed: ", jobs$final[i])
      } else if (!good) {
        write_log("Incomplete download retained: ", jobs$partial[i])
      }
    }
  }
}
write_log("GEO download stage completed in mode=", mode)
