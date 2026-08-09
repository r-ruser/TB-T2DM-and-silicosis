source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table"))

aria2 <- Sys.which("aria2c")
if (!nzchar(aria2)) stop("aria2c is not installed or not on PATH")
manifest_file <- file.path(path_data, "00_manifest", "M06_GEO_selected_downloads.csv")
if (!file.exists(manifest_file)) stop("Selected GEO manifest missing")
m <- data.table::fread(manifest_file, encoding = "UTF-8")

prepare_stage <- function(final, expected) {
  stage <- paste0(final, ".aria")
  if (file.exists(stage)) return(stage)
  candidates <- c(paste0(final, ".download.curltmp"), paste0(final, ".part.curltmp"),
                  paste0(final, ".download"), paste0(final, ".part"),
                  paste0(final, ".aria_test"), paste0(final, ".ftp_test.curltmp"))
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates)) {
    sizes <- file.info(candidates)$size
    candidates <- candidates[sizes > 0 & sizes <= expected]
    if (length(candidates)) {
      source <- candidates[which.max(file.info(candidates)$size)]
      if (!file.rename(source, stage)) stop("Could not stage partial file: ", source)
    }
  }
  stage
}

m[, verified_before := file.exists(final) & file.info(final)$size == bytes]
jobs <- m[verified_before == FALSE]
if (!nrow(jobs)) {
  write_log("All selected GEO files already verified; aria2 stage skipped")
  quit(save = "no", status = 0)
}
jobs[, stage := mapply(prepare_stage, final, bytes, USE.NAMES = FALSE)]
safe_write_csv(jobs, file.path(path_data, "00_manifest", "M09_GEO_aria2_jobs.csv"))

worker <- function(i, aria2, urls, stages, expected) {
  stage <- stages[i]
  dir.create(dirname(stage), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(stage) || file.info(stage)$size != expected[i]) {
    args <- c("--allow-overwrite=true", "--auto-file-renaming=false", "--continue=true",
              "--file-allocation=none", "--max-connection-per-server=1", "--split=1",
              "--max-tries=0", "--retry-wait=10", "--timeout=60", "--connect-timeout=30",
              "--summary-interval=0", "--console-log-level=warn", "--download-result=hide",
              paste0("--dir=", dirname(stage)), paste0("--out=", basename(stage)), urls[i])
    status <- tryCatch(system2(aria2, args, stdout = FALSE, stderr = FALSE),
                       error = function(e) 99L)
  } else status <- 0L
  actual <- if (file.exists(stage)) file.info(stage)$size else NA_real_
  data.frame(i = i, status = status, stage_bytes = actual,
             stage_verified = !is.na(actual) && actual == expected[i])
}

cl <- parallel::makePSOCKcluster(min(4L, nrow(jobs)))
res_list <- parallel::parLapplyLB(cl, seq_len(nrow(jobs)), worker, aria2 = aria2,
                                 urls = jobs$url, stages = jobs$stage, expected = jobs$bytes)
parallel::stopCluster(cl)
res <- cbind(jobs[, .(accession, file, url, final, stage, expected_bytes = bytes)],
             data.table::rbindlist(res_list))
for (i in seq_len(nrow(res))) {
  if (res$stage_verified[i] && !file.exists(res$final[i])) {
    if (!file.rename(res$stage[i], res$final[i])) res$status[i] <- 98L
  }
  res$final_bytes[i] <- if (file.exists(res$final[i])) file.info(res$final[i])$size else NA_real_
  res$verified[i] <- !is.na(res$final_bytes[i]) && res$final_bytes[i] == res$expected_bytes[i]
}
safe_write_csv(res, file.path(path_data, "00_manifest", "M10_GEO_aria2_download_status.csv"))
write_log("aria2 local GEO download finished: verified ", sum(res$verified), "/", nrow(res))
if (!all(res$verified)) quit(save = "no", status = 2)

