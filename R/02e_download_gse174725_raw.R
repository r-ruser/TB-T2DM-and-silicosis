source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table"))

m04 <- data.table::fread(file.path(path_data, "00_manifest", "M04_GEO_supplementary_files.csv"))
row <- m04[accession == "GSE174725" & file == "GSE174725_RAW.tar"]
if (nrow(row) != 1L) stop("GSE174725 RAW manifest row not found")
final <- file.path(path_data, "03_scRNA", "GSE174725", "raw", row$file)
stage <- paste0(final, ".aria")
if (!file.exists(final) || file.info(final)$size != row$bytes) {
  aria2 <- Sys.which("aria2c")
  args <- c("--allow-overwrite=true", "--auto-file-renaming=false", "--continue=true",
            "--file-allocation=none", "--max-connection-per-server=1", "--split=1",
            "--max-tries=0", "--retry-wait=10", "--timeout=60", "--connect-timeout=30",
            "--summary-interval=0", "--console-log-level=warn", "--download-result=hide",
            paste0("--dir=", dirname(stage)), paste0("--out=", basename(stage)), row$url)
  status <- system2(aria2, args, stdout = FALSE, stderr = FALSE)
  if (!identical(status, 0L)) stop("aria2 failed with status ", status)
  if (!file.exists(stage) || file.info(stage)$size != row$bytes) stop("GSE174725 RAW size mismatch")
  if (!file.exists(final) && !file.rename(stage, final)) stop("Could not finalize GSE174725 RAW")
}
if (file.info(final)$size != row$bytes) stop("Final GSE174725 RAW verification failed")
write_log("GSE174725 donor-resolved RAW archive downloaded and byte-verified")

