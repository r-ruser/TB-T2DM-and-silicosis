source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table"))

split_geo_line <- function(line) {
  z <- strsplit(line, "\t", fixed = TRUE)[[1]][-1]
  gsub('^"|"$', "", z)
}

matrix_files <- list.files(path_data, pattern = "_series_matrix[.]txt[.]gz$",
                           recursive = TRUE, full.names = TRUE)
if (!length(matrix_files)) stop("No GEO series matrices found; run R/02_download_geo.R --mode=matrix")

summaries <- list()
for (f in matrix_files) {
  acc <- sub("_series_matrix[.]txt[.]gz$", "", basename(f))
  con <- gzfile(f, "rt")
  lines <- readLines(con, warn = FALSE, encoding = "UTF-8")
  close(con)
  get_first <- function(prefix) {
    hit <- grep(paste0("^", prefix, "\t"), lines, value = TRUE)
    if (length(hit)) split_geo_line(hit[1]) else character()
  }
  gsm <- get_first("!Sample_geo_accession")
  title <- get_first("!Sample_title")
  source <- get_first("!Sample_source_name_ch1")
  organism <- get_first("!Sample_organism_ch1")
  n <- max(length(gsm), length(title), length(source), length(organism))
  if (n == 0) n <- 1L
  pad <- function(x) c(x, rep(NA_character_, n - length(x)))
  meta <- data.table::data.table(accession = acc, sample_index = seq_len(n),
                                 sample_id = pad(gsm), title = pad(title),
                                 source_name = pad(source), organism = pad(organism))
  char_lines <- grep("^!Sample_characteristics_ch1\t", lines, value = TRUE)
  if (length(char_lines)) {
    for (k in seq_along(char_lines)) {
      vals <- pad(split_geo_line(char_lines[k]))
      labels <- trimws(sub(":.*$", "", vals))
      values <- trimws(sub("^[^:]+:[[:space:]]*", "", vals))
      field <- if (length(unique(labels[!is.na(labels)])) == 1L) unique(labels[!is.na(labels)]) else paste0("characteristic_", k)
      field <- make.names(field, unique = TRUE)
      if (field %in% names(meta)) field <- paste0(field, "_", k)
      meta[[field]] <- values
    }
  }
  target_layer <- geo_catalog$layer[match(acc, geo_catalog$accession)]
  target <- file.path(path_result, if (target_layer == "bulk") "02_bulk" else "04_scRNA", acc)
  dir.create(target, recursive = TRUE, showWarnings = FALSE)
  safe_write_csv(meta, file.path(target, paste0("A_", acc, "_sample_metadata_raw.csv")))
  table_begin <- match("!series_matrix_table_begin", lines)
  table_end <- match("!series_matrix_table_end", lines)
  n_features <- if (is.na(table_begin) || is.na(table_end)) NA_integer_ else max(0L, table_end - table_begin - 2L)
  summaries[[length(summaries) + 1L]] <- data.table::data.table(
    accession = acc, layer = target_layer, bytes = file.info(f)$size,
    n_samples = sum(!is.na(meta$sample_id)), n_metadata_fields = ncol(meta) - 3L,
    has_expression_table = isTRUE(n_features > 0L), n_expression_features = n_features,
    title = sub('^!Series_title\t"?|"$', "", grep("^!Series_title\t", lines, value = TRUE)[1]),
    md5 = unname(tools::md5sum(f))
  )
}
summary_table <- data.table::rbindlist(summaries, fill = TRUE)
summary_table <- merge(data.table::as.data.table(geo_catalog), summary_table,
                       by = c("accession", "layer"), all.x = TRUE, sort = FALSE)
safe_write_csv(summary_table, file.path(path_result, "00_audit", "A06_GEO_matrix_audit.csv"))

supp_files <- list.files(path_data, pattern = "_supplementary_manifest[.]csv$",
                         recursive = TRUE, full.names = TRUE)
if (length(supp_files)) {
  supp <- data.table::rbindlist(lapply(supp_files, data.table::fread), fill = TRUE)
  safe_write_csv(supp, file.path(path_data, "00_manifest", "M04_GEO_supplementary_files.csv"))
  size_summary <- supp[, .(n_files = .N, known_size_files = sum(!is.na(bytes)),
                           total_gib_known = sum(bytes, na.rm = TRUE) / 1024^3,
                           largest_gib = max(bytes, na.rm = TRUE) / 1024^3), by = accession]
  safe_write_csv(size_summary, file.path(path_result, "00_audit", "A07_GEO_supplementary_size_summary.csv"))
}
write_log("GEO metadata audit completed: ", nrow(summary_table), " series")
