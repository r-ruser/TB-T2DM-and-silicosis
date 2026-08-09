source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "curl"))

sources <- data.frame(
  id = c("WHO01", "WHO02", "WHO03", "WHO04"),
  title = c("TB systematic screening Q&A", "WHO TB Module 2 screening",
            "WHO TB Module 1 preventive treatment", "WHO TB Module 6 comorbidities"),
  url = c(
    "https://www.who.int/news-room/questions-and-answers/item/systematic-screening-for-tb",
    "https://www.who.int/publications-detail-redirect/9789240022676",
    "https://www.who.int/publications/i/item/9789240096196",
    "https://www.who.int/publications/i/item/9789240103276"
  ), stringsAsFactors = FALSE
)

for (i in seq_len(nrow(sources))) {
  dest <- file.path(path_data, "04_WHO", paste0(sources$id[i], ".html"))
  if (!file.exists(dest) || file.info(dest)$size == 0) {
    tmp <- paste0(dest, ".part")
    ok <- tryCatch({
      curl::curl_download(sources$url[i], tmp, mode = "wb", quiet = FALSE)
      file.rename(tmp, dest)
    }, error = function(e) {
      write_log("WHO download failed: ", sources$id[i], " :: ", conditionMessage(e)); FALSE
    })
    sources$status[i] <- if (isTRUE(ok)) "downloaded" else "failed"
  } else sources$status[i] <- "existing"
}
sources$retrieved_at <- timestamp()
safe_write_csv(sources, file.path(path_data, "00_manifest", "M03_WHO_sources.csv"))
write_log("WHO source download stage completed")

