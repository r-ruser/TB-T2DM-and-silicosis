source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "jsonlite", "digest"))

all_pkgs <- c("data.table", "dplyr", "tidyr", "readr", "ggplot2", "patchwork", "sf",
              "spdep", "rnaturalearth", "rnaturalearthdata", "GEOquery", "Biobase", "limma",
              "edgeR", "DESeq2", "fgsea", "WGCNA", "ComplexHeatmap", "svglite", "ragg")
pkg_table <- data.frame(
  package = all_pkgs,
  installed = vapply(all_pkgs, requireNamespace, logical(1), quietly = TRUE),
  version = vapply(all_pkgs, function(x) if (requireNamespace(x, quietly = TRUE))
    as.character(utils::packageVersion(x)) else NA_character_, character(1))
)
safe_write_csv(pkg_table, file.path(path_result, "00_audit", "A01_R_package_audit.csv"))

files <- list.files(path_data, recursive = TRUE, full.names = TRUE, all.files = FALSE)
file_info <- if (length(files)) {
  info <- file.info(files)
  data.frame(path = substring(normalizePath(files, winslash = "/"), nchar(project_root) + 2L),
             bytes = info$size, modified = format(info$mtime, "%Y-%m-%d %H:%M:%S"),
             md5 = unname(tools::md5sum(files)), stringsAsFactors = FALSE)
} else data.frame(path = character(), bytes = numeric(), modified = character(), md5 = character())
safe_write_csv(file_info, file.path(path_data, "00_manifest", "M01_file_inventory.csv"))

capture.output(sessionInfo(), file = file.path(path_result, "99_logs", "sessionInfo.txt"))
write_log("Environment audit completed: ", nrow(file_info), " data files inventoried")

