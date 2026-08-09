options(stringsAsFactors = FALSE, warn = 1)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
# Protocol file check - skip if not present
if (!file.exists(file.path(project_root, "TB_DM_Silicosis_GBD2023_transcriptomics_research_protocol_V1.md"))) {
  cat("Note: Protocol file not found, continuing anyway\n")
}

path_data <- file.path(project_root, "data")
path_r <- file.path(project_root, "R")
path_result <- file.path(project_root, "result")

project_dirs <- c(
  file.path(path_data, "00_manifest"),
  file.path(path_data, "01_GBD2023", "processed"),
  file.path(path_data, "02_GEO_bulk"),
  file.path(path_data, "03_scRNA"),
  file.path(path_data, "04_WHO"),
  file.path(path_data, "05_reference", "R_library"),
  file.path(path_result, "00_audit"),
  file.path(path_result, "01_GBD2023", "tables"),
  file.path(path_result, "01_GBD2023", "models"),
  file.path(path_result, "01_GBD2023", "figures"),
  file.path(path_result, "01_GBD2023", "source_data"),
  file.path(path_result, "02_bulk"),
  file.path(path_result, "03_cross_disease"),
  file.path(path_result, "04_scRNA"),
  file.path(path_result, "05_GRN_KO"),
  file.path(path_result, "06_final", "tables"),
  file.path(path_result, "06_final", "figures"),
  file.path(path_result, "99_logs")
)
invisible(lapply(project_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

local_library <- file.path(path_data, "05_reference", "R_library")
.libPaths(unique(c(local_library, .libPaths())))

set.seed(20260807)

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")
write_log <- function(..., file = file.path(path_result, "99_logs", "pipeline.log")) {
  line <- paste0("[", timestamp(), "] ", paste(..., collapse = ""))
  cat(line, "\n")
  cat(line, "\n", file = file, append = TRUE)
}

assert_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
  invisible(TRUE)
}

safe_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, path, bom = TRUE, na = "")
}

assert_english_only_svg <- function(path) {
  if (!file.exists(path)) stop("SVG not found: ", path)
  x <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  cp <- utf8ToInt(enc2utf8(x))
  n_cjk <- sum((cp >= 0x3400 & cp <= 0x9FFF) |
                 (cp >= 0xF900 & cp <= 0xFAFF))
  if (n_cjk > 0L) stop("English-only figure QA failed: ", n_cjk, " CJK characters in ", path)
  invisible(TRUE)
}

geo_catalog <- data.frame(
  accession = c("GSE114192", "GSE165489", "GSE181143", "GSE193978", "GSE193979",
                "GSE249102", "GSE283452", "GSE264182", "GSE174725", "GSE326212",
                "GSE192483", "GSE268210"),
  layer = c(rep("bulk", 8), rep("scRNA", 4)),
  role = c("TB-DM discovery", "silicosis discovery", "TB-DM external validation",
           "TB-DM longitudinal validation", "TB outcome validation", "glycaemic gradient",
           "T2DM-Mtb macrophage bridge", "silicosis lung validation",
           "silicosis BALF localization", "TB airway localization",
           "TB lung lesion validation", "T2DM PBMC supplement"),
  priority = c("primary", "primary", "validation", "validation", "validation", "exploratory",
               "mechanistic", "validation", "primary", "primary", "validation", "supplementary"),
  stringsAsFactors = FALSE
)
