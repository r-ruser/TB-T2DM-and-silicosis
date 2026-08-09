source(file.path("R", "00_config.R"), encoding = "UTF-8")

cran_packages <- c("data.table", "dplyr", "tidyr", "readr", "ggplot2", "patchwork",
                   "sf", "spdep", "rnaturalearth", "rnaturalearthdata", "svglite",
                   "ragg", "jsonlite", "yaml", "digest", "curl", "cluster", "broom")
bioc_packages <- c("GEOquery", "Biobase", "limma", "edgeR", "DESeq2", "fgsea",
                   "WGCNA", "ComplexHeatmap")

repos <- c(CRAN = "https://cloud.r-project.org")
missing_cran <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran)) {
  install.packages(missing_cran, lib = local_library, repos = repos, dependencies = TRUE)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = local_library, repos = repos)
}
missing_bioc <- bioc_packages[!vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) {
  BiocManager::install(missing_bioc, lib = local_library, ask = FALSE, update = FALSE)
}

write_log("Bootstrap completed; local R library: ", local_library)

