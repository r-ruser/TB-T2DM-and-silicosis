scripts <- c("01_environment_audit.R", "02_download_geo.R", "03_download_who.R",
             "10_gbd_audit_clean.R", "11_gbd_proxy_analysis.R", "12_gbd_proxy_figures.R")
for (s in scripts) {
  message("\n===== Running ", s, " =====")
  status <- system2(file.path(R.home("bin"), "Rscript.exe"), file.path("R", s))
  if (!identical(status, 0L)) stop("Stage failed: ", s, " (exit ", status, ")")
}

