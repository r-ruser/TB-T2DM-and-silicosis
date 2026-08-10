source(file.path("R", "00_config.R"), encoding = "UTF-8")

cat("=== Running All New Analyses ===\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ============================================================
# Analysis execution order
# ============================================================
analyses <- list(
  list(
    name = "SCENIC Heatmap",
    script = "220_SCENIC_heatmap.R",
    description = "SCENIC regulon activity heatmap visualization"
  ),
  list(
    name = "CellChat",
    script = "221_CellChat_analysis.R",
    description = "Cell-cell communication analysis"
  ),
  list(
    name = "NicheNet",
    script = "222_NicheNet_analysis.R",
    description = "Ligand-receptor analysis"
  ),
  list(
    name = "Monocle3 Pseudotime",
    script = "223_Monocle3_pseudotime.R",
    description = "Pseudotime trajectory analysis"
  ),
  list(
    name = "Circlize UMAP",
    script = "224_circlize_UMAP.R",
    description = "Circular UMAP visualization"
  )
)

# ============================================================
# Execute analyses
# ============================================================
results <- list()
for (i in seq_along(analyses)) {
  analysis <- analyses[[i]]
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("Analysis", i, "/", length(analyses), ":", analysis$name, "\n")
  cat(analysis$description, "\n")
  cat(paste(rep("=", 60), collapse = ""), "\n\n")

  script_path <- file.path(path_r, analysis$script)

  if (!file.exists(script_path)) {
    cat("WARNING: Script not found:", script_path, "\n")
    results[[analysis$name]] <- list(status = "MISSING", error = "Script not found")
    next
  }

  tryCatch({
    start_time <- Sys.time()

    # Source the script
    source(script_path, encoding = "UTF-8")

    end_time <- Sys.time()
    elapsed <- difftime(end_time, start_time, units = "mins")

    cat("\n", analysis$name, "completed in", round(elapsed, 2), "minutes\n")
    results[[analysis$name]] <- list(status = "SUCCESS", time = elapsed)

  }, error = function(e) {
    cat("\nERROR in", analysis$name, ":\n")
    cat(e$message, "\n")
    results[[analysis$name]] <- list(status = "ERROR", error = e$message)
  })
}

# ============================================================
# Summary
# ============================================================
cat("\n\n", paste(rep("=", 60), collapse = ""), "\n")
cat("ANALYSIS SUMMARY\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

for (name in names(results)) {
  r <- results[[name]]
  cat(name, ": ", r$status, sep = "")
  if (!is.null(r$time)) cat(" (", round(r$time, 2), " min)", sep = "")
  if (!is.null(r$error)) cat(" - ", r$error, sep = "")
  cat("\n")
}

# Count successes
n_success <- sum(sapply(results, function(r) r$status == "SUCCESS"))
n_total <- length(results)

cat("\n", n_success, "/", n_total, " analyses completed successfully\n")
cat("\nEnd time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

# Save results summary
results_dt <- rbindlist(lapply(names(results), function(name) {
  data.table(
    analysis = name,
    status = results[[name]]$status,
    time_min = ifelse(!is.null(results[[name]]$time),
      round(results[[name]]$time, 2), NA_real_),
    error = ifelse(!is.null(results[[name]]$error),
      results[[name]]$error, "")
  )
}))

fwrite(results_dt, file.path(path_result, "00_audit", "analysis_execution_summary.csv"))

cat("\nResults summary saved to:", file.path(path_result, "00_audit", "analysis_execution_summary.csv"), "\n")
cat("\n=== All analyses completed ===\n")
