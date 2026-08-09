source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "Matrix", "S4Vectors", "SingleCellExperiment", "SummarizedExperiment", "scuttle",
  "ggplot2", "patchwork", "svglite", "ragg"))

acc <- "GSE326212"
result_dir <- file.path(path_result, "04_scRNA", acc)
processed_dir <- file.path(path_data, "03_scRNA", acc, "processed", "selected_samples")
selection <- data.table::fread(file.path(result_dir, "tables", "T05_GSE326212_primary_sample_selection.csv"))
paths <- file.path(processed_dir, paste0(selection$sample_id, "_sce_qc.rds"))
if (!all(file.exists(paths))) stop("Run R/49_gse326212_selected_sample_qc.R to completion first")

write_log("Assembling GSE326212 QC-passed samples sequentially")
# Do not hold all 33 SCE objects in memory: their combined serialized objects alone
# exceed 30 GB.  Derive the shared Ensembl universe one file at a time, then
# append sparse count blocks one file at a time.
ensembl_key <- function(x) {
  feature <- as.character(SummarizedExperiment::rowData(x)[["gene_id"]])
  key <- regmatches(feature, regexpr("ENSG[0-9]+", feature))
  key[grepl("^ENSG[0-9]+$", key) & !duplicated(key)]
}
common <- NULL
for (pth in paths) {
  x <- readRDS(pth)
  keys <- ensembl_key(x)
  common <- if (is.null(common)) keys else base::intersect(common, keys)
  rm(x, keys); gc(verbose = FALSE)
}
if (length(common) < 10000L) stop("Unexpectedly small shared Ensembl gene universe: ", length(common))
sample_counts <- NULL
sample_meta <- NULL
for (pth in paths) {
  x <- readRDS(pth)
  keys <- ensembl_key(x)
  idx <- match(common, keys)
  keep <- grepl("^ENSG[0-9]+$", regmatches(as.character(SummarizedExperiment::rowData(x)[["gene_id"]]),
    regexpr("ENSG[0-9]+", as.character(SummarizedExperiment::rowData(x)[["gene_id"]]))))
  keep <- which(keep)[!duplicated(keys)]
  block <- SummarizedExperiment::assay(x, "counts")[keep[idx], , drop = FALSE]
  # One selected sample represents one participant.  Sum its QC-passed cells
  # before combining, so this creates a valid donor-level matrix without an
  # infeasible 300k-cell in-memory concatenation.
  sid <- as.character(SummarizedExperiment::colData(x)$sample_id[1])
  sample_counts <- cbind(sample_counts, Matrix::Matrix(Matrix::rowSums(block), sparse = TRUE))
  colnames(sample_counts)[ncol(sample_counts)] <- sid
  sample_meta <- data.table::rbindlist(list(sample_meta, data.table::as.data.table(as.data.frame(
    SummarizedExperiment::colData(x)[1, , drop = FALSE]))), fill = TRUE)
  rm(x, keys, idx, keep, block); gc(verbose = FALSE)
}
rownames(sample_counts) <- common
saveRDS(list(counts = sample_counts, sample_metadata = sample_meta,
  gene_id = common, unit = "QC-passed cells summed within each selected sample"),
  file.path(path_data, "03_scRNA", acc, "processed", paste0(acc, "_selected_sample_pseudobulk.rds")), compress = FALSE)

summary_dt <- data.table::fread(file.path(result_dir, "tables", "T07_GSE326212_selected_QC_progress.csv"))
summary_dt <- merge(selection[, .(sample_id, subject_id, analysis_group, status)], summary_dt,
  by = c("sample_id", "subject_id", "analysis_group", "status"), all.x = TRUE, sort = FALSE)
summary_dt[, `:=`(median_counts = NA_real_, median_detected = NA_real_, median_mito_percent = NA_real_)]
safe_write_csv(summary_dt, file.path(result_dir, "tables", "T09_GSE326212_QC_summary.csv"))

q_long <- summary_dt[, .(sample_id, subject_id, analysis_group, status,
  metric = "Final singlets", value = n_final_singlets)]
pal <- c("Active TB" = "#B84A4A", "Stable controller" = "#416A9A", "Pre-progression" = "#D08A3E")
p <- ggplot2::ggplot(q_long, ggplot2::aes(analysis_group, value, fill = analysis_group)) +
  ggplot2::geom_boxplot(width = 0.62, outlier.size = 0.55, linewidth = 0.25) +
  ggplot2::geom_jitter(width = 0.11, size = 0.6, alpha = 0.75) +
  ggplot2::facet_wrap(~metric, scales = "free_y", nrow = 1) +
  ggplot2::scale_fill_manual(values = pal) +
  ggplot2::labs(title = "GSE326212 selected BAL samples: quality-control summary",
    subtitle = "One sample per subject; sample-level distributions shown for audit only",
    x = NULL, y = NULL, fill = NULL) +
  ggplot2::theme_classic(base_size = 7, base_family = "Arial") +
  ggplot2::theme(axis.line = ggplot2::element_line(linewidth = 0.3),
    axis.ticks = ggplot2::element_line(linewidth = 0.3), axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
    strip.text = ggplot2::element_text(size = 6.5, face = "bold"), plot.title = ggplot2::element_text(size = 8, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 6, colour = "#767676"), legend.position = "none")
out <- file.path(result_dir, "figures", "F18_GSE326212_QC_summary")
w <- 183 / 25.4; h <- 72 / 25.4
svglite::svglite(paste0(out, ".svg"), width = w, height = h); print(p); grDevices::dev.off()
assert_english_only_svg(paste0(out, ".svg"))
grDevices::cairo_pdf(paste0(out, ".pdf"), width = w, height = h, family = "Arial"); print(p); grDevices::dev.off()
ragg::agg_tiff(paste0(out, ".tiff"), width = w, height = h, units = "in", res = 600, compression = "lzw"); print(p); grDevices::dev.off()
ragg::agg_png(paste0(out, ".png"), width = w, height = h, units = "in", res = 300); print(p); grDevices::dev.off()
write_log("GSE326212 donor-level pseudobulk completed: samples=", nrow(summary_dt), "; final singlets=", sum(summary_dt$n_final_singlets),
  "; common genes=", length(common), "; English-only SVG QA passed")
