source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "rhdf5", "Matrix", "SingleCellExperiment", "SummarizedExperiment",
  "S4Vectors", "scDblFinder", "BiocParallel"))

acc <- "GSE326212"
result_dir <- file.path(path_result, "04_scRNA", acc)
processed_dir <- file.path(path_data, "03_scRNA", acc, "processed", "selected_samples")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
selection <- data.table::fread(file.path(result_dir, "tables", "T05_GSE326212_primary_sample_selection.csv"))
if (nrow(selection) != 33L) stop("Run R/48_gse326212_select_extract.R first")

read_10x_h5 <- function(path, sample_id, subject_id, analysis_group, status) {
  shape <- as.integer(rhdf5::h5read(path, "/matrix/shape"))
  p <- as.integer(rhdf5::h5read(path, "/matrix/indptr"))
  i <- as.integer(rhdf5::h5read(path, "/matrix/indices"))
  x <- as.numeric(rhdf5::h5read(path, "/matrix/data"))
  if (length(shape) != 2L || length(p) != shape[2] + 1L) stop("Invalid sparse H5 structure: ", path)
  m <- methods::new("dgCMatrix", Dim = shape, p = p, i = i, x = x)
  gene_id <- as.character(rhdf5::h5read(path, "/matrix/features/id"))
  gene_symbol <- as.character(rhdf5::h5read(path, "/matrix/features/name"))
  barcodes <- as.character(rhdf5::h5read(path, "/matrix/barcodes"))
  rownames(m) <- gene_id
  colnames(m) <- paste(sample_id, barcodes, sep = "__")
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = m),
    rowData = S4Vectors::DataFrame(gene_id = gene_id, gene_symbol = gene_symbol),
    colData = S4Vectors::DataFrame(cell_id = colnames(m), sample_id = sample_id,
      subject_id = as.character(subject_id), analysis_group = analysis_group, status = status))
  sce
}

summary_list <- list(); threshold_list <- list()
for (k in seq_len(nrow(selection))) {
  s <- selection[k]
  out_file <- file.path(processed_dir, paste0(s$sample_id, "_sce_qc.rds"))
  if (file.exists(out_file)) {
    write_log("GSE326212 QC checkpoint exists: ", s$sample_id)
    q <- readRDS(out_file)
    summary_list[[s$sample_id]] <- data.table::data.table(sample_id = s$sample_id,
      subject_id = s$subject_id, analysis_group = s$analysis_group, status = s$status,
      n_raw = q$qc_n_raw[1], n_after_adaptive_qc = q$qc_n_adaptive[1],
      n_doublets = q$qc_n_adaptive[1] - ncol(q), n_final_singlets = ncol(q),
      final_retention_percent = 100 * ncol(q) / q$qc_n_raw[1])
    next
  }
  write_log("GSE326212 QC reading sample ", k, "/", nrow(selection), ": ", s$sample_id)
  sce <- read_10x_h5(s$local_h5, s$sample_id, s$subject_id, s$analysis_group, s$status)
  n_raw <- ncol(sce)
  lib <- Matrix::colSums(SummarizedExperiment::assay(sce, "counts"))
  detected <- Matrix::colSums(SummarizedExperiment::assay(sce, "counts") > 0)
  mito <- grepl("^MT-", SummarizedExperiment::rowData(sce)$gene_symbol, ignore.case = TRUE)
  if (!any(mito)) stop("No mitochondrial genes found: ", s$sample_id)
  mito_percent <- 100 * Matrix::colSums(SummarizedExperiment::assay(sce, "counts")[mito, , drop = FALSE]) /
    pmax(lib, 1)
  lc <- log10(lib + 1); lf <- log10(detected + 1)
  threshold <- data.table::data.table(sample_id = s$sample_id,
    min_log10_counts = max(log10(300), stats::median(lc) - 3 * stats::mad(lc)),
    max_log10_counts = stats::median(lc) + 4 * stats::mad(lc),
    min_log10_features = max(log10(200), stats::median(lf) - 3 * stats::mad(lf)),
    max_log10_features = stats::median(lf) + 4 * stats::mad(lf),
    max_mito_percent = min(30, max(10, stats::median(mito_percent) + 3 * stats::mad(mito_percent))))
  keep <- lc >= threshold$min_log10_counts & lc <= threshold$max_log10_counts &
    lf >= threshold$min_log10_features & lf <= threshold$max_log10_features &
    mito_percent <= threshold$max_mito_percent
  keep[is.na(keep)] <- FALSE
  sce$sum <- lib; sce$detected <- detected; sce$subsets_Mito_percent <- mito_percent
  sce <- sce[, keep]
  n_adaptive <- ncol(sce)
  if (n_adaptive < 500L) stop("Too few cells after adaptive QC: ", s$sample_id)
  set.seed(20260808 + k)
  sce <- scDblFinder::scDblFinder(sce, samples = "sample_id", clusters = FALSE,
    BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
  sce <- sce[, sce$scDblFinder.class == "singlet"]
  sce$qc_n_raw <- n_raw
  sce$qc_n_adaptive <- n_adaptive
  # Save counts and compact QC metadata only; no dense normalized assay is persisted here.
  keep_coldata <- c("cell_id", "sample_id", "subject_id", "analysis_group", "status", "sum",
    "detected", "subsets_Mito_percent", "scDblFinder.score", "scDblFinder.class",
    "qc_n_raw", "qc_n_adaptive")
  SummarizedExperiment::colData(sce) <- SummarizedExperiment::colData(sce)[, keep_coldata, drop = FALSE]
  saveRDS(sce, out_file, compress = FALSE)
  threshold_list[[s$sample_id]] <- threshold
  summary_list[[s$sample_id]] <- data.table::data.table(sample_id = s$sample_id,
    subject_id = s$subject_id, analysis_group = s$analysis_group, status = s$status,
    n_raw = n_raw, n_after_adaptive_qc = n_adaptive, n_doublets = n_adaptive - ncol(sce),
    n_final_singlets = ncol(sce), final_retention_percent = 100 * ncol(sce) / n_raw)
  safe_write_csv(data.table::rbindlist(summary_list, fill = TRUE),
    file.path(result_dir, "tables", "T07_GSE326212_selected_QC_progress.csv"))
  if (length(threshold_list)) safe_write_csv(data.table::rbindlist(threshold_list, fill = TRUE),
    file.path(result_dir, "tables", "T08_GSE326212_selected_QC_thresholds.csv"))
  write_log("GSE326212 QC saved ", s$sample_id, ": raw=", n_raw, "; adaptive=", n_adaptive,
    "; singlets=", ncol(sce))
  rm(sce); gc(verbose = FALSE)
}
summary_dt <- data.table::rbindlist(summary_list, fill = TRUE)
safe_write_csv(summary_dt, file.path(result_dir, "tables", "T07_GSE326212_selected_QC_progress.csv"))
write_log("GSE326212 selected-sample QC completed: samples=", nrow(summary_dt),
  "; raw cells=", sum(summary_dt$n_raw), "; final singlets=", sum(summary_dt$n_final_singlets))
