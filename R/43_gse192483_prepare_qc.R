source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "Matrix", "SingleCellExperiment", "SummarizedExperiment", "S4Vectors",
  "scuttle", "scDblFinder", "BiocParallel", "ggplot2", "patchwork", "svglite", "ragg"))

acc <- "GSE192483"
raw_dir <- file.path(path_data, "03_scRNA", acc, "raw")
processed_dir <- file.path(path_data, "03_scRNA", acc, "processed")
unpack_dir <- file.path(raw_dir, "sample_matrices")
result_dir <- file.path(path_result, "04_scRNA", acc)
invisible(lapply(c(processed_dir, unpack_dir,
  file.path(result_dir, c("tables", "source_data", "models", "figures"))),
  dir.create, recursive = TRUE, showWarnings = FALSE))

tar_file <- file.path(raw_dir, paste0(acc, "_RAW.tar"))
if (!file.exists(tar_file)) stop("RAW archive missing: ", tar_file)
members <- utils::untar(tar_file, list = TRUE)
matrix_members <- members[grepl("[.](matrix[.]mtx|barcodes[.]tsv|features[.]tsv)[.]gz$", members)]
if (length(matrix_members) != 33L) stop("Expected 11 sample triplets, found ", length(matrix_members), " files")
missing_members <- matrix_members[!file.exists(file.path(unpack_dir, matrix_members))]
if (length(missing_members)) utils::untar(tar_file, files = missing_members, exdir = unpack_dir)

feature_members <- matrix_members[grepl("features[.]tsv[.]gz$", matrix_members)]
prefix <- sub("[.]features[.]tsv[.]gz$", "", feature_members)
sample_code <- sub("^GSM[0-9]+_", "", prefix)
patient <- sub("[HL]$", "", sample_code)
region_code <- sub("^.*([HL])$", "\\1", sample_code)
sample_map <- data.table::data.table(
  sample_id = prefix, sample_code = sample_code, patient = patient, region_code = region_code,
  region = ifelse(region_code == "H", "TB lesion", "Less-involved lung")
)
data.table::setorder(sample_map, patient, region_code)
paired_patients <- sample_map[, .N, by = patient][N == 2L, patient]
sample_map[, paired_patient := patient %in% paired_patients]
if (data.table::uniqueN(sample_map$patient) != 6L || nrow(sample_map) != 11L) stop("Unexpected patient/sample structure")
safe_write_csv(sample_map, file.path(result_dir, "source_data", "SD01_GSE192483_sample_manifest.csv"))

read_sample <- function(i) {
  sid <- sample_map$sample_id[i]
  f_feature <- file.path(unpack_dir, paste0(sid, ".features.tsv.gz"))
  f_barcode <- file.path(unpack_dir, paste0(sid, ".barcodes.tsv.gz"))
  f_matrix <- file.path(unpack_dir, paste0(sid, ".matrix.mtx.gz"))
  write_log("Reading GSE192483 sample: ", sid)
  features <- data.table::fread(f_feature, header = FALSE)
  barcodes <- data.table::fread(f_barcode, header = FALSE)[[1]]
  con <- gzfile(f_matrix, open = "rt")
  on.exit(close(con), add = TRUE)
  m <- Matrix::readMM(con)
  if (nrow(m) != nrow(features) || ncol(m) != length(barcodes)) stop("10x dimension mismatch: ", sid)
  gene_id <- as.character(features[[1]])
  gene_symbol <- as.character(features[[2]])
  if (anyDuplicated(gene_id)) stop("Duplicated feature identifiers: ", sid)
  rownames(m) <- gene_id
  cell_id <- paste(sid, barcodes, sep = "__")
  colnames(m) <- cell_id
  lib <- Matrix::colSums(m)
  detected <- Matrix::colSums(m > 0)
  mito <- grepl("^MT-", gene_symbol, ignore.case = TRUE)
  if (!any(mito)) stop("No mitochondrial genes found: ", sid)
  mito_percent <- 100 * Matrix::colSums(m[mito, , drop = FALSE]) / pmax(lib, 1)
  list(counts = methods::as(m, "dgCMatrix"), gene_id = gene_id, gene_symbol = gene_symbol,
    qc = data.table::data.table(cell_id = cell_id, sum = lib, detected = detected,
      subsets_Mito_percent = mito_percent))
}

sample_objects <- lapply(seq_len(nrow(sample_map)), read_sample)
gene_sets <- lapply(sample_objects, `[[`, "gene_id")
common_gene_id <- Reduce(intersect, gene_sets)
if (length(common_gene_id) < 15000L) stop("Unexpectedly small common feature universe")
reference_symbols <- sample_objects[[1]]$gene_symbol[match(common_gene_id, sample_objects[[1]]$gene_id)]
gene_audit <- data.table::rbindlist(lapply(seq_along(sample_objects), function(i) {
  data.table::data.table(sample_id = sample_map$sample_id[i], raw_genes = length(sample_objects[[i]]$gene_id),
    common_genes = length(common_gene_id), raw_cells = ncol(sample_objects[[i]]$counts))
}))
safe_write_csv(gene_audit, file.path(result_dir, "tables", "T01_GSE192483_gene_cell_audit.csv"))
counts <- Reduce(Matrix::cbind2, lapply(sample_objects, function(z) {
  z$counts[match(common_gene_id, z$gene_id), , drop = FALSE]
}))
rownames(counts) <- common_gene_id
native_qc <- data.table::rbindlist(lapply(sample_objects, `[[`, "qc"))
rm(sample_objects); gc(verbose = FALSE)

sid_vec <- sub("__.*$", "", colnames(counts))
sample_idx <- match(sid_vec, sample_map$sample_id)
if (anyNA(sample_idx)) stop("Cell-to-sample mapping failed")
sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts),
  rowData = S4Vectors::DataFrame(gene_id = common_gene_id, gene_symbol = reference_symbols),
  colData = S4Vectors::DataFrame(cell_id = colnames(counts), sample_id = sid_vec,
    sample_code = sample_map$sample_code[sample_idx], patient = sample_map$patient[sample_idx],
    region = sample_map$region[sample_idx], paired_patient = sample_map$paired_patient[sample_idx]))
native_qc <- native_qc[match(colnames(sce), cell_id)]
sce$sum <- native_qc$sum
sce$detected <- native_qc$detected
sce$subsets_Mito_percent <- native_qc$subsets_Mito_percent

thresholds <- data.table::rbindlist(lapply(unique(sce$sample_id), function(sid) {
  j <- which(sce$sample_id == sid)
  lc <- log10(sce$sum[j] + 1); lf <- log10(sce$detected[j] + 1); mt <- sce$subsets_Mito_percent[j]
  data.table::data.table(sample_id = sid,
    min_log10_counts = max(log10(300), stats::median(lc) - 3 * stats::mad(lc)),
    max_log10_counts = stats::median(lc) + 4 * stats::mad(lc),
    min_log10_features = max(log10(200), stats::median(lf) - 3 * stats::mad(lf)),
    max_log10_features = stats::median(lf) + 4 * stats::mad(lf),
    max_mito_percent = min(30, max(10, stats::median(mt) + 3 * stats::mad(mt))))
}))
thr <- thresholds[match(sce$sample_id, sample_id)]
sce$qc_pass_initial <- log10(sce$sum + 1) >= thr$min_log10_counts &
  log10(sce$sum + 1) <= thr$max_log10_counts &
  log10(sce$detected + 1) >= thr$min_log10_features &
  log10(sce$detected + 1) <= thr$max_log10_features &
  sce$subsets_Mito_percent <= thr$max_mito_percent
sce$qc_pass_initial[is.na(sce$qc_pass_initial)] <- FALSE
sce_qc <- sce[, sce$qc_pass_initial]
write_log("GSE192483 adaptive QC retained ", ncol(sce_qc), " of ", ncol(sce), " cells")
set.seed(20260808)
sce_qc <- scDblFinder::scDblFinder(sce_qc, samples = "sample_id", clusters = FALSE,
  BPPARAM = BiocParallel::SerialParam(), verbose = TRUE)
sce_final <- sce_qc[, sce_qc$scDblFinder.class == "singlet"]
sce_final <- scuttle::logNormCounts(sce_final)

safe_write_csv(thresholds, file.path(result_dir, "tables", "T02_GSE192483_QC_thresholds.csv"))
cell_qc <- data.table::as.data.table(as.data.frame(SummarizedExperiment::colData(sce_qc)))
safe_write_csv(cell_qc, file.path(result_dir, "source_data", "SD02_GSE192483_cell_QC.csv"))
qc_summary <- data.table::rbindlist(lapply(sample_map$sample_id, function(sid) {
  n_raw <- sum(sce$sample_id == sid); n_adaptive <- sum(sce_qc$sample_id == sid)
  n_final <- sum(sce_final$sample_id == sid)
  data.table::data.table(sample_id = sid, sample_code = sample_map[sample_id == sid, sample_code],
    patient = sample_map[sample_id == sid, patient], region = sample_map[sample_id == sid, region],
    n_raw = n_raw, n_after_adaptive_qc = n_adaptive, n_doublets = n_adaptive - n_final,
    n_final_singlets = n_final, final_retention_percent = 100 * n_final / n_raw)
}))
safe_write_csv(qc_summary, file.path(result_dir, "tables", "T03_GSE192483_QC_summary.csv"))
saveRDS(sce_final, file.path(processed_dir, paste0(acc, "_sce_qc.rds")), compress = FALSE)

plot_qc <- data.table::as.data.table(as.data.frame(SummarizedExperiment::colData(sce)))
plot_qc[, sample_code := factor(sample_code, levels = sample_map$sample_code)]
region_pal <- c("Less-involved lung" = "#6A9FB5", "TB lesion" = "#B84A4A")
p_a <- ggplot2::ggplot(plot_qc, ggplot2::aes(sample_code, detected, fill = region)) +
  ggplot2::geom_violin(scale = "width", linewidth = 0.2, trim = TRUE) +
  ggplot2::scale_y_log10() + ggplot2::scale_fill_manual(values = region_pal) +
  ggplot2::labs(title = "Detected features before filtering", x = NULL, y = "Detected genes", fill = NULL)
p_b <- ggplot2::ggplot(plot_qc, ggplot2::aes(sample_code, subsets_Mito_percent, fill = region)) +
  ggplot2::geom_violin(scale = "width", linewidth = 0.2, trim = TRUE) +
  ggplot2::scale_fill_manual(values = region_pal) +
  ggplot2::labs(title = "Mitochondrial fraction before filtering", x = NULL,
    y = "Mitochondrial reads (%)", fill = NULL) + ggplot2::theme(legend.position = "none")
q_long <- data.table::melt(qc_summary, id.vars = c("sample_id", "sample_code", "patient", "region"),
  measure.vars = c("n_raw", "n_after_adaptive_qc", "n_final_singlets"),
  variable.name = "stage", value.name = "cells")
q_long[, stage := factor(stage, levels = c("n_raw", "n_after_adaptive_qc", "n_final_singlets"),
  labels = c("Raw", "Adaptive QC", "Final singlets"))]
q_long[, sample_code := factor(sample_code, levels = sample_map$sample_code)]
p_c <- ggplot2::ggplot(q_long, ggplot2::aes(sample_code, cells, fill = stage)) +
  ggplot2::geom_col(position = "dodge", width = 0.75) +
  ggplot2::scale_fill_manual(values = c("Raw" = "#C9C9C9", "Adaptive QC" = "#6AA6A1",
    "Final singlets" = "#416A9A")) +
  ggplot2::labs(title = "Cell retention by sample", x = NULL, y = "Cells", fill = NULL)
theme_sc <- ggplot2::theme_classic(base_size = 7, base_family = "Arial") +
  ggplot2::theme(axis.line = ggplot2::element_line(linewidth = 0.3),
    axis.ticks = ggplot2::element_line(linewidth = 0.3), legend.position = "top",
    axis.text.x = ggplot2::element_text(angle = 38, hjust = 1),
    plot.title = ggplot2::element_text(size = 8, face = "bold"),
    plot.tag = ggplot2::element_text(size = 8, face = "bold"))
fig <- (p_a | p_b | p_c) + patchwork::plot_annotation(tag_levels = "a",
  title = "GSE192483 patient-resolved quality control",
  subtitle = "Adaptive sample-specific thresholds followed by per-sample doublet detection") & theme_sc
out <- file.path(result_dir, "figures", "F14_GSE192483_QC")
w <- 183 / 25.4; h <- 90 / 25.4
svglite::svglite(paste0(out, ".svg"), width = w, height = h); print(fig); grDevices::dev.off()
assert_english_only_svg(paste0(out, ".svg"))
grDevices::cairo_pdf(paste0(out, ".pdf"), width = w, height = h, family = "Arial"); print(fig); grDevices::dev.off()
ragg::agg_tiff(paste0(out, ".tiff"), width = w, height = h, units = "in", res = 600,
  compression = "lzw"); print(fig); grDevices::dev.off()
ragg::agg_png(paste0(out, ".png"), width = w, height = h, units = "in", res = 300)
print(fig); grDevices::dev.off()
write_log("GSE192483 QC completed: raw cells=", ncol(sce), "; final singlets=", ncol(sce_final),
  "; patients=", data.table::uniqueN(sample_map$patient), "; samples=", nrow(sample_map),
  "; English-only SVG QA passed")
