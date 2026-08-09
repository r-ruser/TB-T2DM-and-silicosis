source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "Matrix", "SingleCellExperiment", "SummarizedExperiment", "S4Vectors",
                  "scuttle", "scDblFinder", "BiocParallel", "ggplot2",
                  "patchwork", "svglite", "ragg"))

acc <- "GSE174725"
raw_dir <- file.path(path_data, "03_scRNA", acc, "raw")
processed_dir <- file.path(path_data, "03_scRNA", acc, "processed")
result_dir <- file.path(path_result, "04_scRNA", acc)
invisible(lapply(c(processed_dir, file.path(raw_dir, "donor_matrices"),
                   file.path(result_dir, c("tables", "source_data", "models", "figures"))),
                 dir.create, recursive = TRUE, showWarnings = FALSE))

tar_file <- file.path(raw_dir, paste0(acc, "_RAW.tar"))
unpack_dir <- file.path(raw_dir, "donor_matrices")
if (!file.exists(tar_file)) stop("Donor-level RAW archive missing: ", tar_file)
expected <- utils::untar(tar_file, list = TRUE)
if (length(expected) != 5L) stop("Expected five donor matrices, found ", length(expected))
missing_unpacked <- expected[!file.exists(file.path(unpack_dir, expected))]
if (length(missing_unpacked)) utils::untar(tar_file, files = missing_unpacked, exdir = unpack_dir)
files <- file.path(unpack_dir, expected)

sample_map <- data.table::data.table(
  file = expected,
  sample_id = sub("[.]csv[.]gz$", "", expected),
  donor = sub("^GSM[0-9]+_", "", sub("[.]csv[.]gz$", "", expected)),
  group = ifelse(grepl("^GSM[0-9]+_silicosis", expected), "Silicosis", "Silica-exposed control")
)

# Audit gene universes before loading the dense author-supplied CSV matrices.
genes_by_sample <- lapply(files, function(f) as.character(data.table::fread(f, select = 1)[[1]]))
names(genes_by_sample) <- sample_map$sample_id
date_like <- "^[0-9]{1,2}-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)$"
common_genes <- Reduce(intersect, genes_by_sample)
common_genes <- common_genes[!grepl(date_like, common_genes, ignore.case = TRUE)]
common_genes <- common_genes[!duplicated(common_genes)]
if (length(common_genes) < 10000L) stop("Unexpectedly small common gene universe")

gene_audit <- data.table::rbindlist(lapply(seq_along(files), function(i) {
  g <- genes_by_sample[[i]]
  data.table::data.table(
    sample_id = sample_map$sample_id[i], donor = sample_map$donor[i], group = sample_map$group[i],
    n_genes_raw = length(g), n_duplicate_symbols = sum(duplicated(g)),
    n_date_like_symbols = sum(grepl(date_like, g, ignore.case = TRUE)),
    n_common_unique_genes = length(common_genes), n_missing_counts_replaced_with_zero = NA_integer_
  )
}))
safe_write_csv(sample_map, file.path(result_dir, "source_data", "SD01_GSE174725_donor_manifest.csv"))

read_donor <- function(i) {
  f <- files[i]
  write_log("Reading GSE174725 donor matrix: ", basename(f))
  x <- data.table::fread(f, showProgress = TRUE)
  gene_col <- names(x)[1]
  idx <- match(common_genes, x[[gene_col]])
  if (anyNA(idx)) stop("Common-gene matching failed for ", basename(f))
  barcodes <- names(x)[-1]
  native_genes <- as.character(x[[gene_col]])
  m_native <- as.matrix(x[, -1, with = FALSE])
  storage.mode(m_native) <- "integer"
  n_missing <- sum(is.na(m_native))
  if (n_missing) m_native[is.na(m_native)] <- 0L
  gene_audit[i, n_missing_counts_replaced_with_zero := n_missing]
  write_log("Missing count cells replaced with zero for ", basename(f), ": ", n_missing)
  native_sum <- colSums(m_native)
  native_detected <- colSums(m_native > 0)
  native_mito <- grepl("^MT[-.]", native_genes, ignore.case = TRUE)
  if (!any(native_mito)) stop("No mitochondrial genes detected in native matrix: ", basename(f))
  native_mito_percent <- 100 * colSums(m_native[native_mito, , drop = FALSE]) / pmax(native_sum, 1)
  cell_ids <- paste(sample_map$sample_id[i], barcodes, sep = "__")
  qc_native <- data.table::data.table(cell_id = cell_ids, sum = native_sum,
    detected = native_detected, subsets_Mito_percent = native_mito_percent)
  m <- m_native[idx, , drop = FALSE]
  rm(m_native)
  m <- Matrix::Matrix(m, sparse = TRUE)
  rownames(m) <- common_genes
  colnames(m) <- cell_ids
  rm(x); gc(verbose = FALSE)
  list(counts = m, qc = qc_native)
}

donor_objects <- lapply(seq_along(files), read_donor)
safe_write_csv(gene_audit, file.path(result_dir, "tables", "T01_GSE174725_gene_universe_audit.csv"))
counts <- Reduce(Matrix::cbind2, lapply(donor_objects, `[[`, "counts"))
native_qc <- data.table::rbindlist(lapply(donor_objects, `[[`, "qc"))
rm(donor_objects); gc(verbose = FALSE)

donor_vec <- sub("__.*$", "", colnames(counts))
sample_idx <- match(donor_vec, sample_map$sample_id)
if (anyNA(sample_idx)) stop("Cell-to-donor mapping failed")
sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = counts),
  colData = S4Vectors::DataFrame(
    cell_id = colnames(counts), sample_id = donor_vec,
    donor = sample_map$donor[sample_idx], group = sample_map$group[sample_idx]
  )
)
native_qc <- native_qc[match(colnames(sce), cell_id)]
if (anyNA(native_qc$cell_id)) stop("Native QC metric alignment failed")
sce$sum <- native_qc$sum
sce$detected <- native_qc$detected
sce$subsets_Mito_percent <- native_qc$subsets_Mito_percent

thresholds <- data.table::rbindlist(lapply(unique(sce$sample_id), function(sid) {
  j <- which(sce$sample_id == sid)
  lc <- log10(sce$sum[j] + 1)
  lf <- log10(sce$detected[j] + 1)
  mt <- sce$subsets_Mito_percent[j]
  mt_finite <- mt[is.finite(mt)]
  mt_limit <- if (length(mt_finite)) {
    min(30, max(10, stats::median(mt_finite) + 3 * stats::mad(mt_finite)))
  } else 30
  data.table::data.table(
    sample_id = sid,
    min_log10_counts = stats::median(lc, na.rm = TRUE) - 3 * stats::mad(lc, na.rm = TRUE),
    max_log10_counts = stats::median(lc, na.rm = TRUE) + 4 * stats::mad(lc, na.rm = TRUE),
    min_log10_features = stats::median(lf, na.rm = TRUE) - 3 * stats::mad(lf, na.rm = TRUE),
    max_log10_features = stats::median(lf, na.rm = TRUE) + 4 * stats::mad(lf, na.rm = TRUE),
    max_mito_percent = mt_limit
  )
}))
thr <- thresholds[match(sce$sample_id, sample_id)]
sce$qc_pass_initial <- is.finite(sce$sum) & is.finite(sce$detected) &
  is.finite(sce$subsets_Mito_percent) &
  log10(sce$sum + 1) >= thr$min_log10_counts &
  log10(sce$sum + 1) <= thr$max_log10_counts &
  log10(sce$detected + 1) >= thr$min_log10_features &
  log10(sce$detected + 1) <= thr$max_log10_features &
  sce$subsets_Mito_percent <= thr$max_mito_percent
sce$qc_pass_initial[is.na(sce$qc_pass_initial)] <- FALSE
write_log("GSE174725 adaptive QC retained cells by donor: ", paste(
  names(table(sce$sample_id[sce$qc_pass_initial])),
  as.integer(table(sce$sample_id[sce$qc_pass_initial])), sep = "=", collapse = "; "))

sce_qc <- sce[, sce$qc_pass_initial]
set.seed(20260808)
sce_qc <- scDblFinder::scDblFinder(
  sce_qc, samples = "sample_id", clusters = FALSE,
  BPPARAM = BiocParallel::SerialParam(), verbose = TRUE
)
sce_qc$qc_pass_final <- sce_qc$scDblFinder.class == "singlet"
sce_final <- sce_qc[, sce_qc$qc_pass_final]
sce_final <- scuttle::logNormCounts(sce_final)

cell_meta <- data.table::as.data.table(as.data.frame(SummarizedExperiment::colData(sce_qc)))
safe_write_csv(cell_meta, file.path(result_dir, "source_data", "SD02_GSE174725_cell_QC.csv"))
safe_write_csv(thresholds, file.path(result_dir, "tables", "T02_GSE174725_QC_thresholds.csv"))

qc_summary <- data.table::rbindlist(lapply(sample_map$sample_id, function(sid) {
  n_raw <- sum(sce$sample_id == sid)
  n_initial <- sum(sce_qc$sample_id == sid)
  n_final <- sum(sce_final$sample_id == sid)
  data.table::data.table(
    sample_id = sid, donor = sample_map[sample_id == sid, donor],
    group = sample_map[sample_id == sid, group], n_raw = n_raw,
    n_after_adaptive_qc = n_initial, n_doublets = n_initial - n_final,
    n_final_singlets = n_final, final_retention_percent = 100 * n_final / n_raw
  )
}))
safe_write_csv(qc_summary, file.path(result_dir, "tables", "T03_GSE174725_QC_summary.csv"))
saveRDS(sce_final, file.path(processed_dir, paste0(acc, "_sce_qc.rds")), compress = FALSE)

plot_qc <- data.table::as.data.table(as.data.frame(SummarizedExperiment::colData(sce)))
plot_qc[, donor_label := factor(gsub("_", " ", donor), levels = gsub("_", " ", sample_map$donor))]
p_a <- ggplot2::ggplot(plot_qc, ggplot2::aes(donor_label, detected, fill = group)) +
  ggplot2::geom_violin(scale = "width", linewidth = 0.25, trim = TRUE) +
  ggplot2::scale_y_log10() +
  ggplot2::scale_fill_manual(values = c("Silica-exposed control" = "#416A9A", "Silicosis" = "#B84A4A")) +
  ggplot2::labs(title = "Detected features before filtering", x = NULL, y = "Detected genes", fill = NULL)
p_b <- ggplot2::ggplot(plot_qc, ggplot2::aes(donor_label, subsets_Mito_percent, fill = group)) +
  ggplot2::geom_violin(scale = "width", linewidth = 0.25, trim = TRUE) +
  ggplot2::scale_fill_manual(values = c("Silica-exposed control" = "#416A9A", "Silicosis" = "#B84A4A")) +
  ggplot2::labs(title = "Mitochondrial fraction before filtering", x = NULL,
                y = "Mitochondrial reads (%)", fill = NULL) +
  ggplot2::theme(legend.position = "none")
q_long <- data.table::melt(qc_summary,
  id.vars = c("sample_id", "donor", "group"),
  measure.vars = c("n_raw", "n_after_adaptive_qc", "n_final_singlets"),
  variable.name = "stage", value.name = "cells")
q_long[, stage := factor(stage, levels = c("n_raw", "n_after_adaptive_qc", "n_final_singlets"),
  labels = c("Raw", "Adaptive QC", "Final singlets"))]
q_long[, donor_label := factor(gsub("_", " ", donor), levels = gsub("_", " ", sample_map$donor))]
p_c <- ggplot2::ggplot(q_long, ggplot2::aes(donor_label, cells, fill = stage)) +
  ggplot2::geom_col(position = "dodge", width = 0.75) +
  ggplot2::scale_fill_manual(values = c("Raw" = "#C9C9C9", "Adaptive QC" = "#6AA6A1",
                                        "Final singlets" = "#416A9A")) +
  ggplot2::labs(title = "Cell retention by donor", x = NULL, y = "Cells", fill = NULL)

theme_sc <- ggplot2::theme_classic(base_size = 7, base_family = "Arial") +
  ggplot2::theme(axis.line = ggplot2::element_line(linewidth = 0.3),
    axis.ticks = ggplot2::element_line(linewidth = 0.3), legend.position = "top",
    axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
    plot.title = ggplot2::element_text(size = 8, face = "bold"),
    plot.tag = ggplot2::element_text(size = 8, face = "bold"))
fig <- (p_a | p_b | p_c) + patchwork::plot_annotation(
  tag_levels = "a", title = "GSE174725 donor-resolved quality control",
  subtitle = "Sample-specific MAD thresholds followed by per-donor doublet detection") & theme_sc
out <- file.path(result_dir, "figures", "F06_GSE174725_QC")
w <- 183 / 25.4; h <- 82 / 25.4
svglite::svglite(paste0(out, ".svg"), width = w, height = h); print(fig); grDevices::dev.off()
assert_english_only_svg(paste0(out, ".svg"))
grDevices::cairo_pdf(paste0(out, ".pdf"), width = w, height = h, family = "Arial"); print(fig); grDevices::dev.off()
ragg::agg_tiff(paste0(out, ".tiff"), width = w, height = h, units = "in", res = 600,
               compression = "lzw"); print(fig); grDevices::dev.off()
ragg::agg_png(paste0(out, ".png"), width = w, height = h, units = "in", res = 300);
print(fig); grDevices::dev.off()

write_log("GSE174725 donor QC completed: raw cells=", ncol(sce),
          "; final singlets=", ncol(sce_final), "; common genes=", nrow(sce_final),
          "; English-only SVG QA passed")
