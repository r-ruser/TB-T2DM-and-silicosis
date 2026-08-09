source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "SingleCellExperiment", "SummarizedExperiment", "Matrix",
  "ggplot2", "patchwork", "svglite", "ragg"))

up_file <- file.path(path_result, "03_cross_disease", "RRHO", "tables", "T03_RRHO_shared_up_peak_genes.csv")
down_file <- file.path(path_result, "03_cross_disease", "RRHO", "tables", "T04_RRHO_shared_down_peak_genes.csv")
up_genes <- unique(data.table::fread(up_file)$gene_symbol)[1:100]
down_genes <- unique(data.table::fread(down_file)$gene_symbol)[1:100]
up_genes <- up_genes[!is.na(up_genes)]; down_genes <- down_genes[!is.na(down_genes)]

out_dir <- file.path(path_result, "04_scRNA", "cross_disease")
invisible(lapply(file.path(out_dir, c("figures", "tables", "source_data")), dir.create,
  recursive = TRUE, showWarnings = FALSE))

process_dataset <- function(acc, dataset_label) {
  sce_file <- file.path(path_data, "03_scRNA", acc, "processed", paste0(acc, "_sce_qc.rds"))
  ann_file <- file.path(path_result, "04_scRNA", acc, "source_data",
    paste0("SD03_", acc, "_cell_annotation_UMAP.csv"))
  write_log("Shared signature localization: loading ", acc)
  sce <- readRDS(sce_file)
  ann <- data.table::fread(ann_file)
  idx <- match(colnames(sce), ann$cell_id)
  if (anyNA(idx)) stop("Annotation mismatch for ", acc)
  ann <- ann[idx]
  if (acc == "GSE192483") {
    symbols <- as.character(SummarizedExperiment::rowData(sce)$gene_symbol)
  } else {
    symbols <- rownames(sce)
  }
  first_row <- !duplicated(symbols) & nzchar(symbols) & !is.na(symbols)
  symbol_row <- stats::setNames(which(first_row), symbols[first_row])
  up_rows <- unname(symbol_row[intersect(up_genes, names(symbol_row))])
  down_rows <- unname(symbol_row[intersect(down_genes, names(symbol_row))])
  if (length(up_rows) < 50L || length(down_rows) < 50L) stop("Insufficient signature coverage in ", acc)
  rows <- c(up_rows, down_rows)
  x <- as.matrix(SummarizedExperiment::assay(sce, "logcounts")[rows, , drop = FALSE])
  row_mu <- rowMeans(x)
  row_sd <- apply(x, 1, stats::sd)
  row_sd[!is.finite(row_sd) | row_sd == 0] <- 1
  z <- sweep(sweep(x, 1, row_mu, "-"), 1, row_sd, "/")
  score <- colMeans(z[seq_along(up_rows), , drop = FALSE]) -
    colMeans(z[length(up_rows) + seq_along(down_rows), , drop = FALSE])
  rm(x, z, sce); gc(verbose = FALSE)
  ann[, signature_score := score]
  ann[, dataset := dataset_label]
  cap <- stats::quantile(ann$signature_score, c(0.01, 0.99), na.rm = TRUE)
  ann[, score_plot := pmax(cap[1], pmin(cap[2], signature_score))]
  donor <- ann[, .(donor_mean_score = mean(signature_score), cells = .N,
    percent_positive = 100 * mean(signature_score > 0)),
    by = .(dataset, sample_id, cell_type)]
  coverage <- data.table::data.table(dataset = dataset_label, accession = acc,
    requested_up = length(up_genes), matched_up = length(up_rows),
    requested_down = length(down_genes), matched_down = length(down_rows), cells = nrow(ann))
  list(cells = ann[, .(cell_id, UMAP1, UMAP2, cell_type, sample_id, dataset,
    signature_score, score_plot)], donor = donor, coverage = coverage)
}

x174 <- process_dataset("GSE174725", "Silicosis BALF")
x192 <- process_dataset("GSE192483", "TB lung")
cells <- data.table::rbindlist(list(x174$cells, x192$cells), fill = TRUE)
donor <- data.table::rbindlist(list(x174$donor, x192$donor), fill = TRUE)
coverage <- data.table::rbindlist(list(x174$coverage, x192$coverage))
summary_dt <- donor[, .(donors = .N, cells = sum(cells),
  median_donor_score = stats::median(donor_mean_score),
  q1_donor_score = stats::quantile(donor_mean_score, 0.25),
  q3_donor_score = stats::quantile(donor_mean_score, 0.75),
  median_percent_positive = stats::median(percent_positive)), by = .(dataset, cell_type)]
safe_write_csv(cells, file.path(out_dir, "source_data", "SD01_RRHO_local_signature_cell_scores.csv"))
safe_write_csv(donor, file.path(out_dir, "source_data", "SD02_RRHO_local_signature_donor_scores.csv"))
safe_write_csv(summary_dt, file.path(out_dir, "tables", "T01_RRHO_local_signature_cell_type_summary.csv"))
safe_write_csv(coverage, file.path(out_dir, "tables", "T02_RRHO_local_signature_coverage.csv"))

lim <- max(abs(stats::quantile(cells$score_plot, c(0.01, 0.99), na.rm = TRUE)))
score_scale <- ggplot2::scale_colour_gradient2(low = "#416A9A", mid = "#F2F2F2", high = "#B84A4A",
  midpoint = 0, limits = c(-lim, lim), oob = scales::squish, name = "Signature score")
p_a <- ggplot2::ggplot(cells[dataset == "Silicosis BALF"],
  ggplot2::aes(UMAP1, UMAP2, colour = score_plot)) +
  ggplot2::geom_point(size = 0.12, alpha = 0.62) + score_scale +
  ggplot2::labs(title = "Silicosis BALF", subtitle = "GSE174725",
    x = "UMAP 1", y = "UMAP 2")
p_b <- ggplot2::ggplot(cells[dataset == "TB lung"],
  ggplot2::aes(UMAP1, UMAP2, colour = score_plot)) +
  ggplot2::geom_point(size = 0.07, alpha = 0.58) + score_scale +
  ggplot2::labs(title = "TB lung", subtitle = "GSE192483",
    x = "UMAP 1", y = "UMAP 2")
p_c <- ggplot2::ggplot(summary_dt, ggplot2::aes(cell_type, dataset,
  colour = median_donor_score, size = median_percent_positive)) +
  ggplot2::geom_point(alpha = 0.90) +
  ggplot2::scale_colour_gradient2(low = "#416A9A", mid = "#F2F2F2", high = "#B84A4A",
    midpoint = 0, name = "Median donor score") +
  ggplot2::scale_size_continuous(range = c(1.5, 6), limits = c(0, 100), name = "Cells above zero (%)") +
  ggplot2::labs(title = "Donor-resolved cell-type localization",
    subtitle = "Colour: median donor mean; size: median positive-cell fraction",
    x = NULL, y = NULL)
theme_sc <- ggplot2::theme_classic(base_size = 7, base_family = "Arial") +
  ggplot2::theme(axis.line = ggplot2::element_line(linewidth = 0.3),
    axis.ticks = ggplot2::element_line(linewidth = 0.3),
    legend.key.height = grid::unit(3, "mm"), legend.key.width = grid::unit(4, "mm"),
    legend.text = ggplot2::element_text(size = 5.4), legend.title = ggplot2::element_text(size = 5.8),
    plot.title = ggplot2::element_text(size = 8, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 6, colour = "#767676"),
    plot.tag = ggplot2::element_text(size = 8, face = "bold"))
p_a <- p_a + theme_sc + ggplot2::theme(legend.position = "bottom", aspect.ratio = 1)
p_b <- p_b + theme_sc + ggplot2::theme(legend.position = "bottom", aspect.ratio = 1)
p_c <- p_c + theme_sc + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
  legend.position = "right")
fig <- ((p_a | p_b) / p_c) + patchwork::plot_layout(heights = c(2.15, 1)) +
  patchwork::plot_annotation(tag_levels = "a", title = "Exploratory RRHO-local signature localization",
    subtitle = "Top 100 shared-up minus top 100 shared-down genes; within-dataset standardized expression")
out <- file.path(out_dir, "figures", "F17_RRHO_local_signature_single_cell")
w <- 183 / 25.4; h <- 164 / 25.4
svglite::svglite(paste0(out, ".svg"), width = w, height = h); print(fig); grDevices::dev.off()
assert_english_only_svg(paste0(out, ".svg"))
grDevices::cairo_pdf(paste0(out, ".pdf"), width = w, height = h, family = "Arial"); print(fig); grDevices::dev.off()
ragg::agg_tiff(paste0(out, ".tiff"), width = w, height = h, units = "in", res = 600,
  compression = "lzw"); print(fig); grDevices::dev.off()
ragg::agg_png(paste0(out, ".png"), width = w, height = h, units = "in", res = 300)
print(fig); grDevices::dev.off()
write_log("Shared RRHO-local signature localization completed; English-only SVG QA passed")
