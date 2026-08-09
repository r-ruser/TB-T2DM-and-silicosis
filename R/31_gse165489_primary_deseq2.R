source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "DESeq2", "edgeR", "ggplot2", "patchwork",
                  "svglite", "ragg"))

acc <- "GSE165489"
processed_dir <- file.path(path_data, "02_GEO_bulk", acc, "processed")
txi_file <- file.path(processed_dir, paste0(acc, "_tximport.rds"))
if (!file.exists(txi_file)) stop("Run R/30_prepare_gse165489_tximport.R first")
x <- readRDS(txi_file)
manifest <- data.table::fread(file.path(path_result, "02_bulk", acc,
                                        "A_GSE165489_salmon_quant_manifest.csv"))
library_threshold_log10 <- stats::median(log10(manifest$total_salmon_numreads + 1)) -
  3 * stats::mad(log10(manifest$total_salmon_numreads + 1))
manifest[, `:=`(
  library_threshold_log10 = library_threshold_log10,
  library_qc_pass = log10(total_salmon_numreads + 1) >= library_threshold_log10,
  catastrophic_qc_pass = total_salmon_numreads >= 1e6
)]
meta_all <- merge(data.table::copy(x$metadata), manifest[, .(sample_id, total_salmon_numreads,
  imported_count_fraction, library_threshold_log10, library_qc_pass, catastrophic_qc_pass)],
  by = "sample_id", all.x = TRUE, sort = FALSE)
meta_inclusive <- meta_all[include_primary == TRUE]
meta <- meta_inclusive[library_qc_pass == TRUE]
meta[, group := factor(group, levels = c("exposed_no_silicosis", "silicosis"))]
meta[, exposure_z := as.numeric(scale(exposure_years))]
samples <- meta$sample_id
if (!all(samples %in% colnames(x$raw$counts))) stop("GSE165489 metadata/count mismatch")

scaled_counts <- x$length_scaled$counts[, samples, drop = FALSE]
keep <- edgeR::filterByExpr(scaled_counts, group = meta$group)
if (sum(keep) < 5000L) stop("Unexpectedly few expressed genes")
txi_de <- lapply(x$raw[c("abundance", "counts", "length")], function(z) z[keep, samples, drop = FALSE])
txi_de$countsFromAbundance <- x$raw$countsFromAbundance
rownames(meta) <- meta$sample_id
design <- stats::model.matrix(~ exposure_z + group, data = meta)
if (qr(design)$rank < ncol(design)) stop("Adjusted GSE165489 design is rank deficient")

dds <- DESeq2::DESeqDataSetFromTximport(txi_de, colData = as.data.frame(meta),
                                        design = ~ exposure_z + group)
dds <- DESeq2::DESeq(dds, quiet = TRUE)
coef_name <- grep("^group_silicosis_vs_exposed_no_silicosis$", DESeq2::resultsNames(dds), value = TRUE)
if (length(coef_name) != 1L) stop("Expected adjusted DESeq2 coefficient not found")
res <- DESeq2::results(dds, name = coef_name, alpha = 0.05)
shrunk <- DESeq2::lfcShrink(dds, coef = coef_name, res = res, type = "normal")

dds_unadjusted <- dds
DESeq2::design(dds_unadjusted) <- ~ group
dds_unadjusted <- DESeq2::DESeq(dds_unadjusted, quiet = TRUE)
coef_unadj <- grep("^group_silicosis_vs_exposed_no_silicosis$",
                   DESeq2::resultsNames(dds_unadjusted), value = TRUE)
res_unadj <- DESeq2::results(dds_unadjusted, name = coef_unadj)

fit_adjusted_sensitivity <- function(meta_sensitivity) {
  meta_sensitivity <- data.table::copy(meta_sensitivity)
  meta_sensitivity[, group := factor(group, levels = c("exposed_no_silicosis", "silicosis"))]
  meta_sensitivity[, exposure_z := as.numeric(scale(exposure_years))]
  ids <- meta_sensitivity$sample_id
  rownames(meta_sensitivity) <- ids
  txi_s <- lapply(x$raw[c("abundance", "counts", "length")],
                  function(z) z[keep, ids, drop = FALSE])
  txi_s$countsFromAbundance <- x$raw$countsFromAbundance
  d <- DESeq2::DESeqDataSetFromTximport(txi_s, colData = as.data.frame(meta_sensitivity),
                                       design = ~ exposure_z + group)
  d <- DESeq2::DESeq(d, quiet = TRUE)
  cn <- grep("^group_silicosis_vs_exposed_no_silicosis$", DESeq2::resultsNames(d), value = TRUE)
  if (length(cn) != 1L) stop("Sensitivity model coefficient missing")
  list(dds = d, result = DESeq2::results(d, name = cn))
}
inclusive_fit <- fit_adjusted_sensitivity(meta_inclusive)
catastrophic_fit <- fit_adjusted_sensitivity(meta_inclusive[catastrophic_qc_pass == TRUE])

gene_map <- data.table::fread(file.path(processed_dir, paste0(acc, "_gene_map.csv")))
symbol <- gene_map$gene_symbol[match(rownames(shrunk), gene_map$entrez_id)]
de <- data.table::data.table(
  entrez_id = rownames(shrunk), gene_symbol = symbol,
  baseMean = shrunk$baseMean, log2FoldChange = shrunk$log2FoldChange,
  lfcSE = shrunk$lfcSE, stat = shrunk$stat, pvalue = shrunk$pvalue, padj = shrunk$padj,
  unadjusted_log2FoldChange = res_unadj$log2FoldChange,
  unadjusted_padj = res_unadj$padj,
  inclusive_adjusted_log2FoldChange = inclusive_fit$result$log2FoldChange,
  inclusive_adjusted_padj = inclusive_fit$result$padj,
  catastrophic_adjusted_log2FoldChange = catastrophic_fit$result$log2FoldChange,
  catastrophic_adjusted_padj = catastrophic_fit$result$padj
)
de[, `:=`(ci_low = log2FoldChange - 1.96 * lfcSE,
          ci_high = log2FoldChange + 1.96 * lfcSE)]
de[, evidence_class := data.table::fcase(
  !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1, "FDR<0.05 and |log2FC|>=1",
  !is.na(padj) & padj < 0.05, "FDR<0.05", default = "Not FDR-significant")]
de[, abs_lfc := abs(log2FoldChange)]
data.table::setorder(de, padj, -abs_lfc, na.last = TRUE)

y <- edgeR::DGEList(counts = scaled_counts[keep, , drop = FALSE], samples = as.data.frame(meta))
y <- edgeR::calcNormFactors(y)
y <- edgeR::estimateDisp(y, design, robust = TRUE)
fit_q <- edgeR::glmQLFit(y, design, robust = TRUE)
coef_idx <- grep("^groupsilicosis$", colnames(design))
if (length(coef_idx) != 1L) stop("edgeR silicosis coefficient not found")
qlf <- edgeR::glmQLFTest(fit_q, coef = coef_idx)
edger <- data.table::as.data.table(edgeR::topTags(qlf, n = Inf, adjust.method = "BH")$table,
                                   keep.rownames = "entrez_id")
data.table::setnames(edger, c("FDR", "PValue"), c("edgeR_FDR", "edgeR_PValue"), skip_absent = TRUE)

comparison <- merge(de[, .(entrez_id, adjusted_log2FC = log2FoldChange, adjusted_padj = padj,
                            unadjusted_log2FoldChange, unadjusted_padj)],
                    edger[, .(entrez_id, edgeR_log2FC = logFC, edgeR_PValue, edgeR_FDR)],
                    by = "entrez_id")
summary_tab <- data.table::data.table(
  metric = c("DESeq2 adjusted FDR<0.05", "DESeq2 unadjusted FDR<0.05", "edgeR adjusted FDR<0.05",
             "DESeq2 inclusive adjusted FDR<0.05", "DESeq2 exclude <1M adjusted FDR<0.05",
             "Adjusted vs unadjusted Spearman log2FC", "DESeq2 vs edgeR Spearman log2FC",
             "DESeq2 vs edgeR direction agreement"),
  value = c(sum(comparison$adjusted_padj < 0.05, na.rm = TRUE),
            sum(comparison$unadjusted_padj < 0.05, na.rm = TRUE),
            sum(comparison$edgeR_FDR < 0.05, na.rm = TRUE),
            sum(inclusive_fit$result$padj < 0.05, na.rm = TRUE),
            sum(catastrophic_fit$result$padj < 0.05, na.rm = TRUE),
            stats::cor(comparison$adjusted_log2FC, comparison$unadjusted_log2FoldChange,
                       method = "spearman", use = "complete.obs"),
            stats::cor(comparison$adjusted_log2FC, comparison$edgeR_log2FC,
                       method = "spearman", use = "complete.obs"),
            mean(sign(comparison$adjusted_log2FC) == sign(comparison$edgeR_log2FC), na.rm = TRUE))
)

result_dir <- file.path(path_result, "02_bulk", acc)
invisible(lapply(file.path(result_dir, c("tables", "source_data", "models", "figures")),
                 dir.create, recursive = TRUE, showWarnings = FALSE))
safe_write_csv(de, file.path(result_dir, "tables", "T01_GSE165489_silicosis_vs_exposed_DESeq2.csv"))
safe_write_csv(edger, file.path(result_dir, "tables", "T02_GSE165489_silicosis_vs_exposed_edgeR.csv"))
safe_write_csv(summary_tab, file.path(result_dir, "tables", "T03_GSE165489_sensitivity_summary.csv"))
safe_write_csv(meta_all, file.path(result_dir, "tables", "T04_GSE165489_library_QC.csv"))
safe_write_csv(comparison, file.path(result_dir, "source_data", "SD01_GSE165489_method_concordance.csv"))
saveRDS(list(dds = dds, dds_unadjusted = dds_unadjusted,
             dds_inclusive = inclusive_fit$dds, dds_exclude_lt_1M = catastrophic_fit$dds,
             edgeR_fit = fit_q,
             design = design, filter = keep),
        file.path(result_dir, "models", "M01_GSE165489_primary_models.rds"))

vsd <- DESeq2::vst(dds, blind = FALSE)
pca <- stats::prcomp(t(SummarizedExperiment::assay(vsd)))
ve <- 100 * pca$sdev^2 / sum(pca$sdev^2)
pca_df <- data.table::data.table(sample_id = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2])
pca_df <- merge(pca_df, meta[, .(sample_id, group, exposure_years)], by = "sample_id", sort = FALSE)
safe_write_csv(pca_df, file.path(result_dir, "source_data", "SD02_GSE165489_PCA.csv"))

pal <- c(exposed_no_silicosis = "#416A9A", silicosis = "#B84A4A")
labels <- c(exposed_no_silicosis = "Exposed, no silicosis", silicosis = "Silicosis")
theme_nature <- ggplot2::theme_classic(base_size = 7, base_family = "Arial") +
  ggplot2::theme(axis.line = ggplot2::element_line(linewidth = 0.3),
    axis.ticks = ggplot2::element_line(linewidth = 0.3), legend.position = "top",
    plot.title = ggplot2::element_text(size = 8, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 6, colour = "#767676"),
    plot.tag = ggplot2::element_text(size = 8, face = "bold"))
p_a <- ggplot2::ggplot(meta, ggplot2::aes(group, exposure_years, fill = group)) +
  ggplot2::geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.35) +
  ggplot2::geom_jitter(width = 0.12, size = 0.8, alpha = 0.65) +
  ggplot2::scale_fill_manual(values = pal, labels = labels) +
  ggplot2::scale_x_discrete(labels = labels) +
  ggplot2::labs(title = "Exposure-duration imbalance", subtitle = "Primary model adjusts for exposure duration",
                x = NULL, y = "Silica exposure (years)") + theme_nature +
  ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))
p_b <- ggplot2::ggplot(pca_df, ggplot2::aes(PC1, PC2, colour = group)) +
  ggplot2::geom_point(size = 1.7, alpha = 0.85) +
  ggplot2::scale_colour_manual(values = pal, labels = labels) +
  ggplot2::labs(title = "Variance-stabilized PCA", subtitle = "Two low-depth libraries excluded",
                x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]), colour = NULL) + theme_nature
vol <- data.table::copy(de)
vol[, neglog10_fdr := -log10(pmax(padj, .Machine$double.xmin))]
vol[is.na(neglog10_fdr), neglog10_fdr := 0]
lab <- vol[!is.na(gene_symbol) & evidence_class != "Not FDR-significant"][1:min(.N, 12)]
p_c <- ggplot2::ggplot(vol, ggplot2::aes(log2FoldChange, neglog10_fdr, colour = evidence_class)) +
  ggplot2::geom_point(size = 0.5, alpha = 0.55) +
  ggplot2::geom_vline(xintercept = c(-1, 1), linetype = 2, linewidth = 0.3, colour = "#767676") +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.3, colour = "#767676") +
  ggplot2::geom_text(data = lab, ggplot2::aes(label = gene_symbol), size = 1.8,
                     check_overlap = TRUE, colour = "#272727", vjust = -0.5) +
  ggplot2::scale_colour_manual(values = c("FDR<0.05 and |log2FC|>=1" = "#B84A4A",
    "FDR<0.05" = "#416A9A", "Not FDR-significant" = "#C9C9C9"), drop = FALSE) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.20))) +
  ggplot2::labs(title = "Silicosis versus exposed controls",
                subtitle = "DESeq2; exposure-duration adjusted; BH-FDR",
                x = "Shrunken log2 fold change", y = expression(-log[10](FDR))) +
  theme_nature + ggplot2::theme(legend.position = "none")

fig <- (p_a | p_b | p_c) + patchwork::plot_annotation(
  tag_levels = "a", title = "Primary silicosis host-response contrast",
  subtitle = "GSE165489 whole blood RNA-seq; exposure-duration adjusted; low-depth libraries removed")
out <- file.path(result_dir, "figures", "F07_GSE165489_primary_silicosis_vs_exposed")
w <- 183 / 25.4; h <- 78 / 25.4
svglite::svglite(paste0(out, ".svg"), width = w, height = h); print(fig); grDevices::dev.off()
assert_english_only_svg(paste0(out, ".svg"))
grDevices::cairo_pdf(paste0(out, ".pdf"), width = w, height = h, family = "Arial"); print(fig); grDevices::dev.off()
ragg::agg_tiff(paste0(out, ".tiff"), width = w, height = h, units = "in", res = 600,
               compression = "lzw"); print(fig); grDevices::dev.off()
ragg::agg_png(paste0(out, ".png"), width = w, height = h, units = "in", res = 300);
print(fig); grDevices::dev.off()

write_log("GSE165489 adjusted primary analysis completed: DESeq2 FDR-significant genes=",
          sum(de$padj < 0.05, na.rm = TRUE), "; English-only SVG QA passed")
