source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "DESeq2", "edgeR", "AnnotationDbi", "org.Hs.eg.db",
                  "ggplot2", "patchwork", "svglite", "ragg"))

acc <- "GSE114192"
counts_file <- file.path(path_data, "02_GEO_bulk", acc, "processed", paste0(acc, "_counts.tsv.gz"))
meta_file <- file.path(path_data, "02_GEO_bulk", acc, "processed", paste0(acc, "_metadata_curated.csv"))
if (!file.exists(counts_file) || !file.exists(meta_file)) stop("Prepared GSE114192 inputs missing")

ct <- data.table::fread(counts_file)
gene_id <- ct[[1]]
counts <- as.matrix(ct[, -1, with = FALSE])
rownames(counts) <- gene_id
storage.mode(counts) <- "integer"
meta <- data.table::fread(meta_file, encoding = "UTF-8")
meta <- meta[include_primary == TRUE]
meta[, group := stats::relevel(factor(group), ref = "TB_only")]
meta[, site := factor(site)]
if (!all(meta$sample_id %in% colnames(counts))) stop("Metadata/count mismatch")
counts <- counts[, meta$sample_id, drop = FALSE]
rownames(meta) <- meta$sample_id

design_check <- stats::model.matrix(~ site + group, data = meta)
if (qr(design_check)$rank < ncol(design_check)) stop("Design ~ site + group is not full rank")
keep <- edgeR::filterByExpr(counts, group = meta$group)
counts_f <- counts[keep, , drop = FALSE]

dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts_f, colData = as.data.frame(meta),
                                      design = ~ site + group)
dds <- DESeq2::DESeq(dds, quiet = TRUE)
coef_name <- grep("^group_TB_DM_vs_TB_only$", DESeq2::resultsNames(dds), value = TRUE)
if (length(coef_name) != 1L) stop("Expected DESeq2 coefficient not found: ", paste(DESeq2::resultsNames(dds), collapse = ", "))
res <- DESeq2::results(dds, name = coef_name, alpha = 0.05, independentFiltering = TRUE)
shrunk <- DESeq2::lfcShrink(dds, coef = coef_name, res = res, type = "normal")

ens <- sub("[.][0-9]+$", "", rownames(shrunk))
symbol <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = ens, keytype = "ENSEMBL",
                                column = "SYMBOL", multiVals = "first")
de <- data.table::data.table(
  gene_id = rownames(shrunk), gene_symbol = unname(symbol),
  baseMean = shrunk$baseMean, log2FoldChange = shrunk$log2FoldChange,
  lfcSE = shrunk$lfcSE, stat = shrunk$stat, pvalue = shrunk$pvalue, padj = shrunk$padj
)
de[, `:=`(ci_low = log2FoldChange - 1.96 * lfcSE, ci_high = log2FoldChange + 1.96 * lfcSE)]
de[, evidence_class := data.table::fcase(
  !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1, "FDR<0.05 and |log2FC|>=1",
  !is.na(padj) & padj < 0.05, "FDR<0.05",
  default = "Not FDR-significant")]
de[, abs_lfc := abs(log2FoldChange)]
data.table::setorder(de, padj, -abs_lfc, na.last = TRUE)

# edgeR quasi-likelihood sensitivity using the same design and filter.
y <- edgeR::DGEList(counts = counts_f, samples = as.data.frame(meta))
y <- edgeR::calcNormFactors(y)
y <- edgeR::estimateDisp(y, design_check, robust = TRUE)
fit_q <- edgeR::glmQLFit(y, design_check, robust = TRUE)
coef_idx <- grep("^groupTB_DM$", colnames(design_check))
if (length(coef_idx) != 1L) stop("edgeR group coefficient not found")
qlf <- edgeR::glmQLFTest(fit_q, coef = coef_idx)
edger <- data.table::as.data.table(edgeR::topTags(qlf, n = Inf, adjust.method = "BH")$table,
                                   keep.rownames = "gene_id")
data.table::setnames(edger, c("FDR", "PValue"), c("edgeR_FDR", "edgeR_PValue"), skip_absent = TRUE)
comparison <- merge(de[, .(gene_id, DESeq2_log2FC = log2FoldChange, DESeq2_padj = padj)],
                    edger[, .(gene_id, edgeR_log2FC = logFC, edgeR_PValue, edgeR_FDR)], by = "gene_id")
concordance <- data.frame(
  metric = c("Spearman log2FC", "Direction agreement", "DESeq2 FDR<0.05", "edgeR FDR<0.05"),
  value = c(stats::cor(comparison$DESeq2_log2FC, comparison$edgeR_log2FC, method = "spearman", use = "complete.obs"),
            mean(sign(comparison$DESeq2_log2FC) == sign(comparison$edgeR_log2FC), na.rm = TRUE),
            sum(comparison$DESeq2_padj < 0.05, na.rm = TRUE), sum(comparison$edgeR_FDR < 0.05, na.rm = TRUE))
)

result_dir <- file.path(path_result, "02_bulk", acc)
dirs <- file.path(result_dir, c("tables", "figures", "source_data", "models"))
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
safe_write_csv(de, file.path(result_dir, "tables", "T01_GSE114192_TBDM_vs_TB_DESeq2.csv"))
safe_write_csv(edger, file.path(result_dir, "tables", "T02_GSE114192_TBDM_vs_TB_edgeR.csv"))
safe_write_csv(comparison, file.path(result_dir, "source_data", "SD01_GSE114192_method_concordance.csv"))
safe_write_csv(concordance, file.path(result_dir, "tables", "T03_GSE114192_sensitivity_summary.csv"))
saveRDS(list(dds = dds, edgeR_fit = fit_q, design = design_check),
        file.path(result_dir, "models", "M01_GSE114192_primary_models.rds"))

vsd <- DESeq2::vst(dds, blind = FALSE)
pca <- stats::prcomp(t(SummarizedExperiment::assay(vsd)))
ve <- 100 * pca$sdev^2 / sum(pca$sdev^2)
pca_df <- data.table::data.table(sample_id = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2])
pca_df <- merge(pca_df, meta[, .(sample_id, group, site)], by = "sample_id", sort = FALSE)
pca_df[, site := factor(gsub("_", " ", as.character(site)))]
qc <- data.table::data.table(sample_id = colnames(counts_f), library_size = colSums(counts_f))
qc <- merge(qc, meta[, .(sample_id, group, site)], by = "sample_id", sort = FALSE)

pal <- c(TB_only = "#416A9A", TB_DM = "#B84A4A")
group_labels <- c(TB_only = "TB-only", TB_DM = "TB-DM")
theme_nature <- ggplot2::theme_classic(base_size = 7, base_family = "Arial") +
  ggplot2::theme(axis.line = ggplot2::element_line(linewidth = 0.3),
    axis.ticks = ggplot2::element_line(linewidth = 0.3), legend.position = "top",
    plot.title = ggplot2::element_text(size = 8, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 6, colour = "#767676"),
    plot.tag = ggplot2::element_text(size = 8, face = "bold"))
p_a <- ggplot2::ggplot(qc, ggplot2::aes(group, library_size / 1e6, fill = group)) +
  ggplot2::geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.35) +
  ggplot2::geom_jitter(width = 0.12, size = 0.7, alpha = 0.55) +
  ggplot2::scale_fill_manual(values = pal, labels = group_labels) +
  ggplot2::scale_x_discrete(labels = group_labels) +
  ggplot2::labs(title = "Library-size audit", subtitle = sprintf("%d genes retained by filterByExpr", nrow(counts_f)),
                x = NULL, y = "Library size (million reads)") + theme_nature +
  ggplot2::theme(legend.position = "none")
p_b <- ggplot2::ggplot(pca_df, ggplot2::aes(PC1, PC2, colour = group, shape = site)) +
  ggplot2::geom_point(size = 1.7, alpha = 0.85) +
  ggplot2::scale_colour_manual(values = pal, labels = group_labels) +
  ggplot2::labs(title = "Variance-stabilized PCA", subtitle = "Site retained in the DESeq2 design",
                x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]), colour = NULL, shape = "Site") + theme_nature
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
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.22))) +
  ggplot2::labs(title = "TB-DM versus TB-only", subtitle = "DESeq2; site-adjusted; BH-FDR",
                x = "Shrunken log2 fold change", y = expression(-log[10](FDR))) + theme_nature +
  ggplot2::theme(legend.position = "none")

fig <- (p_a | p_b | p_c) + patchwork::plot_annotation(tag_levels = "a",
  title = "Primary TB-diabetes host-response contrast",
  subtitle = "GSE114192 whole blood RNA-seq; association within a tuberculosis background")
safe_write_csv(pca_df, file.path(result_dir, "source_data", "SD02_GSE114192_PCA.csv"))
safe_write_csv(qc, file.path(result_dir, "source_data", "SD03_GSE114192_library_QC.csv"))
out <- file.path(result_dir, "figures", "F03_GSE114192_primary_TBDM_vs_TB")
w <- 183 / 25.4; h <- 72 / 25.4
svglite::svglite(paste0(out, ".svg"), width = w, height = h); print(fig); grDevices::dev.off()
assert_english_only_svg(paste0(out, ".svg"))
grDevices::cairo_pdf(paste0(out, ".pdf"), width = w, height = h, family = "Arial"); print(fig); grDevices::dev.off()
ragg::agg_tiff(paste0(out, ".tiff"), width = w, height = h, units = "in", res = 600, compression = "lzw"); print(fig); grDevices::dev.off()
ragg::agg_png(paste0(out, ".png"), width = w, height = h, units = "in", res = 300); print(fig); grDevices::dev.off()
write_log("GSE114192 primary DESeq2 analysis completed: ", sum(de$padj < 0.05, na.rm = TRUE),
          " FDR-significant genes; English-only SVG QA passed")
