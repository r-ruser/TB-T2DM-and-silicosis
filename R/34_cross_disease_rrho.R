source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "RRHO", "ggplot2", "patchwork", "svglite", "ragg"))

tb_file <- file.path(path_result, "02_bulk", "GSE114192", "tables",
                     "T01_GSE114192_TBDM_vs_TB_DESeq2.csv")
sil_file <- file.path(path_result, "02_bulk", "GSE165489", "tables",
                      "T01_GSE165489_silicosis_vs_exposed_DESeq2.csv")
tb_gsea_file <- file.path(path_result, "03_cross_disease", "GSE114192", "tables",
                          "T01_GSE114192_Hallmark_Reactome_GSEA.csv")
sil_gsea_file <- file.path(path_result, "03_cross_disease", "GSE165489", "tables",
                           "T01_GSE165489_Hallmark_Reactome_GSEA.csv")
if (!all(file.exists(c(tb_file, sil_file, tb_gsea_file, sil_gsea_file)))) {
  stop("Primary DE or GSEA inputs are missing")
}

collapse_rank <- function(x) {
  x <- x[!is.na(gene_symbol) & gene_symbol != "" & is.finite(stat)]
  x[, abs_stat := abs(stat)]
  data.table::setorder(x, gene_symbol, -abs_stat)
  x[, .SD[1], by = gene_symbol]
}
tb <- collapse_rank(data.table::fread(tb_file, encoding = "UTF-8"))
sil <- collapse_rank(data.table::fread(sil_file, encoding = "UTF-8"))
gene <- merge(tb[, .(gene_symbol, tb_stat = stat, tb_log2FC = log2FoldChange, tb_padj = padj)],
              sil[, .(gene_symbol, sil_stat = stat, sil_log2FC = log2FoldChange, sil_padj = padj)],
              by = "gene_symbol")
data.table::setorder(gene, -tb_stat)
n_common <- nrow(gene)
if (n_common < 10000L) stop("Too few common ranked genes for RRHO")

list_tb <- as.data.frame(gene[, .(gene_symbol, tb_stat)])
list_sil <- as.data.frame(gene[, .(gene_symbol, sil_stat)])
step_size <- RRHO:::defaultStepSize(list_tb, list_sil)
rrho <- RRHO::RRHO(list_tb, list_sil, stepsize = step_size,
                   labels = c("TB-DM", "Silicosis"), alternative = "enrichment",
                   plots = FALSE, BY = FALSE, log10.ind = FALSE)
raw_p <- exp(-rrho$hypermat)
adj_p <- matrix(stats::p.adjust(as.vector(raw_p), method = "BY"),
                nrow = nrow(raw_p), ncol = ncol(raw_p))
log10_by <- -log10(pmax(adj_p, 1e-300))
ranks1 <- seq(1, n_common, by = step_size)[seq_len(nrow(log10_by))]
ranks2 <- seq(1, n_common, by = step_size)[seq_len(ncol(log10_by))]
rrho_dt <- data.table::as.data.table(expand.grid(
  row_index = seq_len(nrow(log10_by)), col_index = seq_len(ncol(log10_by))))
rrho_dt[, by_log10_p := as.vector(log10_by)]
rrho_dt[, `:=`(tb_rank = ranks1[as.integer(row_index)],
               silicosis_rank = ranks2[as.integer(col_index)])]
rrho_dt[, c("row_index", "col_index") := NULL]

half1 <- seq_len(floor(nrow(log10_by) / 2))
half2 <- (floor(nrow(log10_by) / 2) + 1):nrow(log10_by)
up_idx <- which(log10_by[half1, half1, drop = FALSE] ==
                  max(log10_by[half1, half1, drop = FALSE], na.rm = TRUE), arr.ind = TRUE)[1, ]
down_local <- which(log10_by[half2, half2, drop = FALSE] ==
                      max(log10_by[half2, half2, drop = FALSE], na.rm = TRUE), arr.ind = TRUE)[1, ]
down_idx <- c(half2[down_local[1]], half2[down_local[2]])
tb_sorted <- list_tb[order(list_tb[, 2], decreasing = TRUE), ]
sil_sorted <- list_sil[order(list_sil[, 2], decreasing = TRUE), ]
up_genes <- intersect(tb_sorted[seq_len(ranks1[up_idx[1]]), 1],
                      sil_sorted[seq_len(ranks2[up_idx[2]]), 1])
down_genes <- intersect(tb_sorted[ranks1[down_idx[1]]:n_common, 1],
                        sil_sorted[ranks2[down_idx[2]]:n_common, 1])

tb_gsea <- data.table::fread(tb_gsea_file, encoding = "UTF-8")
sil_gsea <- data.table::fread(sil_gsea_file, encoding = "UTF-8")
pathway <- merge(tb_gsea[, .(framework, pathway, pathway_label, tb_NES = NES, tb_padj = padj)],
                 sil_gsea[, .(framework, pathway, sil_NES = NES, sil_padj = padj)],
                 by = c("framework", "pathway"))
pathway[, both_fdr := tb_padj < 0.05 & sil_padj < 0.05]
pathway[, concordant := sign(tb_NES) == sign(sil_NES)]

summary_tab <- data.table::data.table(
  metric = c("Common ranked genes", "Gene-statistic Spearman correlation",
             "Gene-effect direction agreement", "Common pathways",
             "Pathway NES Spearman correlation", "Pathway direction agreement",
             "Pathways FDR<0.05 in both", "Both-FDR pathway direction agreement",
             "RRHO step size", "RRHO shared-up genes at peak", "RRHO shared-down genes at peak"),
  value = c(n_common,
            stats::cor(gene$tb_stat, gene$sil_stat, method = "spearman"),
            mean(sign(gene$tb_log2FC) == sign(gene$sil_log2FC), na.rm = TRUE),
            nrow(pathway), stats::cor(pathway$tb_NES, pathway$sil_NES, method = "spearman"),
            mean(pathway$concordant), sum(pathway$both_fdr),
            mean(pathway[both_fdr == TRUE, concordant]), step_size,
            length(up_genes), length(down_genes))
)

result_dir <- file.path(path_result, "03_cross_disease", "RRHO")
invisible(lapply(file.path(result_dir, c("tables", "source_data", "models", "figures")),
                 dir.create, recursive = TRUE, showWarnings = FALSE))
safe_write_csv(summary_tab, file.path(result_dir, "tables", "T01_cross_disease_concordance_summary.csv"))
safe_write_csv(pathway, file.path(result_dir, "tables", "T02_cross_disease_pathway_concordance.csv"))
safe_write_csv(gene, file.path(result_dir, "source_data", "SD01_cross_disease_gene_statistics.csv"))
safe_write_csv(rrho_dt, file.path(result_dir, "source_data", "SD02_RRHO_BY_map.csv"))
safe_write_csv(data.table::data.table(gene_symbol = up_genes),
               file.path(result_dir, "tables", "T03_RRHO_shared_up_peak_genes.csv"))
safe_write_csv(data.table::data.table(gene_symbol = down_genes),
               file.path(result_dir, "tables", "T04_RRHO_shared_down_peak_genes.csv"))
saveRDS(rrho, file.path(result_dir, "models", "M01_cross_disease_RRHO.rds"))

gene[, both_fdr := tb_padj < 0.05 & sil_padj < 0.05]
p_a <- ggplot2::ggplot(gene, ggplot2::aes(tb_stat, sil_stat, colour = both_fdr)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, colour = "#BEBEBE") +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, colour = "#BEBEBE") +
  ggplot2::geom_point(size = 0.45, alpha = 0.35) +
  ggplot2::scale_colour_manual(values = c(`FALSE` = "#AFAFAF", `TRUE` = "#B84A4A"), guide = "none") +
  ggplot2::labs(title = "Gene-level rank relation", subtitle = sprintf("%s common genes", format(n_common, big.mark = ",")),
                x = "TB-DM Wald statistic", y = "Silicosis Wald statistic")
p_b <- ggplot2::ggplot(rrho_dt, ggplot2::aes(tb_rank, silicosis_rank, fill = by_log10_p)) +
  ggplot2::geom_raster() +
  ggplot2::geom_vline(xintercept = n_common / 2, linetype = 2, linewidth = 0.25, colour = "white") +
  ggplot2::geom_hline(yintercept = n_common / 2, linetype = 2, linewidth = 0.25, colour = "white") +
  ggplot2::annotate("text", x = n_common * 0.16, y = n_common * 0.10,
                    label = "Shared up", size = 2.2, colour = "#272727") +
  ggplot2::annotate("text", x = n_common * 0.84, y = n_common * 0.90,
                    label = "Shared down", size = 2.2, colour = "#272727") +
  ggplot2::scale_y_reverse() +
  ggplot2::scale_fill_gradientn(colours = c("white", "#9AC8C1", "#416A9A", "#272727"),
                                name = expression(-log[10](P[BY]))) +
  ggplot2::labs(title = "Rank-rank hypergeometric overlap", subtitle = "Map-wide Benjamini-Yekutieli adjustment",
                x = "TB-DM rank threshold", y = "Silicosis rank threshold") +
  ggplot2::coord_fixed()
pathway[, evidence := ifelse(both_fdr, "FDR<0.05 in both", "Other")]
lab_path <- pathway[both_fdr == TRUE & concordant == TRUE][order(pmax(tb_padj, sil_padj))][1:min(.N, 6)]
lab_path[, pathway_label_short := ifelse(nchar(pathway_label) > 30,
  paste0(substr(pathway_label, 1, 27), "..."), pathway_label)]
p_c <- ggplot2::ggplot(pathway, ggplot2::aes(tb_NES, sil_NES, colour = evidence)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, colour = "#BEBEBE") +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, colour = "#BEBEBE") +
  ggplot2::geom_point(size = 0.7, alpha = 0.55) +
  ggplot2::geom_text(data = lab_path, ggplot2::aes(label = pathway_label_short),
                     size = 1.7, check_overlap = TRUE, colour = "#272727",
                     hjust = 1, nudge_x = -0.05, vjust = -0.6) +
  ggplot2::scale_colour_manual(values = c("Other" = "#AFAFAF", "FDR<0.05 in both" = "#B84A4A")) +
  ggplot2::labs(title = "Pathway-direction concordance", subtitle = "Hallmark and Reactome",
                x = "TB-DM normalized enrichment score", y = "Silicosis normalized enrichment score",
                colour = NULL)

theme_cross <- ggplot2::theme_classic(base_size = 7, base_family = "Arial") +
  ggplot2::theme(axis.line = ggplot2::element_line(linewidth = 0.3),
    axis.ticks = ggplot2::element_line(linewidth = 0.3), legend.position = "top",
    plot.title = ggplot2::element_text(size = 8, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 6, colour = "#767676"),
    plot.tag = ggplot2::element_text(size = 8, face = "bold"))
fig <- (p_a | p_b | p_c) + patchwork::plot_layout(widths = c(0.9, 1.05, 1.05)) +
  patchwork::plot_annotation(tag_levels = "a", title = "Cross-disease transcriptomic comparison",
    subtitle = "TB-DM versus TB-only and silicosis versus silica-exposed controls") & theme_cross
out <- file.path(result_dir, "figures", "F09_cross_disease_RRHO_pathway_concordance")
w <- 183 / 25.4; h <- 88 / 25.4
svglite::svglite(paste0(out, ".svg"), width = w, height = h); print(fig); grDevices::dev.off()
assert_english_only_svg(paste0(out, ".svg"))
grDevices::cairo_pdf(paste0(out, ".pdf"), width = w, height = h, family = "Arial"); print(fig); grDevices::dev.off()
ragg::agg_tiff(paste0(out, ".tiff"), width = w, height = h, units = "in", res = 600,
               compression = "lzw"); print(fig); grDevices::dev.off()
ragg::agg_png(paste0(out, ".png"), width = w, height = h, units = "in", res = 300);
print(fig); grDevices::dev.off()
write_log("Cross-disease RRHO/pathway concordance completed: common genes=", n_common,
          "; both-FDR pathways=", sum(pathway$both_fdr), "; English-only SVG QA passed")
