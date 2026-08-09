source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "fgsea", "msigdbr", "ggplot2", "svglite", "ragg"))

acc <- "GSE165489"
de_file <- file.path(path_result, "02_bulk", acc, "tables",
                     "T01_GSE165489_silicosis_vs_exposed_DESeq2.csv")
if (!file.exists(de_file)) stop("Run R/31_gse165489_primary_deseq2.R first")
de <- data.table::fread(de_file, encoding = "UTF-8")
rank_dt <- de[!is.na(gene_symbol) & gene_symbol != "" & is.finite(stat)]
rank_dt[, abs_stat := abs(stat)]
data.table::setorder(rank_dt, gene_symbol, -abs_stat)
rank_dt <- rank_dt[, .SD[1], by = gene_symbol]
ranks <- sort(stats::setNames(rank_dt$stat, rank_dt$gene_symbol), decreasing = TRUE)

get_sets <- function(collection, subcollection = NULL) {
  legacy <- "category" %in% names(formals(msigdbr::msigdbr))
  x <- if (legacy && is.null(subcollection)) {
    msigdbr::msigdbr(species = "Homo sapiens", category = collection)
  } else if (legacy) {
    msigdbr::msigdbr(species = "Homo sapiens", category = collection, subcategory = subcollection)
  } else if (is.null(subcollection)) {
    msigdbr::msigdbr(species = "Homo sapiens", collection = collection)
  } else {
    msigdbr::msigdbr(species = "Homo sapiens", collection = collection,
                     subcollection = subcollection)
  }
  split(x$gene_symbol, x$gs_name)
}
sets <- list(Hallmark = get_sets("H"), Reactome = get_sets("C2", "CP:REACTOME"))
run_fgsea <- function(pathways, framework) {
  set.seed(20260808)
  z <- fgsea::fgseaMultilevel(pathways = pathways, stats = ranks, minSize = 15,
                              maxSize = 500, eps = 0, scoreType = "std")
  z <- data.table::as.data.table(z)
  z[, `:=`(framework = framework,
           leadingEdge_text = vapply(leadingEdge, paste, collapse = ";", character(1)))]
  z[, leadingEdge := NULL]
  z
}
gsea <- data.table::rbindlist(Map(run_fgsea, sets, names(sets)), fill = TRUE)
gsea[, pathway_label := gsub("^(HALLMARK_|REACTOME_)", "", pathway)]
gsea[, pathway_label := gsub("_", " ", pathway_label)]
gsea[, abs_NES := abs(NES)]
data.table::setorder(gsea, framework, padj, -abs_NES)

result_dir <- file.path(path_result, "03_cross_disease", acc)
invisible(lapply(file.path(result_dir, c("tables", "source_data", "figures")),
                 dir.create, recursive = TRUE, showWarnings = FALSE))
safe_write_csv(gsea, file.path(result_dir, "tables", "T01_GSE165489_Hallmark_Reactome_GSEA.csv"))
top <- gsea[order(padj, -abs_NES), head(.SD, 10), by = framework]
top[, label := factor(pathway_label, levels = rev(unique(pathway_label)))]
top[, direction := factor(ifelse(NES > 0, "Higher in silicosis", "Lower in silicosis"),
                          levels = c("Lower in silicosis", "Higher in silicosis"))]
safe_write_csv(top, file.path(result_dir, "source_data", "SD01_GSE165489_top_GSEA.csv"))

p <- ggplot2::ggplot(top, ggplot2::aes(NES, label, size = -log10(pmax(padj, 1e-300)),
                                        colour = direction)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "#767676") +
  ggplot2::geom_point(alpha = 0.88) +
  ggplot2::facet_wrap(~ framework, scales = "free_y", ncol = 2) +
  ggplot2::scale_colour_manual(values = c("Lower in silicosis" = "#416A9A",
                                          "Higher in silicosis" = "#B84A4A")) +
  ggplot2::scale_size_continuous(range = c(1.5, 4.5), name = expression(-log[10](FDR))) +
  ggplot2::labs(title = "Pathway-level host-response shifts in silicosis",
                subtitle = "GSE165489; exposure-duration adjusted DESeq2 Wald statistics",
                x = "Normalized enrichment score", y = NULL, colour = NULL) +
  ggplot2::theme_classic(base_size = 7, base_family = "Arial") +
  ggplot2::theme(axis.line = ggplot2::element_line(linewidth = 0.3),
    axis.ticks = ggplot2::element_line(linewidth = 0.3), legend.position = "top",
    strip.text = ggplot2::element_text(size = 7, face = "bold"),
    plot.title = ggplot2::element_text(size = 8, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 6, colour = "#767676"))
out <- file.path(result_dir, "figures", "F08_GSE165489_Hallmark_Reactome_GSEA")
w <- 183 / 25.4; h <- 110 / 25.4
svglite::svglite(paste0(out, ".svg"), width = w, height = h); print(p); grDevices::dev.off()
assert_english_only_svg(paste0(out, ".svg"))
grDevices::cairo_pdf(paste0(out, ".pdf"), width = w, height = h, family = "Arial"); print(p); grDevices::dev.off()
ragg::agg_tiff(paste0(out, ".tiff"), width = w, height = h, units = "in", res = 600,
               compression = "lzw"); print(p); grDevices::dev.off()
ragg::agg_png(paste0(out, ".png"), width = w, height = h, units = "in", res = 300);
print(p); grDevices::dev.off()
write_log("GSE165489 Hallmark/Reactome GSEA completed: ", sum(gsea$padj < 0.05, na.rm = TRUE),
          " FDR-significant pathways; English-only SVG QA passed")
