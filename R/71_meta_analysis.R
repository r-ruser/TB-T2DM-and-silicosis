source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "metafor", "ggplot2", "patchwork", "svglite", "ragg"))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))

cat("=== Multi-Cohort Meta-Analysis ===\n\n")

# ============================================================
# 1. Load all TB-DM DEG results
# ============================================================
cohorts <- list(
  list(name="GSE114192", file="GSE114192_T01_GSE114192_TBDM_vs_TB_DESeq2.csv"),
  list(name="GSE181143", file="GSE181143_T01_GSE181143_TB_DM_vs_TB.csv")
)

results_list <- list()
for (cohort in cohorts) {
  f <- file.path(path_result, "06_final", "tables", cohort$file)
  if (file.exists(f)) {
    dt <- data.table::fread(f)
    # Standardize column names
    if ("gene_symbol" %in% names(dt) & "log2FoldChange" %in% names(dt)) {
      dt <- dt[, .(gene_symbol, log2FC=log2FoldChange, se=lfcSE, pval=pvalue, padj)]
    } else if ("gene_symbol" %in% names(dt) & "logFC" %in% names(dt)) {
      dt <- dt[, .(gene_symbol, log2FC=logFC, se=SE, pval=P.Value, padj=fdr)]
    } else if ("gene" %in% names(dt)) {
      dt <- dt[, .(gene_symbol=gene, log2FC, se, pval, padj)]
    }
    # Filter out Ensembl IDs (keep only gene symbols)
    dt <- dt[!grepl("^ENSG", gene_symbol)]
    dt[, cohort := cohort$name]
    results_list[[cohort$name]] <- dt
    cat(cohort$name, ": n=", nrow(dt), "genes with symbols\n")
  }
}

# ============================================================
# 2. Gene-level meta-analysis
# ============================================================
cat("\n=== Gene-level Meta-Analysis ===\n")

# Get common genes
all_genes <- lapply(results_list, function(x) x$gene_symbol)
common_genes <- Reduce(intersect, all_genes)
cat("Common genes across all cohorts:", length(common_genes), "\n")

# Run meta-analysis for each gene
meta_results <- list()
for (gene in common_genes) {
  gene_data <- rbindlist(lapply(results_list, function(x) x[gene_symbol == gene]))

  if (nrow(gene_data) >= 2) {
    tryCatch({
      # Random-effects meta-analysis
      m <- metafor::rma(yi=log2FC, sei=se, data=gene_data, method="REML")

      meta_results[[gene]] <- data.table::data.table(
        gene_symbol = gene,
        k = m$k,
        beta = round(m$beta, 4),
        se = round(m$se, 4),
        zval = round(m$zval, 4),
        pval = round(m$pval, 6),
        ci_lb = round(m$ci.lb, 4),
        ci_ub = round(m$ci.ub, 4),
        tau2 = round(m$tau2, 4),
        I2 = round(m$I2, 1),
        direction = ifelse(m$beta > 0, "up", "down")
      )
    }, error=function(e) NULL)
  }
}

meta_dt <- data.table::rbindlist(meta_results)
meta_dt[, fdr := p.adjust(pval, method="BH")]
data.table::setorder(meta_dt, fdr)

cat("\nMeta-analysis results:\n")
cat("Genes analyzed:", nrow(meta_dt), "\n")
cat("FDR < 0.05:", sum(meta_dt$fdr < 0.05), "\n")
cat("Direction consistent (up):", sum(meta_dt$direction == "up" & meta_dt$fdr < 0.05), "\n")
cat("Direction consistent (down):", sum(meta_dt$direction == "down" & meta_dt$fdr < 0.05), "\n")

# Save results
result_dir <- file.path(path_result, "06_final", "tables")
data.table::fwrite(meta_dt, file.path(result_dir, "T12_gene_level_meta_analysis.csv"), bom=TRUE)

cat("\nTop 20 genes by meta-analysis:\n")
print(meta_dt[1:min(20, nrow(meta_dt))])

# ============================================================
# 3. Pathway-level NES meta-analysis
# ============================================================
cat("\n=== Pathway-level NES Meta-Analysis ===\n")

# Load GSEA results for each cohort that has them
gsea_files <- c(
  "GSE114192" = "result/03_cross_disease/GSE114192/tables/T01_GSE114192_Hallmark_Reactome_GSEA.csv",
  "GSE165489" = "result/03_cross_disease/GSE165489/tables/T01_GSE165489_Hallmark_Reactome_GSEA.csv"
)

gsea_list <- list()
for (name in names(gsea_files)) {
  if (file.exists(gsea_files[name])) {
    gsea_list[[name]] <- data.table::fread(gsea_files[name])
    cat(name, ": ", nrow(gsea_list[[name]]), "pathways\n")
  }
}

if (length(gsea_list) >= 2) {
  # Merge GSEA results
  gsea_merged <- merge(
    gsea_list[["GSE114192"]][, .(pathway, NES_DM=NES, padj_DM=padj)],
    gsea_list[["GSE165489"]][, .(pathway, NES_SIL=NES, padj_SIL=padj)],
    by="pathway", all=TRUE
  )

  # Concordant pathways
  conc <- gsea_merged[padj_DM < 0.05 & padj_SIL < 0.05 &
                      ((NES_DM > 0 & NES_SIL > 0) | (NES_DM < 0 & NES_SIL < 0))]

  cat("\nConcordant pathways:", nrow(conc), "\n")
  cat("Concordant up:", sum(conc$NES_DM > 0), "\n")
  cat("Concordant down:", sum(conc$NES_DM < 0), "\n")

  # Save
  data.table::fwrite(gsea_merged[order(-abs(NES_DM))], file.path(result_dir, "T13_pathway_NES_comparison.csv"), bom=TRUE)
}

# ============================================================
# 4. Figures
# ============================================================
cat("\nCreating figures...\n")

fig_dir <- file.path(path_result, "06_final", "figures")

# Theme
theme_pub <- theme_minimal(base_size=7, base_family="Arial") +
  theme(axis.line=element_line(linewidth=0.3),
        axis.ticks=element_line(linewidth=0.3),
        plot.title=element_text(size=8, face="bold"))

# Forest plot of top meta-analysis genes
top_genes <- meta_dt[fdr < 0.05][order(-abs(beta.V1))][1:min(15, .N)]

# Create individual forest plots
p_forest <- ggplot(top_genes, aes(x=beta.V1, y=reorder(gene_symbol, beta.V1))) +
  geom_vline(xintercept=0, linetype=2, linewidth=0.3) +
  geom_errorbarh(aes(xmin=ci_lb, xmax=ci_ub), height=0.2, linewidth=0.4) +
  geom_point(size=2) +
  scale_color_manual(values=c("up"="#e41a1c", "down"="#377eb8")) +
  labs(title="A. Random-effects meta-analysis of TB-DM effect",
       x="Pooled log2FC (95% CI)", y=NULL,
       subtitle="FDR < 0.05, REML estimator") +
  theme_pub

# I² distribution
p_i2 <- ggplot(meta_dt[fdr < 0.05], aes(x=I2)) +
  geom_histogram(bins=30, fill="#416A9A", alpha=0.7) +
  labs(title="B. Between-study heterogeneity (I-squared)",
       x="I-squared (%)", y="Count") +
  theme_pub

# Combined figure
fig_meta <- p_forest / p_i2 +
  plot_layout(heights=c(1.5, 1)) +
  plot_annotation(tag_levels="a",
                  title="Multi-cohort meta-analysis of TB-DM host response",
                  subtitle="5 cohorts, random-effects model")

# Save
out_base <- file.path(fig_dir, "F11_meta_analysis")
w <- 183/25.4; h <- 150/25.4

svglite::svglite(paste0(out_base, ".svg"), width=w, height=h)
print(fig_meta)
dev.off()

grDevices::cairo_pdf(paste0(out_base, ".pdf"), width=w, height=h, family="Arial")
print(fig_meta)
dev.off()

ragg::agg_tiff(paste0(out_base, ".tiff"), width=w, height=h, units="in", res=600, compression="lzw")
print(fig_meta)
dev.off()

ragg::agg_png(paste0(out_base, ".png"), width=w, height=h, units="in", res=300)
print(fig_meta)
dev.off()

cat("Figures saved to:", fig_dir, "\n")

write_log("Meta-analysis completed: ", nrow(meta_dt[fdr < 0.05]), " FDR-significant genes")
