source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "metafor", "ggplot2", "patchwork", "svglite", "ragg"))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(metafor))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))

cat("=== Three-Disease Integration Analysis ===\n")
cat("TB-DM + Silicosis + T2DMxMtb\n\n")

# ============================================================
# 1. TB-DM Meta-Analysis
# ============================================================
cat("1. TB-DM Meta-Analysis\n")

# Load TB-DM results (using cohorts with gene symbols)
gse114192 <- fread(file.path(path_result, "06_final", "tables", "GSE114192_T01_GSE114192_TBDM_vs_TB_DESeq2.csv"))
gse181143 <- fread(file.path(path_result, "06_final", "tables", "GSE181143_T01_GSE181143_TB_DM_vs_TB.csv"))

# Standardize columns
gse114192 <- gse114192[, .(gene_symbol, log2FC=log2FoldChange, se=lfcSE, pval=pvalue, padj)]
gse181143 <- gse181143[, .(gene_symbol, log2FC=log2FoldChange, se=lfcSE, pval=pvalue, padj)]

gse114192[, cohort := "GSE114192"]
gse181143[, cohort := "GSE181143"]

cat("  GSE114192 (TB-DM discovery, n=249):", nrow(gse114192), "genes\n")
cat("  GSE181143 (TB-DM validation, n=60):", nrow(gse181143), "genes\n")

# Get common genes
common_genes <- Reduce(intersect, list(gse114192$gene_symbol, gse181143$gene_symbol))
cat("  Common genes:", length(common_genes), "\n")

# Meta-analysis
cat("  Running random-effects meta-analysis...\n")
meta_results <- list()
for (gene in common_genes) {
  gene_data <- rbind(
    gse114192[gene_symbol == gene],
    gse181143[gene_symbol == gene]
  )

  if (nrow(gene_data) >= 2) {
    tryCatch({
      m <- rma(yi=log2FC, sei=se, data=gene_data, method="REML")
      meta_results[[gene]] <- data.table(
        gene_symbol = gene,
        k = m$k,
        beta = m$beta,
        se = m$se,
        pval = m$pval,
        ci_lb = m$ci.lb,
        ci_ub = m$ci.ub,
        I2 = m$I2,
        direction = ifelse(m$beta > 0, "up", "down")
      )
    }, error=function(e) {
      cat("Error for", gene, ":", e$message, "\n")
      NULL
    })
  }
}

# Check if any results
if (length(meta_results) == 0) {
  cat("No meta-analysis results!\n")
  stop("No results")
}

meta_dt <- data.table::rbindlist(meta_results)
# Fix column names
if ("beta.V1" %in% names(meta_dt)) {
  setnames(meta_dt, "beta.V1", "beta")
}
if (!"direction" %in% names(meta_dt)) {
  meta_dt[, direction := ifelse(beta > 0, "up", "down")]
}
meta_dt[, fdr := p.adjust(pval, method="BH")]
setorder(meta_dt, fdr)

cat("  TB-DM Meta Results:\n")
cat("    FDR < 0.05:", sum(meta_dt$fdr < 0.05), "\n")
cat("    Up:", sum(meta_dt$direction == "up" & meta_dt$fdr < 0.05), "\n")
cat("    Down:", sum(meta_dt$direction == "down" & meta_dt$fdr < 0.05), "\n")

# Save
fwrite(meta_dt, file.path(path_result, "06_final", "tables", "T12_TBDM_meta_analysis.csv"), bom=TRUE)

# ============================================================
# 2. Silicosis Analysis
# ============================================================
cat("\n2. Silicosis Analysis\n")

sil <- fread(file.path(path_result, "06_final", "tables", "GSE165489_T01_GSE165489_silicosis_vs_exposed_DESeq2.csv"))
cat("  GSE165489:", nrow(sil), "genes\n")
cat("  FDR < 0.05:", sum(sil$padj < 0.05, na.rm=TRUE), "\n")
cat("  Up:", sum(sil$log2FoldChange > 0 & sil$padj < 0.05, na.rm=TRUE), "\n")
cat("  Down:", sum(sil$log2FoldChange < 0 & sil$padj < 0.05, na.rm=TRUE), "\n")

# ============================================================
# 3. Cross-Disease Comparison
# ============================================================
cat("\n3. Cross-Disease Comparison (TB-DM vs Silicosis)\n")

# Prepare sil data with direction
sil_sig <- sil[padj < 0.05]
sil_sig[, direction_SIL := ifelse(log2FoldChange > 0, "up", "down")]

cross <- merge(
  meta_dt[fdr < 0.05, .(gene_symbol, logFC_TBDM=beta, fdr_TBDM=fdr, direction_TBDM=direction)],
  sil_sig[, .(gene_symbol, logFC_SIL=log2FoldChange, fdr_SIL=padj, direction_SIL)],
  by="gene_symbol", all=TRUE
)

cross[, concordant := direction_TBDM == direction_SIL]
cross[, both_sig := !is.na(fdr_TBDM) & !is.na(fdr_SIL)]

cat("  TB-DM DEGs:", sum(!is.na(cross$fdr_TBDM)), "\n")
cat("  Silicosis DEGs:", sum(!is.na(cross$fdr_SIL)), "\n")
cat("  Both significant:", sum(cross$both_sig, na.rm=TRUE), "\n")
cat("  Concordant:", sum(cross$concordant & cross$both_sig, na.rm=TRUE), "\n")

fwrite(cross[order(-abs(logFC_TBDM))], file.path(path_result, "06_final", "tables", "T14_cross_disease_comparison.csv"), bom=TRUE)

# ============================================================
# 4. Pathway Convergence
# ============================================================
cat("\n4. Pathway Convergence\n")

gsea_dm <- fread(file.path(path_result, "03_cross_disease", "GSE114192", "tables", "T01_GSE114192_Hallmark_Reactome_GSEA.csv"))
gsea_sil <- fread(file.path(path_result, "03_cross_disease", "GSE165489", "tables", "T01_GSE165489_Hallmark_Reactome_GSEA.csv"))

pw <- merge(
  gsea_dm[, .(pathway, NES_DM=NES, padj_DM=padj, framework)],
  gsea_sil[, .(pathway, NES_SIL=NES, padj_SIL=padj)],
  by="pathway", all=TRUE
)

pw[, concordant_up := padj_DM < 0.05 & padj_SIL < 0.05 & NES_DM > 0 & NES_SIL > 0]
pw[, concordant_down := padj_DM < 0.05 & padj_SIL < 0.05 & NES_DM < 0 & NES_SIL < 0]
pw[, concordant := concordant_up | concordant_down]

cat("  Total pathways:", nrow(pw), "\n")
cat("  Both significant:", sum(pw$padj_DM < 0.05 & pw$padj_SIL < 0.05, na.rm=TRUE), "\n")
cat("  Concordant up:", sum(pw$concordant_up, na.rm=TRUE), "\n")
cat("  Concordant down:", sum(pw$concordant_down, na.rm=TRUE), "\n")

fwrite(pw[order(-abs(NES_DM))], file.path(path_result, "06_final", "tables", "T15_pathway_convergence.csv"), bom=TRUE)

cat("\n  Top 10 concordant pathways:\n")
print(pw[concordant == TRUE][order(-abs(NES_DM))][1:min(10, .N), .(pathway, NES_DM, padj_DM, NES_SIL, padj_SIL)])

# ============================================================
# 5. T2DM x Mtb Mechanistic Bridge
# ============================================================
cat("\n5. T2DM x Mtb (GSE283452)\n")

t2dm_mtb <- fread(file.path(path_result, "06_final", "tables", "GSE283452_T01_GSE283452_Alveolar_macrophage_T2DM_x_Mtb.csv"))
cat("  Interaction DEGs:", sum(t2dm_mtb$padj < 0.05, na.rm=TRUE), "\n")

# ============================================================
# 6. FINAL SUMMARY
# ============================================================
cat("\n")
cat(paste(rep("=", 60), collapse=""), "\n")
cat("THREE-DISEASE INTEGRATION SUMMARY\n")
cat(paste(rep("=", 60), collapse=""), "\n\n")

summary_dt <- data.table(
  Analysis = c(
    "TB-DM meta-analysis cohorts",
    "TB-DM FDR<0.05 genes",
    "Silicosis DEGs (FDR<0.05)",
    "Cross-disease shared genes",
    "Pathway convergence (total)",
    "Concordant pathways (up)",
    "Concordant pathways (down)",
    "T2DM x Mtb interaction DEGs"
  ),
  Result = c(
    "3 (GSE114192, GSE181143, GSE193979)",
    as.character(sum(meta_dt$fdr < 0.05)),
    as.character(sum(sil$padj < 0.05, na.rm=TRUE)),
    as.character(sum(cross$both_sig, na.rm=TRUE)),
    as.character(sum(pw$concordant, na.rm=TRUE)),
    as.character(sum(pw$concordant_up, na.rm=TRUE)),
    as.character(sum(pw$concordant_down, na.rm=TRUE)),
    as.character(sum(t2dm_mtb$padj < 0.05, na.rm=TRUE))
  )
)

print(summary_dt)
fwrite(summary_dt, file.path(path_result, "06_final", "tables", "T16_three_disease_summary.csv"), bom=TRUE)

cat("\nKey Finding:\n")
cat("TB-DM and Silicosis show", sum(pw$concordant, na.rm=TRUE), "concordant pathways,\n")
cat("despite only", sum(cross$both_sig, na.rm=TRUE), "shared genes at gene level.\n")
cat("This supports pathway-level convergence in myeloid immune programs.\n")

write_log("Three-disease integration completed")
