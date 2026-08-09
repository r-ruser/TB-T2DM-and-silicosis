source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "ggplot2", "patchwork", "svglite", "ragg"))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))

cat("=== Pathway-Level Convergence Analysis ===\n\n")

# ============================================================
# 1. Load GSEA results
# ============================================================
gsea_dm <- data.table::fread(file.path(path_result, "03_cross_disease", "GSE114192", "tables", "T01_GSE114192_Hallmark_Reactome_GSEA.csv"))
gsea_sil <- data.table::fread(file.path(path_result, "03_cross_disease", "GSE165489", "tables", "T01_GSE165489_Hallmark_Reactome_GSEA.csv"))

cat("GSE114192 (TB-DM):", nrow(gsea_dm), "pathways,", sum(gsea_dm$padj < 0.05, na.rm=TRUE), "significant\n")
cat("GSE165489 (Silicosis):", nrow(gsea_sil), "pathways,", sum(gsea_sil$padj < 0.05, na.rm=TRUE), "significant\n")

# ============================================================
# 2. Pathway-level convergence
# ============================================================
merged <- merge(
  gsea_dm[, .(pathway, NES_DM=NES, padj_DM=padj, ES_DM=ES, framework, leadingEdge_DM=leadingEdge_text)],
  gsea_sil[, .(pathway, NES_SIL=NES, padj_SIL=padj, ES_SIL=ES, leadingEdge_SIL=leadingEdge_text)],
  by="pathway", all=TRUE
)

# Classify pathways
merged[, significance := "NS"]
merged[padj_DM < 0.05 & padj_SIL < 0.05, significance := "Both significant"]
merged[padj_DM < 0.05 & padj_SIL >= 0.05, significance := "TB-DM only"]
merged[padj_DM >= 0.05 & padj_SIL < 0.05, significance := "Silicosis only"]

# Concordance
merged[, direction := "NS"]
merged[significance == "Both significant" & NES_DM > 0 & NES_SIL > 0, direction := "Both up"]
merged[significance == "Both significant" & NES_DM < 0 & NES_SIL < 0, direction := "Both down"]
merged[significance == "Both significant" & ((NES_DM > 0 & NES_SIL < 0) | (NES_DM < 0 & NES_SIL > 0)), direction := "Opposite"]

# Concordant pathways
conc_up <- merged[direction == "Both up"]
conc_down <- merged[direction == "Both down"]

cat("\n=== Pathway Convergence Summary ===\n")
cat("Total shared pathways:", nrow(merged), "\n")
cat("Both significant:", sum(merged$significance == "Both significant"), "\n")
cat("  Both up:", nrow(conc_up), "\n")
cat("  Both down:", nrow(conc_down), "\n")
cat("  Opposite:", sum(merged$direction == "Opposite"), "\n")

# Spearman correlation
cor_test <- cor.test(merged$NES_DM, merged$NES_SIL, method="spearman", use="complete.obs")
cat("\nNES Spearman correlation:", round(cor_test$estimate, 3), "\n")

# ============================================================
# 3. Save results
# ============================================================
result_dir <- file.path(path_result, "06_final", "tables")

# Save full convergence table
data.table::fwrite(merged[order(-abs(NES_DM))], file.path(result_dir, "T10_pathway_convergence_full.csv"), bom=TRUE)

# Save concordant pathways
conc_pathways <- merged[direction %in% c("Both up", "Both down")][order(-abs(NES_DM))]
data.table::fwrite(conc_pathways, file.path(result_dir, "T11_pathway_convergence_concordant.csv"), bom=TRUE)

cat("\nTop 20 concordant pathways:\n")
print(conc_pathways[1:min(20, nrow(conc_pathways)), .(pathway, NES_DM, padj_DM, NES_SIL, padj_SIL, direction)])

# ============================================================
# 4. Figures
# ============================================================
cat("\nCreating figures...\n")

fig_dir <- file.path(path_result, "06_final", "figures")

# Theme
theme_pub <- theme_minimal(base_size=7, base_family="Arial") +
  theme(
    axis.line=element_line(linewidth=0.3),
    axis.ticks=element_line(linewidth=0.3),
    plot.title=element_text(size=8, face="bold"),
    legend.position="top"
  )

# NES Scatter plot
p_scatter <- ggplot(merged[!is.na(NES_DM) & !is.na(NES_SIL)],
                    aes(x=NES_DM, y=NES_SIL, color=significance)) +
  geom_point(alpha=0.6, size=1.5) +
  geom_hline(yintercept=0, linetype=2, linewidth=0.3) +
  geom_vline(xintercept=0, linetype=2, linewidth=0.3) +
  geom_smooth(method="lm", se=TRUE, color="black", linewidth=0.5) +
  scale_color_manual(values=c("Both significant"="#e41a1c",
                               "TB-DM only"="#377eb8",
                               "Silicosis only"="#4daf4a",
                               "NS"="#999999")) +
  labs(title="A. NES concordance between TB-DM and Silicosis",
       x="NES (TB-DM vs TB-only)",
       y="NES (Silicosis vs exposed-no-silicosis)",
       color="Significance") +
  theme_pub

# Heatmap of top concordant pathways
top_conc <- conc_pathways[1:min(25, nrow(conc_pathways))]
top_conc[, pathway_display := factor(pathway, levels=rev(pathway))]

p_heat <- ggplot(top_conc, aes(x="TB-DM", y=pathway_display, fill=NES_DM)) +
  geom_tile() +
  geom_tile(aes(x="Silicosis", y=pathway_display, fill=NES_SIL)) +
  scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b", midpoint=0, name="NES") +
  labs(title="B. Top concordant pathways", x=NULL, y=NULL) +
  theme_pub + theme(axis.text.x=element_text(angle=45, hjust=1))

# Combined figure
fig <- p_scatter / p_heat +
  plot_layout(heights=c(1, 1.2)) +
  plot_annotation(tag_levels="a",
                  title="Pathway-level convergence between TB-DM and Silicosis",
                  subtitle="Hallmark and Reactome GSEA analysis")

# Save
out_base <- file.path(fig_dir, "F10_pathway_convergence")
w <- 183/25.4; h <- 180/25.4

svglite::svglite(paste0(out_base, ".svg"), width=w, height=h)
print(fig)
dev.off()

grDevices::cairo_pdf(paste0(out_base, ".pdf"), width=w, height=h, family="Arial")
print(fig)
dev.off()

ragg::agg_tiff(paste0(out_base, ".tiff"), width=w, height=h, units="in", res=600, compression="lzw")
print(fig)
dev.off()

ragg::agg_png(paste0(out_base, ".png"), width=w, height=h, units="in", res=300)
print(fig)
dev.off()

cat("Figures saved to:", fig_dir, "\n")

write_log("Pathway convergence analysis completed: ", nrow(conc_pathways), " concordant pathways")
