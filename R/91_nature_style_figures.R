source(file.path("R", "00_config.R"), encoding = "UTF-8")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(ragg))

cat("=== Nature-Style Figure Generation ===\n\n")

# ============================================================
# Nature figure style theme
# ============================================================
theme_nature <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      # Axis
      axis.line = element_line(linewidth = 0.3, color = "black"),
      axis.ticks = element_line(linewidth = 0.3, color = "black"),
      axis.ticks.length = unit(1.5, "mm"),
      axis.text = element_text(size = 6, color = "black"),
      axis.title = element_text(size = 7, color = "black"),
      # Legend
      legend.position = "right",
      legend.key.size = unit(2, "mm"),
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 7, face = "bold"),
      legend.box.margin = margin(0, 0, 0, 0),
      # Plot
      plot.title = element_text(size = 8, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 6, color = "#666666", hjust = 0),
      plot.tag = element_text(size = 8, face = "bold"),
      plot.tag.position = "topleft",
      plot.margin = margin(3, 3, 3, 3),
      # Panel
      strip.background = element_blank(),
      strip.text = element_text(size = 7, face = "bold")
    )
}

# Color palettes
colors_disease <- c(
  "TB" = "#D62728",
  "T2DM" = "#FF7F0E",
  "Silicosis" = "#1F77B4",
  "TB-DM" = "#9467BD",
  "Control" = "#2CA02C",
  "Exposed" = "#17BECF"
)

colors_celltype <- c(
  "Monocyte" = "#E41A1C",
  "Macrophage" = "#377EB8",
  "Neutrophil" = "#4DAF4A",
  "CD4 T" = "#984EA3",
  "CD8 T" = "#FF7F00",
  "NK" = "#A65628",
  "B" = "#F781BF",
  "Mast" = "#999999",
  "Epithelial" = "#66C2A5",
  "Other" = "#BBBBBB"
)

fig_dir <- file.path(path_result, "06_final", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Figure 1: GBD Global Co-burden
# ============================================================
cat("1. Generating Figure 1: GBD Global Co-burden...\n")

# Load GBD data
gbd_cbi <- fread(file.path(path_result, "06_final", "tables", "GBD_T01_GBD2023_prevalence_CBI.csv"))

# Panel A-C: Disease prevalence maps (already generated)
# Panel D: Combined burden index bar chart

top20 <- gbd_cbi[order(-cbi)][1:20]
top20[, location_name := factor(location_name, levels = rev(location_name))]

p_bar <- ggplot(top20, aes(x = cbi, y = location_name)) +
  geom_col(aes(fill = burden_class), width = 0.7) +
  scale_fill_manual(values = c(
    "TB + DM + Silicosis high" = "#E41A1C",
    "TB + DM high" = "#FF7F00",
    "TB + Silicosis high" = "#984EA3",
    "DM + Silicosis high" = "#A65628",
    "TB high only" = "#377EB8",
    "DM high only" = "#4DAF4A",
    "Silicosis high only" = "#999999",
    "Low-Low-Low" = "#F0F0F0"
  ), name = "Burden class") +
  labs(x = "Combined Burden Index", y = NULL,
       title = "Top 20 locations by co-burden index",
       subtitle = "GBD 2023, age-standardized rates") +
  theme_nature() +
  theme(legend.position = "right")

# Panel E: Sensitivity analysis
sens <- fread(file.path(path_result, "06_final", "tables", "GBD_T03_threshold_sensitivity_prevalence.csv"))

p_sens <- ggplot(sens, aes(x = factor(threshold), y = n_countries, group = 1)) +
  geom_line(linewidth = 0.5, color = "#416A9A") +
  geom_point(size = 2, shape = 21, fill = "#B84A4A", color = "black") +
  labs(x = "Percentile threshold", y = "Number of locations",
       title = "Sensitivity analysis",
       subtitle = "Triple-high burden count by threshold") +
  theme_nature()

# Save Figure 1
fig1 <- p_bar | p_sens +
  plot_layout(widths = c(1.2, 1)) +
  plot_annotation(tag_levels = "a")

w <- 183/25.4; h <- 90/25.4
out <- file.path(fig_dir, "Fig1_GBD_global_coburden")

svglite::svglite(paste0(out, ".svg"), width = w, height = h); print(fig1); dev.off()
grDevices::cairo_pdf(paste0(out, ".pdf"), width = w, height = h, family = "Arial"); print(fig1); dev.off()
ragg::agg_tiff(paste0(out, ".tiff"), width = w, height = h, units = "in", res = 600, compression = "lzw"); print(fig1); dev.off()
ragg::agg_png(paste0(out, ".png"), width = w, height = h, units = "in", res = 300); print(fig1); dev.off()
cat("  Fig1 saved\n")

# ============================================================
# Figure 2: Cross-disease Transcriptomic Convergence
# ============================================================
cat("\n2. Generating Figure 2: Pathway Convergence...\n")

# Load pathway convergence data
pw <- fread(file.path(path_result, "06_final", "tables", "T15_pathway_convergence.csv"))

# Panel A: NES scatter
pw_plot <- pw[!is.na(NES_DM) & !is.na(NES_SIL)]
pw_plot[, significance := "NS"]
pw_plot[padj_DM < 0.05 & padj_SIL < 0.05, significance := "Both significant"]
pw_plot[padj_DM < 0.05 & padj_SIL >= 0.05, significance := "TB-DM only"]
pw_plot[padj_DM >= 0.05 & padj_SIL < 0.05, significance := "Silicosis only"]

p_nes <- ggplot(pw_plot, aes(x = NES_DM, y = NES_SIL, color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3, color = "grey50") +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3, color = "grey50") +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.5) +
  scale_color_manual(values = c(
    "Both significant" = "#E41A1C",
    "TB-DM only" = "#377EB8",
    "Silicosis only" = "#4DAF4A",
    "NS" = "#999999"
  ), name = "Significance") +
  labs(x = "NES (TB-DM vs TB-only)", y = "NES (Silicosis vs exposed)",
       title = "Pathway NES concordance",
       subtitle = "Hallmark + Reactome GSEA") +
  theme_nature()

# Panel B: Top concordant pathways bar
conc <- pw[concordant == TRUE][order(-abs(NES_DM))][1:min(15, .N)]
conc[, pathway_short := substr(pathway, 1, 40)]
conc[, pathway_short := factor(pathway_short, levels = rev(pathway_short))]

p_conc <- ggplot(conc, aes(x = NES_DM, y = pathway_short)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3, color = "grey50") +
  geom_segment(aes(xend = 0, yend = pathway_short), linewidth = 0.4, color = "#666666") +
  geom_point(size = 2, color = "#E41A1C") +
  labs(x = "NES (TB-DM)", y = NULL,
       title = "Top concordant pathways",
       subtitle = paste("n =", nrow(conc), "pathways, both FDR < 0.05")) +
  theme_nature() +
  theme(axis.text.y = element_text(size = 5))

# Save Figure 2
fig2 <- p_nes | p_conc +
  plot_layout(widths = c(1, 1.2)) +
  plot_annotation(tag_levels = "a")

out <- file.path(fig_dir, "Fig2_pathway_convergence")
svglite::svglite(paste0(out, ".svg"), width = w, height = h); print(fig2); dev.off()
grDevices::cairo_pdf(paste0(out, ".pdf"), width = w, height = h, family = "Arial"); print(fig2); dev.off()
ragg::agg_tiff(paste0(out, ".tiff"), width = w, height = h, units = "in", res = 600, compression = "lzw"); print(fig2); dev.off()
ragg::agg_png(paste0(out, ".png"), width = w, height = h, units = "in", res = 300); print(fig2); dev.off()
cat("  Fig2 saved\n")

# ============================================================
# Figure 3: TB-DM Meta-analysis
# ============================================================
cat("\n3. Generating Figure 3: Meta-analysis...\n")

meta <- fread(file.path(path_result, "06_final", "tables", "T12_TBDM_meta_analysis.csv"))

# Panel A: Forest plot of top genes
top_genes <- meta[fdr < 0.05][order(-abs(beta))][1:min(15, .N)]

p_forest <- ggplot(top_genes, aes(x = beta, y = reorder(gene_symbol, beta))) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3, color = "grey50") +
  geom_errorbarh(aes(xmin = ci_lb, xmax = ci_ub), height = 0.2, linewidth = 0.4, color = "#666666") +
  geom_point(size = 2, color = "#E41A1C") +
  labs(x = "Pooled log2FC (95% CI)", y = NULL,
       title = "Top reproducible TB-DM genes",
       subtitle = "Random-effects meta-analysis, FDR < 0.05") +
  theme_nature()

# Panel B: I² distribution
p_i2 <- ggplot(meta[fdr < 0.05], aes(x = I2)) +
  geom_histogram(bins = 25, fill = "#416A9A", alpha = 0.7, color = "white", linewidth = 0.2) +
  geom_vline(xintercept = 50, linetype = 2, linewidth = 0.3, color = "grey50") +
  labs(x = "I² heterogeneity (%)", y = "Count",
       title = "Between-study heterogeneity",
       subtitle = "Among FDR-significant genes") +
  theme_nature()

# Save Figure 3
fig3 <- p_forest | p_i2 +
  plot_layout(widths = c(1.2, 1)) +
  plot_annotation(tag_levels = "a")

out <- file.path(fig_dir, "Fig3_TBDM_meta_analysis")
svglite::svglite(paste0(out, ".svg"), width = w, height = h); print(fig3); dev.off()
grDevices::cairo_pdf(paste0(out, ".pdf"), width = w, height = h, family = "Arial"); print(fig3); dev.off()
ragg::agg_tiff(paste0(out, ".tiff"), width = w, height = h, units = "in", res = 600, compression = "lzw"); print(fig3); dev.off()
ragg::agg_png(paste0(out, ".png"), width = w, height = h, units = "in", res = 300); print(fig3); dev.off()
cat("  Fig3 saved\n")

# ============================================================
# Figure 4: Single-cell Localization
# ============================================================
cat("\n4. Generating Figure 4: Single-cell localization...\n")

# Load cell type scores
sc_summary <- fread(file.path(path_result, "04_scRNA", "cross_disease", "tables", "T01_RRHO_local_signature_cell_type_summary.csv"))

# Panel A: Silicosis BALF scores
sil_scores <- sc_summary[dataset == "Silicosis BALF"]
sil_scores[, cell_type := factor(cell_type, levels = rev(cell_type[order(median_donor_score)]))]

p_sil <- ggplot(sil_scores, aes(x = median_donor_score, y = cell_type)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3, color = "grey50") +
  geom_point(size = 2, color = "#E41A1C") +
  geom_errorbarh(aes(xmin = q1_donor_score, xmax = q3_donor_score), height = 0.2, linewidth = 0.4) +
  labs(x = "Shared module score", y = NULL,
       title = "Silicosis BALF",
       subtitle = "GSE174725, n = 5 donors") +
  theme_nature()

# Panel B: TB lung scores
tb_scores <- sc_summary[dataset == "TB lung"]
tb_scores[, cell_type := factor(cell_type, levels = rev(cell_type[order(median_donor_score)]))]

p_tb <- ggplot(tb_scores, aes(x = median_donor_score, y = cell_type)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3, color = "grey50") +
  geom_point(size = 2, color = "#377EB8") +
  geom_errorbarh(aes(xmin = q1_donor_score, xmax = q3_donor_score), height = 0.2, linewidth = 0.4) +
  labs(x = "Shared module score", y = NULL,
       title = "TB lung lesions",
       subtitle = "GSE192483, n = 11 patients") +
  theme_nature()

# Save Figure 4
fig4 <- p_sil | p_tb +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(tag_levels = "a")

out <- file.path(fig_dir, "Fig4_single_cell_localization")
svglite::svglite(paste0(out, ".svg"), width = w, height = h); print(fig4); dev.off()
grDevices::cairo_pdf(paste0(out, ".pdf"), width = w, height = h, family = "Arial"); print(fig4); dev.off()
ragg::agg_tiff(paste0(out, ".tiff"), width = w, height = h, units = "in", res = 600, compression = "lzw"); print(fig4); dev.off()
ragg::agg_png(paste0(out, ".png"), width = w, height = h, units = "in", res = 300); print(fig4); dev.off()
cat("  Fig4 saved\n")

# ============================================================
# Figure 5: T2DM x Mtb Macrophage Interaction
# ============================================================
cat("\n5. Generating Figure 5: T2DM x Mtb interaction...\n")

# Load GSE283452 results
t2dm_mtb <- fread(file.path(path_result, "06_final", "tables", "GSE283452_T01_GSE283452_Alveolar_macrophage_T2DM_x_Mtb.csv"))
t2dm_mtb[, fdr := padj]
t2dm_mtb[, significance := "NS"]
t2dm_mtb[fdr < 0.05 & log2FoldChange > 0, significance := "Up in T2DM+Mtb"]
t2dm_mtb[fdr < 0.05 & log2FoldChange < 0, significance := "Down in T2DM+Mtb"]

# Panel A: Volcano plot
p_volcano <- ggplot(t2dm_mtb[!is.na(fdr)], aes(x = log2FoldChange, y = -log10(fdr), color = significance)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.3, color = "grey50") +
  geom_vline(xintercept = c(-1, 1), linetype = 2, linewidth = 0.3, color = "grey50") +
  scale_color_manual(values = c(
    "Up in T2DM+Mtb" = "#E41A1C",
    "Down in T2DM+Mtb" = "#377EB8",
    "NS" = "#999999"
  ), name = "Significance") +
  labs(x = "log2FC (T2DM × Mtb interaction)", y = "-log10(FDR)",
       title = "T2DM reprograms macrophage Mtb response",
       subtitle = "GSE283452, alveolar macrophages") +
  theme_nature()

# Panel B: Top interaction genes
top_int <- t2dm_mtb[fdr < 0.05][order(fdr)][1:min(10, .N)]

p_top <- ggplot(top_int, aes(x = reorder(gene_symbol, -log10(fdr)), y = -log10(fdr))) +
  geom_col(fill = "#E41A1C", width = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.3, color = "grey50") +
  labs(x = NULL, y = "-log10(FDR)",
       title = "Top interaction genes",
       subtitle = "T2DM × Mtb in alveolar macrophages") +
  theme_nature() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5))

# Save Figure 5
fig5 <- p_volcano | p_top +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(tag_levels = "a")

out <- file.path(fig_dir, "Fig5_T2DM_Mtb_interaction")
svglite::svglite(paste0(out, ".svg"), width = w, height = h); print(fig5); dev.off()
grDevices::cairo_pdf(paste0(out, ".pdf"), width = w, height = h, family = "Arial"); print(fig5); dev.off()
ragg::agg_tiff(paste0(out, ".tiff"), width = w, height = h, units = "in", res = 600, compression = "lzw"); print(fig5); dev.off()
ragg::agg_png(paste0(out, ".png"), width = w, height = h, units = "in", res = 300); print(fig5); dev.off()
cat("  Fig5 saved\n")

# ============================================================
# Summary
# ============================================================
cat("\n=== All figures generated ===\n")
cat("Figure 1: GBD Global Co-burden\n")
cat("Figure 2: Cross-disease Transcriptomic Convergence\n")
cat("Figure 3: Cross-cohort TB-DM Meta-analysis\n")
cat("Figure 4: Single-cell Localization\n")
cat("Figure 5: T2DM × Mtb Macrophage Interaction\n")

cat("\nFiles in", fig_dir, ":\n")
main_figs <- list.files(fig_dir, pattern = "^Fig[1-5]")
print(main_figs)

write_log("Nature-style figures generated")
