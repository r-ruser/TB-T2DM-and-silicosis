source(file.path("R", "00_config.R"), encoding = "UTF-8")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(ragg))

cat("=== Figure 6: Integrated Mechanism Model ===\n\n")

# Create a conceptual diagram using ggplot2
# This is a schematic, not data-driven

# Define positions for components
components <- data.table(
  x = c(0.5, 0.5, 0.15, 0.5, 0.85, 0.5, 0.5),
  y = c(0.95, 0.75, 0.5, 0.5, 0.5, 0.25, 0.05),
  label = c(
    "GBD 2023: Global Co-burden\n5 Pacific Island locations",
    "Bulk Transcriptomics\nPathway-level convergence\n36 concordant pathways",
    "GSE174725\nSilicosis BALF\nMyeloid localization",
    "Mechanism Axes\n1. Innate immune activation\n2. IFN & antigen presentation\n3. Complement & coagulation\n4. Macrophage metabolism\n5. ECM remodeling",
    "GSE192483\nTB lung lesions\nPatient-level validation",
    "GSE283452\nT2DM × Mtb\nMacrophage reprogramming",
    "Conclusion\nDivergent transcriptional programs\nconverge on myeloid immune biology"
  ),
  color = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#4DAF4A", "#FF7F0E", "#66C2A5")
)

# Create the diagram
p_model <- ggplot(components, aes(x = x, y = y)) +
  # Draw boxes
  geom_label(aes(label = label, fill = color), size = 2.5, fontface = "plain",
             label.padding = unit(0.3, "lines"), label.size = 0.3, alpha = 0.8) +
  # Draw arrows
  annotate("segment", x = 0.5, xend = 0.5, y = 0.88, yend = 0.82,
           arrow = arrow(length = unit(0.15, "cm")), linewidth = 0.5) +
  annotate("segment", x = 0.5, xend = 0.15, y = 0.68, yend = 0.58,
           arrow = arrow(length = unit(0.15, "cm")), linewidth = 0.5) +
  annotate("segment", x = 0.5, xend = 0.5, y = 0.68, yend = 0.58,
           arrow = arrow(length = unit(0.15, "cm")), linewidth = 0.5) +
  annotate("segment", x = 0.5, xend = 0.85, y = 0.68, yend = 0.58,
           arrow = arrow(length = unit(0.15, "cm")), linewidth = 0.5) +
  annotate("segment", x = 0.5, xend = 0.5, y = 0.42, yend = 0.32,
           arrow = arrow(length = unit(0.15, "cm")), linewidth = 0.5) +
  annotate("segment", x = 0.5, xend = 0.5, y = 0.18, yend = 0.12,
           arrow = arrow(length = unit(0.15, "cm")), linewidth = 0.5) +
  # Scale and theme
  scale_fill_identity() +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "Integrated Mechanism Model",
       subtitle = "TB-T2DM-Silicosis convergent immune program") +
  theme_void(base_size = 7, base_family = "Arial") +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 8, color = "#666666", hjust = 0.5),
    plot.margin = margin(10, 10, 10, 10)
  )

# Save
fig_dir <- file.path(path_result, "06_final", "figures")
out <- file.path(fig_dir, "Figure_6_integrated_mechanism_model")
w <- 183/25.4; h <- 150/25.4

svglite::svglite(paste0(out, ".svg"), width = w, height = h); print(p_model); dev.off()
grDevices::cairo_pdf(paste0(out, ".pdf"), width = w, height = h, family = "Arial"); print(p_model); dev.off()
ragg::agg_tiff(paste0(out, ".tiff"), width = w, height = h, units = "in", res = 600, compression = "lzw"); print(p_model); dev.off()
ragg::agg_png(paste0(out, ".png"), width = w, height = h, units = "in", res = 300); print(p_model); dev.off()

cat("Figure 6 saved\n")
