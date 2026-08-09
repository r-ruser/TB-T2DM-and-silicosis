source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "ggplot2", "patchwork", "svglite", "ragg"))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))

cat("=== GBD 2023 Prevalence Figures ===\n\n")

# ============================================================
# 1. Load data
# ============================================================
result_dir <- file.path(path_result, "01_GBD2023")
fig_dir <- file.path(result_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

prev_wide <- data.table::fread(file.path(result_dir, "tables", "T01_GBD2023_prevalence_CBI.csv"))

# ============================================================
# 2. Create non-map figures (bar charts, scatter plots)
# ============================================================
cat("Creating figures...\n")

# Theme
theme_pub <- theme_minimal(base_size = 7, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3),
    legend.position = "top",
    plot.title = element_text(size = 8, face = "bold"),
    plot.subtitle = element_text(size = 6, colour = "#767676"),
    strip.text = element_text(size = 7, face = "bold")
  )

# Bar chart: Top 20 countries by CBI
top20 <- prev_wide[order(-cbi)][1:min(20, nrow(prev_wide))]
top20[, location_name := factor(location_name, levels = rev(location_name))]

p_bar <- ggplot(top20, aes(x = cbi, y = location_name)) +
  geom_col(aes(fill = burden_class), width = 0.7) +
  scale_fill_manual(values = c(
    "TB + DM + Silicosis high" = "#f781bf",
    "TB + DM high" = "#ff7f00",
    "TB + Silicosis high" = "#984ea3",
    "DM + Silicosis high" = "#a65628",
    "TB high only" = "#e41a1c",
    "DM high only" = "#377eb8",
    "Silicosis high only" = "#4daf4a",
    "Low-Low-Low" = "#f0f0f0"
  )) +
  labs(title = "A. Top 20 countries by combined burden index",
       x = "Combined burden index", y = NULL, fill = "Burden class") +
  theme_pub

# Scatter: TB vs T2DM prevalence
p_scatter <- ggplot(prev_wide, aes(x = log10(val_T2DM), y = log10(val_TB))) +
  geom_point(aes(color = burden_class), size = 1.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.5) +
  scale_color_manual(values = c(
    "TB + DM + Silicosis high" = "#f781bf",
    "TB + DM high" = "#ff7f00",
    "TB + Silicosis high" = "#984ea3",
    "DM + Silicosis high" = "#a65628",
    "TB high only" = "#e41a1c",
    "DM high only" = "#377eb8",
    "Silicosis high only" = "#4daf4a",
    "Low-Low-Low" = "#999999"
  )) +
  labs(title = "B. TB vs T2DM prevalence (log10 scale)",
       x = "log10(T2DM prevalence per 100k)",
       y = "log10(TB prevalence per 100k)",
       color = "Burden class") +
  theme_pub

# Boxplot: Prevalence by burden class
df_long <- data.table::melt(prev_wide,
                            id.vars = c("location_name", "burden_class"),
                            measure.vars = c("val_TB", "val_T2DM", "val_Silicosis"),
                            variable.name = "disease",
                            value.name = "prevalence")
df_long[, disease := gsub("val_", "", disease)]

p_box <- ggplot(df_long, aes(x = disease, y = log10(prevalence + 1), fill = disease)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  facet_wrap(~ burden_class, ncol = 4) +
  scale_fill_manual(values = c("TB" = "#e41a1c", "T2DM" = "#377eb8", "Silicosis" = "#4daf4a")) +
  labs(title = "C. Prevalence distribution by burden class",
       x = NULL, y = "log10(prevalence + 1)") +
  theme_pub +
  theme(legend.position = "none")

# Sensitivity analysis plot
sens_dt <- data.table::fread(file.path(result_dir, "tables", "T03_threshold_sensitivity_prevalence.csv"))

p_sens <- ggplot(sens_dt, aes(x = factor(threshold), y = n_countries, group = 1)) +
  geom_line(linewidth = 0.5, color = "#416A9A") +
  geom_point(size = 2, shape = 21, fill = "#B84A4A") +
  labs(title = "D. Triple-high burden counts by threshold",
       x = "Percentile threshold", y = "Number of countries") +
  theme_pub

# ============================================================
# 3. Combine and save
# ============================================================
cat("\nSaving figures...\n")

fig_combined <- (p_bar | p_scatter) / (p_box | p_sens) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "a",
                  title = "GBD 2023 global prevalence of TB, T2DM, and silicosis",
                  subtitle = "Age-standardized rates, both sexes, 204 locations")

out_base <- file.path(fig_dir, "F01_GBD2023_prevalence_analysis")
w <- 183 / 25.4; h <- 180 / 25.4

svglite::svglite(paste0(out_base, ".svg"), width = w, height = h)
print(fig_combined)
dev.off()

grDevices::cairo_pdf(paste0(out_base, ".pdf"), width = w, height = h, family = "Arial")
print(fig_combined)
dev.off()

ragg::agg_tiff(paste0(out_base, ".tiff"), width = w, height = h, units = "in", res = 600, compression = "lzw")
print(fig_combined)
dev.off()

ragg::agg_png(paste0(out_base, ".png"), width = w, height = h, units = "in", res = 300)
print(fig_combined)
dev.off()

cat("Figures saved to:", fig_dir, "\n")

# ============================================================
# 4. Summary table
# ============================================================
cat("\nCreating summary table...\n")

summary_dt <- data.table::data.table(
  disease = c("TB", "T2DM", "Silicosis"),
  prevalence_mean = c(mean(prev_wide$val_TB, na.rm = TRUE), mean(prev_wide$val_T2DM, na.rm = TRUE), mean(prev_wide$val_Silicosis, na.rm = TRUE)),
  prevalence_median = c(median(prev_wide$val_TB, na.rm = TRUE), median(prev_wide$val_T2DM, na.rm = TRUE), median(prev_wide$val_Silicosis, na.rm = TRUE)),
  prevalence_iqr = c(
    paste(round(quantile(prev_wide$val_TB, 0.25, na.rm = TRUE), 1), "-", round(quantile(prev_wide$val_TB, 0.75, na.rm = TRUE), 1)),
    paste(round(quantile(prev_wide$val_T2DM, 0.25, na.rm = TRUE), 1), "-", round(quantile(prev_wide$val_T2DM, 0.75, na.rm = TRUE), 1)),
    paste(round(quantile(prev_wide$val_Silicosis, 0.25, na.rm = TRUE), 2), "-", round(quantile(prev_wide$val_Silicosis, 0.75, na.rm = TRUE), 2))
  )
)

data.table::fwrite(summary_dt, file.path(result_dir, "tables", "T04_prevalence_summary.csv"), bom = TRUE)

cat("\n=== Summary ===\n")
print(summary_dt)

write_log("GBD 2023 prevalence figures completed")
