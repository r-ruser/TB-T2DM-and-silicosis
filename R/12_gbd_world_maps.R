source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "ggplot2", "patchwork", "sf", "svglite", "ragg"))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))

cat("=== GBD 2023 World Maps ===\n\n")

# ============================================================
# 1. Load data with English names
# ============================================================
result_dir <- file.path(path_result, "01_GBD2023")
fig_dir <- file.path(result_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

prev_wide <- data.table::fread(file.path(result_dir, "tables", "T01_GBD2023_prevalence_CBI_en.csv"))
cat("Loaded:", nrow(prev_wide), "locations\n")

# Filter to complete cases
prev_complete <- prev_wide[!is.na(location_name_en) & !is.na(val_TB) & !is.na(val_T2DM) & !is.na(val_Silicosis)]
cat("Complete cases:", nrow(prev_complete), "\n")

# ============================================================
# 2. Load world map
# ============================================================
cat("Loading world map...\n")
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
cat("World map countries:", nrow(world), "\n")

# Create matching column for merging
prev_complete[, country_match := location_name_en]

# Merge with world map
world_data <- merge(world, prev_complete, by.x = "name", by.y = "country_match", all.x = TRUE)
cat("Matched countries:", sum(!is.na(world_data$val_TB)), "\n")

# ============================================================
# 3. Create theme
# ============================================================
theme_map <- theme_minimal(base_size = 7, base_family = "Arial") +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(size = 8, face = "bold"),
    plot.subtitle = element_text(size = 6, colour = "#767676"),
    legend.position = "bottom",
    legend.key.width = unit(1.5, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 6)
  )

# ============================================================
# 4. Create maps
# ============================================================
cat("\nCreating maps...\n")

# TB prevalence map
p_tb <- ggplot(world_data) +
  geom_sf(aes(fill = val_TB), color = "white", linewidth = 0.1) +
  scale_fill_viridis_c(option = "inferno", name = "TB prevalence\n(per 100k)", na.value = "grey90") +
  labs(title = "A. TB prevalence rate, 2023",
       subtitle = "Age-standardized, both sexes") +
  theme_map

# T2DM prevalence map
p_dm <- ggplot(world_data) +
  geom_sf(aes(fill = val_T2DM), color = "white", linewidth = 0.1) +
  scale_fill_viridis_c(option = "inferno", name = "T2DM prevalence\n(per 100k)", na.value = "grey90") +
  labs(title = "B. T2DM prevalence rate, 2023",
       subtitle = "Age-standardized, both sexes") +
  theme_map

# Silicosis prevalence map
p_sil <- ggplot(world_data) +
  geom_sf(aes(fill = val_Silicosis), color = "white", linewidth = 0.1) +
  scale_fill_viridis_c(option = "inferno", name = "Silicosis prevalence\n(per 100k)", na.value = "grey90") +
  labs(title = "C. Silicosis prevalence rate, 2023",
       subtitle = "Age-standardized, both sexes") +
  theme_map

# Combined burden index map
p_cbi <- ggplot(world_data) +
  geom_sf(aes(fill = cbi), color = "white", linewidth = 0.1) +
  scale_fill_viridis_c(option = "inferno", name = "Combined\nburden index", na.value = "grey90") +
  labs(title = "D. Combined burden index, 2023",
       subtitle = "Geometric mean of percentiles") +
  theme_map

# Burden class map
burden_colors <- c(
  "Low-Low-Low" = "#f0f0f0",
  "TB high only" = "#e41a1c",
  "DM high only" = "#377eb8",
  "Silicosis high only" = "#4daf4a",
  "TB + DM high" = "#ff7f00",
  "TB + Silicosis high" = "#984ea3",
  "DM + Silicosis high" = "#a65628",
  "TB + DM + Silicosis high" = "#f781bf"
)

p_class <- ggplot(world_data) +
  geom_sf(aes(fill = burden_class), color = "white", linewidth = 0.1) +
  scale_fill_manual(values = burden_colors, name = "Burden class", na.value = "grey90") +
  labs(title = "E. Triple-burden classification, 2023",
       subtitle = "75th percentile threshold") +
  theme_map

# ============================================================
# 5. Combine and save
# ============================================================
cat("\nSaving maps...\n")

# Combined figure (2x2 + 1)
fig_maps <- (p_tb | p_dm) / (p_sil | p_cbi) / p_class +
  plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(tag_levels = "a",
                  title = "GBD 2023 global prevalence of TB, T2DM, and silicosis",
                  subtitle = "Age-standardized rates, both sexes, 204 locations")

out_base <- file.path(fig_dir, "F01_GBD2023_world_maps")
w <- 183 / 25.4; h <- 240 / 25.4

svglite::svglite(paste0(out_base, ".svg"), width = w, height = h)
print(fig_maps)
dev.off()

grDevices::cairo_pdf(paste0(out_base, ".pdf"), width = w, height = h, family = "Arial")
print(fig_maps)
dev.off()

ragg::agg_tiff(paste0(out_base, ".tiff"), width = w, height = h, units = "in", res = 600, compression = "lzw")
print(fig_maps)
dev.off()

ragg::agg_png(paste0(out_base, ".png"), width = w, height = h, units = "in", res = 300)
print(fig_maps)
dev.off()

cat("Maps saved to:", fig_dir, "\n")

# ============================================================
# 6. Additional figures
# ============================================================
cat("\nCreating additional figures...\n")

# Top 20 countries bar chart
top20 <- prev_complete[order(-cbi)][1:min(20, nrow(prev_complete))]
top20[, location_name_en := factor(location_name_en, levels = rev(location_name_en))]

p_bar <- ggplot(top20, aes(x = cbi, y = location_name_en)) +
  geom_col(aes(fill = burden_class), width = 0.7) +
  scale_fill_manual(values = burden_colors) +
  labs(title = "F. Top 20 countries by combined burden index",
       x = "Combined burden index", y = NULL, fill = "Burden class") +
  theme_minimal(base_size = 7, base_family = "Arial") +
  theme(legend.position = "right")

# Scatter: TB vs T2DM
p_scatter <- ggplot(prev_complete, aes(x = log10(val_T2DM), y = log10(val_TB))) +
  geom_point(aes(color = burden_class), size = 1.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.5) +
  scale_color_manual(values = burden_colors) +
  labs(title = "G. TB vs T2DM prevalence (log10 scale)",
       x = "log10(T2DM prevalence per 100k)",
       y = "log10(TB prevalence per 100k)",
       color = "Burden class") +
  theme_minimal(base_size = 7, base_family = "Arial")

# Combine additional figures
fig_additional <- p_bar / p_scatter +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "a")

out_additional <- file.path(fig_dir, "F02_GBD2023_additional_analysis")
w2 <- 183 / 25.4; h2 <- 180 / 25.4

svglite::svglite(paste0(out_additional, ".svg"), width = w2, height = h2)
print(fig_additional)
dev.off()

grDevices::cairo_pdf(paste0(out_additional, ".pdf"), width = w2, height = h2, family = "Arial")
print(fig_additional)
dev.off()

ragg::agg_tiff(paste0(out_additional, ".tiff"), width = w2, height = h2, units = "in", res = 600, compression = "lzw")
print(fig_additional)
dev.off()

ragg::agg_png(paste0(out_additional, ".png"), width = w2, height = h2, units = "in", res = 300)
print(fig_additional)
dev.off()

cat("Additional figures saved\n")

write_log("GBD 2023 world maps completed")
