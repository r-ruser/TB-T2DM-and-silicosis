source(file.path("R", "00_config.R"), encoding = "UTF-8")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(sf))
suppressPackageStartupMessages(library(rnaturalearth))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(ragg))

cat("=== GBD 2023 Figure 1 Rebuild (TB Incidence) ===\n\n")

# ============================================================
# 1. Load audited data
# ============================================================
cat("1. Loading audited data...\n")

# TB Incidence (204 locations)
tb_inc <- fread(file.path(path_result, "01_GBD2023", "tables", "T03_TB_incidence_audit.csv"))
cat("  TB Incidence:", nrow(tb_inc), "locations\n")

# T2DM Prevalence
t2dm_prev <- fread(file.path(path_result, "01_GBD2023", "tables", "T04_T2DM_audit.csv"))
cat("  T2DM Prevalence:", nrow(t2dm_prev), "locations\n")

# Silicosis Prevalence
sil_prev <- fread(file.path(path_result, "01_GBD2023", "tables", "T05_silicosis_audit.csv"))
cat("  Silicosis Prevalence:", nrow(sil_prev), "locations\n")

# ============================================================
# 2. Location crosswalk
# ============================================================
cat("\n2. Creating location crosswalk...\n")

# Get common locations
common_locs <- Reduce(intersect, list(
  tb_inc$location_id,
  t2dm_prev$location_id,
  sil_prev$location_id
))
cat("  Common locations:", length(common_locs), "\n")

# Merge data
tb_sub <- tb_inc[location_id %in% common_locs, .(location_id, location, tb_value = value, tb_lower = lower, tb_upper = upper)]
t2dm_sub <- t2dm_prev[location_id %in% common_locs, .(location_id, t2dm_value = value, t2dm_lower = lower, t2dm_upper = upper)]
sil_sub <- sil_prev[location_id %in% common_locs, .(location_id, sil_value = value, sil_lower = lower, sil_upper = upper)]

merged <- merge(tb_sub, t2dm_sub, by = "location_id", all = TRUE)
merged <- merge(merged, sil_sub, by = "location_id", all = TRUE)

# Save crosswalk
fwrite(merged, file.path(path_result, "01_GBD2023", "tables", "T06_location_crosswalk_final.csv"), bom = TRUE)
cat("  Crosswalk saved:", nrow(merged), "locations\n")

# ============================================================
# 3. Percentile ranks and co-burden scores
# ============================================================
cat("\n3. Calculating co-burden scores...\n")

# Percentile ranks (empirical)
calc_pct <- function(x) (rank(x, na.last = "keep") - 0.5) / sum(!is.na(x))

merged[, p_tb := calc_pct(tb_value)]
merged[, p_dm := calc_pct(t2dm_value)]
merged[, p_sil := calc_pct(sil_value)]

# Geometric mean (primary)
merged[, geo_score := (p_tb * p_dm * p_sil)^(1/3)]

# Arithmetic mean (sensitivity)
merged[, arith_score := (p_tb + p_dm + p_sil) / 3]

# Z-score composite (sensitivity)
merged[, z_tb := scale(log1p(tb_value))]
merged[, z_dm := scale(log1p(t2dm_value))]
merged[, z_sil := scale(log1p(sil_value))]
merged[, z_score := (z_tb + z_dm + z_sil) / 3]

# Rank
merged[, rank_geo := rank(-geo_score)]

# Save co-burden scores
fwrite(merged, file.path(path_result, "01_GBD2023", "tables", "T07_GBD2023_coburden_scores.csv"), bom = TRUE)

cat("  Co-burden scores calculated for", nrow(merged), "locations\n")

# ============================================================
# 4. Triple-high classification
# ============================================================
cat("\n4. Triple-high classification...\n")

# Primary: P75
merged[, triple_high_p75 := p_tb >= 0.75 & p_dm >= 0.75 & p_sil >= 0.75]
cat("  P75 triple-high:", sum(merged$triple_high_p75, na.rm = TRUE), "locations\n")

# Sensitivity thresholds
for (thr in c(0.70, 0.80, 0.90)) {
  col <- paste0("triple_high_p", thr * 100)
  merged[, (col) := p_tb >= thr & p_dm >= thr & p_sil >= thr]
  n <- sum(merged[[col]], na.rm = TRUE)
  cat("  P", thr * 100, " triple-high:", n, "locations\n")
}

# Eight-type classification (based on P75)
merged[, burden_class := "Low-Low-Low"]
merged[p_tb >= 0.75 & p_dm < 0.75 & p_sil < 0.75, burden_class := "TB high only"]
merged[p_tb < 0.75 & p_dm >= 0.75 & p_sil < 0.75, burden_class := "T2DM high only"]
merged[p_tb < 0.75 & p_dm < 0.75 & p_sil >= 0.75, burden_class := "Silicosis high only"]
merged[p_tb >= 0.75 & p_dm >= 0.75 & p_sil < 0.75, burden_class := "TB + T2DM high"]
merged[p_tb >= 0.75 & p_dm < 0.75 & p_sil >= 0.75, burden_class := "TB + Silicosis high"]
merged[p_tb < 0.75 & p_dm >= 0.75 & p_sil >= 0.75, burden_class := "T2DM + Silicosis high"]
merged[p_tb >= 0.75 & p_dm >= 0.75 & p_sil >= 0.75, burden_class := "TB + T2DM + Silicosis high"]

# Save triple-high results
fwrite(merged, file.path(path_result, "01_GBD2023", "tables", "T08_triple_high_threshold_sensitivity.csv"), bom = TRUE)

# ============================================================
# 5. Threshold sensitivity
# ============================================================
cat("\n5. Threshold sensitivity analysis...\n")

sens_results <- list()
for (thr in c(0.70, 0.75, 0.80, 0.90)) {
  n_triple <- sum(merged$p_tb >= thr & merged$p_dm >= thr & merged$p_sil >= thr, na.rm = TRUE)
  triple_locs <- merged[p_tb >= thr & p_dm >= thr & p_sil >= thr]$location
  sens_results[[as.character(thr)]] <- data.table(
    threshold = thr,
    n_locations = n_triple,
    locations = paste(triple_locs, collapse = "; ")
  )
}
sens_dt <- rbindlist(sens_results)
fwrite(sens_dt, file.path(path_result, "01_GBD2023", "tables", "T09_threshold_sensitivity.csv"), bom = TRUE)

# ============================================================
# 6. Regional distribution
# ============================================================
cat("\n6. Regional distribution...\n")

# Add region classification (simplified)
merged[, region := "Other"]
merged[grepl("Pacific|Marshall|Fiji|Samoa|Tonga|Kiribati|Tuvalu|Nauru|Papua", location, ignore.case = TRUE), region := "Pacific Islands"]
merged[grepl("China|Japan|Korea|India|Indonesia|Thailand|Vietnam|Philippines|Malaysia|Myanmar|Cambodia|Laos|Bangladesh|Pakistan|Nepal|Sri Lanka", location, ignore.case = TRUE), region := "Asia"]
merged[grepl("United States|Canada|Mexico|Brazil|Argentina|Colombia|Chile|Peru|Venezuela|Ecuador|Bolivia|Paraguay|Uruguay|Guatemala|Honduras|El Salvador|Nicaragua|Costa Rica|Panama|Cuba|Dominican|Haiti|Jamaica|Trinidad|Guyana|Suriname|Belize", location, ignore.case = TRUE), region := "Americas"]
merged[grepl("United Kingdom|Germany|France|Italy|Spain|Portugal|Netherlands|Belgium|Switzerland|Austria|Sweden|Norway|Denmark|Finland|Ireland|Poland|Czech|Hungary|Romania|Bulgaria|Greece|Croatia|Slovenia|Slovakia|Estonia|Latvia|Lithuania|Serbia|Montenegro|Albania|North Macedonia|Bosnia|Moldova|Ukraine|Belarus|Russia|Cyprus|Malta|Luxembourg|Iceland", location, ignore.case = TRUE), region := "Europe"]
merged[grepl("Nigeria|Ethiopia|DR Congo|Tanzania|South Africa|Kenya|Uganda|Ghana|Mozambique|Madagascar|Cameroon|Angola|Cote|Senegal|Mali|Burkina|Mozambique|Zambia|Zimbabwe|Malawi|Somalia|Rwanda|Burundi|Guinea|Benin|Togo|Niger|Chad|Central African|Sierra Leone|Liberia|Gambia|Cape Verde|Comoros|Mauritius|Seychelles|Equatorial Guinea|Gabon", location, ignore.case = TRUE), region := "Africa"]

# Regional summary
region_summary <- merged[, .(
  locations_total = .N,
  triple_high_n = sum(triple_high_p75, na.rm = TRUE),
  triple_high_pct = round(100 * sum(triple_high_p75, na.rm = TRUE) / .N, 1)
), by = region]

fwrite(region_summary, file.path(path_result, "01_GBD2023", "tables", "T10_regional_distribution_triple_high.csv"), bom = TRUE)

cat("  Regional distribution:\n")
print(region_summary)

# ============================================================
# 7. Co-burden score robustness
# ============================================================
cat("\n7. Co-burden score robustness...\n")

# Spearman correlations
cor_geo_arith <- cor.test(merged$geo_score, merged$arith_score, method = "spearman")
cor_geo_z <- cor.test(merged$geo_score, merged$z_score, method = "spearman")

cat("  Geo vs Arithmetic:", round(cor_geo_arith$estimate, 3), "\n")
cat("  Geo vs Z-score:", round(cor_geo_z$estimate, 3), "\n")

# Top-location overlap
top20_geo <- merged[order(-geo_score)][1:20]$location
top20_arith <- merged[order(-arith_score)][1:20]$location
top20_z <- merged[order(-z_score)][1:20]$location

overlap_geo_arith <- length(intersect(top20_geo, top20_arith))
overlap_geo_z <- length(intersect(top20_geo, top20_z))

cat("  Top20 overlap (Geo vs Arithmetic):", overlap_geo_arith, "/20\n")
cat("  Top20 overlap (Geo vs Z-score):", overlap_geo_z, "/20\n")

# Save robustness results
robustness <- data.table(
  comparison = c("Geo vs Arithmetic", "Geo vs Z-score"),
  spearman_rho = c(cor_geo_arith$estimate, cor_geo_z$estimate),
  top20_overlap = c(overlap_geo_arith, overlap_geo_z)
)
fwrite(robustness, file.path(path_result, "01_GBD2023", "tables", "T11_score_robustness.csv"), bom = TRUE)

# ============================================================
# 8. Generate Figure 1
# ============================================================
cat("\n8. Generating Figure 1...\n")

# Nature theme
theme_nature <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.3, color = "black"),
      axis.ticks = element_line(linewidth = 0.3, color = "black"),
      axis.ticks.length = unit(1.5, "mm"),
      axis.text = element_text(size = 6, color = "black"),
      axis.title = element_text(size = 7, color = "black"),
      legend.position = "right",
      legend.key.size = unit(2, "mm"),
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 7, face = "bold"),
      plot.title = element_text(size = 8, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 6, color = "#666666", hjust = 0),
      plot.tag = element_text(size = 8, face = "bold"),
      plot.margin = margin(3, 3, 3, 3)
    )
}

# Panel A: TB Incidence (bar chart of top locations)
top30_tb <- merged[order(-tb_value)][1:min(30, nrow(merged))]
top30_tb[, location := factor(location, levels = rev(location))]

p_a <- ggplot(top30_tb, aes(x = tb_value, y = location)) +
  geom_col(fill = "#E41A1C", width = 0.7) +
  labs(x = "Age-standardized TB incidence rate\n(per 100,000)", y = NULL,
       title = "A. TB incidence, 2023") +
  theme_nature()

# Panel B: T2DM Prevalence (bar chart)
top30_dm <- merged[order(-t2dm_value)][1:min(30, nrow(merged))]
top30_dm[, location := factor(location, levels = rev(location))]

p_b <- ggplot(top30_dm, aes(x = t2dm_value, y = location)) +
  geom_col(fill = "#FF7F0E", width = 0.7) +
  labs(x = "Age-standardized T2DM prevalence rate\n(per 100,000)", y = NULL,
       title = "B. T2DM prevalence, 2023") +
  theme_nature()

# Panel C: Silicosis Prevalence (bar chart)
top30_sil <- merged[order(-sil_value)][1:min(30, nrow(merged))]
top30_sil[, location := factor(location, levels = rev(location))]

p_c <- ggplot(top30_sil, aes(x = sil_value, y = location)) +
  geom_col(fill = "#1F77B4", width = 0.7) +
  labs(x = "Age-standardized silicosis prevalence rate\n(per 100,000)", y = NULL,
       title = "C. Silicosis prevalence, 2023") +
  theme_nature()

# Panel D: Co-burden score (histogram)
p_d <- ggplot(merged, aes(x = geo_score)) +
  geom_histogram(bins = 30, fill = "#9467BD", alpha = 0.7, color = "white", linewidth = 0.2) +
  geom_vline(xintercept = quantile(merged$geo_score, 0.75, na.rm = TRUE),
             linetype = 2, linewidth = 0.3, color = "red") +
  labs(x = "Descriptive co-burden score (geometric mean)", y = "Count",
       title = "D. Distribution of co-burden scores",
       subtitle = "Dashed line = 75th percentile") +
  theme_nature()

# Panel E: Triple-high classification (bar)
burden_counts <- merged[, .N, by = burden_class]
burden_counts[, burden_class := factor(burden_class, levels = c(
  "TB + T2DM + Silicosis high", "TB + T2DM high", "TB + Silicosis high",
  "T2DM + Silicosis high", "TB high only", "T2DM high only",
  "Silicosis high only", "Low-Low-Low"))]

p_e <- ggplot(burden_counts, aes(x = burden_class, y = N, fill = burden_class)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c(
    "TB + T2DM + Silicosis high" = "#E41A1C",
    "TB + T2DM high" = "#FF7F00",
    "TB + Silicosis high" = "#984EA3",
    "T2DM + Silicosis high" = "#A65628",
    "TB high only" = "#377EB8",
    "T2DM high only" = "#4DAF4A",
    "Silicosis high only" = "#999999",
    "Low-Low-Low" = "#F0F0F0"
  ), guide = "none") +
  labs(x = NULL, y = "Number of locations",
       title = "E. Triple-burden classification (P75)") +
  theme_nature() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5))

# Combine
fig1 <- (p_a / p_b / p_c) | (p_d / p_e) +
  plot_layout(widths = c(1, 1), heights = c(1, 1, 1)) +
  plot_annotation(tag_levels = "a",
                  title = "Global co-burden of TB, T2DM, and silicosis, 2023",
                  subtitle = "GBD 2023, age-standardized rates, 204 locations")

# Save
fig_dir <- file.path(path_result, "06_final", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(fig_dir, "Figure_1_GBD2023_global_coburden")
w <- 183/25.4; h <- 240/25.4

svglite::svglite(paste0(out, ".svg"), width = w, height = h); print(fig1); dev.off()
grDevices::cairo_pdf(paste0(out, ".pdf"), width = w, height = h, family = "Arial"); print(fig1); dev.off()
ragg::agg_tiff(paste0(out, ".tiff"), width = w, height = h, units = "in", res = 600, compression = "lzw"); print(fig1); dev.off()
ragg::agg_png(paste0(out, ".png"), width = w, height = h, units = "in", res = 300); print(fig1); dev.off()

cat("  Figure 1 saved\n")

# ============================================================
# 9. Final audit
# ============================================================
cat("\n=== FINAL GBD AUDIT ===\n")
cat("\nTB:\n")
cat("  Cause = Tuberculosis\n")
cat("  Measure = Incidence\n")
cat("  Metric = Rate\n")
cat("  Age = Age-standardized\n")
cat("  Sex = Both\n")
cat("  Year = 2023\n")
cat("  Unit = per 100,000\n")
cat("  N locations =", nrow(tb_inc), "\n")
cat("  Min =", min(tb_inc$value, na.rm = TRUE), "\n")
cat("  Median =", median(tb_inc$value, na.rm = TRUE), "\n")
cat("  Max =", max(tb_inc$value, na.rm = TRUE), "\n")

cat("\nT2DM:\n")
cat("  Cause = Type 2 diabetes mellitus\n")
cat("  Measure = Prevalence\n")
cat("  Metric = Rate\n")
cat("  Age = Age-standardized\n")
cat("  Sex = Both\n")
cat("  Year = 2023\n")
cat("  Unit = per 100,000\n")
cat("  N locations =", nrow(t2dm_prev), "\n")
cat("  Min =", min(t2dm_prev$value, na.rm = TRUE), "\n")
cat("  Median =", median(t2dm_prev$value, na.rm = TRUE), "\n")
cat("  Max =", max(t2dm_prev$value, na.rm = TRUE), "\n")

cat("\nSilicosis:\n")
cat("  Cause = Silicosis\n")
cat("  Measure = Prevalence\n")
cat("  Metric = Rate\n")
cat("  Age = Age-standardized\n")
cat("  Sex = Both\n")
cat("  Year = 2023\n")
cat("  Unit = per 100,000\n")
cat("  N locations =", nrow(sil_prev), "\n")
cat("  Min =", min(sil_prev$value, na.rm = TRUE), "\n")
cat("  Median =", median(sil_prev$value, na.rm = TRUE), "\n")
cat("  Max =", max(sil_prev$value, na.rm = TRUE), "\n")

cat("\nTriple-high P75:\n")
cat("  N =", sum(merged$triple_high_p75, na.rm = TRUE), "\n")
cat("  Locations:", paste(merged[triple_high_p75 == TRUE]$location, collapse = ", "), "\n")

cat("\nSensitivity:\n")
cat("  P70 N =", sum(merged$p_tb >= 0.70 & merged$p_dm >= 0.70 & merged$p_sil >= 0.70, na.rm = TRUE), "\n")
cat("  P80 N =", sum(merged$p_tb >= 0.80 & merged$p_dm >= 0.80 & merged$p_sil >= 0.80, na.rm = TRUE), "\n")
cat("  P90 N =", sum(merged$p_tb >= 0.90 & merged$p_dm >= 0.90 & merged$p_sil >= 0.90, na.rm = TRUE), "\n")

cat("\n=== AUDIT COMPLETE ===\n")

write_log("GBD Figure 1 rebuild completed: TB INCIDENCE used")
