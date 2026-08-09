source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "ggplot2", "patchwork", "sf", "rnaturalearth",
                  "rnaturalearthdata", "svglite", "ragg"))

f <- file.path(path_data, "01_GBD2023", "processed", "GBD2023_proxy_derived.csv")
if (!file.exists(f)) stop("Run R/11_gbd_proxy_analysis.R first")
d <- data.table::fread(f, encoding = "UTF-8")
d[, triple_high_proxy := factor(as.logical(triple_high_proxy), levels = c(FALSE, TRUE),
                                labels = c("No", "Yes"))]
d[, pam_cluster_proxy := factor(pam_cluster_proxy)]

pal <- c(neutral_dark = "#272727", neutral_mid = "#767676", neutral_light = "#D8D8D8",
         blue = "#416A9A", blue_light = "#A8BED8", rose = "#C9828D", rose_light = "#E7C1C7",
         ochre = "#C49A55", teal = "#4F9A91", red = "#B84A4A")
theme_nature <- function(base_size = 6.5) {
  ggplot2::theme_classic(base_size = base_size, base_family = "Arial") +
    ggplot2::theme(axis.line = ggplot2::element_line(linewidth = 0.3),
      axis.ticks = ggplot2::element_line(linewidth = 0.3),
      plot.title = ggplot2::element_text(size = 7, face = "bold", margin = ggplot2::margin(b = 2)),
      plot.subtitle = ggplot2::element_text(size = 5.8, colour = pal["neutral_mid"]),
      legend.title = ggplot2::element_text(size = 5.8), legend.text = ggplot2::element_text(size = 5.3),
      strip.text = ggplot2::element_text(size = 6, face = "bold"),
      plot.tag = ggplot2::element_text(size = 8, face = "bold"), panel.grid = ggplot2::element_blank())
}
ggplot2::theme_set(theme_nature())

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
normalize_zh <- function(x) {
  x <- gsub("[[:space:]]", "", x)
  x <- gsub("共和国|民主共和国|联合共和国|联邦|公国|国$|岛$|群岛$", "", x)
  x
}
world$join_zh <- normalize_zh(world$name_zh)
d$join_zh <- normalize_zh(d$location_name)

aliases <- c(
  "中华人民共和国" = "中国", "中华民国" = "台湾", "美国" = "美利坚合众",
  "老挝" = "老挝人民", "马耳他" = "马尔他", "伯利兹" = "伯利兹城",
  "特立尼达和多巴哥" = "特立尼达拉岛和多巴哥", "坦桑尼亚" = "坦桑尼亚联合",
  "科特迪瓦" = "科特廸亚", "韩国" = "大韩民", "朝鲜" = "朝鲜人民",
  "俄罗斯" = "俄罗斯", "越南" = "越南", "老挝" = "老挝人民民主",
  "坦桑尼亚" = "坦桑尼亚联合", "玻利维亚" = "玻利维亚",
  "委内瑞拉" = "委内瑞拉玻利瓦尔", "伊朗" = "伊朗伊斯兰",
  "叙利亚" = "阿拉伯叙利亚", "摩尔多瓦" = "摩尔多瓦",
  "文莱" = "文莱达鲁萨兰", "巴勒斯坦" = "巴勒斯坦"
)
for (nm in names(aliases)) world$join_zh[world$join_zh == nm] <- aliases[[nm]]

map <- merge(world, d, by = "join_zh", all.x = TRUE, sort = FALSE)
matched_names <- unique(map$location_name[!is.na(map$location_id)])
unmatched <- d[!location_name %in% matched_names, .(location_id, location_name)]
safe_write_csv(unmatched, file.path(path_result, "00_audit", "A05_GBD_map_unmatched_locations.csv"))

map_panel <- function(variable, title, fill_name, trans = "log10") {
  ggplot2::ggplot(map) +
    ggplot2::geom_sf(ggplot2::aes(fill = .data[[variable]]), colour = "white", linewidth = 0.05) +
    ggplot2::scale_fill_gradientn(colours = c("#F1F3F5", pal["blue_light"], pal["blue"], "#243A5A"),
      trans = trans, na.value = "#ECECEC", name = fill_name) +
    ggplot2::coord_sf(crs = "+proj=robin", expand = FALSE) +
    ggplot2::labs(title = title, subtitle = "2023, age-standardized rate; proxy analysis") +
    ggplot2::theme_void(base_family = "Arial", base_size = 6.5) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 7, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 5.5, colour = pal["neutral_mid"]),
      legend.title = ggplot2::element_text(size = 5.5), legend.text = ggplot2::element_text(size = 5),
      legend.key.height = grid::unit(7, "mm"), plot.tag = ggplot2::element_text(size = 8, face = "bold"))
}

p_a <- map_panel("tb_asir", "Tuberculosis incidence", "Rate")
p_b <- map_panel("t2dm_asir_proxy", "Type 2 diabetes incidence proxy", "Rate")
p_c <- map_panel("pneumoconiosis_asir_proxy", "Pneumoconiosis incidence proxy", "Rate")

class_levels <- c("Low-Low-Low", "TB high only", "DM high only", "Pneumoconiosis high only",
                  "TB + DM high", "TB + Pneumoconiosis high", "DM + Pneumoconiosis high",
                  "TB + DM + Pneumoconiosis high")
map$burden_class_proxy <- factor(map$burden_class_proxy, levels = class_levels)
class_cols <- c("#E3E3E3", "#6F8FB7", "#D59BA4", "#C6A563", "#967AA1", "#4F8A88", "#B37B5B", "#B84A4A")
p_d <- ggplot2::ggplot(map) +
  ggplot2::geom_sf(ggplot2::aes(fill = burden_class_proxy), colour = "white", linewidth = 0.05) +
  ggplot2::scale_fill_manual(values = stats::setNames(class_cols, class_levels), drop = FALSE,
                             na.value = "#ECECEC", name = "P75 class") +
  ggplot2::coord_sf(crs = "+proj=robin", expand = FALSE) +
  ggplot2::labs(title = "Eight-class overlap", subtitle = "Disease-total incidence proxies") +
  ggplot2::theme_void(base_family = "Arial", base_size = 6.5) +
  ggplot2::theme(plot.title = ggplot2::element_text(size = 7, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 5.5, colour = pal["neutral_mid"]),
    legend.title = ggplot2::element_text(size = 5.5), legend.text = ggplot2::element_text(size = 4.8),
    legend.key.height = grid::unit(3, "mm"), plot.tag = ggplot2::element_text(size = 8, face = "bold"))

ranked <- d[order(cbi_proxy)]
ranked[, rank := seq_len(.N)]
name_lookup <- data.table::as.data.table(sf::st_drop_geometry(world))[, .(join_zh, name_en)]
ranked <- merge(ranked, name_lookup, by = "join_zh", all.x = TRUE, sort = FALSE)
ranked[, label_en := data.table::fifelse(is.na(name_en) | name_en == "",
                                         paste0("Location ", location_id), name_en)]
data.table::setorder(ranked, rank)
lab <- ranked[rank > .N - 12]
p_e <- ggplot2::ggplot(ranked, ggplot2::aes(rank, cbi_proxy)) +
  ggplot2::geom_point(ggplot2::aes(colour = triple_high_proxy), size = 0.7, alpha = 0.85) +
  ggplot2::geom_text(data = lab, ggplot2::aes(label = label_en), size = 1.65,
                     hjust = -0.05, check_overlap = TRUE, family = "Arial") +
  ggplot2::scale_colour_manual(values = c(No = unname(pal["neutral_mid"]), Yes = unname(pal["red"])),
                               name = "Triple-high proxy") +
  ggplot2::coord_cartesian(xlim = c(0, nrow(ranked) * 1.12), clip = "off") +
  ggplot2::labs(title = "Combined burden index", subtitle = "Geometric mean of three percentiles",
                x = "Country rank", y = "CBI proxy") +
  ggplot2::theme(legend.position = "top")

long <- data.table::melt(d, id.vars = c("location_name", "pam_cluster_proxy"),
  measure.vars = c("z_tb", "z_dm_proxy", "z_pneu_proxy"),
  variable.name = "metric", value.name = "z")
long[, pam_cluster_proxy := factor(pam_cluster_proxy)]
long[, metric := factor(metric, levels = c("z_tb", "z_dm_proxy", "z_pneu_proxy"),
                        labels = c("TB", "T2DM proxy", "Pneumoconiosis proxy"))]
p_f <- ggplot2::ggplot(long, ggplot2::aes(metric, z, colour = pam_cluster_proxy)) +
  ggplot2::stat_summary(ggplot2::aes(group = pam_cluster_proxy), fun = mean, geom = "line", linewidth = 0.55) +
  ggplot2::stat_summary(fun = mean, geom = "point", size = 1.4) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, linetype = 2, colour = pal["neutral_mid"]) +
  ggplot2::scale_colour_brewer(palette = "Dark2", name = "PAM cluster") +
  ggplot2::labs(title = "Unsupervised burden profiles", subtitle = "Selected by mean silhouette",
                x = NULL, y = "Mean log-rate z-score") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1), legend.position = "top")

fig <- ((p_a | p_b) / (p_c | p_d) / (p_e | p_f)) +
  patchwork::plot_annotation(tag_levels = "a", title = "2023 global burden overlap — proxy sensitivity analysis",
    subtitle = "T2DM and pneumoconiosis incidence are proxies; this is not the preregistered T2DM/silicosis prevalence analysis") &
  ggplot2::theme(plot.tag = ggplot2::element_text(size = 8, face = "bold"),
                 plot.title = ggplot2::element_text(family = "Arial"))

out <- file.path(path_result, "01_GBD2023", "figures", "F01_GBD2023_burden_proxy")
w <- 183 / 25.4; h <- 220 / 25.4
svglite::svglite(paste0(out, ".svg"), width = w, height = h); print(fig); grDevices::dev.off()
assert_english_only_svg(paste0(out, ".svg"))
grDevices::cairo_pdf(paste0(out, ".pdf"), width = w, height = h, family = "Arial"); print(fig); grDevices::dev.off()
ragg::agg_tiff(paste0(out, ".tiff"), width = w, height = h, units = "in", res = 600,
               compression = "lzw"); print(fig); grDevices::dev.off()
ragg::agg_png(paste0(out, ".png"), width = w, height = h, units = "in", res = 300); print(fig); grDevices::dev.off()
write_log("GBD proxy figure exported; map match=", length(matched_names), "/", nrow(d))
