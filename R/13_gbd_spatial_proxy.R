source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "sf", "spdep", "rnaturalearth", "rnaturalearthdata"))

f <- file.path(path_data, "01_GBD2023", "processed", "GBD2023_proxy_derived.csv")
if (!file.exists(f)) stop("Run GBD clean and proxy analysis first")
d <- data.table::fread(f, encoding = "UTF-8")
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
normalize_zh <- function(x) {
  x <- gsub("[[:space:]]", "", x)
  gsub("共和国|民主共和国|联合共和国|联邦|公国|国$|岛$|群岛$", "", x)
}
world$join_zh <- normalize_zh(world$name_zh)
d$join_zh <- normalize_zh(d$location_name)
aliases <- c("中华人民共和国" = "中国", "中华民国" = "台湾", "美国" = "美利坚合众",
  "老挝" = "老挝人民", "马耳他" = "马尔他", "伯利兹" = "伯利兹城",
  "特立尼达和多巴哥" = "特立尼达拉岛和多巴哥", "坦桑尼亚" = "坦桑尼亚联合",
  "科特迪瓦" = "科特廸亚", "韩国" = "大韩民", "朝鲜" = "朝鲜人民",
  "玻利维亚" = "玻利维亚", "委内瑞拉" = "委内瑞拉玻利瓦尔",
  "伊朗" = "伊朗伊斯兰", "叙利亚" = "阿拉伯叙利亚")
for (nm in names(aliases)) world$join_zh[world$join_zh == nm] <- aliases[[nm]]
sp <- merge(world, d, by = "join_zh", all = FALSE, sort = FALSE)
sp <- sf::st_make_valid(sp)

# Keep one Natural Earth feature per GBD location; avoid topology-unsafe cross-feature unions.
metric_cols <- c("tb_asir", "t2dm_asir_proxy", "pneumoconiosis_asir_proxy", "cbi_proxy")
sp <- sp[, c("location_id", "location_name", metric_cols, "geometry")]
sp <- sp[!duplicated(sp$location_id), ]
sp <- sf::st_make_valid(sp)

nb_queen <- spdep::poly2nb(sp, queen = TRUE, row.names = as.character(sp$location_id))
lw_queen <- spdep::nb2listw(nb_queen, style = "W", zero.policy = TRUE)
cent <- sf::st_transform(sf::st_point_on_surface(sf::st_transform(sp, 8857)), 8857)
coords <- sf::st_coordinates(cent)
dup <- duplicated(data.frame(x = coords[, 1], y = coords[, 2])) |
       duplicated(data.frame(x = coords[, 1], y = coords[, 2]), fromLast = TRUE)
if (any(dup)) {
  idx <- which(dup)
  coords[idx, 1] <- coords[idx, 1] + seq_along(idx) * 1e-4
  coords[idx, 2] <- coords[idx, 2] + seq_along(idx) * 1e-4
}
make_knn <- function(k) spdep::nb2listw(spdep::knn2nb(spdep::knearneigh(coords, k = k),
  row.names = as.character(sp$location_id)), style = "W", zero.policy = TRUE)
weights <- list(Queen = lw_queen, KNN4 = make_knn(4), KNN6 = make_knn(6))

global <- list(); local <- list()
for (metric in metric_cols) {
  x <- sp[[metric]]
  for (wn in names(weights)) {
    mt <- spdep::moran.test(x, weights[[wn]], zero.policy = TRUE, na.action = na.fail,
                            randomisation = TRUE)
    global[[length(global) + 1L]] <- data.table::data.table(
      metric = metric, weight = wn, n_locations = length(x),
      moran_i = unname(mt$estimate["Moran I statistic"]),
      expected_i = unname(mt$estimate["Expectation"]), variance = unname(mt$estimate["Variance"]),
      p_value = mt$p.value, alternative = mt$alternative)
    lm <- spdep::localmoran(x, weights[[wn]], zero.policy = TRUE, na.action = na.fail)
    z <- as.numeric(scale(x))
    lag_z <- spdep::lag.listw(weights[[wn]], z, zero.policy = TRUE)
    cluster <- data.table::fcase(
      z >= 0 & lag_z >= 0, "High-High", z < 0 & lag_z < 0, "Low-Low",
      z >= 0 & lag_z < 0, "High-Low", z < 0 & lag_z >= 0, "Low-High")
    local[[length(local) + 1L]] <- data.table::data.table(
      location_id = sp$location_id, location_name = sp$location_name,
      metric = metric, weight = wn, value = x, z = z, spatial_lag_z = lag_z,
      local_moran_i = lm[, "Ii"], z_score = lm[, "Z.Ii"], p_value = lm[, "Pr(z != E(Ii))"],
      p_fdr = stats::p.adjust(lm[, "Pr(z != E(Ii))"], method = "BH"), cluster = cluster)
  }
}
global <- data.table::rbindlist(global)
local <- data.table::rbindlist(local)
global[, p_fdr := stats::p.adjust(p_value, method = "BH")]
queen_isolates <- data.table::data.table(
  location_id = sp$location_id[spdep::card(nb_queen) == 0L],
  location_name = sp$location_name[spdep::card(nb_queen) == 0L])

safe_write_csv(global, file.path(path_result, "01_GBD2023", "models", "M04_global_Moran_proxy.csv"))
safe_write_csv(local, file.path(path_result, "01_GBD2023", "source_data", "SD03_local_Moran_proxy.csv"))
safe_write_csv(queen_isolates, file.path(path_result, "00_audit", "A10_Queen_island_locations.csv"))
saveRDS(list(spatial = sp, queen = nb_queen, weights = weights),
        file.path(path_result, "01_GBD2023", "models", "M05_spatial_objects_proxy.rds"))
write_log("GBD spatial proxy analysis completed: ", nrow(sp), " mapped locations; ",
          nrow(queen_isolates), " Queen isolates; KNN4/KNN6 sensitivity retained")
