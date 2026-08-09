source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "cluster", "broom"))

f <- file.path(path_data, "01_GBD2023", "processed", "GBD2023_proxy_rates.csv")
if (!file.exists(f)) stop("Run R/10_gbd_audit_clean.R first")
d <- data.table::fread(f, encoding = "UTF-8")

percentile_midrank <- function(x) (rank(x, ties.method = "average") - 0.5) / length(x)
d[, p_tb := percentile_midrank(tb_asir)]
d[, p_dm_proxy := percentile_midrank(t2dm_asir_proxy)]
d[, p_pneu_proxy := percentile_midrank(pneumoconiosis_asir_proxy)]
d[, cbi_proxy := (p_tb * p_dm_proxy * p_pneu_proxy)^(1 / 3)]
d[, cbi_arithmetic_proxy := (p_tb + p_dm_proxy + p_pneu_proxy) / 3]

z <- scale(log(d[, .(tb_asir, t2dm_asir_proxy, pneumoconiosis_asir_proxy)]))
d[, `:=`(z_tb = z[, 1], z_dm_proxy = z[, 2], z_pneu_proxy = z[, 3],
         z_combined_proxy = rowMeans(z))]

class_at <- function(dt, threshold) {
  q_tb <- as.numeric(quantile(dt$tb_asir, threshold, type = 7))
  q_dm <- as.numeric(quantile(dt$t2dm_asir_proxy, threshold, type = 7))
  q_pn <- as.numeric(quantile(dt$pneumoconiosis_asir_proxy, threshold, type = 7))
  h_tb <- dt$tb_asir >= q_tb
  h_dm <- dt$t2dm_asir_proxy >= q_dm
  h_pn <- dt$pneumoconiosis_asir_proxy >= q_pn
  cls <- data.table::fifelse(h_tb & h_dm & h_pn, "TB + DM + Pneumoconiosis high",
          data.table::fifelse(h_tb & h_dm, "TB + DM high",
          data.table::fifelse(h_tb & h_pn, "TB + Pneumoconiosis high",
          data.table::fifelse(h_dm & h_pn, "DM + Pneumoconiosis high",
          data.table::fifelse(h_tb, "TB high only", data.table::fifelse(h_dm, "DM high only",
          data.table::fifelse(h_pn, "Pneumoconiosis high only", "Low-Low-Low")))))))
  data.table::data.table(location_id = dt$location_id, threshold = threshold,
                         q_tb = q_tb, q_dm_proxy = q_dm, q_pneu_proxy = q_pn,
                         high_tb = h_tb, high_dm_proxy = h_dm, high_pneu_proxy = h_pn,
                         burden_class_proxy = cls,
                         triple_high_proxy = h_tb & h_dm & h_pn)
}
threshold_results <- data.table::rbindlist(lapply(c(0.67, 0.75, 0.80, 0.90), class_at, dt = d))
primary_class <- threshold_results[threshold == 0.75,
  .(location_id, burden_class_proxy, triple_high_proxy)]
d <- merge(d, primary_class, by = "location_id", all.x = TRUE, sort = FALSE)

scaled <- scale(d[, .(tb_asir, t2dm_asir_proxy, pneumoconiosis_asir_proxy)])
diss <- dist(scaled)
k_grid <- 2:min(8, nrow(d) - 1L)
sil <- vapply(k_grid, function(k) {
  fit <- cluster::pam(scaled, k = k)
  mean(cluster::silhouette(fit$clustering, diss)[, "sil_width"])
}, numeric(1))
best_k <- k_grid[which.max(sil)]
pam_fit <- cluster::pam(scaled, k = best_k)
d[, pam_cluster_proxy := factor(pam_fit$clustering)]

model1 <- lm(log(tb_asir) ~ log(t2dm_asir_proxy) + log(pneumoconiosis_asir_proxy), data = d)
model2 <- lm(log(tb_asir) ~ log(t2dm_asir_proxy) * log(pneumoconiosis_asir_proxy), data = d)
model_table <- data.table::rbindlist(list(
  data.table::as.data.table(broom::tidy(model1, conf.int = TRUE))[, model := "M1 additive"],
  data.table::as.data.table(broom::tidy(model2, conf.int = TRUE))[, model := "M2 interaction"]
), fill = TRUE)
model_table[, interpretation_boundary := "Country-level ecological association; not an individual or causal effect"]

ranked <- data.table::copy(d)[order(-cbi_proxy)]
ranked[, cbi_rank := seq_len(.N)]
top_table <- ranked[1:min(30L, .N), .(cbi_rank, location_id, location_name, tb_asir,
  t2dm_asir_proxy, pneumoconiosis_asir_proxy, p_tb, p_dm_proxy, p_pneu_proxy, cbi_proxy,
  burden_class_proxy, pam_cluster_proxy)]

class_counts <- threshold_results[, .(n_countries = .N, n_triple_high = sum(triple_high_proxy)),
                                  by = threshold]
sil_table <- data.table::data.table(k = k_grid, mean_silhouette = sil,
                                    selected = k_grid == best_k)

safe_write_csv(d, file.path(path_data, "01_GBD2023", "processed", "GBD2023_proxy_derived.csv"))
safe_write_csv(d, file.path(path_result, "01_GBD2023", "source_data", "SD02_GBD2023_proxy_derived.csv"))
safe_write_csv(top_table, file.path(path_result, "01_GBD2023", "tables", "T01_top30_combined_burden_proxy.csv"))
safe_write_csv(threshold_results, file.path(path_result, "01_GBD2023", "tables", "T02_threshold_sensitivity_proxy.csv"))
safe_write_csv(class_counts, file.path(path_result, "01_GBD2023", "tables", "T03_triple_high_counts_proxy.csv"))
safe_write_csv(sil_table, file.path(path_result, "01_GBD2023", "models", "M01_PAM_silhouette_proxy.csv"))
safe_write_csv(model_table, file.path(path_result, "01_GBD2023", "models", "M02_ecological_models_proxy.csv"))
saveRDS(list(model1 = model1, model2 = model2, pam = pam_fit),
        file.path(path_result, "01_GBD2023", "models", "M03_proxy_model_objects.rds"))
write_log("GBD proxy analysis completed: k=", best_k, "; triple-high at p75=",
          sum(d$triple_high_proxy))
