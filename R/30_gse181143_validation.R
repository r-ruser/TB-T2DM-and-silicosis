source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "DESeq2", "edgeR", "ggplot2", "patchwork", "svglite", "ragg"))

acc <- "GSE181143"

# ============================================================
# 1. Load count data
# ============================================================
count_file <- file.path(path_data, "02_GEO_bulk", acc, "raw", paste0(acc, "_count_data.csv.gz"))
counts_raw <- data.table::fread(count_file)
cat("Loaded counts:", nrow(counts_raw), "genes x", ncol(counts_raw) - 1, "samples\n")

counts <- as.matrix(counts_raw[, -1])
rownames(counts) <- counts_raw$V1
# Remove leading X from column names
colnames(counts) <- sub("^X", "", colnames(counts))

# ============================================================
# 2. Load metadata from pre-existing curated CSV
# ============================================================
meta_file <- file.path(path_result, "02_bulk", acc, paste0("A_", acc, "_sample_metadata_raw.csv"))
meta_raw <- data.table::fread(meta_file)
cat("Metadata rows:", nrow(meta_raw), "\n")

# Build matching key: subject_id + "_" + timepoint + "m"
meta_raw[, count_key := paste0(characteristic_7, "_", timepoint, "m")]

# Standardize group labels
meta_raw[, group := NA_character_]
meta_raw[disease.state == "Diabetes" & tb.status == "TB", group := "TB-DM"]
meta_raw[disease.state == "non Diabetes" & tb.status == "TB", group := "TB-only"]
meta_raw[disease.state == "Diabetes" & tb.status == "non TB", group := "DM-only"]
meta_raw[disease.state == "non Diabetes" & tb.status == "non TB", group := "HC"]

# ============================================================
# 3. Match counts to metadata
# ============================================================
matched <- intersect(colnames(counts), meta_raw$count_key)
cat("Matched samples:", length(matched), "/", ncol(counts), "\n")

meta_matched <- meta_raw[match(matched, count_key)]
counts_matched <- counts[, matched, drop = FALSE]

cat("\nGroup distribution:\n")
print(table(meta_matched$group, useNA = "ifany"))
cat("\nSite distribution:\n")
print(table(meta_matched$site, useNA = "ifany"))
cat("\nTimepoint distribution:\n")
print(table(meta_matched$timepoint, useNA = "ifany"))

# ============================================================
# 4. PRIMARY ANALYSIS: TB-DM vs TB-only at baseline (timepoint 0)
# ============================================================
meta_base <- meta_matched[group %in% c("TB-DM", "TB-only") & timepoint == "0"]
cat("\n=== PRIMARY: TB-DM vs TB-only at baseline ===\n")
cat("Samples:", nrow(meta_base), "\n")
print(table(meta_base$group, meta_base$site))

if (nrow(meta_base) < 10) {
  stop("Too few baseline TB samples for primary contrast")
}

counts_base <- counts_matched[, meta_base$count_key, drop = FALSE]

# Low-expression filter
dge <- edgeR::DGEList(counts = counts_base)
keep <- edgeR::filterByExpr(dge, group = meta_base$group, min.count = 10, min.total.count = 50)
dge <- dge[keep, , keep.lib.sizes = FALSE]
cat("Genes after filtering:", nrow(dge), "\n")

# DESeq2
group_factor <- factor(meta_base$group, levels = c("TB-only", "TB-DM"))
site_factor <- factor(meta_base$site)
coldata <- data.frame(group = group_factor, site = site_factor, row.names = colnames(dge))

# Only include site in design if there are multiple sites
if (length(levels(site_factor)) > 1) {
  design_formula <- ~ site + group
} else {
  design_formula <- ~ group
}
dds <- DESeq2::DESeqDataSetFromMatrix(countData = dge$counts, colData = coldata, design = design_formula)
dds <- DESeq2::DESeq(dds)
res <- DESeq2::results(dds, contrast = c("group", "TB-DM", "TB-only"), alpha = 0.05)

res_dt <- data.table::as.data.table(as.data.frame(res), keep.rownames = "gene_symbol")
if ("rn" %in% colnames(res_dt)) data.table::setnames(res_dt, "rn", "gene_symbol")
res_dt[, fdr := padj]
res_dt[, is_deg := fdr < 0.05]
res_dt[, is_strong := fdr < 0.05 & abs(log2FoldChange) > 0.5]

cat("\nDEG (FDR<0.05):", sum(res_dt$is_deg, na.rm = TRUE), "\n")
cat("Strong DEG (FDR<0.05 & |logFC|>0.5):", sum(res_dt$is_strong, na.rm = TRUE), "\n")

# ============================================================
# 5. Save results
# ============================================================
result_dir <- file.path(path_result, "02_bulk", acc)
invisible(lapply(file.path(result_dir, c("tables", "source_data", "figures")),
                 dir.create, recursive = TRUE, showWarnings = FALSE))
safe_write_csv(res_dt, file.path(result_dir, "tables", paste0("T01_", acc, "_TB_DM_vs_TB.csv")))

# ============================================================
# 6. Compare with discovery (GSE114192)
# ============================================================
discovery_file <- file.path(path_result, "02_bulk", "GSE114192", "tables", "T02_GSE114192_DESeq2_results.csv")
if (file.exists(discovery_file)) {
  discovery <- data.table::fread(discovery_file)
  disc_sig <- discovery[gene_symbol != "" & fdr < 0.05]
  val_sig <- res_dt[gene_symbol != "" & fdr < 0.05]

  merged <- merge(
    disc_sig[, .(gene_symbol, logFC_discovery = log2FC)],
    val_sig[, .(gene_symbol, logFC_validation = log2FoldChange)],
    by = "gene_symbol"
  )

  if (nrow(merged) > 0) {
    concordant <- merged[sign(logFC_discovery) == sign(logFC_validation)]
    cat("\nDirection concordance with discovery:\n")
    cat("Shared significant genes:", nrow(merged), "\n")
    cat("Concordant direction:", nrow(concordant), "(", round(100 * nrow(concordant) / nrow(merged), 1), "%)\n")

    cor_test <- cor.test(merged$logFC_discovery, merged$logFC_validation, method = "spearman")
    cat("Spearman correlation:", round(cor_test$estimate, 3), "; P =", format.pval(cor_test$p.value, digits = 3), "\n")

    # Save concordance
    safe_write_csv(merged, file.path(result_dir, "tables", paste0("T02_", acc, "_discovery_concordance.csv")))
  }
}

# ============================================================
# 7. Site-stratified analysis
# ============================================================
cat("\n=== Site-stratified analysis ===\n")
for (site_name in unique(meta_base$site)) {
  meta_site <- meta_base[site == site_name]
  if (length(unique(meta_site$group)) < 2) {
    cat(site_name, ": only one group, skipping\n")
    next
  }

  counts_site <- counts_base[, meta_site$count_key, drop = FALSE]
  dge_site <- edgeR::DGEList(counts = counts_site)
  keep_site <- edgeR::filterByExpr(dge_site, group = meta_site$group, min.count = 5)
  dge_site <- dge_site[keep_site, , keep.lib.sizes = FALSE]

  group_site <- factor(meta_site$group, levels = c("TB-only", "TB-DM"))
  coldata_site <- data.frame(group = group_site, row.names = colnames(dge_site))

  dds_site <- DESeq2::DESeqDataSetFromMatrix(countData = dge_site$counts, colData = coldata_site, design = ~ group)
  dds_site <- DESeq2::DESeq(dds_site)
  res_site <- DESeq2::results(dds_site, contrast = c("group", "TB-DM", "TB-only"), alpha = 0.05)

  cat(site_name, ": samples =", nrow(meta_site),
      "; DEG(FDR<0.05) =", sum(res_site$padj < 0.05, na.rm = TRUE), "\n")
}

# ============================================================
# 8. Longitudinal: module score trajectories
# ============================================================
cat("\n=== Longitudinal module score trajectories ===\n")
# Load discovery shared module genes
if (file.exists(discovery_file)) {
  discovery <- data.table::fread(discovery_file)
  shared_up <- discovery[gene_symbol != "" & fdr < 0.05 & log2FC > 0, gene_symbol]
  shared_down <- discovery[gene_symbol != "" & fdr < 0.05 & log2FC < 0, gene_symbol]

  shared_up_match <- intersect(shared_up, rownames(counts_matched))
  shared_down_match <- intersect(shared_down, rownames(counts_matched))

  cat("Shared up genes matched:", length(shared_up_match), "\n")
  cat("Shared down genes matched:", length(shared_down_match), "\n")

  if (length(shared_up_match) > 5) {
    # Calculate module scores using log-CPM
    dge_all <- edgeR::DGEList(counts = counts_matched)
    dge_all <- edgeR::calcNormFactors(dge_all)
    logcpm <- edgeR::cpm(dge_all, log = TRUE)

    # Module scores per sample
    meta_matched$module_up <- colMeans(logcpm[shared_up_match, , drop = FALSE], na.rm = TRUE)
    meta_matched$module_down <- colMeans(logcpm[shared_down_match, , drop = FALSE], na.rm = TRUE)

    # Trajectory plot for TB-DM and TB-only
    meta_traj <- meta_matched[group %in% c("TB-DM", "TB-only")]
    meta_traj[, timepoint_num := as.numeric(timepoint)]

    for (module_name in c("module_up", "module_down")) {
      # Compute mean ± SE per group per timepoint
      traj_summary <- meta_traj[, .(
        mean = mean(.SD[[module_name]], na.rm = TRUE),
        se = stats::sd(.SD[[module_name]], na.rm = TRUE) / sqrt(.N),
        n = .N
      ), by = .(group, timepoint, timepoint_num)]

      cat(module_name, "trajectory:\n")
      print(traj_summary)
    }

    safe_write_csv(meta_traj, file.path(result_dir, "source_data", paste0("SD01_", acc, "_module_scores.csv")))
  }
}

# ============================================================
# 9. Summary
# ============================================================
summary_dt <- data.table::data.table(
  metric = c("Total samples", "TB-DM baseline", "TB-only baseline", "Genes tested",
             "DEG (FDR<0.05)", "Strong DEG (FDR<0.05 & |logFC|>0.5)",
             "India samples", "Brazil samples"),
  value = c(nrow(meta_matched), sum(meta_base$group == "TB-DM"), sum(meta_base$group == "TB-only"),
            nrow(res_dt), sum(res_dt$is_deg, na.rm = TRUE), sum(res_dt$is_strong, na.rm = TRUE),
            sum(meta_matched$site == "Chennai, India"), sum(meta_matched$site == "Salvador, Brazil"))
)
safe_write_csv(summary_dt, file.path(result_dir, "tables", paste0("T00_", acc, "_summary.csv")))

write_log(acc, " TB-DM validation completed: samples=", nrow(meta_base),
          "; DEG=", sum(res_dt$is_deg, na.rm = TRUE))
