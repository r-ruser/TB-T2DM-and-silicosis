source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "limma", "GEOquery", "ggplot2", "patchwork", "svglite", "ragg"))
suppressPackageStartupMessages(library(GEOquery))

acc <- "GSE249102"

# Load ExpressionSet
geo_file <- file.path(path_data, "02_GEO_bulk", acc, "metadata", paste0(acc, "_series_matrix.txt.gz"))
eset <- getGEO(filename = geo_file, GSEMatrix = TRUE, getGPL = TRUE)
if (is.list(eset)) eset <- eset[[1]]
cat("Loaded:", acc, "; samples:", ncol(eset), "; genes:", nrow(eset), "\n")

# Metadata
meta <- data.table::as.data.table(Biobase::pData(eset))
meta[, group := `group:ch1`]
cat("\nGroup distribution:\n")
print(table(meta$group))

# Expression matrix (already log2 transformed for microarray)
expr <- Biobase::exprs(eset)

# ============================================================
# 1. Glycaemic gradient: CTRL -> DM2 -> PDM2
# ============================================================
cat("\n=== Glycaemic gradient: CTRL -> DM2 -> PDM2 ===\n")
meta_grad <- meta[group %in% c("CTRL", "DM2", "PDM2")]
cat("Gradient samples:", nrow(meta_grad), "\n")

expr_grad <- expr[, meta_grad$geo_accession, drop = FALSE]
group_grad <- factor(meta_grad$group, levels = c("CTRL", "DM2", "PDM2"))
group_num <- as.numeric(group_grad)

design <- model.matrix(~ group_num)
fit <- limma::lmFit(expr_grad, design)
fit <- limma::eBayes(fit, robust = TRUE, trend = TRUE)

res_grad <- limma::topTable(fit, coef = "group_num", number = Inf, sort.by = "none")
res_grad_dt <- data.table::as.data.table(res_grad, keep.rownames = "probe_id")
res_grad_dt[, fdr := stats::p.adjust(P.Value, method = "BH")]

result_dir <- file.path(path_result, "02_bulk", acc)
dir.create(file.path(result_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
safe_write_csv(res_grad_dt, file.path(result_dir, "tables", paste0("T01_", acc, "_glycaemic_trend.csv")))

cat("Glycaemic trend DEG (FDR<0.05):", sum(res_grad_dt$fdr < 0.05, na.rm = TRUE), "\n")

# ============================================================
# 2. TB-DM vs TB
# ============================================================
cat("\n=== TB-DM vs TB ===\n")
meta_tbdm <- meta[group %in% c("TB", "TBDM")]
cat("TB-DM vs TB samples:", nrow(meta_tbdm), "\n")

expr_tbdm <- expr[, meta_tbdm$geo_accession, drop = FALSE]
group_tbdm <- factor(meta_tbdm$group, levels = c("TB", "TBDM"))

design2 <- model.matrix(~ group_tbdm)
fit2 <- limma::lmFit(expr_tbdm, design2)
fit2 <- limma::eBayes(fit2, robust = TRUE, trend = TRUE)

res_tbdm <- limma::topTable(fit2, coef = "group_tbdmTBDM", number = Inf, sort.by = "none")
res_tbdm_dt <- data.table::as.data.table(res_tbdm, keep.rownames = "probe_id")
res_tbdm_dt[, fdr := stats::p.adjust(P.Value, method = "BH")]

safe_write_csv(res_tbdm_dt, file.path(result_dir, "tables", paste0("T02_", acc, "_TBDM_vs_TB.csv")))

cat("TB-DM vs TB DEG (FDR<0.05):", sum(res_tbdm_dt$fdr < 0.05, na.rm = TRUE), "\n")

# ============================================================
# 3. Module score analysis
# ============================================================
cat("\n=== Module score analysis ===\n")
discovery_file <- file.path(path_result, "02_bulk", "GSE114192", "tables", "T02_GSE114192_DESeq2_results.csv")
if (file.exists(discovery_file)) {
  discovery <- data.table::fread(discovery_file)
  shared_up <- discovery[gene_symbol != "" & fdr < 0.05 & log2FC > 0, gene_symbol]
  shared_down <- discovery[gene_symbol != "" & fdr < 0.05 & log2FC < 0, gene_symbol]

  # Match to platform probes
  fdata <- data.table::as.data.table(Biobase::fData(eset))
  gene_col <- grep("Gene.Symbol|gene_symbol|Symbol|GENE_SYMBOL", colnames(fdata), value = TRUE, ignore.case = TRUE)
  if (length(gene_col) > 0) {
    fdata[, matched_gene := fdata[[gene_col[1]]]]
  } else {
    fdata[, matched_gene := rownames(fdata)]
  }

  rownames(expr) <- make.names(rownames(expr))
  shared_up_match <- intersect(shared_up, fdata$matched_gene)
  shared_down_match <- intersect(shared_down, fdata$matched_gene)

  cat("Shared up genes matched:", length(shared_up_match), "\n")
  cat("Shared down genes matched:", length(shared_down_match), "\n")

  if (length(shared_up_match) > 3) {
    for (grp in c("CTRL", "DM2", "PDM2", "TB", "TBDM")) {
      meta_sub <- meta[group == grp]
      expr_sub <- expr[, meta_sub$geo_accession, drop = FALSE]
      up_score <- colMeans(expr_sub[shared_up_match, , drop = FALSE], na.rm = TRUE)
      down_score <- colMeans(expr_sub[shared_down_match, , drop = FALSE], na.rm = TRUE)
      cat(grp, ": up module =", round(mean(up_score), 2), "±", round(sd(up_score), 2),
          "; down module =", round(mean(down_score), 2), "±", round(sd(down_score), 2), "\n")
    }
  }
}

# ============================================================
# 4. Summary
# ============================================================
summary_dt <- data.table::data.table(
  metric = c("Total samples", "CTRL", "DM2", "PDM2", "TB", "TBDM",
             "Glycaemic trend DEG (FDR<0.05)", "TB-DM vs TB DEG (FDR<0.05)"),
  value = c(nrow(meta), sum(meta$group == "CTRL"), sum(meta$group == "DM2"),
            sum(meta$group == "PDM2"), sum(meta$group == "TB"), sum(meta$group == "TBDM"),
            sum(res_grad_dt$fdr < 0.05, na.rm = TRUE), sum(res_tbdm_dt$fdr < 0.05, na.rm = TRUE))
)
safe_write_csv(summary_dt, file.path(result_dir, "tables", paste0("T00_", acc, "_summary.csv")))

write_log(acc, " glycaemic gradient completed: gradient_DEG=", sum(res_grad_dt$fdr < 0.05, na.rm = TRUE),
          "; TBDM_vs_TB_DEG=", sum(res_tbdm_dt$fdr < 0.05, na.rm = TRUE))
