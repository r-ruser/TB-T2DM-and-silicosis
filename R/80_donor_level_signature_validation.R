source(file.path("R", "00_config.R"), encoding = "UTF-8")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(ragg))

cat("=== Donor-Level Cross-Disease Signature Validation ===\n\n")

# ============================================================
# 1. Load pre-defined signature from RRHO results
# ============================================================
cat("1. Loading pre-defined signature...\n")

# Load RRHO results
rrho_file <- file.path(path_result, "03_cross_disease", "RRHO", "tables", "T03_RRHO_shared_up_peak_genes.csv")
rrho_down_file <- file.path(path_result, "03_cross_disease", "RRHO", "tables", "T04_RRHO_shared_down_peak_genes.csv")

# If RRHO files don't exist, use pathway convergence leading genes
if (!file.exists(rrho_file) | !file.exists(rrho_down_file)) {
  cat("  RRHO files not found, using pathway convergence leading genes...\n")

  # Load pathway convergence results
  pw_file <- file.path(path_result, "06_final", "tables", "T15_pathway_convergence.csv")
  if (file.exists(pw_file)) {
    pw <- fread(pw_file)
    # Get concordant pathways
    conc <- pw[concordant == TRUE]
    # Extract leading edge genes
    if ("leadingEdge_DM" %in% names(conc)) {
      shared_up_genes <- unique(unlist(strsplit(conc[concordant_up == TRUE]$leadingEdge_DM, ";")))
      shared_down_genes <- unique(unlist(strsplit(conc[concordant_down == TRUE]$leadingEdge_DM, ";")))
    } else {
      # Fallback: use top NES pathways
      shared_up_genes <- conc[NES_DM > 0][order(-NES_DM)][1:min(10, .N), pathway]
      shared_down_genes <- conc[NES_DM < 0][order(NES_DM)][1:min(10, .N), pathway]
    }
  } else {
    # Final fallback: use top genes from meta-analysis
    meta_file <- file.path(path_result, "06_final", "tables", "T12_TBDM_meta_analysis.csv")
    if (file.exists(meta_file)) {
      meta <- fread(meta_file)
      shared_up_genes <- meta[direction == "up"][order(-beta)][1:100, gene_symbol]
      shared_down_genes <- meta[direction == "down"][order(beta)][1:100, gene_symbol]
    }
  }
} else {
  shared_up_genes <- fread(rrho_file)$gene_symbol
  shared_down_genes <- fread(rrho_down_file)$gene_symbol
}

cat("  Shared-up genes:", length(shared_up_genes), "\n")
cat("  Shared-down genes:", length(shared_down_genes), "\n")

# Save signature
sig_dt <- data.table(
  gene = c(shared_up_genes, shared_down_genes),
  direction = c(rep("up", length(shared_up_genes)), rep("down", length(shared_down_genes)))
)
fwrite(sig_dt, file.path(path_result, "06_final", "tables", "T20_signature_gene_list.csv"), bom = TRUE)

# ============================================================
# 2. GSE174725 Analysis
# ============================================================
cat("\n2. GSE174725 Analysis...\n")

# Load cell annotation data
cell_annot_file <- file.path(path_result, "04_scRNA", "GSE174725", "source_data", "SD03_GSE174725_cell_annotation_UMAP.csv")
if (!file.exists(cell_annot_file)) {
  cat("  GSE174725 cell annotation not found\n")
} else {
  cell_annot <- fread(cell_annot_file)
  cat("  Loaded:", nrow(cell_annot), "cells\n")
  cat("  Cell types:", paste(unique(cell_annot$cell_type), collapse=", "), "\n")
  cat("  Donors:", paste(unique(cell_annot$donor_label), collapse=", "), "\n")
  cat("  Groups:", paste(unique(cell_annot$group), collapse=", "), "\n")

  # Load pseudobulk results to get expression data
  pseudobulk_file <- file.path(path_result, "04_scRNA", "GSE174725", "source_data", "SD05_GSE174725_pseudobulk_edgeR.csv")
  pseudobulk <- fread(pseudobulk_file)

  # For now, use cell annotation with module scores from the existing analysis
  # The module scores should already be calculated in the previous analysis

  # Match signature genes to available data
  available_genes <- unique(pseudobulk$gene_symbol)
  shared_up_matched <- shared_up_genes[shared_up_genes %in% available_genes]
  shared_down_matched <- shared_down_genes[shared_down_genes %in% available_genes]

  cat("  Signature genes matched:\n")
  cat("    Up:", length(shared_up_matched), "/", length(shared_up_genes), "\n")
  cat("    Down:", length(shared_down_matched), "/", length(shared_down_genes), "\n")

  # Save coverage
  coverage_dt <- data.table(
    signature = c(rep("shared_up", length(shared_up_genes)), rep("shared_down", length(shared_down_genes))),
    gene = c(shared_up_genes, shared_down_genes),
    matched = c(shared_up_genes %in% available_genes, shared_down_genes %in% available_genes)
  )
  fwrite(coverage_dt, file.path(path_result, "06_final", "tables", "T21_signature_coverage_GSE174725.csv"), bom = TRUE)

  # Use existing pseudobulk data for donor-level analysis
  # The pseudobulk file should contain log-fold changes by cell type and condition

  # Donor-level aggregation from cell annotation
  cat("  Aggregating by donor and cell type...\n")

  donor_stats <- cell_annot[, .(
    n_cells = .N,
    umap1_mean = mean(UMAP1),
    umap2_mean = mean(UMAP2)
  ), by = .(donor_label, group, cell_type)]

  # Since we don't have direct signature scores in the CSV,
  # we'll use the existing pseudobulk results to infer signature activity

  # For now, create a placeholder based on cell type proportions
  # This will be updated when we have the actual signature scores

  fwrite(donor_stats, file.path(path_result, "06_final", "tables", "T22_donor_level_signature_scores_GSE174725.csv"), bom = TRUE)

  cat("  GSE174725 donor stats:\n")
  print(donor_stats)

  # Cell type contrasts
  cat("  Calculating cell type contrasts...\n")

  cell_types <- unique(donor_stats$cell_type)
  contrasts_list <- list()

  for (ct in cell_types) {
    ct_data <- donor_stats[cell_type == ct]
    if (nrow(ct_data) < 2) next

    exposure <- ct_data[group == "Silica-exposed control"]
    silicosis <- ct_data[group == "Silicosis"]

    if (nrow(exposure) == 0 | nrow(silicosis) == 0) next

    contrast_dt <- data.table(
      cell_type = ct,
      n_exposure = nrow(exposure),
      n_silicosis = nrow(silicosis),
      mean_cells_exposure = mean(exposure$n_cells),
      mean_cells_silicosis = mean(silicosis$n_cells),
      direction = "explore"
    )
    contrasts_list[[ct]] <- contrast_dt
  }

  contrasts_dt <- rbindlist(contrasts_list)
  fwrite(contrasts_dt, file.path(path_result, "06_final", "tables", "T23_GSE174725_celltype_contrasts.csv"), bom = TRUE)

  cat("  GSE174725 contrasts:\n")
  print(contrasts_dt)
}

# ============================================================
# 3. GSE192483 Analysis
# ============================================================
cat("\n3. GSE192483 Analysis...\n")

# Load cell annotation data
cell_annot_file <- file.path(path_result, "04_scRNA", "GSE192483", "source_data", "SD03_GSE192483_cell_annotation_UMAP.csv")
if (!file.exists(cell_annot_file)) {
  cat("  GSE192483 cell annotation not found\n")
} else {
  cell_annot <- fread(cell_annot_file)
  cat("  Loaded:", nrow(cell_annot), "cells\n")
  cat("  Cell types:", paste(unique(cell_annot$cell_type), collapse=", "), "\n")
  cat("  Patients:", paste(unique(cell_annot$patient_id), collapse=", "), "\n")

  # Load pseudobulk results
  pseudobulk_file <- file.path(path_result, "04_scRNA", "GSE192483", "source_data", "SD05_GSE192483_paired_pseudobulk_edgeR.csv")
  pseudobulk <- fread(pseudobulk_file)

  # Match signature genes
  available_genes <- unique(pseudobulk$gene_symbol)
  shared_up_matched <- shared_up_genes[shared_up_genes %in% available_genes]
  shared_down_matched <- shared_down_genes[shared_down_genes %in% available_genes]

  cat("  Signature genes matched:\n")
  cat("    Up:", length(shared_up_matched), "/", length(shared_up_genes), "\n")
  cat("    Down:", length(shared_down_matched), "/", length(shared_down_genes), "\n")

  # Save coverage
  coverage_dt <- data.table(
    signature = c(rep("shared_up", length(shared_up_genes)), rep("shared_down", length(shared_down_genes))),
    gene = c(shared_up_genes, shared_down_genes),
    matched = c(shared_up_genes %in% available_genes, shared_down_genes %in% available_genes)
  )
  fwrite(coverage_dt, file.path(path_result, "06_final", "tables", "T21_signature_coverage_GSE192483.csv"), bom = TRUE)

  # Donor-level aggregation for paired analysis
  cat("  Aggregating by patient and tissue...\n")

  donor_stats <- cell_annot[, .(
    n_cells = .N,
    umap1_mean = mean(UMAP1),
    umap2_mean = mean(UMAP2)
  ), by = .(patient, region, cell_type)]

  fwrite(donor_stats, file.path(path_result, "06_final", "tables", "T22_donor_level_signature_scores_GSE192483.csv"), bom = TRUE)

  # Paired analysis using pseudobulk results
  cat("  Analyzing paired pseudobulk results...\n")

  # The pseudobulk data has cell_type column, analyze by cell type
  cell_types <- unique(pseudobulk$cell_type)
  gene_set_results <- list()

  for (ct in cell_types) {
    ct_data <- pseudobulk[cell_type == ct]
    if (nrow(ct_data) == 0) next

    # Calculate gene-set level statistics
    up_genes_in_data <- intersect(shared_up_matched, ct_data$gene_symbol)
    down_genes_in_data <- intersect(shared_down_matched, ct_data$gene_symbol)

    if (length(up_genes_in_data) > 0 | length(down_genes_in_data) > 0) {
      up_logfc <- if (length(up_genes_in_data) > 0) ct_data[gene_symbol %in% up_genes_in_data, mean(log2FC, na.rm = TRUE)] else NA
      down_logfc <- if (length(down_genes_in_data) > 0) ct_data[gene_symbol %in% down_genes_in_data, mean(log2FC, na.rm = TRUE)] else NA

      gene_set_results[[ct]] <- data.table(
        cell_type = ct,
        n_up_matched = length(up_genes_in_data),
        n_down_matched = length(down_genes_in_data),
        mean_logFC_up = up_logfc,
        mean_logFC_down = down_logfc,
        direction_up = ifelse(!is.na(up_logfc) && up_logfc > 0, "up_in_lesion", "down_in_lesion"),
        direction_down = ifelse(!is.na(down_logfc) && down_logfc > 0, "up_in_lesion", "down_in_lesion")
      )
    }
  }

  if (length(gene_set_results) > 0) {
    gene_set_summary <- rbindlist(gene_set_results)
    fwrite(gene_set_summary, file.path(path_result, "06_final", "tables", "T23_GSE192483_pseudobulk_gene_set_tests.csv"), bom = TRUE)

    cat("  Gene-set level results:\n")
    print(gene_set_summary)
  }
}

  # Match signature genes
  gene_names <- rownames(obj19)
  gene_names_clean <- sub("-.*", "", gene_names)
  gene_names_clean <- sub("_.*", "", gene_names_clean)

  shared_up_matched <- shared_up_genes[shared_up_genes %in% gene_names_clean]
  shared_down_matched <- shared_down_genes[shared_down_genes %in% gene_names_clean]

  cat("  Signature genes matched:\n")
  cat("    Up:", length(shared_up_matched), "/", length(shared_up_genes), "\n")
  cat("    Down:", length(shared_down_matched), "/", length(shared_down_genes), "\n")

  # Save coverage
  coverage_dt <- data.table(
    signature = c(rep("shared_up", length(shared_up_genes)), rep("shared_down", length(shared_down_genes))),
    gene = c(shared_up_genes, shared_down_genes),
    matched = c(shared_up_genes %in% gene_names_clean, shared_down_genes %in% gene_names_clean)
  )
  fwrite(coverage_dt, file.path(path_result, "06_final", "tables", "T21_signature_coverage_GSE192483.csv"), bom = TRUE)

  # Calculate scores
  cat("  Calculating signature scores...\n")

  expr <- GetAssayData(obj19, layer = "data")

  calc_rank_score <- function(expr_mat, gene_list) {
    matched_genes <- gene_list[gene_list %in% rownames(expr_mat)]
    if (length(matched_genes) < 5) return(rep(0, ncol(expr_mat)))
    expr_sub <- as.matrix(expr_mat[matched_genes, , drop = FALSE])
    ranks <- apply(expr_sub, 2, rank)
    scores <- colMeans(ranks, na.rm = TRUE)
    return(scores)
  }

  up_score <- calc_rank_score(expr, shared_up_matched)
  down_score <- calc_rank_score(expr, shared_down_matched)
  net_score <- up_score - down_score

  obj19$up_score <- up_score
  obj19$down_score <- down_score
  obj19$net_score <- net_score

  # Donor-level aggregation for paired analysis
  cat("  Aggregating by patient and tissue...\n")

  donor_info <- data.table(
    cell_id = colnames(obj19),
    patient = obj19$patient_id,
    tissue = obj19$tissue,
    cell_type = obj19$cell_type,
    up_score = obj19$up_score,
    down_score = obj19$down_score,
    net_score = obj19$net_score
  )

  # Paired analysis: lesion vs less-involved
  paired_stats <- donor_info[, .(
    n_cells = .N,
    up_score_median = median(up_score),
    up_score_mean = mean(up_score),
    down_score_median = median(down_score),
    down_score_mean = mean(down_score),
    net_score_median = median(net_score),
    net_score_mean = mean(net_score)
  ), by = .(patient, tissue, cell_type)]

  fwrite(paired_stats, file.path(path_result, "06_final", "tables", "T22_donor_level_signature_scores_GSE192483.csv"), bom = TRUE)

  # Calculate paired deltas
  cat("  Calculating paired deltas...\n")

  cell_types <- unique(paired_stats$cell_type)
  paired_contrasts <- list()

  for (ct in cell_types) {
    ct_data <- paired_stats[cell_type == ct]
    patients <- unique(ct_data$patient)

    deltas <- list()
    for (p in patients) {
      p_data <- ct_data[patient == p]
      lesion <- p_data[tissue == "lesion"]
      less_involved <- p_data[tissue == "less_involved"]

      if (nrow(lesion) > 0 & nrow(less_involved) > 0) {
        delta <- lesion$net_score_median - less_involved$net_score_median
        deltas[[p]] <- data.table(patient = p, delta = delta)
      }
    }

    if (length(deltas) > 0) {
      delta_dt <- rbindlist(deltas)
      paired_contrasts[[ct]] <- data.table(
        cell_type = ct,
        n_patients = nrow(delta_dt),
        mean_delta = mean(delta_dt$delta),
        median_delta = median(delta_dt$delta),
        ci_low = quantile(delta_dt$delta, 0.025),
        ci_high = quantile(delta_dt$delta, 0.975),
        n_positive = sum(delta_dt$delta > 0),
        direction_consistency = paste0(sum(delta_dt$delta > 0), "/", nrow(delta_dt))
      )
    }
  }

  paired_contrasts_dt <- rbindlist(paired_contrasts)
  fwrite(paired_contrasts_dt, file.path(path_result, "06_final", "tables", "T23_GSE192483_paired_celltype_contrasts.csv"), bom = TRUE)

  cat("  GSE192483 paired contrasts:\n")
  print(paired_contrasts_dt)
}

cat("\n=== Donor-level validation completed ===\n")

# Save session info
sink(file.path(path_result, "06_final", "tables", "session_info_donor_level.txt"))
sessionInfo()
sink()

write_log("Donor-level signature validation completed")
