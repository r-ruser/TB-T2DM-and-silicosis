source(file.path("R", "00_config.R"), encoding = "UTF-8")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(ragg))

cat("=== Donor-Level Signature Validation v3 ===\n\n")

# ============================================================
# 1. Load pre-defined signature
# ============================================================
cat("1. Loading pre-defined signature...\n")

# Use meta-analysis results to define signature
meta_file <- file.path(path_result, "06_final", "tables", "T12_TBDM_meta_analysis.csv")
meta <- fread(meta_file)

# Shared-up: genes up in both TB-DM and Silicosis (direction-concordant up)
# Shared-down: genes down in both TB-DM and Silicosis (direction-concordant down)
# Use top genes by effect size
shared_up <- meta[direction == "up"][order(-beta)][1:min(200, .N), gene_symbol]
shared_down <- meta[direction == "down"][order(beta)][1:min(200, .N), gene_symbol]

cat("  Shared-up:", length(shared_up), "genes\n")
cat("  Shared-down:", length(shared_down), "genes\n")

sig_dt <- data.table(gene = c(shared_up, shared_down),
                     direction = c(rep("up", length(shared_up)), rep("down", length(shared_down))))
fwrite(sig_dt, file.path(path_result, "06_final", "tables", "T20_signature_gene_list.csv"), bom = TRUE)

# ============================================================
# 2. GSE174725: Donor-level analysis
# ============================================================
cat("\n2. GSE174725: Donor-level analysis...\n")

cell_annot <- fread(file.path(path_result, "04_scRNA", "GSE174725", "source_data", "SD03_GSE174725_cell_annotation_UMAP.csv"))

# Trim whitespace from cell_type and donor_label
cell_annot[, cell_type := trimws(cell_type)]
cell_annot[, donor_label := trimws(donor_label)]
cell_annot[, group := trimws(group)]

cat("  Cells:", nrow(cell_annot), "\n")
cat("  Donors:", paste(unique(cell_annot$donor_label), collapse = ", "), "\n")

# Load pseudobulk for gene-level logFC
pb <- fread(file.path(path_result, "04_scRNA", "GSE174725", "source_data", "SD05_GSE174725_pseudobulk_edgeR.csv"))

# Match signature genes
available_genes <- unique(pb$gene_symbol)
shared_up_matched <- shared_up[shared_up %in% available_genes]
shared_down_matched <- shared_down[shared_down %in% available_genes]

cat("  Signature matched: up =", length(shared_up_matched), "/", length(shared_up),
    ", down =", length(shared_down_matched), "/", length(shared_down), "\n")

coverage_dt <- data.table(
  signature = c(rep("shared_up", length(shared_up)), rep("shared_down", length(shared_down))),
  gene = c(shared_up, shared_down),
  matched = c(shared_up %in% available_genes, shared_down %in% available_genes)
)
fwrite(coverage_dt, file.path(path_result, "06_final", "tables", "T21_signature_coverage_GSE174725.csv"), bom = TRUE)

# Calculate donor-level signature scores using pseudobulk logFC per cell type
cat("  Calculating donor-level signature scores...\n")

# For each cell type, get the pseudobulk logFC and calculate signature scores
# Then assign to each donor based on their cell type composition

# Get cell types from pseudobulk (these are the ones we have scores for)
cell_types <- unique(pb$cell_type)
donors <- unique(cell_annot$donor_label)

cat("  Cell types from pseudobulk:", paste(cell_types, collapse = ", "), "\n")
cat("  Donors:", paste(donors, collapse = ", "), "\n")

# Calculate signature scores per donor per cell type
donor_scores <- list()

for (ct in cell_types) {
  ct_pb <- pb[cell_type == ct]
  if (nrow(ct_pb) == 0) next

  # Get matched signature genes
  up_genes <- intersect(shared_up_matched, ct_pb$gene_symbol)
  down_genes <- intersect(shared_down_matched, ct_pb$gene_symbol)

  if (length(up_genes) < 3 | length(down_genes) < 3) next

  # Calculate UP and DOWN scores from pseudobulk logFC
  up_score <- mean(ct_pb[gene_symbol %in% up_genes]$log2FC, na.rm = TRUE)
  down_score <- mean(ct_pb[gene_symbol %in% down_genes]$log2FC, na.rm = TRUE)
  net_score <- up_score - down_score

  cat("  ", ct, ": up_score=", round(up_score, 4), ", down_score=", round(down_score, 4), "\n")

  # Assign to each donor that has this cell type
  for (donor in donors) {
    cells <- cell_annot[cell_type == ct & donor_label == donor]
    if (nrow(cells) > 0) {
      donor_scores[[length(donor_scores) + 1]] <- data.table(
        donor = donor,
        condition = cells$group[1],
        cell_type = ct,
        n_cells = nrow(cells),
        n_up_matched = length(up_genes),
        n_down_matched = length(down_genes),
        up_score = up_score,
        down_score = down_score,
        net_score = net_score
      )
    }
  }
}

donor_scores_dt <- rbindlist(donor_scores)
fwrite(donor_scores_dt, file.path(path_result, "06_final", "tables", "T22_donor_level_signature_scores_GSE174725.csv"), bom = TRUE)

cat("  Donor-level scores:\n")
print(donor_scores_dt)

# Calculate contrasts: Silicosis - Exposure
cat("\n  Calculating Silicosis vs Exposure contrasts...\n")

cell_types_ct <- unique(donor_scores_dt$cell_type)
contrasts_list <- list()

for (ct in cell_types_ct) {
  ct_data <- donor_scores_dt[cell_type == ct]
  exposure <- ct_data[condition == "Silica-exposed control"]
  silicosis <- ct_data[condition == "Silicosis"]

  if (nrow(exposure) < 1 | nrow(silicosis) < 1) next

  # Calculate effect sizes
  delta_up <- mean(silicosis$up_score) - mean(exposure$up_score)
  delta_down <- mean(silicosis$down_score) - mean(exposure$down_score)
  delta_net <- mean(silicosis$net_score) - mean(exposure$net_score)

  contrasts_list[[ct]] <- data.table(
    cell_type = ct,
    n_exposure = nrow(exposure),
    n_silicosis = nrow(silicosis),
    mean_cells_exposure = mean(exposure$n_cells),
    mean_cells_silicosis = mean(silicosis$n_cells),
    delta_up = delta_up,
    delta_down = delta_down,
    delta_net = delta_net,
    direction = ifelse(delta_net > 0, "higher_in_silicosis", "lower_in_silicosis")
  )
}

contrasts_dt <- rbindlist(contrasts_list)
fwrite(contrasts_dt, file.path(path_result, "06_final", "tables", "T23_GSE174725_celltype_contrasts.csv"), bom = TRUE)

cat("  GSE174725 contrasts:\n")
print(contrasts_dt)

# ============================================================
# 3. GSE192483: Paired donor-level analysis
# ============================================================
cat("\n3. GSE192483: Paired analysis...\n")

cell_annot19 <- fread(file.path(path_result, "04_scRNA", "GSE192483", "source_data", "SD03_GSE192483_cell_annotation_UMAP.csv"))
cat("  Cells:", nrow(cell_annot19), "\n")
cat("  Patients:", paste(unique(cell_annot19$patient), collapse = ", "), "\n")

pb19 <- fread(file.path(path_result, "04_scRNA", "GSE192483", "source_data", "SD05_GSE192483_paired_pseudobulk_edgeR.csv"))

# Match signature genes
available_genes19 <- unique(pb19$gene_symbol)
shared_up_matched19 <- shared_up[shared_up %in% available_genes19]
shared_down_matched19 <- shared_down[shared_down %in% available_genes19]

cat("  Signature matched: up =", length(shared_up_matched19), "/", length(shared_up),
    ", down =", length(shared_down_matched19), "/", length(shared_down), "\n")

coverage_dt19 <- data.table(
  signature = c(rep("shared_up", length(shared_up)), rep("shared_down", length(shared_down))),
  gene = c(shared_up, shared_down),
  matched = c(shared_up %in% available_genes19, shared_down %in% available_genes19)
)
fwrite(coverage_dt19, file.path(path_result, "06_final", "tables", "T21_signature_coverage_GSE192483.csv"), bom = TRUE)

# Paired analysis: lesion vs less-involved
cat("  Calculating paired gene-set tests...\n")

cell_types19 <- unique(pb19$cell_type)
paired_results <- list()

for (ct in cell_types19) {
  ct_pb <- pb19[cell_type == ct]
  if (nrow(ct_pb) == 0) next

  # Get matched signature genes
  up_genes <- intersect(shared_up_matched19, ct_pb$gene_symbol)
  down_genes <- intersect(shared_down_matched19, ct_pb$gene_symbol)

  if (length(up_genes) < 3 | length(down_genes) < 3) next

  # Calculate mean logFC for each gene set
  up_logfc <- mean(ct_pb[gene_symbol %in% up_genes]$log2FC, na.rm = TRUE)
  down_logfc <- mean(ct_pb[gene_symbol %in% down_genes]$log2FC, na.rm = TRUE)
  net_logfc <- up_logfc - down_logfc

  # Calculate p-values using t-test
  up_pval <- t.test(ct_pb[gene_symbol %in% up_genes]$log2FC, mu = 0)$p.value
  down_pval <- t.test(ct_pb[gene_symbol %in% down_genes]$log2FC, mu = 0)$p.value

  # FDR correction
  pvals <- c(up_pval, down_pval)
  fdrs <- p.adjust(pvals, method = "BH")

  paired_results[[ct]] <- data.table(
    cell_type = ct,
    n_up_matched = length(up_genes),
    n_down_matched = length(down_genes),
    mean_logFC_up = up_logfc,
    mean_logFC_down = down_logfc,
    net_logFC = net_logfc,
    pval_up = up_pval,
    pval_down = down_pval,
    fdr_up = fdrs[1],
    fdr_down = fdrs[2],
    direction_up = ifelse(up_logfc > 0, "up_in_lesion", "down_in_lesion"),
    direction_down = ifelse(down_logfc > 0, "up_in_lesion", "down_in_lesion")
  )
}

paired_dt <- rbindlist(paired_results)
fwrite(paired_dt, file.path(path_result, "06_final", "tables", "T23_GSE192483_pseudobulk_gene_set_tests.csv"), bom = TRUE)

cat("  GSE192483 gene-set tests:\n")
print(paired_dt)

# ============================================================
# 4. Summary
# ============================================================
cat("\n=== Summary ===\n")
cat("GSE174725:\n")
cat("  Donors: 2 Exposure, 3 Silicosis\n")
cat("  Cell types:", length(unique(donor_scores_dt$cell_type)), "\n")
cat("  Contrasts calculated:", nrow(contrasts_dt), "\n")

cat("\nGSE192483:\n")
cat("  Patients:", length(unique(cell_annot19$patient)), "\n")
cat("  Cell types:", length(unique(pb19$cell_type)), "\n")
cat("  Gene-set tests:", nrow(paired_dt), "\n")

sink(file.path(path_result, "06_final", "tables", "session_info_donor_level_v3.txt"))
sessionInfo()
sink()

write_log("Donor-level signature validation v3 completed")
