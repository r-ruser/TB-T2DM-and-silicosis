source(file.path("R", "00_config.R"), encoding = "UTF-8")
suppressPackageStartupMessages(library(data.table))

cat("=== Simple Donor-Level Validation ===\n\n")

# Load data
cell_annot <- fread(file.path(path_result, "04_scRNA", "GSE174725", "source_data", "SD03_GSE174725_cell_annotation_UMAP.csv"))
cell_annot[, cell_type := trimws(cell_type)]
cell_annot[, donor_label := trimws(donor_label)]
cell_annot[, group := trimws(group)]

pb <- fread(file.path(path_result, "04_scRNA", "GSE174725", "source_data", "SD05_GSE174725_pseudobulk_edgeR.csv"))

# Load signature from meta-analysis
meta <- fread(file.path(path_result, "06_final", "tables", "T12_TBDM_meta_analysis.csv"))
shared_up <- meta[direction == "up"][order(-beta)][1:200, gene_symbol]
shared_down <- meta[direction == "down"][order(beta)][1:200, gene_symbol]

# Match genes
available_genes <- unique(pb$gene_symbol)
shared_up_matched <- shared_up[shared_up %in% available_genes]
shared_down_matched <- shared_down[shared_down %in% available_genes]

cat("Signature: up=", length(shared_up_matched), ", down=", length(shared_down_matched), "\n\n")

# Calculate scores per cell type
cell_types <- unique(pb$cell_type)
donors <- unique(cell_annot$donor_label)

results <- list()

for (ct in cell_types) {
  ct_pb <- pb[cell_type == ct]
  up_genes <- intersect(shared_up_matched, ct_pb$gene_symbol)
  down_genes <- intersect(shared_down_matched, ct_pb$gene_symbol)

  if (length(up_genes) < 3 | length(down_genes) < 3) next

  up_score <- mean(ct_pb[gene_symbol %in% up_genes]$log2FC, na.rm = TRUE)
  down_score <- mean(ct_pb[gene_symbol %in% down_genes]$log2FC, na.rm = TRUE)
  net_score <- up_score - down_score

  for (donor in donors) {
    cells <- cell_annot[cell_type == ct & donor_label == donor]
    n_cells <- nrow(cells)
    if (n_cells > 0) {
      condition <- cells$group[1]
      results[[length(results) + 1]] <- data.table(
        donor = donor,
        condition = condition,
        cell_type = ct,
        n_cells = n_cells,
        n_up = length(up_genes),
        n_down = length(down_genes),
        up_score = up_score,
        down_score = down_score,
        net_score = net_score
      )
    }
  }
}

dt <- rbindlist(results)
fwrite(dt, file.path(path_result, "06_final", "tables", "T22_donor_level_signature_scores_GSE174725.csv"), bom = TRUE)

cat("Results:\n")
print(dt)

# Contrasts
cat("\nSilicosis vs Exposure contrasts:\n")
for (ct in unique(dt$cell_type)) {
  exp <- dt[cell_type == ct & condition == "Silica-exposed control"]
  sil <- dt[cell_type == ct & condition == "Silicosis"]
  if (nrow(exp) > 0 & nrow(sil) > 0) {
    delta <- mean(sil$net_score) - mean(exp$net_score)
    cat(ct, ": delta_net =", round(delta, 4), "\n")
  }
}

write_log("Simple donor-level validation completed")
