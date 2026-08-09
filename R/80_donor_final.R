source(file.path("R", "00_config.R"), encoding = "UTF-8")
suppressPackageStartupMessages(library(data.table))

cat("=== Donor-Level Signature Validation (Final) ===\n\n")

# ── 1. Signature ──
meta <- fread(file.path(path_result, "06_final", "tables", "T12_TBDM_meta_analysis.csv"))
shared_up <- meta[direction == "up"][order(-beta)][1:200, gene_symbol]
shared_down <- meta[direction == "down"][order(beta)][1:200, gene_symbol]
sig_dt <- data.table(gene = c(shared_up, shared_down), direction = c(rep("up", 200), rep("down", 200)))
fwrite(sig_dt, file.path(path_result, "06_final", "tables", "T20_signature_gene_list.csv"), bom = TRUE)
cat("Signature: 200 up + 200 down\n")

# ── 2. GSE174725 ──
cat("\n--- GSE174725 ---\n")
annot <- fread(file.path(path_result, "04_scRNA", "GSE174725", "source_data", "SD03_GSE174725_cell_annotation_UMAP.csv"))
annot[, cell_type := trimws(cell_type)]
annot[, donor_label := trimws(donor_label)]
annot[, group := trimws(group)]

pb <- fread(file.path(path_result, "04_scRNA", "GSE174725", "source_data", "SD05_GSE174725_pseudobulk_edgeR.csv"))
avail <- unique(pb$gene_symbol)
up_m <- shared_up[shared_up %in% avail]
dn_m <- shared_down[shared_down %in% avail]
cat("Matched: up=", length(up_m), "down=", length(dn_m), "\n")

# Coverage
cov <- data.table(signature = c(rep("up", 200), rep("down", 200)),
                  gene = c(shared_up, shared_down),
                  matched = c(shared_up %in% avail, shared_down %in% avail))
fwrite(cov, file.path(path_result, "06_final", "tables", "T21_signature_coverage_GSE174725.csv"), bom = TRUE)

# Build donor×celltype score table using merge
ct_list <- unique(pb$cell_type)
donor_list <- unique(annot$donor_label)

score_grid <- CJ(donor = donor_list, cell_type = ct_list)

# Calculate per-celltype scores from pseudobulk
ct_scores <- pb[, .(
  up_score = mean(log2FC[gene_symbol %in% up_m], na.rm = TRUE),
  down_score = mean(log2FC[gene_symbol %in% dn_m], na.rm = TRUE),
  n_up = sum(gene_symbol %in% up_m),
  n_down = sum(gene_symbol %in% dn_m)
), by = cell_type]

# Count cells per donor×celltype
cell_counts <- annot[, .(n_cells = .N), by = .(donor_label, cell_type, group)]
setnames(cell_counts, "donor_label", "donor")

# Merge
score_grid <- merge(score_grid, ct_scores, by = "cell_type", all.x = TRUE)
score_grid <- merge(score_grid, cell_counts, by = c("donor", "cell_type"), all.x = TRUE)
score_grid[is.na(n_cells), n_cells := 0L]
score_grid[, net_score := up_score - down_score]

fwrite(score_grid, file.path(path_result, "06_final", "tables", "T22_donor_level_signature_scores_GSE174725.csv"), bom = TRUE)
cat("Donor scores:\n")
print(score_grid)

# Contrasts
cat("\nSilicosis vs Exposure contrasts:\n")
contrast_list <- list()
for (ct in ct_list) {
  sub <- score_grid[cell_type == ct & n_cells > 0]
  exp <- sub[group == "Silica-exposed control"]
  sil <- sub[group == "Silicosis"]
  if (nrow(exp) >= 1 & nrow(sil) >= 1) {
    contrast_list[[ct]] <- data.table(
      cell_type = ct,
      n_exposure = nrow(exp), n_silicosis = nrow(sil),
      cells_exposure = mean(exp$n_cells), cells_silicosis = mean(sil$n_cells),
      delta_up = mean(sil$up_score) - mean(exp$up_score),
      delta_down = mean(sil$down_score) - mean(exp$down_score),
      delta_net = mean(sil$net_score) - mean(exp$net_score),
      direction = ifelse(mean(sil$net_score) > mean(exp$net_score), "higher_silicosis", "lower_silicosis")
    )
  }
}
ct_contrasts <- rbindlist(contrast_list)
fwrite(ct_contrasts, file.path(path_result, "06_final", "tables", "T23_GSE174725_celltype_contrasts.csv"), bom = TRUE)
cat("Contrasts:\n"); print(ct_contrasts)

# ── 3. GSE192483 paired gene-set tests ──
cat("\n--- GSE192483 ---\n")
pb19 <- fread(file.path(path_result, "04_scRNA", "GSE192483", "source_data", "SD05_GSE192483_paired_pseudobulk_edgeR.csv"))
avail19 <- unique(pb19$gene_symbol)
up19 <- shared_up[shared_up %in% avail19]
dn19 <- shared_down[shared_down %in% avail19]
cat("Matched: up=", length(up19), "down=", length(dn19), "\n")

# Coverage
cov19 <- data.table(signature = c(rep("up", 200), rep("down", 200)),
                    gene = c(shared_up, shared_down),
                    matched = c(shared_up %in% avail19, shared_down %in% avail19))
fwrite(cov19, file.path(path_result, "06_final", "tables", "T21_signature_coverage_GSE192483.csv"), bom = TRUE)

# Gene-set tests per cell type
gs_list <- list()
for (ct in unique(pb19$cell_type)) {
  ctp <- pb19[cell_type == ct]
  ug <- intersect(up19, ctp$gene_symbol)
  dg <- intersect(dn19, ctp$gene_symbol)
  if (length(ug) < 5 | length(dg) < 5) next

  up_fc <- ctp[gene_symbol %in% ug, log2FC]
  dn_fc <- ctp[gene_symbol %in% dg, log2FC]

  up_t <- t.test(up_fc, mu = 0)
  dn_t <- t.test(dn_fc, mu = 0)
  ps <- p.adjust(c(up_t$p.value, dn_t$p.value), method = "BH")

  gs_list[[ct]] <- data.table(
    cell_type = ct,
    n_up = length(ug), n_down = length(dg),
    mean_logFC_up = mean(up_fc), mean_logFC_down = mean(dn_fc),
    net_logFC = mean(up_fc) - mean(dn_fc),
    pval_up = up_t$p.value, pval_down = dn_t$p.value,
    fdr_up = ps[1], fdr_down = ps[2],
    direction_up = ifelse(mean(up_fc) > 0, "up_in_lesion", "down_in_lesion"),
    direction_down = ifelse(mean(dn_fc) > 0, "up_in_lesion", "down_in_lesion")
  )
}
gs_dt <- rbindlist(gs_list)
fwrite(gs_dt, file.path(path_result, "06_final", "tables", "T23_GSE192483_pseudobulk_gene_set_tests.csv"), bom = TRUE)
cat("Gene-set tests:\n"); print(gs_dt)

cat("\n=== Done ===\n")
write_log("Donor-level validation final completed")
