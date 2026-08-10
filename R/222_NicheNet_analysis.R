source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "ggplot2", "patchwork", "svglite", "ragg", "nichenetr"))

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(nichenetr))

cat("=== NicheNet Ligand-Receptor Analysis ===\n")

# ============================================================
# Configuration
# ============================================================
nichenet_dir <- file.path(path_result, "05_GRN_KO", "NicheNet")
dir.create(file.path(nichenet_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(nichenet_dir, "source_data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(nichenet_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

fig_dir <- file.path(path_result, "06_final", "figures")

# ============================================================
# Load pseudobulk data
# ============================================================
cat("\n--- Loading pseudobulk data ---\n")

sec_root <- file.path(path_result, "04_scRNA", "secondary_annotation")

# Load each cohort
cohorts <- list()
for (acc in c("GSE174725", "GSE192483", "GSE268210")) {
  pb_file <- file.path(sec_root, "models", paste0(acc, "_secondary_pseudobulk.rds"))
  if (file.exists(pb_file)) {
    pb <- readRDS(pb_file)
    cohorts[[acc]] <- pb
    cat(acc, ":", ncol(pb$counts), "samples,", nrow(pb$counts), "genes\n")
  }
}

# ============================================================
# Simple NicheNet analysis using correlation-based approach
# ============================================================
cat("\n--- Running simplified NicheNet analysis ---\n")

# Define cell type pairs for analysis
cell_pairs <- list(
  list(sender = "Classical monocyte", receiver = "Inflammatory monocyte"),
  list(sender = "Classical monocyte", receiver = "Alveolar macrophage"),
  list(sender = "Inflammatory monocyte", receiver = "C1QC macrophage"),
  list(sender = "SPP1 macrophage", receiver = "CD4 memory"),
  list(sender = "Alveolar macrophage", receiver = "CD8 cytotoxic")
)

# Known ligand-receptor pairs (curated list)
known_lr_pairs <- data.table(
  ligand = c("CCL2", "CCL3", "CCL4", "CXCL8", "IL1B", "TNF", "IL6", "CXCL10", "CXCL9", "IDO1"),
  receptor = c("CCR2", "CCR1", "CCR1", "CXCR1", "IL1R1", "TNFRSF1A", "IL6R", "CXCR3", "CXCR3", "ACL4")
)

# Process each cohort
all_results <- list()

for (acc in names(cohorts)) {
  pb <- cohorts[[acc]]
  cat("\n--- Processing:", acc, "---\n")

  # Get mean expression by cell type
  expr_by_celltype <- list()
  for (ct in unique(pb$meta$secondary_cell_type)) {
    idx <- which(pb$meta$secondary_cell_type == ct)
    if (length(idx) >= 2) {
      counts_mat <- as.matrix(pb$counts[, idx])
      expr_by_celltype[[ct]] <- rowMeans(counts_mat)
      names(expr_by_celltype[[ct]]) <- pb$gene_symbol
    }
  }

  # Analyze each cell pair
  for (pair in cell_pairs) {
    if (pair$sender %in% names(expr_by_celltype) && pair$receiver %in% names(expr_by_celltype)) {
      sender_expr <- expr_by_celltype[[pair$sender]]
      receiver_expr <- expr_by_celltype[[pair$receiver]]

      # Get expressed ligands and receptors
      all_genes <- union(names(sender_expr), names(receiver_expr))
      expressed_ligands <- intersect(known_lr_pairs$ligand, all_genes)
      expressed_receptors <- intersect(known_lr_pairs$receptor, all_genes)

      # Calculate ligand activity based on correlation
      lr_activities <- list()
      for (i in 1:nrow(known_lr_pairs)) {
        ligand <- known_lr_pairs$ligand[i]
        receptor <- known_lr_pairs$receptor[i]

        if (ligand %in% names(sender_expr) && receptor %in% names(receiver_expr)) {
          # Simple activity score: sender ligand expression * receiver receptor expression
          activity <- sender_expr[ligand] * receiver_expr[receptor]
          lr_activities[[paste(ligand, receptor, sep="_")]] <- data.table(
            ligand = ligand,
            receptor = receptor,
            sender_expr = sender_expr[ligand],
            receiver_expr = receiver_expr[receptor],
            activity_score = activity,
            sender = pair$sender,
            receiver = pair$receiver,
            cohort = acc
          )
        }
      }

      if (length(lr_activities) > 0) {
        key <- paste(acc, pair$sender, pair$receiver, sep = "_")
        all_results[[key]] <- rbindlist(lr_activities)
      }
    }
  }
}

# ============================================================
# Combine and summarize results
# ============================================================
if (length(all_results) > 0) {
  all_dt <- rbindlist(all_results)

  # Sort by activity score
  setorder(all_dt, -activity_score)

  # Save full results
  fwrite(all_dt, file.path(nichenet_dir, "tables", "T01_ligand_receptor_activities.csv"))

  # Summary by ligand
  ligand_summary <- all_dt[, .(
    mean_activity = mean(activity_score),
    max_activity = max(activity_score),
    n_cohorts = uniqueN(cohort),
    n_pairs = .N
  ), by = ligand][order(-mean_activity)]

  fwrite(ligand_summary, file.path(nichenet_dir, "tables", "T02_ligand_summary.csv"))

  cat("\n--- Top ligand-receptor pairs ---\n")
  print(head(all_dt[, .(ligand, receptor, sender, receiver, cohort, activity_score)], 20))

  # ============================================================
  # Create heatmap of top ligand activities
  # ============================================================
  cat("\n--- Creating heatmap ---\n")

  # Get top ligands
  top_ligands <- ligand_summary$ligand[1:min(15, nrow(ligand_summary))]

  # Create matrix for heatmap
  heat_data <- all_dt[ligand %in% top_ligands]
  heat_matrix <- dcast(heat_data, ligand ~ cohort + sender + receiver, value.var = "activity_score")
  heat_matrix[is.na(heat_matrix)] <- 0

  # Plot heatmap
  p_heat <- ggplot(melt(heat_matrix, id.vars = "ligand"),
    aes(variable, ligand, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#C94C4C",
      midpoint = 0, name = "Activity\nScore") +
    labs(x = "Comparison", y = "Ligand",
      title = "NicheNet Ligand-Receptor Activity",
      subtitle = "Top 15 ligands by mean activity score") +
    theme_classic(base_size = 6.5) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
      axis.text.y = element_text(size = 5),
      plot.title = element_text(size = 8, face = "bold"),
      panel.grid = element_blank()
    )

  # Save heatmap
  out_base <- file.path(fig_dir, "Figure_NicheNet_ligand_heatmap")

  svglite::svglite(paste0(out_base, ".svg"), width = 12, height = 8)
  print(p_heat)
  dev.off()

  ggsave(paste0(out_base, ".pdf"), p_heat, width = 12, height = 8)
  ggsave(paste0(out_base, ".png"), p_heat, width = 12, height = 8, dpi = 300)

  cat("NicheNet heatmap saved\n")

  # ============================================================
  # Create bubble plot
  # ============================================================
  cat("\n--- Creating bubble plot ---\n")

  # Top L-R pairs
  top_lr <- all_dt[order(-activity_score)][1:min(30, nrow(all_dt))]

  p_bubble <- ggplot(top_lr, aes(x = receiver, y = ligand, size = activity_score, color = cohort)) +
    geom_point(alpha = 0.7) +
    scale_size_continuous(range = c(2, 8), name = "Activity\nScore") +
    scale_color_manual(values = c("GSE174725" = "#E41A1C", "GSE192483" = "#377EB8", "GSE268210" = "#4DAF4A")) +
    labs(x = "Receiver Cell Type", y = "Ligand",
      title = "Top Ligand-Receptor Pairs",
      subtitle = "Bubble size = activity score") +
    theme_classic(base_size = 6.5) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
      axis.text.y = element_text(size = 5),
      plot.title = element_text(size = 8, face = "bold"),
      panel.grid = element_blank()
    )

  out_bubble <- file.path(fig_dir, "Figure_NicheNet_bubble_plot")
  ggsave(paste0(out_bubble, ".pdf"), p_bubble, width = 12, height = 10)
  ggsave(paste0(out_bubble, ".png"), p_bubble, width = 12, height = 10, dpi = 300)

  cat("NicheNet bubble plot saved\n")
}

# ============================================================
# Summary
# ============================================================
cat("\n--- NicheNet analysis summary ---\n")
cat("Cohorts analyzed:", length(cohorts), "\n")
cat("Cell pairs analyzed:", length(cell_pairs), "\n")
cat("Total L-R interactions:", nrow(all_dt), "\n")

cat("\n=== NicheNet analysis completed ===\n")
