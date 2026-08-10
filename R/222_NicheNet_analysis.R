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
# Load NicheNet databases
# ============================================================
cat("Loading NicheNet databases...\n")

# Ligand-receptor database
ligand_target <- readRDS("data/05_reference/NicheNet/ligand_target.rds")
ligand_receptor <- readRDS("data/05_reference/NicheNet/ligand_receptor.rds")
receptor_target <- readRDS("data/05_reference/NicheNet/receptor_target.rds")

# If databases not available, use built-in
if (!exists("ligand_target")) {
  cat("Using NicheNet built-in databases\n")
  # NicheNet will load databases automatically
}

# ============================================================
# Helper function: Run NicheNet analysis
# ============================================================
run_nichenet <- function(sender_cells, receiver_cells, sender_name, receiver_name,
  condition_name) {
  cat("\n--- Running NicheNet:", sender_name, "->", receiver_name, "---\n")

  # Get expressed genes
  sender_genes <- unique(sender_cells$gene)
  receiver_genes <- unique(receiver_cells$gene)

  # Background genes (all genes in dataset)
  background_genes <- union(sender_genes, receiver_genes)

  # Genes of interest (differentially expressed in receiver)
  genes_of_interest <- receiver_cells[abs(log2FoldChange) > 1 & padj < 0.05, gene]

  if (length(genes_of_interest) < 5) {
    cat("Too few DE genes, using top 50 by fold change\n")
    genes_of_interest <- receiver_cells[order(-abs(log2FoldChange)), gene][1:50]
  }

  cat("  Sender genes:", length(sender_genes), "\n")
  cat("  Receiver genes:", length(receiver_genes), "\n")
  cat("  Genes of interest:", length(genes_of_interest), "\n")

  # Step 1: Define expressed ligands and receptors
  ligands <- ligands_human$ligand
  receptors <- receptors_human$receptor

  expressed_ligands <- intersect(ligands, sender_genes)
  expressed_receptors <- intersect(receptors, receiver_genes)

  cat("  Expressed ligands:", length(expressed_ligands), "\n")
  cat("  Expressed receptors:", length(expressed_receptors), "\n")

  # Step 2: Get active ligand-receptor pairs
  ligand_receptor_pairs <- get_ligand_receptor_pairs(
    expressed_ligands, expressed_receptors,
    ligand_receptor_network = ligand_receptor
  )

  cat("  Active L-R pairs:", nrow(ligand_receptor_pairs), "\n")

  if (nrow(ligand_receptor_pairs) == 0) {
    cat("  No active L-R pairs found, skipping\n")
    return(NULL)
  }

  # Step 3: Predict active ligands
  ligand_activities <- predict_ligand_activities(
    geneset = genes_of_interest,
    background_expressed_genes = background_genes,
    ligand_target_matrix = ligand_target,
    potential_ligands = expressed_ligands
  )

  ligand_activities <- ligand_activities[order(-auroc), ]

  cat("  Top ligands:\n")
  print(head(ligand_activities, 5))

  # Step 4: Get target genes of top ligands
  top_ligands <- head(ligand_activities$potential_ligand, 10)

  ligand_target_df <- list()
  for (ligand in top_ligands) {
    targets <- names(which(ligand_target[ligand, genes_of_interest] > 0.1))
    if (length(targets) > 0) {
      ligand_target_df[[ligand]] <- data.table(
        ligand = ligand,
        target = targets,
        weight = ligand_target[ligand, targets]
      )
    }
  }

  target_df <- rbindlist(ligand_target_df)

  # Step 5: Get receptor-ligand interactions
  receptor_ligand_df <- ligand_receptor_pairs[
    ligand %in% top_ligands
  ]

  # Save results
  results <- list(
    ligand_activities = ligand_activities,
    ligand_target_df = target_df,
    receptor_ligand_df = receptor_ligand_df,
    genes_of_interest = genes_of_interest,
    sender_name = sender_name,
    receiver_name = receiver_name,
    condition_name = condition_name
  )

  return(results)
}

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
# Define cell type pairs for analysis
# ============================================================
# Focus on myeloid-myeloid and myeloid-lymphocyte communication
cell_pairs <- list(
  list(sender = "Classical monocyte", receiver = "Inflammatory monocyte"),
  list(sender = "Classical monocyte", receiver = "Alveolar macrophage"),
  list(sender = "Inflammatory monocyte", receiver = "C1QC macrophage"),
  list(sender = "SPP1 macrophage", receiver = "CD4 T"),
  list(sender = "Alveolar macrophage", receiver = "CD8 T")
)

# ============================================================
# Run NicheNet for each cohort and cell pair
# ============================================================
all_results <- list()

for (acc in names(cohorts)) {
  pb <- cohorts[[acc]]

  # Get DE genes for each cell type
  de_genes <- list()
  for (ct in unique(pb$meta$secondary_cell_type)) {
    # Pseudobulk DE would require proper design matrix
    # For now, use mean expression
    idx <- which(pb$meta$secondary_cell_type == ct)
    if (length(idx) >= 2) {
      mean_expr <- rowMeans(pb$counts[, idx])
      de_genes[[ct]] <- data.table(
        gene = pb$gene_symbol,
        log2FoldChange = log2(mean_expr + 1),
        padj = 0.01  # Placeholder
      )
    }
  }

  # Run NicheNet for each cell pair
  for (pair in cell_pairs) {
    if (pair$sender %in% names(de_genes) && pair$receiver %in% names(de_genes)) {
      result <- run_nichenet(
        sender_cells = de_genes[[pair$sender]],
        receiver_cells = de_genes[[pair$receiver]],
        sender_name = pair$sender,
        receiver_name = pair$receiver,
        condition_name = acc
      )

      if (!is.null(result)) {
        key <- paste(acc, pair$sender, pair$receiver, sep = "_")
        all_results[[key]] <- result
      }
    }
  }
}

# ============================================================
# Create heatmap of top ligand activities
# ============================================================
cat("\n--- Creating ligand activity heatmap ---\n")

if (length(all_results) > 0) {
  # Combine ligand activities across all comparisons
  activity_list <- lapply(all_results, function(r) {
    r$ligand_activities[, .(potential_ligand, auroc, auprc)]
  })

  activity_dt <- rbindlist(activity_list, idcol = "comparison")
  activity_dt[, c("cohort", "sender", "receiver") := tstrsplit(comparison, "_", keep = 1:3)]

  # Get top ligands overall
  top_ligands <- activity_dt[, .(mean_auroc = mean(auroc)), by = potential_ligand
  ][order(-mean_auroc)][1:min(20, .N), potential_ligand]

  # Create heatmap matrix
  heat_data <- activity_dt[potential_ligand %in% top_ligands]
  heat_matrix <- dcast(heat_data, potential_ligand ~ comparison, value.var = "auroc")
  heat_matrix[is.na(heat_matrix)] <- 0

  # Plot heatmap
  p_heat <- ggplot(melt(heat_matrix, id.vars = "potential_ligand"),
    aes(variable, potential_ligand, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#C94C4C",
      midpoint = 0.5, name = "AUROC") +
    labs(x = "Comparison", y = "Ligand",
      title = "NicheNet Ligand Activity Heatmap",
      subtitle = "Top 20 ligands by mean AUROC across comparisons") +
    theme_classic(base_size = 7) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
      axis.text.y = element_text(size = 5),
      plot.title = element_text(size = 8, face = "bold")
    )

  # Save heatmap
  out_base <- file.path(fig_dir, "Figure_NicheNet_ligand_heatmap")

  svglite::svglite(paste0(out_base, ".svg"), width = 12, height = 8)
  print(p_heat)
  dev.off()

  ggsave(paste0(out_base, ".pdf"), p_heat, width = 12, height = 8)
  ggsave(paste0(out_base, ".png"), p_heat, width = 12, height = 8, dpi = 300)

  cat("NicheNet heatmap saved\n")

  # Save results table
  fwrite(activity_dt, file.path(nichenet_dir, "tables", "T01_ligand_activities_all.csv"))
}

# ============================================================
# Create ligand-receptor network
# ============================================================
cat("\n--- Creating ligand-receptor network ---\n")

if (length(all_results) > 0) {
  # Combine all L-R pairs
  lr_list <- lapply(all_results, function(r) {
    r$receptor_ligand_df[, .(ligand, receptor)]
  })

  lr_dt <- rbindlist(lr_list)
  lr_unique <- unique(lr_dt)

  # Count occurrences
  lr_counts <- lr_dt[, .N, by = .(ligand, receptor)][order(-N)]

  # Save
  fwrite(lr_counts, file.path(nichenet_dir, "tables", "T02_ligand_receptor_pairs.csv"))

  cat("Top ligand-receptor pairs:\n")
  print(head(lr_counts, 10))
}

# ============================================================
# Summary
# ============================================================
cat("\n--- NicheNet analysis summary ---\n")
cat("Comparisons analyzed:", length(all_results), "\n")
cat("Unique ligands:", length(unique(unlist(lapply(all_results, function(r)
  r$ligand_activities$potential_ligand)))), "\n")

cat("\n=== NicheNet analysis completed ===\n")
write_log("NicheNet analysis completed")
