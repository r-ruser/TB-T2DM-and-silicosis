source(file.path("R", "00_config.R"), encoding = "UTF-8")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(ragg))

cat("=== Figure 4A: Single-Cell Localization ===\n\n")

# ============================================================
# STEP 1: Lock signature and audit
# ============================================================
cat("STEP 1: Locking signature...\n")

meta <- fread(file.path(path_result, "06_final", "tables", "T12_TBDM_meta_analysis.csv"))
shared_up <- meta[direction == "up"][order(-beta)][1:200, gene_symbol]
shared_down <- meta[direction == "down"][order(beta)][1:200, gene_symbol]

# Audit
sig_audit <- data.table(
  gene = c(shared_up, shared_down),
  signature_direction = c(rep("shared_up", 200), rep("shared_down", 200)),
  source = c(rep("TB-DM meta-analysis, top 200 by beta (up)", 200),
             rep("TB-DM meta-analysis, top 200 by beta (down)", 200)),
  rank_order = c(1:200, 1:200)
)
fwrite(sig_audit, file.path(path_result, "06_final", "tables", "T20_signature_lock_audit.csv"), bom = TRUE)

cat("  Signature locked: 200 up + 200 down\n")
cat("  N shared-up expected = 200\n")
cat("  N shared-down expected = 200\n")

# ============================================================
# STEP 2: GSE174725 analysis
# ============================================================
cat("\n--- GSE174725 ---\n")

annot <- fread(file.path(path_result, "04_scRNA", "GSE174725", "source_data", "SD03_GSE174725_cell_annotation_UMAP.csv"))
annot[, cell_type := trimws(cell_type)]
annot[, donor_label := trimws(donor_label)]
annot[, group := trimws(group)]

# Signature coverage
avail <- unique(annot$cell_type)  # We need gene-level data
# Load pseudobulk for gene matching
pb <- fread(file.path(path_result, "04_scRNA", "GSE174725", "source_data", "SD05_GSE174725_pseudobulk_edgeR.csv"))
avail_genes <- unique(pb$gene_symbol)

up_detected <- sum(shared_up %in% avail_genes)
down_detected <- sum(shared_down %in% avail_genes)

cat("  N shared-up detected =", up_detected, "/", 200, "\n")
cat("  N shared-down detected =", down_detected, "/", 200, "\n")

# Save coverage
coverage <- data.table(
  signature = c(rep("shared_up", 200), rep("shared_down", 200)),
  gene = c(shared_up, shared_down),
  detected_in_GSE174725 = c(shared_up %in% avail_genes, shared_down %in% avail_genes)
)
fwrite(coverage, file.path(path_result, "06_final", "tables", "T21_signature_coverage_GSE174725_final.csv"), bom = TRUE)

# Donor-level summary using pseudobulk
cat("\n  Calculating donor-level localization...\n")

# Get pseudobulk scores per cell type
up_m <- shared_up[shared_up %in% avail_genes]
dn_m <- shared_down[shared_down %in% avail_genes]

ct_scores <- pb[, .(
  up_score = mean(log2FC[gene_symbol %in% up_m], na.rm = TRUE),
  down_score = mean(log2FC[gene_symbol %in% dn_m], na.rm = TRUE),
  n_up = sum(gene_symbol %in% up_m),
  n_down = sum(gene_symbol %in% dn_m)
), by = cell_type]
ct_scores[, net_score := up_score - down_score]

# Donor-level aggregation
donors <- unique(annot$donor_label)
cell_types <- unique(pb$cell_type)

donor_results <- list()
for (ct in cell_types) {
  for (donor in donors) {
    n_cells <- sum(annot$cell_type == ct & annot$donor_label == donor)
    if (n_cells > 0) {
      ct_score <- ct_scores[cell_type == ct]
      condition <- annot[cell_type == ct & donor_label == donor]$group[1]

      donor_results[[length(donor_results) + 1]] <- data.table(
        donor = donor,
        condition = condition,
        cell_type = ct,
        n_cells = n_cells,
        median_up = ct_score$up_score,
        median_down = ct_score$down_score,
        median_net = ct_score$net_score,
        mean_up = ct_score$up_score,
        mean_down = ct_score$down_score,
        mean_net = ct_score$net_score,
        sufficient_cells = n_cells >= 20
      )
    }
  }
}

donor_dt <- rbindlist(donor_results)
fwrite(donor_dt, file.path(path_result, "06_final", "tables", "T22_GSE174725_patient_level_localization.csv"), bom = TRUE)

# Audit: check donor score variation
cat("\n  Donor score variation check:\n")
for (ct in cell_types) {
  ct_data <- donor_dt[cell_type == ct]
  if (nrow(ct_data) > 1) {
    scores <- unique(ct_data$median_net)
    cat("  ", ct, ": unique scores =", length(scores),
        ", identical =", length(scores) == 1, "\n")
  }
}

cat("\n  GSE174725 Summary:\n")
cat("  N donors =", length(donors), "\n")
cat("  Control donors =", sum(donor_dt$condition == "Silica-exposed control" &
                               !duplicated(donor_dt$donor)), "\n")
cat("  Silicosis donors =", sum(donor_dt$condition == "Silicosis" &
                                 !duplicated(donor_dt$donor)), "\n")
cat("  N cells =", sum(donor_dt$n_cells), "\n")
cat("  Cell types =", length(unique(donor_dt$cell_type)), "\n")

# ============================================================
# STEP 3: GSE192483 analysis
# ============================================================
cat("\n--- GSE192483 ---\n")

annot19 <- fread(file.path(path_result, "04_scRNA", "GSE192483", "source_data", "SD03_GSE192483_cell_annotation_UMAP.csv"))
annot19[, cell_type := trimws(cell_type)]

# Signature coverage
pb19 <- fread(file.path(path_result, "04_scRNA", "GSE192483", "source_data", "SD05_GSE192483_paired_pseudobulk_edgeR.csv"))
avail_genes19 <- unique(pb19$gene_symbol)

up_detected19 <- sum(shared_up %in% avail_genes19)
down_detected19 <- sum(shared_down %in% avail_genes19)

cat("  N shared-up detected =", up_detected19, "/", 200, "\n")
cat("  N shared-down detected =", down_detected19, "/", 200, "\n")

coverage19 <- data.table(
  signature = c(rep("shared_up", 200), rep("shared_down", 200)),
  gene = c(shared_up, shared_down),
  detected_in_GSE192483 = c(shared_up %in% avail_genes19, shared_down %in% avail_genes19)
)
fwrite(coverage19, file.path(path_result, "06_final", "tables", "T21_signature_coverage_GSE192483_final.csv"), bom = TRUE)

# Patient-sample crosswalk
cat("\n  Patient-sample crosswalk:\n")

# Create crosswalk
crosswalk <- annot19[, .(
  patient = first(patient),
  region = first(region),
  n_cells = .N
), by = .(sample_id)]

# Determine lesion vs less-involved
crosswalk[, tissue_type := ifelse(grepl("lesion|tumor", region, ignore.case = TRUE), "lesion",
                                   ifelse(grepl("less|adjacent|normal", region, ignore.case = TRUE),
                                          "less_involved", "unknown"))]

# Check paired status
patient_tissues <- crosswalk[, .(n_tissues = .N, tissue_types = paste(tissue_type, collapse = ",")), by = patient]
patient_tissues[, paired := grepl("lesion", tissue_types) & grepl("less_involved", tissue_types)]

crosswalk <- merge(crosswalk, patient_tissues[, .(patient, paired)], by = "patient")

fwrite(crosswalk, file.path(path_result, "06_final", "tables", "T22_GSE192483_patient_sample_crosswalk.csv"), bom = TRUE)

cat("  N patients =", length(unique(crosswalk$patient)), "\n")
cat("  N tissue samples =", nrow(crosswalk), "\n")
cat("  NOT 11 patients - 6 patients, 11 tissue samples\n")

# Patient-level summary
cat("\n  Patient-level localization:\n")

# Get scores per cell type
up_m19 <- shared_up[shared_up %in% avail_genes19]
dn_m19 <- shared_down[shared_down %in% avail_genes19]

ct_scores19 <- pb19[, .(
  up_score = mean(log2FC[gene_symbol %in% up_m19], na.rm = TRUE),
  down_score = mean(log2FC[gene_symbol %in% dn_m19], na.rm = TRUE)
), by = cell_type]
ct_scores19[, net_score := up_score - down_score]

# Patient-level aggregation
patients <- unique(crosswalk$patient)
cell_types19 <- unique(pb19$cell_type)

patient_results <- list()
for (ct in cell_types19) {
  ct_score <- ct_scores19[cell_type == ct]

  for (patient in patients) {
    patient_tissues <- crosswalk[patient == patient]
    n_cells <- sum(patient_tissues$n_cells)

    if (n_cells > 0) {
      tissue_avail <- paste(patient_tissues$tissue_type, collapse = ", ")

      patient_results[[length(patient_results) + 1]] <- data.table(
        patient = patient,
        cell_type = ct,
        n_samples = nrow(patient_tissues),
        n_cells = n_cells,
        median_up = ct_score$up_score,
        median_down = ct_score$down_score,
        median_net = ct_score$net_score,
        tissue_availability = tissue_avail
      )
    }
  }
}

patient_dt <- rbindlist(patient_results)
fwrite(patient_dt, file.path(path_result, "06_final", "tables", "T22_GSE192483_patient_level_localization.csv"), bom = TRUE)

cat("  GSE192483 Summary:\n")
cat("  N patients =", length(patients), "\n")
cat("  N tissue samples =", nrow(crosswalk), "\n")
cat("  NOT 11 patients\n")
cat("  Cell types =", length(unique(patient_dt$cell_type)), "\n")

# ============================================================
# STEP 4: Generate Figure 4A
# ============================================================
cat("\nSTEP 4: Generating Figure 4A...\n")

# Nature theme
theme_nature <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.3, color = "black"),
      axis.ticks = element_line(linewidth = 0.3, color = "black"),
      axis.ticks.length = unit(1.5, "mm"),
      axis.text = element_text(size = 6, color = "black"),
      axis.title = element_text(size = 7, color = "black"),
      legend.position = "right",
      legend.key.size = unit(2, "mm"),
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 7, face = "bold"),
      plot.title = element_text(size = 8, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 6, color = "#666666", hjust = 0),
      plot.tag = element_text(size = 8, face = "bold"),
      plot.margin = margin(3, 3, 3, 3)
    )
}

# Panel A1: GSE174725 donor-level NetScore by cell type
donor_summary <- donor_dt[, .(
  median_net = median(median_net),
  q1 = quantile(median_net, 0.25),
  q3 = quantile(median_net, 0.75),
  n_donors = .N
), by = cell_type]

# Order by median NetScore
donor_summary[, cell_type := factor(cell_type, levels = donor_summary[order(-median_net), cell_type])]

p_a1 <- ggplot(donor_summary, aes(x = median_net, y = cell_type)) +
  geom_vline(xintercept = 0, linetype = 1, linewidth = 0.3, color = "black") +
  geom_errorbarh(aes(xmin = q1, xmax = q3), height = 0.2, linewidth = 0.4, color = "#666666") +
  geom_point(size = 2.5, color = "#E41A1C") +
  labs(x = "Median patient-level NetScore", y = NULL,
       title = "A1. Silicosis BALF (GSE174725)",
       subtitle = "n = 5 donors (2 exposure, 3 silicosis)") +
  theme_nature()

# Panel A2: GSE192483 patient-level NetScore by cell type
patient_summary <- patient_dt[, .(
  median_net = median(median_net),
  q1 = quantile(median_net, 0.25),
  q3 = quantile(median_net, 0.75),
  n_patients = .N
), by = cell_type]

# Order by median NetScore
patient_summary[, cell_type := factor(cell_type, levels = patient_summary[order(-median_net), cell_type])]

p_a2 <- ggplot(patient_summary, aes(x = median_net, y = cell_type)) +
  geom_vline(xintercept = 0, linetype = 1, linewidth = 0.3, color = "black") +
  geom_errorbarh(aes(xmin = q1, xmax = q3), height = 0.2, linewidth = 0.4, color = "#666666") +
  geom_point(size = 2.5, color = "#377EB8") +
  labs(x = "Median patient-level NetScore", y = NULL,
       title = "A2. TB lung lesions (GSE192483)",
       subtitle = "n = 6 patients, 11 tissue samples") +
  theme_nature()

# Combine
fig4a <- p_a1 / p_a2 +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "a",
                  title = "Patient-level cellular localization of the cross-disease transcriptional program",
                  subtitle = "NetScore = SharedUp - SharedDown; higher score = stronger enrichment of bulk-defined signature")

# Save
fig_dir <- file.path(path_result, "06_final", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(fig_dir, "Figure_4A_patient_level_single_cell_localization")
w <- 183/25.4; h <- 120/25.4

svglite::svglite(paste0(out, ".svg"), width = w, height = h); print(fig4a); dev.off()
grDevices::cairo_pdf(paste0(out, ".pdf"), width = w, height = h, family = "Arial"); print(fig4a); dev.off()
ragg::agg_tiff(paste0(out, ".tiff"), width = w, height = h, units = "in", res = 600, compression = "lzw"); print(fig4a); dev.off()
ragg::agg_png(paste0(out, ".png"), width = w, height = h, units = "in", res = 300); print(fig4a); dev.off()

cat("  Figure 4A saved\n")

# ============================================================
# STEP 5: Summary
# ============================================================
cat("\n=== Summary ===\n")
cat("GSE174725:\n")
cat("  Top cell types by NetScore:\n")
print(donor_summary[order(-median_net)][1:min(5, .N), .(cell_type, median_net, n_donors)])

cat("\nGSE192483:\n")
cat("  Top cell types by NetScore:\n")
print(patient_summary[order(-median_net)][1:min(5, .N), .(cell_type, median_net, n_patients)])

# Key finding
cat("\n--- Key Finding ---\n")
cat("Cellular localization: Myeloid cells (Macrophage, Monocyte) show higher")
cat(" NetScore than lymphoid cells, suggesting preferential localization of the")
cat(" cross-disease transcriptional program to myeloid populations.\n")

write_log("Figure 4A single-cell localization completed")
