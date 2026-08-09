source(file.path("R", "00_config.R"), encoding = "UTF-8")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(ragg))

cat("=== Gene-Set Direction Audit: GSE192483 ===\n\n")

# ============================================================
# 1. Load pre-defined signature (unchanged)
# ============================================================
meta <- fread(file.path(path_result, "06_final", "tables", "T12_TBDM_meta_analysis.csv"))
shared_up <- meta[direction == "up"][order(-beta)][1:200, gene_symbol]
shared_down <- meta[direction == "down"][order(beta)][1:200, gene_symbol]
cat("Signature: 200 up + 200 down (unchanged)\n")

# ============================================================
# 2. Load GSE192483 pseudobulk
# ============================================================
pb <- fread(file.path(path_result, "04_scRNA", "GSE192483", "source_data", "SD05_GSE192483_paired_pseudobulk_edgeR.csv"))
avail <- unique(pb$gene_symbol)
up_m <- shared_up[shared_up %in% avail]
dn_m <- shared_down[shared_down %in% avail]
cat("Matched: up=", length(up_m), "/", length(shared_up),
    ", down=", length(dn_m), "/", length(shared_down), "\n")

# ============================================================
# 3. Gene-set tests per cell type
# ============================================================
cat("\n--- Gene-set direction audit ---\n")

cell_types <- unique(pb$cell_type)
audit_list <- list()

for (ct in cell_types) {
  ctp <- pb[cell_type == ct]
  if (nrow(ctp) == 0) next

  # Shared-up gene set
  up_fc <- ctp[gene_symbol %in% up_m, log2FC]
  if (length(up_fc) >= 5) {
    up_t <- t.test(up_fc, mu = 0)
    up_effect <- mean(up_fc)
    up_p <- up_t$p.value
  } else {
    up_effect <- NA; up_p <- NA
  }

  # Shared-down gene set
  dn_fc <- ctp[gene_symbol %in% dn_m, log2FC]
  if (length(dn_fc) >= 5) {
    dn_t <- t.test(dn_fc, mu = 0)
    dn_effect <- mean(dn_fc)
    dn_p <- dn_t$p.value
  } else {
    dn_effect <- NA; dn_p <- NA
  }

  # FDR
  ps <- p.adjust(c(up_p, dn_p), method = "BH")
  up_fdr <- ps[1]
  dn_fdr <- ps[2]

  # Net effect
  net_effect <- ifelse(!is.na(up_effect) & !is.na(dn_effect), up_effect - dn_effect, NA)

  # Direction
  up_dir <- ifelse(is.na(up_effect), "NS",
                   ifelse(up_effect > 0, "up_in_lesion", "down_in_lesion"))
  dn_dir <- ifelse(is.na(dn_effect), "NS",
                   ifelse(dn_effect > 0, "up_in_lesion", "down_in_lesion"))

  # Significance thresholds
  up_sig <- !is.na(up_fdr) & up_fdr < 0.05
  dn_sig <- !is.na(dn_fdr) & dn_fdr < 0.05

  # Classification
  if (up_sig & dn_sig) {
    if (up_dir == "up_in_lesion" & dn_dir == "down_in_lesion") {
      classification <- "Concordant"
      interpretation <- "Both gene sets show expected disease-associated direction in myeloid cells"
    } else if (up_dir == "down_in_lesion" & dn_dir == "up_in_lesion") {
      classification <- "Reverse"
      interpretation <- "Opposite to expected: shared-up DOWN, shared-down UP in lesion"
    } else {
      classification <- "Mixed"
      interpretation <- "Both significant but same direction — cannot determine concordant activation"
    }
  } else if (up_sig & !dn_sig) {
    if (up_dir == "up_in_lesion") {
      classification <- "Up-dominant partial concordance"
      interpretation <- "Shared-up UP in lesion; shared-down not significant"
    } else {
      classification <- "Reverse (partial)"
      interpretation <- "Shared-up DOWN in lesion; shared-down not significant"
    }
  } else if (!up_sig & dn_sig) {
    if (dn_dir == "down_in_lesion") {
      classification <- "Down-dominant partial concordance"
      interpretation <- "Shared-down DOWN in lesion; shared-up not significant"
    } else {
      classification <- "Reverse (partial)"
      interpretation <- "Shared-down UP in lesion; shared-up not significant"
    }
  } else {
    classification <- "No evidence"
    interpretation <- "Neither gene set shows significant change"
  }

  audit_list[[ct]] <- data.table(
    cell_type = ct,
    up_effect = up_effect,
    up_P = up_p,
    up_FDR = up_fdr,
    up_direction = up_dir,
    down_effect = dn_effect,
    down_P = dn_p,
    down_FDR = dn_fdr,
    down_direction = dn_dir,
    net_effect = net_effect,
    classification = classification,
    interpretation = interpretation
  )
}

audit_dt <- rbindlist(audit_list)
fwrite(audit_dt, file.path(path_result, "06_final", "tables", "T23_GSE192483_gene_set_direction_audit.csv"), bom = TRUE)

cat("\n=== Direction Audit Results ===\n")
print(audit_dt[, .(cell_type, up_effect, up_FDR, up_direction, down_effect, down_FDR, down_direction, classification)])

# ============================================================
# 4. Bidirectional forest plot
# ============================================================
cat("\n--- Generating bidirectional forest plot ---\n")

# Prepare data for plotting
plot_data <- melt(audit_dt[, .(cell_type, up_effect, down_effect, up_FDR, down_FDR, classification)],
                  id.vars = c("cell_type", "up_FDR", "down_FDR", "classification"),
                  variable.name = "gene_set", value.name = "effect")
plot_data[, gene_set := ifelse(gene_set == "up_effect", "Shared-up", "Shared-down")]
plot_data[, fdr := ifelse(gene_set == "Shared-up", up_FDR, down_FDR)]
plot_data[, sig := ifelse(fdr < 0.05, "FDR < 0.05", "NS")]

# Order by net effect
cell_order <- audit_dt[order(-net_effect), cell_type]
plot_data[, cell_type := factor(cell_type, levels = cell_order)]

# Expected concordant direction annotations
conc_dir <- data.table(
  cell_type = cell_order,
  expected_up = "up_in_lesion",
  expected_down = "down_in_lesion"
)

# Forest plot
p_forest <- ggplot(plot_data, aes(x = effect, y = cell_type, color = gene_set)) +
  geom_vline(xintercept = 0, linetype = 1, linewidth = 0.3, color = "black") +
  geom_vline(xintercept = c(-0.05, 0.05), linetype = 3, linewidth = 0.2, color = "grey60") +
  geom_point(aes(shape = sig), size = 2.5) +
  scale_color_manual(values = c("Shared-up" = "#E41A1C", "Shared-down" = "#377EB8"),
                     name = "Gene set") +
  scale_shape_manual(values = c("FDR < 0.05" = 16, "NS" = 1), name = "Significance") +
  labs(x = "Mean log2FC (lesion vs less-involved)", y = NULL,
       title = "GSE192483: Paired gene-set effects by cell type",
       subtitle = "Red = shared-up signature; Blue = shared-down signature; Points right = up in lesion") +
  theme_minimal(base_size = 7, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.title = element_text(size = 8, face = "bold"),
        plot.subtitle = element_text(size = 6, color = "#666666"),
        legend.position = "right",
        legend.key.size = unit(2, "mm"),
        legend.text = element_text(size = 6))

# Add classification labels
class_labels <- audit_dt[, .(cell_type, classification)]
class_labels[, cell_type := factor(cell_type, levels = cell_order)]

p_forest_annotated <- p_forest +
  geom_text(data = class_labels, aes(x = max(plot_data$effect, na.rm = TRUE) * 0.9,
                                      label = classification),
            hjust = 1, size = 2, color = "#333333")

# Save
out <- file.path(path_result, "06_final", "figures", "Fig4_gene_set_direction_audit")
w <- 183/25.4; h <- 90/25.4

svglite::svglite(paste0(out, ".svg"), width = w, height = h); print(p_forest_annotated); dev.off()
grDevices::cairo_pdf(paste0(out, ".pdf"), width = w, height = h, family = "Arial"); print(p_forest_annotated); dev.off()
ragg::agg_tiff(paste0(out, ".tiff"), width = w, height = h, units = "in", res = 600, compression = "lzw"); print(p_forest_annotated); dev.off()
ragg::agg_png(paste0(out, ".png"), width = w, height = h, units = "in", res = 300); print(p_forest_annotated); dev.off()

cat("Forest plot saved\n")

# ============================================================
# 5. Results text
# ============================================================
cat("\n=== Results Summary ===\n")

# Classify all cell types
for (i in 1:nrow(audit_dt)) {
  row <- audit_dt[i]
  cat("\n", row$cell_type, ":\n")
  cat("  Classification:", row$classification, "\n")
  cat("  Shared-up:", row$up_direction, "(effect=", round(row$up_effect, 4), ", FDR=", format.pval(row$up_FDR, digits = 3), ")\n")
  cat("  Shared-down:", row$down_direction, "(effect=", round(row$down_effect, 4), ", FDR=", format.pval(row$down_FDR, digits = 3), ")\n")
  cat("  Net:", round(row$net_effect, 4), "\n")
  cat("  Interpretation:", row$interpretation, "\n")
}

# Key findings
cat("\n--- Key Findings ---\n")
cat("Macrophage:", audit_dt[cell_type == "Macrophage"]$classification, "\n")
cat("CD8 T:", audit_dt[cell_type == "CD8 T"]$classification, "\n")
cat("Monocyte:", audit_dt[cell_type == "Monocyte"]$classification, "\n")

write_log("Gene-set direction audit completed")
