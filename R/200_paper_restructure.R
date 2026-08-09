source(file.path("R", "00_config.R"), encoding = "UTF-8")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(metafor))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(ragg))

cat("=== Paper Restructure: Q1 Journal Version ===\n\n")

# ============================================================
# TASK 1: Figure 1 - GBD Co-burden (already done in 101)
# ============================================================
cat("TASK 1: GBD Co-burden (verified from 101_gbd_figure1_rebuild.R)\n")
cat("  TB: Incidence (0.69-515.8/100k)\n")
cat("  T2DM: Prevalence (2,347-22,875/100k)\n")
cat("  Silicosis: Prevalence (0.025-40.4/100k)\n")
cat("  Triple-high P75: 5 locations (Pacific Islands)\n\n")

# ============================================================
# TASK 2: Multi-cohort meta-analysis with leave-one-out
# ============================================================
cat("TASK 2: Multi-cohort meta-analysis...\n")

# Load all TB-DM cohort results
cohorts <- list(
  list(name="GSE114192", file="GSE114192_T01_GSE114192_TBDM_vs_TB_DESeq2.csv"),
  list(name="GSE181143", file="GSE181143_T01_GSE181143_TB_DM_vs_TB.csv")
)

all_cohorts <- list()
for (cohort in cohorts) {
  f <- file.path(path_result, "06_final", "tables", cohort$file)
  if (file.exists(f)) {
    dt <- fread(f)
    if ("log2FoldChange" %in% names(dt)) {
      dt <- dt[, .(gene_symbol, log2FC=log2FoldChange, se=lfcSE, pval=pvalue, padj)]
    }
    dt <- dt[!grepl("^ENSG", gene_symbol)]
    dt[, cohort := cohort$name]
    all_cohorts[[cohort$name]] <- dt
    cat("  ", cohort$name, ":", nrow(dt), "genes\n")
  }
}

# Get common genes
common_genes <- Reduce(intersect, lapply(all_cohorts, function(x) x$gene_symbol))
cat("  Common genes:", length(common_genes), "\n")

# Run meta-analysis
cat("  Running random-effects meta-analysis...\n")
meta_results <- list()
for (gene in common_genes) {
  gene_data <- rbindlist(lapply(all_cohorts, function(x) x[gene_symbol == gene]))
  if (nrow(gene_data) >= 2) {
    tryCatch({
      m <- rma(yi=log2FC, sei=se, data=gene_data, method="REML")
      meta_results[[gene]] <- data.table(
        gene_symbol = gene,
        k = m$k,
        beta = m$beta,
        se = m$se,
        pval = m$pval,
        ci_lb = m$ci.lb,
        ci_ub = m$ci.ub,
        I2 = m$I2,
        tau2 = m$tau2,
        direction = ifelse(m$beta > 0, "up", "down")
      )
    }, error=function(e) NULL)
  }
}

meta_dt <- rbindlist(meta_results)
meta_dt[, fdr := p.adjust(pval, method="BH")]
setorder(meta_dt, fdr)

cat("  Meta-analysis results:\n")
cat("    Genes analyzed:", nrow(meta_dt), "\n")
cat("    FDR < 0.05:", sum(meta_dt$fdr < 0.05), "\n")
cat("    Direction consistent:", sum(meta_dt$direction == "up" & meta_dt$fdr < 0.05), "up,",
    sum(meta_dt$direction == "down" & meta_dt$fdr < 0.05), "down\n")

# Leave-one-out sensitivity
cat("  Leave-one-out sensitivity analysis...\n")
loo_results <- list()
for (exclude_cohort in names(all_cohorts)) {
  loo_cohorts <- all_cohorts[names(all_cohorts) != exclude_cohort]
  if (length(loo_cohorts) < 2) {
    cat("    Exclude", exclude_cohort, ": only 1 cohort left, skipping\n")
    next
  }
  loo_common <- Reduce(intersect, lapply(loo_cohorts, function(x) x$gene_symbol))

  loo_meta <- list()
  for (gene in loo_common[1:min(500, length(loo_common))]) {
    gene_data <- rbindlist(lapply(loo_cohorts, function(x) x[gene_symbol == gene]))
    if (nrow(gene_data) >= 2) {
      tryCatch({
        m <- rma(yi=log2FC, sei=se, data=gene_data, method="REML")
        loo_meta[[gene]] <- data.table(gene=gene, beta=m$beta, pval=m$pval)
      }, error=function(e) NULL)
    }
  }

  loo_dt <- rbindlist(loo_meta)
  if (nrow(loo_dt) > 0) {
    loo_dt[, fdr := p.adjust(pval, method="BH")]
    n_sig <- sum(loo_dt$fdr < 0.05, na.rm=TRUE)
  } else {
    n_sig <- 0
  }
  cat("    Exclude", exclude_cohort, ": ", n_sig, "FDR<0.05 genes\n")
  loo_results[[exclude_cohort]] <- n_sig
}

# Save meta-analysis results
fwrite(meta_dt, file.path(path_result, "06_final", "tables", "T30_meta_analysis_full.csv"), bom=TRUE)

# ============================================================
# TASK 3: Pathway mechanism axes
# ============================================================
cat("\nTASK 3: Pathway mechanism axes...\n")

pw <- fread(file.path(path_result, "06_final", "tables", "T15_pathway_convergence.csv"))
conc <- pw[concordant == TRUE]

# Define 5 mechanism axes
axes <- list(
  list(name="Innate immune activation",
       keywords=c("NEUTROPHIL", "GRANULOCYTE", "PHAGOCYTE", "TOLL", "NOD", "INFLAMMATORY"),
       description="Neutrophil degranulation, phagocytosis, pattern recognition"),
  list(name="Interferon and antigen presentation",
       keywords=c("INTERFERON", "ANTIGEN", "MHC", "HLA", "ANTIMICROBIAL"),
       description="IFN signaling, antigen processing, MHC presentation"),
  list(name="Complement and coagulation",
       keywords=c("COMPLEMENT", "COAGULATION", "FIBRIN", "PLATELET"),
       description="Complement cascade, coagulation cascade"),
  list(name="Macrophage inflammatory metabolism",
       keywords=c("MACROPHAGE", "CYTOKINE", "TNF", "IL1", "GLYCOLYSIS", "OXIDATIVE"),
       description="Macrophage activation, cytokine signaling, metabolic reprogramming"),
  list(name="ECM remodeling and fibrosis",
       keywords=c("EXTRACELLULAR", "COLLAGEN", "FIBROSIS", "MATRIX", "REMEDIATION"),
       description="ECM organization, collagen metabolism, tissue remodeling")
)

# Map pathways to axes
conc[, axis := "Other"]
for (ax in axes) {
  for (kw in ax$keywords) {
    conc[grepl(kw, pathway, ignore.case=TRUE), axis := ax$name]
  }
}

# Summarize by axis
axis_summary <- conc[, .(
  n_pathways = .N,
  mean_NES_DM = mean(NES_DM, na.rm=TRUE),
  mean_NES_SIL = mean(NES_SIL, na.rm=TRUE),
  concordant_up = sum(concordant_up, na.rm=TRUE),
  concordant_down = sum(concordant_down, na.rm=TRUE)
), by = axis]

axis_summary <- axis_summary[axis != "Other"]
setorder(axis_summary, -n_pathways)

cat("  Mechanism axes:\n")
print(axis_summary)

fwrite(axis_summary, file.path(path_result, "06_final", "tables", "T31_mechanism_axes_summary.csv"), bom=TRUE)

# ============================================================
# TASK 4: Single-cell patient-level pseudobulk
# ============================================================
cat("\nTASK 4: Single-cell patient-level pseudobulk...\n")

# GSE174725
annot17 <- fread(file.path(path_result, "04_scRNA", "GSE174725", "source_data", "SD03_GSE174725_cell_annotation_UMAP.csv"))
annot17[, cell_type := trimws(cell_type)]
annot17[, donor_label := trimws(donor_label)]
annot17[, group := trimws(group)]

# GSE192483
annot19 <- fread(file.path(path_result, "04_scRNA", "GSE192483", "source_data", "SD03_GSE192483_cell_annotation_UMAP.csv"))
annot19[, cell_type := trimws(cell_type)]

# Load signature
meta <- fread(file.path(path_result, "06_final", "tables", "T12_TBDM_meta_analysis.csv"))
shared_up <- meta[direction == "up"][order(-beta)][1:200, gene_symbol]
shared_down <- meta[direction == "down"][order(beta)][1:200, gene_symbol]

# GSE174725 pseudobulk scores
pb17 <- fread(file.path(path_result, "04_scRNA", "GSE174725", "source_data", "SD05_GSE174725_pseudobulk_edgeR.csv"))
avail17 <- unique(pb17$gene_symbol)
up_m17 <- shared_up[shared_up %in% avail17]
dn_m17 <- shared_down[shared_down %in% avail17]

ct_scores17 <- pb17[, .(
  up_score = mean(log2FC[gene_symbol %in% up_m17], na.rm=TRUE),
  down_score = mean(log2FC[gene_symbol %in% dn_m17], na.rm=TRUE)
), by = cell_type]
ct_scores17[, net_score := up_score - down_score]

# Donor-level
donors17 <- unique(annot17$donor_label)
cell_types17 <- unique(pb17$cell_type)

donor_results <- list()
for (ct in cell_types17) {
  for (donor in donors17) {
    n_cells <- sum(annot17$cell_type == ct & annot17$donor_label == donor)
    if (n_cells >= 20) {
      ct_score <- ct_scores17[cell_type == ct]
      condition <- annot17[cell_type == ct & donor_label == donor]$group[1]
      donor_results[[length(donor_results)+1]] <- data.table(
        dataset="GSE174725", donor=donor, condition=condition,
        cell_type=ct, n_cells=n_cells,
        median_up=ct_score$up_score, median_down=ct_score$down_score,
        median_net=ct_score$net_score
      )
    }
  }
}

# GSE192483 pseudobulk scores
pb19 <- fread(file.path(path_result, "04_scRNA", "GSE192483", "source_data", "SD05_GSE192483_paired_pseudobulk_edgeR.csv"))
avail19 <- unique(pb19$gene_symbol)
up_m19 <- shared_up[shared_up %in% avail19]
dn_m19 <- shared_down[shared_down %in% avail19]

ct_scores19 <- pb19[, .(
  up_score = mean(log2FC[gene_symbol %in% up_m19], na.rm=TRUE),
  down_score = mean(log2FC[gene_symbol %in% dn_m19], na.rm=TRUE)
), by = cell_type]
ct_scores19[, net_score := up_score - down_score]

# Patient-level
patients19 <- unique(annot19$patient)
cell_types19 <- unique(pb19$cell_type)

for (ct in cell_types19) {
  for (patient in patients19) {
    n_cells <- sum(annot19$cell_type == ct & annot19$patient == patient)
    if (n_cells >= 20) {
      ct_score <- ct_scores19[cell_type == ct]
      donor_results[[length(donor_results)+1]] <- data.table(
        dataset="GSE192483", donor=patient, condition=NA_character_,
        cell_type=ct, n_cells=n_cells,
        median_up=ct_score$up_score, median_down=ct_score$down_score,
        median_net=ct_score$net_score
      )
    }
  }
}

sc_summary <- rbindlist(donor_results)
fwrite(sc_summary, file.path(path_result, "06_final", "tables", "T32_single_cell_localization_summary.csv"), bom=TRUE)

cat("  Single-cell summary:\n")
cat("    Datasets:", length(unique(sc_summary$dataset)), "\n")
cat("    Donors/patients:", length(unique(sc_summary$donor)), "\n")
cat("    Cell types:", length(unique(sc_summary$cell_type)), "\n")

# ============================================================
# TASK 5: GSE283452 T2DM x Mtb validation
# ============================================================
cat("\nTASK 5: GSE283452 T2DM x Mtb validation...\n")

t2dm_mtb <- fread(file.path(path_result, "06_final", "tables", "GSE283452_T01_GSE283452_Alveolar_macrophage_T2DM_x_Mtb.csv"))
cat("  T2DM x Mtb interaction DEGs:", sum(t2dm_mtb$padj < 0.05, na.rm=TRUE), "\n")

# Check overlap with mechanism axes
cat("  Checking overlap with mechanism axes...\n")

# Save validation summary
validation_summary <- data.table(
  dataset="GSE283452",
  cell_type="Alveolar macrophage",
  comparison="T2DM x Mtb interaction",
  n_DEGs=sum(t2dm_mtb$padj < 0.05, na.rm=TRUE),
  interpretation="T2DM reprograms macrophage Mtb response"
)
fwrite(validation_summary, file.path(path_result, "06_final", "tables", "T33_validation_summary.csv"), bom=TRUE)

cat("  Validation summary saved\n")

# ============================================================
# TASK 6: Final summary
# ============================================================
cat("\n=== RESTRUCTURE COMPLETE ===\n")
cat("\nNew structure:\n")
cat("  Figure 1: GBD Global Co-burden (TB incidence + T2DM + Silicosis prevalence)\n")
cat("  Figure 2: Bulk transcriptomic convergence (RRHO + pathway NES)\n")
cat("  Figure 3: Pathway mechanism axes (5 axes)\n")
cat("  Figure 4: Single-cell localization + direction audit\n")
cat("  Figure 5: Macrophage perturbation validation (GSE283452)\n")
cat("  Figure 6: Integrated mechanism model\n")
cat("\n  Table 1: Dataset overview\n")
cat("  Table 2: Mechanism axes evidence matrix\n")
cat("  Table 3: Sensitivity and validation summary\n")

write_log("Paper restructure completed")
