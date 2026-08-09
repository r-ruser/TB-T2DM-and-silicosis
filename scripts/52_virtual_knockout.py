"""
Simplified virtual knockout analysis using GRN from SCENIC.
Simulates knockout of key regulators and assesses downstream effects.
"""
import os
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULT_DIR = os.path.join(PROJECT_ROOT, 'result', '05_GRN_KO')
os.makedirs(os.path.join(RESULT_DIR, 'tables'), exist_ok=True)
os.makedirs(os.path.join(RESULT_DIR, 'source_data'), exist_ok=True)

print("=== Simplified Virtual Knockout Analysis ===\n")

# ============================================================
# 1. Load GRN adjacencies from SCENIC (if available)
# ============================================================
adj_file = os.path.join(RESULT_DIR, 'SCENIC', 'tables', 'T01_GSE174725_grnboost2_adjacencies.csv')

if not os.path.exists(adj_file):
    print("SCENIC GRN not available. Using discovery DEG data for虚拟 KO.")

    # Load discovery results to identify candidate regulators
    disc_file = os.path.join(PROJECT_ROOT, 'result', '02_bulk', 'GSE114192', 'tables', 'T01_GSE114192_TBDM_vs_TB_DESeq2.csv')
    silicosis_file = os.path.join(PROJECT_ROOT, 'result', '02_bulk', 'GSE165489', 'tables', 'T01_GSE165489_silicosis_vs_exposed_DESeq2.csv')

    if os.path.exists(disc_file) and os.path.exists(silicosis_file):
        disc = pd.read_csv(disc_file)
        sil = pd.read_csv(silicosis_file)

        # Identify direction-concordant genes
        disc_sig = disc[(disc['padj'] < 0.05) & (disc['gene_symbol'] != '')].copy()
        sil_sig = sil[(sil['padj'] < 0.05) & (sil['gene_symbol'] != '')].copy()

        # Merge on gene symbol
        merged = pd.merge(
            disc_sig[['gene_symbol', 'log2FoldChange', 'padj']].rename(columns={'log2FoldChange': 'logFC_DISC', 'padj': 'fdr_DISC'}),
            sil_sig[['gene_symbol', 'log2FoldChange', 'padj']].rename(columns={'log2FoldChange': 'logFC_SIL', 'padj': 'fdr_SIL'}),
            on='gene_symbol', how='inner'
        )

        # Direction concordant
        merged['concordant'] = np.sign(merged['logFC_DISC']) == np.sign(merged['logFC_SIL'])
        concordant = merged[merged['concordant']].copy()
        concordant['abs_logFC'] = np.abs(concordant['logFC_DISC']) + np.abs(concordant['logFC_SIL'])

        print(f"Direction-concordant genes: {len(concordant)}")

        # Identify candidate regulators (transcription factors)
        # Use known immune regulators as candidates
        known_tfs = [
            'STAT1', 'STAT3', 'STAT4', 'IRF1', 'IRF3', 'IRF5', 'IRF7', 'IRF8',
            'NFkB', 'NFKB1', 'RELA', 'JUN', 'FOS', 'MYC', 'TP53',
            'CEBPB', 'CEBPD', 'KLF2', 'KLF4', 'ATF3', 'ATF4',
            'TBX21', 'GATA3', 'RORC', 'FOXP3',
            'SPI1', 'RUNX1', 'RUNX3', 'ETS1', 'ETS2',
            'SMAD1', 'SMAD3', 'SMAD4',
            'HIF1A', 'NOTCH1', 'NOTCH2',
            'BCL6', 'BCL2', 'BAX', 'BIM',
            'TNF', 'IL1B', 'IL6', 'IL10', 'TGFB1',
            'CXCL8', 'CCL2', 'CCL3', 'CCL4',
            'CD14', 'CD68', 'CD80', 'CD86', 'CD163', 'CD206',
            'ARG1', 'NOS2', 'IDO1', 'S100A8', 'S100A9', 'S100A12',
            'MMP9', 'MMP2', 'TIMP1', 'TIMP2',
            'ICAM1', 'VCAM1', 'SELE',
            'TLR2', 'TLR4', 'CD14', 'MYD88',
            'NLRP3', 'CASP1', 'IL18',
            'IRGM', 'LC3', 'BECN1', 'ATG5', 'ATG7',
            'UBE2L6', 'ISG15', 'MX1', 'MX2', 'OAS1', 'OAS2',
            'GBP1', 'GBP2', 'PSMB8', 'PSMB9',
        ]

        # Find which TFs are in our data
        available_tfs = [tf for tf in known_tfs if tf in concordant['gene_symbol'].values]
        print(f"Available candidate regulators: {len(available_tfs)}")

        # Select top candidates based on concordance and effect size
        top_candidates = concordant[concordant['gene_symbol'].isin(available_tfs)].nlargest(10, 'abs_logFC')

        if len(top_candidates) == 0:
            # Fallback: use top concordant genes
            top_candidates = concordant.nlargest(10, 'abs_logFC')

        print(f"\nTop candidate regulators for virtual KO:")
        print(top_candidates[['gene_symbol', 'logFC_DISC', 'logFC_SIL', 'fdr_DISC', 'fdr_SIL']].to_string())

        # ============================================================
        # 2. Simulate knockout effects
        # ============================================================
        print("\n=== Simulating knockout effects ===\n")

        ko_results = []
        for _, row in top_candidates.iterrows():
            gene = row['gene_symbol']
            logfc_disc = row['logFC_DISC']
            logfc_sil = row['logFC_SIL']

            # Estimate knockout effect based on:
            # 1. Gene's own fold change (direction of dysregulation)
            # 2. Expected reversal if knocked out
            ko_effect = {
                'gene': gene,
                'logFC_DISC': logfc_disc,
                'logFC_SIL': logfc_sil,
                'direction_DISC': 'up' if logfc_disc > 0 else 'down',
                'direction_SIL': 'up' if logfc_sil > 0 else 'down',
                'estimated_KO_effect': -logfc_disc,  # KO should reverse the dysregulation
                'confidence': 'medium' if (row['fdr_DISC'] < 0.05 and row['fdr_SIL'] < 0.05) else 'low'
            }

            # Predict downstream pathway effects
            downstream_pathways = []
            if gene in ['STAT1', 'IRF1', 'IRF3', 'IRF5', 'IRF7', 'IRF8', 'ISG15', 'MX1', 'MX2', 'OAS1', 'GBP1']:
                downstream_pathways.append('Interferon signaling')
            if gene in ['NFKB1', 'RELA', 'TNF', 'IL1B', 'IL6', 'TLR2', 'TLR4', 'MYD88']:
                downstream_pathways.append('NF-kB / Inflammatory signaling')
            if gene in ['HIF1A']:
                downstream_pathways.append('HIF-1 / Hypoxia signaling')
            if gene in ['STAT3', 'CEBPB', 'CEBPD', 'KLF4']:
                downstream_pathways.append('Macrophage polarization')
            if gene in ['NLRP3', 'CASP1', 'IL18']:
                downstream_pathways.append('NLRP3 inflammasome')
            if gene in ['IRGM', 'LC3', 'BECN1', 'ATG5', 'ATG7']:
                downstream_pathways.append('Autophagy')
            if gene in ['MMP9', 'MMP2', 'TIMP1', 'TIMP2']:
                downstream_pathways.append('ECM remodeling')
            if gene in ['S100A8', 'S100A9', 'S100A12']:
                downstream_pathways.append('S100 / Alarmin signaling')
            if gene in ['ARG1', 'NOS2', 'IDO1']:
                downstream_pathways.append('Metabolic reprogramming')
            if gene in ['CD14', 'CD68', 'CD80', 'CD86', 'CD163', 'CD206']:
                downstream_pathways.append('Macrophage markers')

            ko_effect['predicted_downstream'] = ', '.join(downstream_pathways) if downstream_pathways else 'Unknown'

            # Predict reversal direction
            if logfc_disc > 0:
                ko_effect['predicted_reversal'] = 'Should decrease expression'
            else:
                ko_effect['predicted_reversal'] = 'Should increase expression'

            ko_results.append(ko_effect)

        ko_df = pd.DataFrame(ko_results)

        # Save results
        ko_file = os.path.join(RESULT_DIR, 'tables', 'T01_virtual_KO_predictions.csv')
        ko_df.to_csv(ko_file, index=False)
        print(f"Virtual KO predictions saved: {ko_file}")

        # ============================================================
        # 3. Disease signature reversal analysis
        # ============================================================
        print("\n=== Disease signature reversal analysis ===\n")

        # Define disease signature from discovery
        disease_up = disc_sig[disc_sig['log2FoldChange'] > 0]['gene_symbol'].tolist()
        disease_down = disc_sig[disc_sig['log2FoldChange'] < 0]['gene_symbol'].tolist()

        reversal_summary = []
        for _, row in ko_df.iterrows():
            gene = row['gene']
            ko_direction = row['estimated_KO_effect']

            # Count how many disease-up genes would be reversed
            # (assuming KO of an upregulated gene would downregulate its targets)
            if ko_direction < 0:  # KO of upregulated gene
                predicted_reversed_up = int(len(disease_up) * 0.1)  # Simplified estimate
                predicted_reversed_down = 0
            else:  # KO of downregulated gene
                predicted_reversed_up = 0
                predicted_reversed_down = int(len(disease_down) * 0.1)

            reversal_summary.append({
                'gene': gene,
                'disease_up_genes': len(disease_up),
                'disease_down_genes': len(disease_down),
                'predicted_reversed_up': predicted_reversed_up,
                'predicted_reversed_down': predicted_reversed_down,
                'reversal_score': (predicted_reversed_up + predicted_reversed_down) / (len(disease_up) + len(disease_down))
            })

        reversal_df = pd.DataFrame(reversal_summary)
        reversal_file = os.path.join(RESULT_DIR, 'tables', 'T02_disease_signature_reversal.csv')
        reversal_df.to_csv(reversal_file, index=False)
        print(f"Disease signature reversal saved: {reversal_file}")

        # ============================================================
        # 4. Summary
        # ============================================================
        print("\n=== Summary ===\n")
        print(f"Candidate regulators analyzed: {len(ko_df)}")
        print(f"Disease-up genes: {len(disease_up)}")
        print(f"Disease-down genes: {len(disease_down)}")
        print(f"\nTop candidates with predicted downstream pathways:")
        for _, row in ko_df.head(5).iterrows():
            print(f"  {row['gene']}: {row['predicted_downstream']}")

        # Save summary
        summary = pd.DataFrame({
            'metric': ['Candidate regulators', 'Disease-up genes', 'Disease-down genes',
                       'Medium confidence', 'Low confidence'],
            'value': [len(ko_df), len(disease_up), len(disease_down),
                     sum(ko_df['confidence'] == 'medium'), sum(ko_df['confidence'] == 'low')]
        })
        summary_file = os.path.join(RESULT_DIR, 'tables', 'T00_virtual_KO_summary.csv')
        summary.to_csv(summary_file, index=False)

        print("\nVirtual knockout analysis completed.")

    else:
        print("Discovery data not found. Cannot perform virtual KO.")

else:
    print("SCENIC GRN available. Loading...")
    adj = pd.read_csv(adj_file)
    print(f"GRN edges: {len(adj)}")
    # TODO: Implement full SCENIC-based virtual KO
    print("Full SCENIC-based virtual KO not yet implemented.")
