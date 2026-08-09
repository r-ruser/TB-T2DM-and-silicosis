"""
SCENIC analysis for GSE174725 (silicosis BALF) macrophage/monocyte cells.
Identifies regulons active in silicosis vs exposed-no-silicosis.
Windows-compatible: uses synchronous dask scheduler.
"""
import os
import sys
import pandas as pd
import numpy as np
import gzip
import glob
import warnings
warnings.filterwarnings('ignore')

# Fix Windows multiprocessing
if __name__ == '__main__':
    import dask
    dask.config.set(scheduler='synchronous')

    # Config
    PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    DATA_DIR = os.path.join(PROJECT_ROOT, 'data', '03_scRNA', 'GSE174725', 'raw', 'donor_matrices')
    RESULT_DIR = os.path.join(PROJECT_ROOT, 'result', '05_GRN_KO', 'SCENIC')
    os.makedirs(RESULT_DIR, exist_ok=True)
    os.makedirs(os.path.join(RESULT_DIR, 'tables'), exist_ok=True)
    os.makedirs(os.path.join(RESULT_DIR, 'source_data'), exist_ok=True)

    # Load cell annotation
    ann_file = os.path.join(PROJECT_ROOT, 'result', '04_scRNA', 'GSE174725', 'source_data', 'SD03_GSE174725_cell_annotation_UMAP.csv')
    ann = pd.read_csv(ann_file)
    print(f"Cell annotations: {len(ann)} cells")
    print(f"Cell types:\n{ann['cell_type'].value_counts()}")

    # Load all donor matrices
    print("\nLoading donor matrices...")
    all_data = {}
    for f in sorted(glob.glob(os.path.join(DATA_DIR, '*.csv.gz'))):
        fname = os.path.basename(f)
        with gzip.open(f, 'rt') as fh:
            df = pd.read_csv(fh, index_col=0)
        df = df[~df.index.duplicated(keep='first')]
        all_data[fname] = df
        print(f"  {fname}: {df.shape[1]} cells")

    # Combine all data
    print("\nCombining matrices...")
    all_counts = pd.concat(list(all_data.values()), axis=1, join='inner')
    print(f"Combined matrix: {all_counts.shape[0]} genes x {all_counts.shape[1]} cells")

    # Filter to myeloid cells
    myeloid_types = ['Macrophage', 'Monocyte']
    myeloid_cells = ann[ann['cell_type'].isin(myeloid_types)].copy()
    myeloid_cells['barcode'] = myeloid_cells['cell_id'].str.split('__').str[-1]
    print(f"\nMyeloid cells: {len(myeloid_cells)}")

    common_cells = list(set(all_counts.columns) & set(myeloid_cells['barcode']))
    counts_myeloid = all_counts[common_cells]
    print(f"Myeloid count matrix: {counts_myeloid.shape}")

    myeloid_meta = myeloid_cells.set_index('barcode').loc[common_cells]

    # Run GRNBoost2 with custom dask client
    print("\nRunning GRNBoost2...")
    from arboreto.algo import grnboost2
    from dask.distributed import Client, LocalCluster

    # Create a local cluster with threads only (no multiprocessing)
    cluster = LocalCluster(n_workers=1, threads_per_worker=4, processes=False)
    client = Client(cluster)
    print(f"Dask client: {client}")

    # Use top expressed genes as candidate TFs
    gene_means = counts_myeloid.mean(axis=1)
    top_genes = gene_means.nlargest(200).index.tolist()

    adjacencies = grnboost2(
        expression_data=counts_myeloid.T,
        tf_names=top_genes,
        verbose=True,
        client_or_address=client
    )
    client.close()
    cluster.close()
    print(f"GRN edges: {len(adjacencies)}")

    # Save adjacencies
    adj_file = os.path.join(RESULT_DIR, 'tables', 'T01_GSE174725_grnboost2_adjacencies.csv')
    adjacencies.to_csv(adj_file, index=False)
    print(f"Adjacencies saved: {adj_file}")

    # Create modules from adjacencies
    print("\nCreating modules...")
    from pyscenic.utils import modules_from_adjacencies
    modules_list = list(modules_from_adjacencies(adjacencies, counts_myeloid.T))
    print(f"Modules identified: {len(modules_list)}")

    # AUCell
    print("\nRunning AUCell...")
    from pyscenic.aucell import aucell
    auc_mtx = aucell(counts_myeloid.T, modules_list)
    print(f"AUC matrix: {auc_mtx.shape}")

    # Save AUC matrix
    auc_file = os.path.join(RESULT_DIR, 'source_data', 'SD01_GSE174725_regulon_auc.csv')
    auc_mtx.to_csv(auc_file)
    print(f"AUC matrix saved: {auc_file}")

    # Differential regulon activity
    if 'group' in myeloid_meta.columns:
        print("\nDifferential regulon activity...")
        from scipy import stats
        results = []
        for regulon in auc_mtx.columns:
            for cond in myeloid_meta['group'].unique():
                cells_cond = myeloid_meta[myeloid_meta['group'] == cond].index
                cells_other = myeloid_meta[myeloid_meta['group'] != cond].index
                if len(cells_cond) < 3 or len(cells_other) < 3:
                    continue
                auc_cond = auc_mtx.loc[cells_cond, regulon]
                auc_other = auc_mtx.loc[cells_other, regulon]
                stat, pval = stats.mannwhitneyu(auc_cond, auc_other, alternative='two-sided')
                results.append({
                    'regulon': regulon,
                    'condition': cond,
                    'n_cells_condition': len(cells_cond),
                    'n_cells_other': len(cells_other),
                    'mean_auc_condition': auc_cond.mean(),
                    'mean_auc_other': auc_other.mean(),
                    'delta_auc': auc_cond.mean() - auc_other.mean(),
                    'statistic': stat,
                    'pvalue': pval
                })

        results_df = pd.DataFrame(results)
        if len(results_df) > 0:
            results_df['fdr'] = stats.false_discovery_control(results_df['pvalue'])
            results_df = results_df.sort_values('fdr')

            results_file = os.path.join(RESULT_DIR, 'tables', 'T02_GSE174725_regulon_differential.csv')
            results_df.to_csv(results_file, index=False)
            print(f"Differential regulon results saved: {results_file}")
            print(f"\nTop differential regulons:")
            print(results_df.head(10).to_string())

    print("\nSCENIC analysis completed.")
