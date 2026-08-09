import pandas as pd
import gzip

# Check GSE174725 data structure
with gzip.open('data/03_scRNA/GSE174725/raw/donor_matrices/GSM5324866_exposure_patient_1.csv.gz', 'rt') as f:
    df = pd.read_csv(f, nrows=3, index_col=0)
    print('Shape:', df.shape)
    print('Columns (first 5):', list(df.columns[:5]))
    print('Index (first 5):', list(df.index[:5]))
    print('Values range:', df.values.min(), '-', df.values.max())
