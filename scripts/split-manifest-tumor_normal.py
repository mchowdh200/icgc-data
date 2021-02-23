import pandas as pd

manifest = 'prostate_cancer_manifest.tsv'

# read the manifest

data:pd.DataFrame
data = pd.read_csv(manifest, sep='\t')

# drop the entries with mini in the filename
data = data[~data['file_name'].str.contains('mini')]
data.sort_values('donor_id/donor_count', inplace=True)

# TESTING
data.to_csv('test.txt',
            index=False,
            sep='\t')

