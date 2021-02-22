import pandas as pd

"""
"""

donor_table = pd.read_csv('object_donor_specimen.tsv', sep='\t')
donor_table = donor_table.loc[~donor_table['file_name'].str.contains('mini')]
donor_table.sort_values(by='donor_id', inplace=True)

donors = set(donor_table.donor_id.values)
sample_types = ['normal', 'tumour']

### Get Tumour Normal Pairs
# Group entries in the donor table by donor_id and if there are multiple
# tumour samples, keep only the first one.
X = []
for donor in donors:
    for sample_type in sample_types:
        X.append(donor_table.loc[
            (donor_table['donor_id'] == donor) &
            donor_table['specimen_type']
                        .str.lower()
                        .str.contains(sample_type)
        ].iloc[0].to_dict())
filtered_table = pd.DataFrame(X)
print(len(filtered_table))
filtered_table.to_csv('tumour-normal-table.tsv', sep='\t', index=False)

### Create Manifests from Tumour Normal pairs
manifest = pd.read_csv('prostate_cancer_manifest.tsv', sep='\t')

for donor in donors:
    file_ids = filtered_table[filtered_table['donor_id'] == donor]['file_id'].values
    donor_pair = manifest[manifest['file_id'].isin(file_ids)]
    donor_pair.to_csv(
        f"tumour-normal-manifests/{donor}-tumour-normal.tsv",
        sep='\t', index=False)
