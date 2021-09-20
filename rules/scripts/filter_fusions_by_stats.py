import sys
import pandas as pd

input = sys.argv[1]
output = sys.argv[2]
onekg_support_threshold = 1

fusion_table = pd.read_csv(input, sep='\t')
fusion_table.sort_values('PCA_index_num_reads',
                         ascending=False,
                         inplace=True)
fusion_table = fusion_table[
    (fusion_table.PCA_index_num_reads > 0) &
    (fusion_table.onekg_index_num_samples <= onekg_support_threshold)
]

fusion_table.to_csv(output, sep='\t', index=False, header=True)
