import sys
import numpy as np
import pandas as pd
from scipy import sparse
import matplotlib.pyplot as plt
import networkx as nx


FEATURE_COLUMN = int(sys.argv[1]) # zero based index
COUNT_COLUMN = int(sys.argv[2])
output_file = sys.argv[3]
input_files = sys.argv[4:]

### Get list of all features from one of the bed files
# we are assuming that all bed files will contain the same set of features
# ie (regions, genes, etc.)
features = []
with open(input_files[0]) as f:
    for line in f:
        A = line.rstrip().split()
        features.append(A[FEATURE_COLUMN])


### Find occurrence counts at each feature for each sample
occ = np.zeros((len(input_files), len(features)), dtype=np.float32)
for i, bed in enumerate(input_files):
    # occ[i,:] = np.where(np.loadtxt(bed, delimiter='\t', usecols=COUNT_COLUMN) > 0, 1, 0)
    occ[i,:] = np.loadtxt(bed, delimiter='\t', usecols=COUNT_COLUMN, dtype=np.float32)
print(occ.shape)
occ = sparse.csr_matrix(occ)
assert(sparse.issparse(occ))
### Get cooccurrence counts accross all features
print('co_occ')
co_occ = sparse.triu(occ.T @ occ, k=1)
assert(sparse.issparse(co_occ))
# np.fill_diagonal(co_occ, 0)
# print('set_diag')
# co_occ.setdiag(0)
# co_occ.eliminate_zeros()


### Calculate the pointwise Mutual Information
# PMI(X,Y) = log2(P(x,y)/(P(x)P(y)))
# TODO check this
# col_totals = np.sum(co_occ, axis=0)
# total = np.sum(col_totals, axis=0)
# row_totals = np.sum(co_occ, axis=1)
# expected = np.outer(row_totals, col_totals) / total
print('col_totals')
col_totals = sparse.csr_matrix(co_occ.sum(axis=0))
print('total')
total = col_totals.sum(axis=0)
print('row_totals')
row_totals = sparse.csr_matrix(co_occ.sum(axis=1))
print('expected')
print(f'{row_totals.shape=}')
print(f'{col_totals.shape=}')
expected = (row_totals @ col_totals)/total
assert(sparse.issparse(expected))
exit(1)

with np.errstate(divide='ignore'):
    # ppmi = np.absolute(np.log(co_occ/expected))
    ppmi.data = np.log(co_occ/expected)
    ppmi.data[np.isnan(ppmi)] = 0.0

ppmi[np.tril(np.ones(co_occ.shape)).astype(bool)] = 0

### reshape and filter ppmi > threshold
# TODO make threshold a param
# TODO add output file for tsv of the original matrix
ppmi = pd.DataFrame.sparse.from_spmatrix(data=ppmi, index=features, columns=features)
# ppmi.to_csv('ppmi.tsv', sep='\t')
ppmi = ppmi.stack()
ppmi = ppmi.rename_axis(
    ('source', 'target')).reset_index(name='weight')
ppmi = ppmi[ppmi.weight > 1]
ppmi = ppmi.sparse.to_dense()

### write to a graphml for later visualization
G = nx.from_pandas_edgelist(ppmi, edge_attr=True)
nx.write_graphml(G, 'ppmi.graphml')
