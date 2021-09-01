import sys
import numpy as np
import pandas as pd
from scipy import sparse
from joblib import Parallel, delayed
import matplotlib.pyplot as plt
import networkx as nx


N_JOBS = int(sys.argv[1])
FEATURE_COLUMN = int(sys.argv[2]) # zero based index
COUNT_COLUMN = int(sys.argv[3])
OUTPUT_FILE = sys.argv[4]
INPUT_FILES = sys.argv[5:]

### Get list of all features from one of the bed files
# we are assuming that all bed files will contain the same set of features
# ie (regions, genes, etc.)
features = []
with open(INPUT_FILES[0]) as f:
    for line in f:
        A = line.rstrip().split()
        features.append(A[FEATURE_COLUMN])


### Find occurrence counts at each feature for each sample
occ = np.zeros((len(input_files), len(features)), dtype=np.float32)
for i, bed in enumerate(input_files):
    # occ[i,:] = np.where(np.loadtxt(bed, delimiter='\t', usecols=COUNT_COLUMN) > 0, 1, 0)
    occ[i,:] = np.loadtxt(bed, delimiter='\t', usecols=COUNT_COLUMN, dtype=np.float32)
occ = sparse.csr_matrix(occ)


### Get cooccurrence counts accross all features
print('calculating cooccurrence matrix')
co_occ = sparse.triu(occ.T @ occ, k=0) # outputs a sparse coo_matrix
print(co_occ.shape)
assert(sparse.issparse(co_occ))


### Calculate the pointwise Mutual Information
# PMI(xi, xj) = log2(P(xi, xj)/(P(xi)P(xj)))
#
# The diag will have the single var counts.
# From that we can get the total # of events and normalize
# the matrix to create calculate P(xi, xj) and P(xi), ..., P(xn).
# Then we will iterate over each nonzero element and compute PMI(xi, xj).

# single counts and probabilities
print('getting diag')
single_counts = co_occ.diagonal()
co_occ.setdiag(0)
co_occ.eliminate_zeros()
total = np.sum(single_counts)
P = single_counts*(1/total) # multiplying preservese float32 dtype

# convert co_occ counts into joint probabilities
print('normalizing matrix')
co_occ *= (1/total)

# divide by independent probability and take log to complete
print('iteration over matrix')
co_occ.data = np.array(
    Parallel(n_jobs=N_JOBS, verbose=10)(
        delayed(lambda i, j, v: np.log2(v * (1/(P[i]*P[j]))))
        for i, j, v in zip(co_occ.row, co_occ.col, co_occ.data)),
    dtype=np.float32
)
exit(1)


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
nx.write_graphml(G, OUTPUT_FILE)
