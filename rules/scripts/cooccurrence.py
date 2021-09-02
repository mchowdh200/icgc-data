import os
import sys
import numpy as np
from scipy import sparse
import networkx as nx
from numba import jit, prange

### Get list of all features from one of the bed files
# we are assuming that all bed files will contain the same set of features
# ie (regions, genes, etc.)
def get_features(input_files, feature_column):
    """
    returns list of features under consideration
    """
    features = []
    with open(input_files[0]) as f:
        for line in f:
            A = line.rstrip().split()
            features.append(A[feature_column])
    return features


### Find occurrence counts at each feature for each sample
def occurrence_counts(input_files, features, count_column):
    """
    Returns sparse csr_matrix of occurrence counts.
    rows = samples
    columns = features
    """
    occ = np.zeros((len(input_files), len(features)), dtype=np.float32)
    for i, bed in enumerate(input_files):
        occ[i,:] = np.where(
            np.loadtxt(bed, delimiter='\t', usecols=count_column) > 0, 1, 0)
    return sparse.csr_matrix(occ)


### Get cooccurrence counts accross all features
def cooccurrence_counts(occ):
    """
    Returns sparse coo_matrix of cooccurence counts, as well as the diag
    of the matrix.  The diagonal of the coo_matrix will be set to zero.
    """
    co_occ = (occ.T @ occ).tocoo()
    single_counts = co_occ.diagonal()
    co_occ.setdiag(0)
    co_occ.eliminate_zeros()
    return co_occ, single_counts


### Calculate the pointwise Mutual Information
@jit(nopython=True, parallel=True)
def compute_pmi(row, col, data, single_counts):
    """
    PMI(xi, xj) = log2(P(xi, xj)/(P(xi)P(xj)))
    row: array of row indices
    col: array of column indices
    data: underlying array for coo_matrix
    single_counts: counts vector of individual events

    returns pmi
     * pmi is a 1d array which will replace the data array of the coo matrix
    """

    # compute pmi
    P = single_counts/np.sum(single_counts)
    Ptotal = np.sum(data)
    pmi = np.array(
        [np.log2(data[i]/(Ptotal*P[row[i]]*P[col[i]]))
         for i in prange(len(data))],
        dtype=np.float32)
    # np.nan_to_num(pmi, copy=False, nan=0.0)
    # # get stats and filter
    # mean, std = np.mean(pmi), np.std(pmi)
    # pmi[pmi < (mean + 4*std)] = 0.0
    pmi[pmi < 0] = 0.0
    return pmi
        
### TODO make upper triangular


if __name__ == '__main__':
    feature_column = int(sys.argv[1]) # zero based index
    count_column = int(sys.argv[2])
    output_file = sys.argv[3]
    input_files = sys.argv[4:]

    ## get the ppmi from cooccurence data
    print("getting co_occ")
    features = get_features(input_files, feature_column)

    occ = occurrence_counts(input_files, features, count_column)
    co_occ, single_counts = cooccurrence_counts(occ)

    print("computing pmi")
    co_occ.data = compute_pmi(co_occ.row, co_occ.col, co_occ.data, single_counts)
    np.nan_to_num(co_occ.data, copy=False, nan=0.0)
    co_occ.eliminate_zeros()
    co_occ = sparse.triu(co_occ)

    ## look at stats
    mean = np.mean(co_occ.data)
    std = np.std(co_occ.data)
    max = np.max(co_occ.data)
    min = np.min(co_occ.data)
    print(f'{min = } {max = } {mean = } {std = }')
    print(len(co_occ.data))

    ## filter
    co_occ.data[co_occ.data < (mean + 6*std)] = 0.0
    co_occ.eliminate_zeros()
    mean = np.mean(co_occ.data)
    std = np.std(co_occ.data)
    max = np.max(co_occ.data)
    min = np.min(co_occ.data)
    print(f'{min = } {max = } {mean = } {std = }')
    print(len(co_occ.data))

    ## Write results to disk
    sparse.save_npz(output_file, co_occ)

