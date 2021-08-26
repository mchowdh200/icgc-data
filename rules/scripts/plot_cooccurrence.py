import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import pyvis as pv

output_file = sys.argv[1]
input_files = sys.argv[2:]

### get list of all regions from one of the bed files
regions = []
with open(input_files[0]) as f:
    for line in f:
        A = line.rstrip().split()
        chrom = A[0]
        cytoband = A[3]
        regions.append(''.join([chrom, cytoband]))
    
# get the file ids from the file name {fid}.bins.bed
fid = []

### Find occurrence counts of variants
### at each region for each sample
# rows = fid, columns = genomic region
# used to denote whether a variant occured at a region for a given sample

CYTOBAND_COUNT_COL = 5
occ = np.zeros((len(input_files), len(regions)))
for i, bed in enumerate(input_files):
    fid.append(bed.split('.')[0])
    # occ[i,:] = np.where(np.loadtxt(bed, delimiter='\t', usecols=CYTOBAND_COUNT_COL) > 0, 1, 0)
    occ[i,:] = np.loadtxt(bed, delimiter='\t', usecols=CYTOBAND_COUNT_COL)

### Get co occurrence counts accross all regions
co_occ = occ.T @ occ
np.fill_diagonal(co_occ, 0)
# plt.matshow(co_occ)
# plt.show()
# co_occ = pd.DataFrame(data=co_occ, index=regions, columns=regions)


### Calculate the pointwise Mutual Information
# try to do this with vector/matrix ops
# PMI(X,Y) = log2(P(x,y)/(P(x)P(y)))

col_totals = np.sum(co_occ, axis=0)
total = np.sum(col_totals, axis=0)
row_totals = np.sum(co_occ, axis=1)
expected = np.outer(row_totals, col_totals) / total

with np.errstate(divide='ignore'):
    # ppmi = np.absolute(np.log(co_occ/expected))
    ppmi = np.log(co_occ/expected)
    ppmi[np.isinf(ppmi)] = 0.0
    ppmi[np.isnan(ppmi)] = 0.0
    # ppmi[ppmi > (np.mean(ppmi) + 2*np.std(ppmi))] = 20

# plt.matshow(ppmi.values)
# plt.show()

# set 0 to lower triangular matrix
# ppmi.values[np.tril(np.ones(co_occ.shape)).astype(np.bool)] = 0
ppmi[np.tril(np.ones(co_occ.shape)).astype(bool)] = 0

# reshape and filter only count > 0
################################
# plt.matshow(ppmi)
# plt.colorbar()
# plt.savefig(output_file)
# plt.show()
# exit()
################################
ppmi = pd.DataFrame(data=ppmi, index=regions, columns=regions)
ppmi.to_csv('ppmi.tsv', sep='\t')
# f = sns.heatmap(ppmi)
# f.set_xticks(range(len(regions)))
# f.set_yticks(range(len(regions)))
# f.set_xticklabels(regions)
# f.set_yticklabels(regions)
# plt.savefig(output_file)
# plt.show()
# exit()

ppmi = ppmi.stack()
ppmi = ppmi.rename_axis(
    ('source', 'target')).reset_index(name='weight')
ppmi = ppmi[ppmi.weight > 1]
# ppmi = ppmi[ppmi.weight > (np.mean(ppmi.weight)+2*np.std(ppmi.weight))]
print(ppmi.head())

# print(np.mean(ppmi.weight))
# print(np.median(ppmi.weight))
# print(np.std(ppmi.weight))
# print(np.min(ppmi.weight))
# print(np.max(ppmi.weight))
# ppmi['weight'] /= np.max(ppmi.weight)
# co_occ = co_occ[co_occ.weight > 0.05]
# print(len(co_occ))

G = nx.from_pandas_edgelist(ppmi, edge_attr=True)
nx.write_graphml(G, 'ppmi.graphml')
plt.subplots(figsize=(28,28))
# pos = nx.spring_layout(G, iterations=50)
pos = nx.spring_layout(G, k=0.1, iterations=50,)

degrees = dict(nx.degree(G))
# print([G[u][v] for u, v in G.edges()])
nx.draw_circular(G, with_labels=True, font_size=8,
        node_size=[v for v in degrees.values()],
        width=[0.1*G[u][v]['weight'] for u, v in G.edges()])
plt.savefig(output_file)
plt.show()
## TODO
# instead of plotting the weighted co-occurrence network, let the weights
# be some p-value so we can see which pairs are most significant
