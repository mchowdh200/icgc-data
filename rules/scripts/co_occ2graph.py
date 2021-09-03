import sys
import pandas as pd
import networkx as nx
from scipy import sparse

matrix_file = sys.argv[1]
features_file = sys.argv[2]
graph = sys.argv[3]

print('loading matrix')
co_occ = sparse.load_npz(matrix_file)

print('loading features')
features = [x.rstrip() for x in open(features_file).readlines()]

### convert to networkx and write as graphml
co_occ = pd.DataFrame.sparse.from_spmatrix(data=co_occ, index=features, columns=features)
G = nx.from_pandas_adjacency(co_occ, edge_attr=True)
nx.write_graphml(G, output_graph)
