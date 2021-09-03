import sys
import pandas as pd
import networkx as nx
from scipy import sparse

matrix_file = sys.argv[1]
features_file = sys.argv[2]
output_graph = sys.argv[3]

print('loading matrix')
co_occ = sparse.load_npz(matrix_file)

print('loading features')
features = [x.rstrip() for x in open(features_file).readlines()]
node_mapping = {i: f for i, f in enumerate(features)}

### convert to networkx and write as graphml
G = nx.convert_matrix.from_scipy_sparse_matrix(co_occ)
relabel_nodes(G, node_mapping, copy=False)
nx.write_graphml(G, output_graph)
