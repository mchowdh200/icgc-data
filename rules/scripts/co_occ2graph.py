import pandas as pd
import networkx as nx
from scipy import sparse

print('loading matrix')
co_occ = sparse.load_npz(output_matrix)

### reshape for networkx format
co_occ = pd.DataFrame.sparse.from_spmatrix(data=co_occ, index=features, columns=features)
co_occ = co_occ.stack()
co_occ = co_occ.rename_axis(
    ('source', 'target')).reset_index(name='weight')
co_occ = co_occ[co_occ.weight > 0]
co_occ = co_occ.sparse.to_dense()

### write to a graphml for later visualization
G = nx.from_pandas_edgelist(co_occ, edge_attr=True)
nx.write_graphml(G, output_graph)
