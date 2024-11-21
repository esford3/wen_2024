

## run clustering and visualization 

# cd 
# ipython 

import pandas as pd
from tcrdist.repertoire import TCRrep
from tcrdist.public import _neighbors_fixed_radius
import os
import networkx as nx
import numpy as np
from tcrdist.html_colors import get_html_colors
from matplotlib import pyplot as plt
from community import community_louvain

## cluster wen_2024 paired sequences 
## file prep - see 'database prep 062124.R' 

## INPUT 
'''
path = INPUT FILEPATH HERE
'''
file = 'data/wen_for_tcrdist_consoltoaa_081924.csv' 

try : 
	df = pd.read_csv(os.path.join(path, file))
except FileNotFoundError: 
	print('remember to exchange out file & path')

## OUTPUT (ADJUST ACCORDINGLY)
'''
figure = "FIGURE FILENAME"
CSV = "CSV FILENAME"
'''

## DATA DICTIONARY 
df.columns
# 'tet' - tetramer used to identify TCR
# 'ptid' - subject ID
# 'v_a_gene' - v_alpha
# 'v_b_gene' - v_beta
# 'j_a_gene' - j_alpha
# 'j_b_gene' - j_beta
# 'cdr3_a_aa' - alpha cdr3 amino acid sequence
# 'cdr3_b_aa' - beta cdr3 amino acid sequence
# 'count' - number of times the unique paired sequence was observed

## CHECK 
df.shape #745, 9 

TCRrep.cpus = 4
tr = TCRrep(cell_df = df, 
            organism = 'human', 
            chains = ['alpha','beta'],
            compute_distances = True)

edge_threshold = 24
# Compute Network
# <tr.pw_alpha_beta> is paired chain TCRdist.
tr.pw_alpha_beta = tr.pw_beta + tr.pw_alpha
# <network> initialize a list to populate with edges between TCRs.
network = list()

from tcrdist.public import _neighbors_fixed_radius
row = list()

for i,n in enumerate(_neighbors_fixed_radius(tr.pw_alpha_beta, edge_threshold)):
    for j in n:
        if i != j:
            row = list()
            row.append(i)
            row.append(j)
            row.append(len(n)-1)
            row.append((tr.pw_alpha_beta )[i,j])
            for col in tr.clone_df.columns:
                row.append(tr.clone_df[col].iloc[i])
                row.append(tr.clone_df[col].iloc[j])
        network.append(row)                      
        # 'K_neighbors' - number of neighbors

import numpy as np
cols = ['node_1', 'node_2','K_neighbors','dist'] + np.concatenate([(f"{x}_1", f"{x}_2") for x in tr.clone_df.columns.to_list()]).tolist()
df_net = pd.DataFrame(network, columns = cols)
df_net['weight'] = (edge_threshold - df_net['dist'])/edge_threshold
df_net = df_net.dropna()

from community import community_louvain
# Draw Graph
# <G> Initialize a networkx Graph instance from the columns of df_net.
G = nx.from_pandas_edgelist(pd.DataFrame({'source' : df_net['node_1'],
                                          'target' : df_net['node_2'],
                                          'weight' : df_net['weight'],
                                          'count1' : df_net['count_1'],
                                          'count2' : df_net['count_2'],
                                          'tet1' : df_net['tet_1'],
                                          'tet2' : df_net['tet_2'],
                                          'ptid1' : df_net['ptid_1'],
                                          'ptid2' : df_net['ptid_2']}))

partition= community_louvain.best_partition(G, random_state = 1)

# Change partition such that cluster id is in descending order based on community size

partitions_by_cluster_size = list(pd.Series(partition.values()).value_counts().index)
partition_reorder = {id:rank for id,rank in zip(partitions_by_cluster_size,
                                                range(len(partitions_by_cluster_size)))}

partition = {k:partition_reorder.get(v) for k,v in partition.items()}

tr.clone_df['cluster_id'] = [str(partition.get(i)) if partition.get(i) is not None else None for i in tr.clone_df.index]

df_net['cluster_id_1'] = df_net['node_1'].apply(lambda i: partition.get(i))
df_net['cluster_id_2'] = df_net['node_2'].apply(lambda i: partition.get(i))

color_var = 'tet'
print(G)
print(nx.number_connected_components(G), "connected components")
# layout graphs with positions using graphviz neato
pos = nx.nx_agraph.graphviz_layout(G, prog="neato")
# color nodes the same in each connected subgraph
C = (G.subgraph(c) for c in nx.connected_components(G))

# Nodes
ns = df_net['node_1'].to_list() + df_net['node_2'].to_list()
cts = df_net['count_1'].to_list() + df_net['count_2'].to_list()

color_vars = pd.Series(df_net[f'{color_var}_1'].to_list() + df_net[f'{color_var}_2'].to_list())
node_cts = {node: x for node, x in zip(ns,cts)}
node_to_color_var = {node:c for node, c in zip(ns,color_vars)}

from tcrdist.html_colors import get_html_colors
colors   = get_html_colors(color_vars.nunique())
var_to_color = {v:color for v,color, in zip(color_vars.value_counts().index,colors)}
#by tet
var_to_color = {'none': 'gray', 'DR0101': 'red', 'DR0401': 'green', 'DR0404': 'blue',
				'DR0701': 'orange'}

#by person
#var_to_color = {8.0: 'red', 6.0: 'green', 1.0: 'blue', 7.0: 'cyan', 9.0: 'magenta',
#                 10.0: 'black', 2.0: 'lime'}

from matplotlib.pyplot import figure
plt.figure(1, figsize=(5, 5), dpi = 1200)
#figure(figsize=(5, 5), dpi=1200)

for g in C:
    my_labels = {x:partition.get(x) for x in g.nodes}
    my_ks = list(my_labels.keys())
    my_vs = list(my_labels.values())
    my_vs2 = ["" for x in my_vs]
    my_vs2[0] = my_vs[0]
    my_labels = {k:v for k,v in zip(my_ks, my_vs2)}
    nx.draw(g, pos,
            edge_color = "gray",
            alpha= .5,
            node_color = [var_to_color.get(node_to_color_var.get(x)) for x in g.nodes],
            node_size = [max(5, 5*node_cts.get(x)) for x in g.nodes],
            vmin=0.0, vmax=1.0, with_labels=True, labels = my_labels,
            font_color = "black", font_size = 9)
plt.savefig(figure)
var_to_color
plt.clf()

tr.clone_df.sort_values('cluster_id').to_csv(CSV, index = False)
