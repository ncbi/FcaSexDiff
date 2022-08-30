import pandas as pd
import anndata as ad
import networkx as nx
from scipy.stats import fisher_exact
import os

dsx_occupancys_file = "DSX_targets.tsv"
putative_tfs_file = "../data/FlyMine_FlyTF_putative.tsv"
trusted_tfs_file = "../data/FlyMine_FlyTF_trusted.tsv"
grn_file = "resources/GRN.tsv"

grn  = (
    pd.read_table(grn_file)
    .query("not ((Gene == 'dsx') & (TF == 'dsx'))")
)
print(grn)

G = nx.from_pandas_edgelist(grn,
        source = "TF",
        target = "Gene",
        create_using = nx.DiGraph)
print(G)

dsx_one_hop = set([x for x in nx.all_neighbors(G, "dsx")])

dsx_two_hop = set()
dsx_intermediate = set()
dsx_ff = set()
dsx_one_hop_type_x = {}
neighbors = {}

dsx_one_hop_type_y = {x:set([]) for x in dsx_one_hop}

for x in dsx_one_hop:
    neighbors[x] = set([y for y in nx.all_neighbors(G, x)])
    neighbors_one_hop = neighbors[x].intersection(dsx_one_hop)
    dsx_one_hop_type_x[x] = neighbors_one_hop
    for y in neighbors_one_hop:
        dsx_one_hop_type_y[y] = dsx_one_hop_type_y[y].union(set(x))
    if len(neighbors_one_hop) > 0:
        dsx_intermediate = dsx_intermediate.union(set(x))
    dsx_ff = dsx_ff.union(neighbors[x].intersection(dsx_one_hop))
    dsx_two_hop = dsx_two_hop.union(neighbors[x])

dsx_two_hop = dsx_two_hop - dsx_one_hop

assert(len(set("dsx").intersection(dsx_one_hop)) == 0)
assert(len(dsx_one_hop - dsx_ff) > 0)


dsx_ff = (
    pd.DataFrame({"symbol": list(dsx_ff)})
    .assign(dsx_ff = 1)
)
print(dsx_intermediate)


keep_tfs = ["dsx"] + list(dsx_one_hop)
keep_genes = list(dsx_one_hop)

grn = (
    grn.query("TF in @keep_tfs")
    .query("Gene in @keep_genes")
    .query("TF != Gene")
)
print(grn)

nodes = (
    pd.DataFrame({"node":keep_tfs})
    .assign(hop = lambda tdf: tdf.node.apply(lambda x: 0 if x == "dsx" else 1 if
        x in dsx_ff else 2))
    .assign(symbol = lambda tdf: tdf.node)
)
print(nodes)

dsx_one_hop = list(dsx_one_hop)

dsx_one_hop_type_x_df = (
    pd.DataFrame({
        'x': dsx_one_hop,
        'as_x_num_ys': [
            len(dsx_one_hop_type_x[x])
            for x in dsx_one_hop
        ],
        'num_targets_x': [
            len(neighbors[x])
            for x in dsx_one_hop
        ],
        'as_x_ys': [
            ','.join(sorted(dsx_one_hop_type_x[x]))
            for x in dsx_one_hop
        ],
    })
    .assign(symbol = lambda df: df.x)
    .merge(bias_info, how = "left")
    .fillna("Unbiased")
    .assign(frac_y_among_x_targets = lambda df: df.num_ys/df.num_targets_x)
    [["x", "bias_group", "num_targets_x",  "num_ys",
    "frac_y_among_x_targets", "ys"]]
    .sort_values(by = ["bias_group", "frac_y_among_x_targets", "num_ys"],
        ascending = [True, False, False])
)

print(dsx_one_hop_type_x_df)

dsx_one_hop_type_x_df.to_csv("dsx_ff_coregulators.csv", index=False)
nodes.to_csv("dsx_one_hop_targets_nodes.csv", index=False)
grn.to_csv("dsx_one_hop_targets_edges.csv", index=False)

