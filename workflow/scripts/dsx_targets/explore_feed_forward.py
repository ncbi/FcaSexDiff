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

print(len(dsx_one_hop))


dsx_two_hop = set()
dsx_intermediate = set()
dsx_ff = set()
dsx_one_hop_type_x = {}
neighbors = {}

for x in dsx_one_hop:
    neighbors[x] = set([y for y in nx.all_neighbors(G, x)])
    neighbors_one_hop = neighbors[x].intersection(dsx_one_hop)
    if len(neighbors_one_hop) > 0:
        dsx_one_hop_type_x[x] = neighbors_one_hop
        dsx_intermediate = dsx_intermediate.union(set(x))
    dsx_ff = dsx_ff.union(neighbors[x].intersection(dsx_one_hop))
    dsx_two_hop = dsx_two_hop.union(neighbors[x])

dsx_two_hop = dsx_two_hop - dsx_one_hop

assert(len(set("dsx").intersection(dsx_one_hop)) == 0)
assert(len(dsx_one_hop - dsx_ff) > 0)

keep_tfs = ["dsx"] + list(dsx_one_hop)
keep_genes = list(dsx_one_hop)

grn = (
    grn.query("TF in @keep_tfs")
    .query("Gene in @keep_genes")
    .query("TF != Gene")
)
print(grn)

biased_genes_file = "exports/biased_gene_groups_tsv/cellfilt~NoSexspecArtef/expr~LogNorm/fcaver~stringent/resol~L4.0/tissue~body/biased_gene_groups.tsv"

bias_info = (
    pd.read_table(biased_genes_file)
    [["symbol", "bias_group"]]
)

nodes = (
    pd.DataFrame({"node":keep_tfs})
    .assign(hop = lambda tdf: tdf.node.apply(lambda x: 0 if x == "dsx" else 1 if
        x in dsx_ff else 2))
    .assign(symbol = lambda tdf: tdf.node)
    .merge(bias_info, how="left")
    .fillna("Unbiased")
)
print(nodes)

dsx_one_hop_type_x_df = (
    pd.DataFrame({
        'x': dsx_one_hop_type_x.keys(),
        'num_ys': [
            len(dsx_one_hop_type_x[x])
            for x in dsx_one_hop_type_x.keys()
        ],
        'num_targets_x': [
            len(neighbors[x])
            for x in dsx_one_hop_type_x.keys()
        ],
        'ys': [
            ','.join(sorted(dsx_one_hop_type_x[x]))
            for x in dsx_one_hop_type_x.keys()
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

dsx_traingle_motif = pd.DataFrame()

for v in dsx_one_hop_type_x.keys():
    dsx_traingle_motif = pd.concat([
       dsx_traingle_motif,
       pd.DataFrame(dict(y = list(dsx_one_hop_type_x[v]))).assign(x = v)
    ])

print(dsx_traingle_motif)

dsx_traingle_motif.to_csv("dsx_ff_motifs.csv", index=False)
dsx_one_hop_type_x_df.to_csv("dsx_ff_coregulators.csv", index=False)
nodes.to_csv("dsx_one_hop_targets_nodes.csv", index=False)
grn.to_csv("dsx_one_hop_targets_edges.csv", index=False)

