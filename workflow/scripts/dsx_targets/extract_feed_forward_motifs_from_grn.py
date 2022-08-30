import pandas as pd
import anndata as ad
import networkx as nx
from scipy.stats import fisher_exact
import os

grn_file = snakemake.input[0]

grn  = (
    pd.read_table(grn_file)
    .query("not ((Gene == 'dsx') & (TF == 'dsx'))")
)
print(grn)

G = nx.from_pandas_edgelist(
    grn,
    source = "TF",
    target = "Gene",
    create_using = nx.DiGraph
)
print(G)

dsx_one_hop = set([x for x in G.neighbors("dsx")])

print(len(dsx_one_hop))

dsx_two_hop = set()
dsx_one_hop_type_x = {}
neighbors = {}

for x in dsx_one_hop:
    neighbors[x] = set([y for y in G.neighbors(x)])
    neighbors_one_hop = neighbors[x].intersection(dsx_one_hop)
    if len(neighbors_one_hop) > 0:
        dsx_one_hop_type_x[x] = neighbors_one_hop
    dsx_two_hop = dsx_two_hop.union(neighbors[x])

# keep in two hop only those that are not in one hop
dsx_two_hop = dsx_two_hop - dsx_one_hop

assert(len(set("dsx").intersection(dsx_one_hop)) == 0)

dsx_ff_motif = (
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
    .assign(frac_targets_y = lambda df: df.num_ys/df.num_targets_x)
    [["x", "num_ys", "num_targets_x", "frac_targets_y", "ys"]]
    .sort_values(
        ["frac_targets_y", "num_targets_x"],
        ascending = [False, False]
    )
    .set_index("x")
)

dsx_ff_motif.to_csv(snakemake.output["ff_motif"], sep = "\t")

with open(snakemake.output["one_hop"], "w") as f:
    for x in sorted(dsx_one_hop):
        f.write(f"{x}\n")

with open(snakemake.output["two_hop"], "w") as f:
    for x in sorted(dsx_two_hop):
        f.write(f"{x}\n")
