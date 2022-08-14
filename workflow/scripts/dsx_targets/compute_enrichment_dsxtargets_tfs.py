import pandas as pd
import anndata as ad
import networkx as nx
from scipy.stats import fisher_exact
import os

dsx_occupancys_file = "DSX_targets.tsv"
putative_tfs_file = "../data/FlyMine_FlyTF_putative.tsv"
trusted_tfs_file = "../data/FlyMine_FlyTF_trusted.tsv"
grn_file = "resources/GRN.tsv"

#expr_file = "scraps/expr_norm_h5ad/expr~LogNorm/fcaver~stringent/tissue~head/expr_norm.h5ad"
#biased_genes_file = "head_biased_gene_groups.tsv"
#out_file = "head_dsxtargets_tfs_enrich_biased_groups.tsv"

expr_file = snakemake.input["expr"]

biased_genes_file = snakemake.input["bias_grps"]
out_file = snakemake.output[0]
out_file2 = os.path.splitext(out_file)[0] + "_counts.tsv"
out_file3 = os.path.splitext(out_file)[0] + "_pvals.tsv"

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


dsx_occupancys = (
    pd.read_table(dsx_occupancys_file)[["Gene_Name"]]
    .rename(columns={"Gene_Name": "symbol"})
    .assign(dsx_occupancy = 1)
)
print(dsx_occupancys)

dsx_one_hop = (
    pd.DataFrame({"symbol": list(dsx_one_hop)})
    .assign(dsx_one_hop = 1)
)
print(dsx_one_hop)

dsx_two_hop = (
    pd.DataFrame({"symbol": list(dsx_two_hop)})
    .assign(dsx_two_hop = 1)
)
print(dsx_two_hop)

dsx_intermediate = (
    pd.DataFrame({"symbol": list(dsx_intermediate)})
    .assign(dsx_intermediate = 1)
)
print(dsx_intermediate)

dsx_ff = (
    pd.DataFrame({"symbol": list(dsx_ff)})
    .assign(dsx_ff = 1)
)
print(dsx_intermediate)



putative_tfs = (
    pd.read_table(putative_tfs_file, header = None)[[1]]
    .rename(columns={1: "symbol"})
    .assign(putative_tf = 1)
)
print(putative_tfs)

trusted_tfs = (
    pd.read_table(trusted_tfs_file, header = None)[[1]]
    .rename(columns={1: "symbol"})
    .assign(trusted_tf = 1)
)
print(putative_tfs)


#adata = ad.read_h5ad(expr_file)
#expressed = (
#    adata.var.reset_index()
#    [["symbol"]]
#    .assign(expressed = 1)
#)
#expressed = (
#    pd.read_table(biased_genes_file)
#    [["symbol", "bias_group"]]
#    .merge(expressed, how="right")
#    .fillna("Unbiased")
#    .merge(dsx_occupancys, how="left")
#    .merge(putative_tfs, how="left")
#    .merge(trusted_tfs, how="left")
#    .merge(dsx_one_hop, how="left")
#    .merge(dsx_two_hop, how="left")
#    .merge(dsx_intermediate, how="left")
#    .merge(dsx_ff, how="left")
#    .fillna(0)
#)
#bias = expressed.query("bias_group != 'Unbiased'")
#bias_copy = bias.assign(bias_group = "AllBiased")
#expressed_copy = expressed.assign(bias_group = "AllExpressed")
#bias = pd.concat([bias, bias_copy, expressed, expressed_copy])
#
#print(bias)
#
#enrich_frac = bias.groupby("bias_group").mean()
#enrich_total = bias.groupby("bias_group").sum().astype(int)
#
#print(enrich_total)
#print(enrich_total.index.get_level_values(0))
#print(enrich_total.columns.get_level_values(0))
#
#def do_fisher(row, col):
#    print("row", row)
#    print("col", col)
#    mat = enrich_total.loc[
#            [row, "AllExpressed"],
#            [col, "expressed"]
#    ].values
#    print(mat)
#    table = [[mat[0,0],          mat[0,1]-mat[0,0]],
#             [mat[1,0]-mat[0,0], mat[1,1]+mat[0,0]-mat[0,1]-mat[1,0]]]
#    print(table)
#    return fisher_exact(table, alternative='greater')[1]
#
#rows = ["AllBiased", "AllExpressed", "FemaleOnlyBiased",
#            "MaleOnlyBiased", "MixedBiased", "Unbiased"]
#cols = ["dsx_occupancy", "dsx_one_hop", "dsx_two_hop", "dsx_intermediate", "expressed", "putative_tf", "trusted_tf"]
#
#rows = enrich_total.index.get_level_values(0)
#cols = enrich_total.columns.get_level_values(0)
#
#pvals = pd.DataFrame(
#        [[do_fisher(row, col) for col in cols]
#         for row in rows],
#         index = rows,
#         columns = cols,
#)
#print(pvals)
#
#
#enrich_frac.to_csv(out_file, sep="\t", index=True)
#enrich_total.to_csv(out_file2, sep="\t", index=True)
#pvals.to_csv(out_file3, sep="\t", index=True)


keep_tfs = ["dsx"] + list(dsx_one_hop)
keep_genes = list(dsx_one_hop)

grn = (
    grn.query("TF in @keep_tfs")
    .query("Gene in @keep_genes")
    .query("TF != Gene")
)
print(grn)

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
        'dsx_target': dsx_one_hop,
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

