import pandas as pd
import anndata as ad
import networkx as nx
from scipy.stats import fisher_exact
import os

dsx_occupancys_file = "DSX_targets.tsv"
putative_tfs_file = "../data/FlyMine_FlyTF_putative.tsv"
trusted_tfs_file = "../data/FlyMine_FlyTF_trusted.tsv"
grn_file = "resources/GRN.tsv"

dsx_occupancys_file =  snakemake.input["dsx_occupancy"]
putative_tfs_file = snakemake.input["putative_tfs"]
trusted_tfs_file = snakemake.input["trusted_tfs"]
dsx_ff_motif_file = snakemake.input["ff_motif"]
dsx_one_hop_file = snakemake.input["one_hop"]
dsx_two_hop_file = snakemake.input["two_hop"]
expr_file = snakemake.input["expr"]
biased_genes_file = snakemake.input["bias_grps"]

out_file = snakemake.output[0]
out_file2 = os.path.splitext(out_file)[0] + "_counts.tsv"
out_file3 = os.path.splitext(out_file)[0] + "_pvals.tsv"

ff_motif = (
    pd.read_table(dsx_ff_motif_file)
)
print(ff_motif)

dsx_occupancys = (
    pd.read_table(dsx_occupancys_file)[["Gene_Name"]]
    .rename(columns={"Gene_Name": "symbol"})
    .assign(dsx_occupancy = 1)
)
print(dsx_occupancys)

dsx_one_hop = (
    pd.read_table(dsx_one_hop_file, header = None, names = ["symbol"])
    .assign(dsx_one_hop = 1)
)
print(dsx_one_hop)

dsx_two_hop = (
    pd.read_table(dsx_two_hop_file, header = None, names = ["symbol"])
    .assign(dsx_two_hop = 1)
)
print(dsx_two_hop)

dsx_ff_x = (
    pd.DataFrame(dict(symbol = ff_motif["x"]))
    .assign(dsx_ff_x = 1)
)
print(dsx_ff_x)

dsx_ff_y = (
    pd.DataFrame(dict(symbol = ff_motif["ys"].str.split(",")))
    .explode("symbol")
    .drop_duplicates("symbol")
    .assign(dsx_ff_y = 1)
)
print(dsx_ff_y)

dsx_high_coreg = (
    ff_motif
    .query("frac_targets_y > 0.5")
    [["x"]]
    .rename(columns = {"x": "symbol"})
    .assign(dsx_ff_x_high = 1)
)

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
print(trusted_tfs)


adata = ad.read_h5ad(expr_file)
expressed = (
    adata.var.reset_index()
    [["symbol"]]
    .assign(expressed = 1)
)
expressed = (
    pd.read_table(biased_genes_file)
    [["symbol", "bias_group"]]
    .merge(expressed, how="right")
    .fillna("Unbiased")
    .merge(dsx_occupancys, how="left")
    .merge(putative_tfs, how="left")
    .merge(trusted_tfs, how="left")
    .merge(dsx_one_hop, how="left")
    .merge(dsx_two_hop, how="left")
    .merge(dsx_ff_y, how="left")
    .merge(dsx_ff_x, how="left")
    .merge(dsx_high_coreg, how="left")
    .fillna(0)
)
bias = expressed.query("bias_group != 'Unbiased'")
bias_copy = bias.assign(bias_group = "AllBiased")
expressed_copy = expressed.assign(bias_group = "AllExpressed")
bias = pd.concat([bias, bias_copy, expressed, expressed_copy])

print(bias)

enrich_frac = bias.groupby("bias_group").mean()
enrich_total = bias.groupby("bias_group").sum().astype(int)

print(enrich_total)
print(enrich_total.index.get_level_values(0))
print(enrich_total.columns.get_level_values(0))

def do_fisher(row, col):
    print("row", row)
    print("col", col)
    mat = enrich_total.loc[
            [row, "AllExpressed"],
            [col, "expressed"]
    ].values
    print(mat)
    table = [[mat[0,0],          mat[0,1]-mat[0,0]],
             [mat[1,0]-mat[0,0], mat[1,1]+mat[0,0]-mat[0,1]-mat[1,0]]]
    print(table)
    return fisher_exact(table, alternative='greater')[1]


rows = enrich_total.index.get_level_values(0)
cols = enrich_total.columns.get_level_values(0)

pvals = pd.DataFrame(
        [[do_fisher(row, col) for col in cols]
         for row in rows],
         index = rows,
         columns = cols,
)
print(pvals)


enrich_frac.to_csv(out_file, sep="\t", index=True)
enrich_total.to_csv(out_file2, sep="\t", index=True)
pvals.to_csv(out_file3, sep="\t", index=True)

