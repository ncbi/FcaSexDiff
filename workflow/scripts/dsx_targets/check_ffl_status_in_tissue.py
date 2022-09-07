import pandas as pd
import anndata as ad

resol = snakemake.wildcards["resol"]
expr_file = snakemake.input["expr"]
bias_file = snakemake.input["bias"]
biased_genes_file = snakemake.input["biased_groups"]
ff_motif_file = snakemake.input["ff_motif"]
out_summary_file = snakemake.output["summary"]
out_detailed_file = snakemake.output["detailed"]
regulon_activity_file = snakemake.input["regulon_activity"]
regulon_targets_file = snakemake.input["regulon_targets"]


df = pd.read_table(ff_motif_file)

regulon_activity = (
    pd.read_table(regulon_activity_file)
)

regulons = (
    regulon_activity.regulon.str.split(
        r'(.*)_\(([+-])\)-motif', expand = True
    )
    .rename(columns = {1: "x", 2:"dir"})
    [["x", "dir"]]
)

regulon_activity = pd.concat(
    [regulon_activity, regulons],
    axis = 1
)
regulon_activity = (
    regulon_activity
    .query("dir == '+'")
    [["x", "cluster", "avg_binarized_activity", "binarized_avg_activity"]]
    .set_index(["cluster", "x"])
)
print(regulon_activity)

binarized_avg_activity = (
    regulon_activity[["binarized_avg_activity"]]
    .unstack("x")
    ["binarized_avg_activity"]
)
print(binarized_avg_activity)


avg_binarized_activity = (
    regulon_activity[["avg_binarized_activity"]]
    .unstack("x")
    ["avg_binarized_activity"]
)
print(avg_binarized_activity)



regulon_targets = (
    pd.read_table(regulon_targets_file)
    .rename(columns = {"x": "regulon"})
)
regulons = (
    regulon_targets.regulon.str.split(
        r'(.*)_\(([+-])\)-motif', expand = True
    )
    .rename(columns = {1: "x", 2:"dir"})
    [["x", "dir"]]
)
regulon_targets = pd.concat(
    [regulon_targets, regulons],
    axis = 1
)
regulon_targets = (
    regulon_targets
    .query("dir == '+'")
    [["x", "y"]]
    .groupby("x")
    .agg({"y": list})
)
regulon_targets = {
    x: row['y']
    for x, row in regulon_targets.iterrows()
}


expr = ad.read_h5ad(expr_file)
bias = ad.read_h5ad(bias_file)


bias_info = (
    pd.read_table(biased_genes_file)
    [["symbol", "bias_group"]]
    .rename(columns = {"symbol": "x"})
)

summary = (
    df
    .merge(bias_info, how = "left")
    .fillna("Unbiased")
    .assign(frac_y_among_x_targets = lambda df: df.num_ys/df.num_targets_x)
    [["x", "bias_group", "num_targets_x",  "num_ys", "frac_targets_y", "ys"]]
    .sort_values(
        by = ["bias_group", "frac_targets_y", "num_ys"],
        ascending = [True, False, False]
    )
    .set_index("x")
)
print(summary)

summary.to_csv(out_summary_file, sep = "\t")

df = (
    df
    .assign(y = lambda tdf: tdf["ys"].str.split(","))
    [["x", "y"]]
    .explode("y")
)

print(df)


useful_genes = list(
    set(df.x).union(df.y).union(set(['dsx']))
    .intersection(set(expr.var_names))
    .intersection(set(bias.obs_names))
)

df = df.query("x in @useful_genes").query("y in @useful_genes")


df["y_in_x_regulon"] = df.apply(lambda row: (row["x"] in
    regulon_targets.keys()) and (row["y"] in regulon_targets[row["x"]]), axis=1)

print(df)

expr = expr[:, expr.var_names.isin(useful_genes)]
expr = expr[:, expr.X.todense().sum(axis = 0) > 0]

expr = (
    expr.to_df()
    .assign(cluster = expr.obs[resol])
    .groupby("cluster").median()
)

print(expr)


bias = bias[bias.obs_names.isin(useful_genes), :]
bias = bias.to_df().T

print(bias)

res = pd.DataFrame([
    ['M' if bias.loc[cls,gene] < 0 else 'F' if bias.loc[cls,gene] > 0 else 'E' if
        expr.loc[cls,gene] > 0 else 'N' for gene in bias.columns.tolist()]
    for cls in bias.index.tolist()],
    index = bias.index,
    columns = bias.columns)
print(res)

print(df)

def myfunc_dsx(row):
    #ret = paste(res['dsx'], res[row["x"]], res[row["y"]])
    ret = res['dsx']
    return pd.Series(ret, index = res.index)

def myfunc_x(row):
    #ret = paste(res['dsx'], res[row["x"]], res[row["y"]])
    ret = res[row["x"]]
    return pd.Series(ret, index = res.index)

def myfunc_y(row):
    #ret = paste(res['dsx'], res[row["x"]], res[row["y"]])
    ret = res[row["y"]]
    return pd.Series(ret, index = res.index)

def myfunc_baa(row):
    if row["x"] in binarized_avg_activity.columns:
        ret = binarized_avg_activity[row["x"]]
        return pd.Series(ret, index = binarized_avg_activity.index)
    else:
        return pd.Series(-1, index = binarized_avg_activity.index)

def myfunc_aba(row):
    if row["x"] in avg_binarized_activity.columns:
        ret = avg_binarized_activity[row["x"]]
        return pd.Series(ret, index = avg_binarized_activity.index)
    else:
        return pd.Series(-1, index = avg_binarized_activity.index)

tmp_dsx = (
    df.apply(myfunc_dsx, axis = 1)
)
print(tmp_dsx)

tmp_x = (
    df.apply(myfunc_x, axis = 1)
)
print(tmp_x)

tmp_y = (
    df.apply(myfunc_y, axis = 1)
)
print(tmp_y)

tmp_baa = (
    df.apply(myfunc_baa, axis = 1)
)
print(tmp_baa)

tmp_aba = (
    df.apply(myfunc_aba, axis = 1)
)
print(tmp_aba)

tmp = pd.concat(
    [tmp_dsx, tmp_x, tmp_y, tmp_baa, tmp_aba],
    keys = ["1dsx", "2x", "3y", "4baa", "5aba"],
    names = ["ffelement", "cluster"],
    axis = 1
)

print(tmp)
tmp = tmp.sort_index(level="cluster", axis = 1)
print(tmp)
#tmp = tmp.sort_index(['dsx', 'x', 'y', 'baa', 'aba'], level='ffelement', axis =
#        1)
#print(tmp)
#if 1: quit()

tmp.columns = ["_".join(x)[1:] for x in tmp.columns]

#if 1: quit()
#
#not_interesting = (
#    (tmp == "FNN") | (tmp == "MNN") | (tmp == "ENN") | (tmp == "NNN") |
#    (tmp == "FEN") | (tmp == "MEN") | (tmp == "EEN") | (tmp == "NEN") |
#    (tmp == "FMN") | (tmp == "MMN") | (tmp == "EMN") | (tmp == "NMN") |
#    (tmp == "FFN") | (tmp == "MFN") | (tmp == "EFN") | (tmp == "NFN") |
#    (tmp == "FNE") | (tmp == "MNE") | (tmp == "ENE") | (tmp == "NNE") |
#    (tmp == "FNM") | (tmp == "MNM") | (tmp == "ENM") | (tmp == "NNM") |
#    (tmp == "FNF") | (tmp == "MNF") | (tmp == "ENF") | (tmp == "NNF")
#)
#
#df["interesting"] = ~(not_interesting.sum(axis = 1)  ==  tmp.shape[1])
#print(df)

df = (
    pd.concat([df, tmp], axis = 1)
    .sort_values(["x", "y_in_x_regulon"], ascending = [True, False])
)
print(df)

df.to_csv(out_detailed_file, index = False, sep = "\t")
