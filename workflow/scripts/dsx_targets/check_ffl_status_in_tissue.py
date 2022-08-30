import pandas as pd
import anndata as ad

import itertools

def paste(*args, sep = '', collapse = None):
    """
    Port of paste from R
    Args:
        *args: lists to be combined
        sep: a string to separate the terms
        collapse: an optional string to separate the results
    Returns:
        A list of combined results or a string of combined results if collapse is not None
    """
    combs = list(zip(*args))
    out = [sep.join(str(j) for j in i) for i in combs]
    if collapse is not None:
        out = collapse.join(out)
    return out


print(paste(['a','b'], ['x', 'y'], ['p', 'q']))

expr_file = "imports/imported_from_loom_h5ad/fcaver~stringent/tissue~body/imported_from_loom.h5ad"
bias_file = "exports/sexdiff_h5ad/cellfilt~NoSexspecArtef/expr~LogNorm/fcaver~stringent/resol~L4.0/tissue~body/sexdiff.h5ad"

resol = snakemake.wildcards["resol"]
expr_file = snakemake.input["expr"]
bias_file = snakemake.input["bias"]
out_file = snakemake.output[0]


expr = ad.read_h5ad(expr_file)
bias = ad.read_h5ad(bias_file)

df = pd.read_csv("dsx_ff_motifs.csv")
print(df)
print(df.drop_duplicates())

useful_genes = list(
    set(df.x).union(df.y).union(set(['dsx']))
    .intersection(set(expr.var_names))
    .intersection(set(bias.obs_names))
)

df = df.query("x in @useful_genes").query("y in @useful_genes")



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

def myfunc(row):
    ret = paste(res['dsx'], res[row["x"]], res[row["y"]])
    return pd.Series(ret, index = res.index)


tmp = (
    df.apply(myfunc, axis = 1)
)
print(tmp)

df = pd.concat([df, tmp], axis = 1)
print(df)

df.to_csv(out_file, index = False, sep = "\t")
