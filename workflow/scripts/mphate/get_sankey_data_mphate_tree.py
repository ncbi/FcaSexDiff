import pandas as pd
import numpy as np
import anndata as ad
import scprep
import pickle
import os

nlevel_use = 7
ad_file = "imports/expr_h5ad/expr_fat_body_stringent.h5ad"
pkl_file = "scraps/mphate_pkl/mphate_fat_body_stringent.pkl"


adata = ad.read_h5ad(ad_file)
print(adata)

mp_op = pickle.load(open(pkl_file, "rb"))
print(mp_op)

print(mp_op.levels)

# mphate saves resolution levels fine to coarse, we want reverse
levels = mp_op.levels[::-1]
print(levels)
for i, l in enumerate(levels):
    print(f'Level {i}, nclusters={len(set(mp_op.NxTs[l]))}')


if nlevel_use <= 0:
    # the highest level has singleton cell groups, not interesting
    nlevel_use = len(levels)-1

print(nlevel_use)

# make a dataframe of it
mp_clusters = pd.concat([
    pd.DataFrame(
        mp_op.NxTs[levels[i]],
        columns = [f"MP{i}"],
        index = adata.obs_names,
    )
    for i in range(nlevel_use)
], axis = 1)

print(mp_clusters)

def rename_clusters(df, cols = None):
    if cols is None:
        cols = df.columns.tolist()

    for col in cols[::-1]:
        df = df.sort_values(col)
    print(df)

    df["_tmp"] = ""
    cols = ["_tmp"] + cols
    
    def rename(x):
        old2new = {v:(i+1) for i, v in enumerate(sorted(x.unique()))}
        return [old2new[old] for old in x]

    for i, (prev, curr) in enumerate(zip(cols, cols[1:])):
        df[curr] = df.groupby(prev)[[curr]].transform(rename)
        sep = '.' if i > 0 else 'C'
        df[curr] = df[[prev, curr]].apply(lambda x: sep.join([str(xx) for xx in x]), axis=1)
        print(df)

    return df.drop("_tmp", axis = 1)
    
mp_clusters = rename_clusters(mp_clusters)
print(mp_clusters)

filter_cells = True

if filter_cells:
    print(adata.shape)
    adata = adata[adata.obs.sex.isin(["male", "female"])]
    print(adata.shape)
    print(pd.read_csv("resources/sex_specific_annotations.csv")
        .query("sex_specific != 'no'")
        .annotation.tolist()
        + ["artefact"]
    )
    adata = adata[~adata.obs.annotation.isin(
        pd.read_csv("resources/sex_specific_annotations.csv")
        .query("sex_specific != 'no'")
        .annotation.tolist()
        + ["artefact"]
    )]
    mp_clusters = mp_clusters.loc[adata.obs_names,:]

print(adata.shape)

print(mp_clusters)

links = pd.DataFrame()

for i in range(nlevel_use-1):
    df = pd.DataFrame(dict(
        src = mp_clusters.iloc[:,i], #.apply(lambda x: f"MP{i}.C{x}"),
        tgt = mp_clusters.iloc[:,i+1], #.apply(lambda x: f"MP{i+1}.C{x}"),
    ))
    df['count'] = 1
    df = df.groupby(['src', 'tgt']).sum().reset_index()
    links = links.append(df, ignore_index=True)
    print(links)


# add one level containing the annotated clusters
i = nlevel_use - 1

df = pd.DataFrame(dict(
    src = mp_clusters.iloc[:,i], #.apply(lambda x: f"MP{i}.C{x}"),
    tgt = adata.obs.annotation,
))
df['count'] = 1
df = df.groupby(['src', 'tgt']).sum().reset_index()
links = links.append(df, ignore_index=True)
print(links)

links.to_csv('links_filtered.tsv', index=False, sep='\t')
