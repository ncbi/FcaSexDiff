#!/usr/bin/env python
# coding: utf-8

import multiscale_phate as mp
import numpy as np
import pandas as pd
import anndata as ad

import scprep
import os


# ## 2. Loading and filter data
# use raw UMI counts

infile = snakemake.input[0] #"imports/expr_h5ad/expr_fat_body_stringent.h5ad"
outfile = snakemake.output[0] #"scraps/mphate_pkl/mphate_fat_body_stringent.pkl"
mphate_cluster_file = snakemake.output[1] #"scraps/mphate_clusters/mphate_clusters_fat_body_stringent.tsv"

adata = ad.read_h5ad(infile)
print(adata)

data = pd.DataFrame(
    adata.X.todense(),
    index = adata.obs_names,
    columns = adata.var_names
)

print(data)

# Now that we have loaded the data, we will remove cells with low transcript counts and unexpressed genes:

# FCA folks have already filtered cells
# data = scprep.filter.filter_library_size(data, cutoff=1000, keep_cells='above')

# However, we can filter genes if that helps m-phate
data = scprep.filter.filter_rare_genes(data)


# Finally, we will library size normalize and square root transform the expression data.

data_norm, libsize = scprep.normalize.library_size_normalize(data, return_library_size=True)

data_sqrt = np.sqrt(data_norm)

print(data_sqrt)


# ## 3. Creating multi-resolution embeddings and clusters with Multiscale PHATE

# Computing Multiscale PHATE tree involves three successive steps:
#
# 1. Building you Multiscale PHATE operator
# 2. Fitting your data with your operator to construct a diffusion condensation tree and identify salient resolutions for clustering and visualization
# 3. Embedding and clustering your data at selected resolutions


mp_op = mp.Multiscale_PHATE()

mp_op.fit(data_sqrt)


import pickle
pickle.dump(mp_op, open(outfile,"wb"))


# mphate saves resolution levels fine to coarse, we want reverse
levels = mp_op.levels[::-1]
print(levels)
for i, l in enumerate(levels):
    print(f'Level {i}, nclusters={len(set(mp_op.NxTs[l]))}')


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

mp_clusters.to_csv(mphate_cluster_file, sep="\t")
