#!/usr/bin/env python3
import re
import json
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
from scipy import sparse
from collections import OrderedDict


from define_sexbias_cutoffs import *
from get_bias_info import *

exprfile = snakemake.input["expr"]
cellsfile = snakemake.input["cells"]
count_bias_file = snakemake.input["cntbias"]
expr_bias_file = snakemake.input["exprbias"]
marker_file = snakemake.input["mark"]
resol = snakemake.wildcards['resol']
cell_filter = snakemake.wildcards['cellfilt']
outfile = snakemake.output[0]

print(exprfile)
print(cellsfile)
print(count_bias_file)
print(expr_bias_file)
print(marker_file)
print(resol)
print(cell_filter)
print(outfile)


def get_tissue_stats(adata, name):
    resols = ["annotation"] + [
        x for x in adata.obs.columns if re.match("L\d+.\d+", x)
    ]
    print(resols)
    num_clusters = {x:len(adata.obs[x].unique()) for x in resols}
    print(num_clusters)

    return pd.DataFrame(dict(
        n_gene               = [adata.shape[1]],
        n_cluster            = [json.dumps(num_clusters)],
        tissue_count         = [adata.shape[0]],
        tissue_count_female  = [adata.obs.query("sex == 'female'").shape[0]],
        tissue_count_male    = [adata.obs.query("sex == 'male'").shape[0]],
        tissue_rep_counts_female = [json.dumps(
            adata.obs.query("sex == 'female'")["sample"].astype(str)
            .value_counts(sort=False).sort_index().to_dict()
        )],
        tissue_rep_counts_male = [json.dumps(
            adata.obs.query("sex == 'male'")["sample"].astype(str)
            .value_counts(sort=False).sort_index().to_dict()
        )],
    ), index=[name])


def do_all(exprfile, cellsfile, resol, count_bias_file, expr_bias_file, marker_file, outfile):

    adata = ad.read_h5ad(exprfile)

    # first discard cells that have 'mix' sex
    adata = adata[adata.obs.sex.isin(["female", "male"])]
    tissue_stats = get_tissue_stats(adata, "AllCells")

    cells = pd.read_csv(cellsfile)
    adata = adata[cells.CellID]

    if cell_filter != "AllCells":
        tissue_stats_filtered = get_tissue_stats(adata, cell_filter)
        tissue_stats = pd.concat([tissue_stats, tissue_stats_filtered])
    print(tissue_stats)

    count_bias = pd.read_hdf(count_bias_file)
    print(count_bias)

    expr_bias = pd.read_hdf(expr_bias_file)
    expr_bias.columns = expr_bias.columns.swaplevel()
    # there is some issue with saving dataframe with nan in index to hdf
    # so dmel_ was prefixed to gene names in the index, now remove that prefix
    expr_bias.index = expr_bias.index.map(lambda x: x[len("dmel_"):])
    print(expr_bias)

    padj = expr_bias['padj']
    log2fc_scanpy = expr_bias['log2fc_scanpy']
    log2fc = expr_bias['log2fc']

    gene_bias = np.sign(log2fc)
    gene_bias[~((padj < PADJ_CUTOFF_EXPR) & (abs(log2fc) > LFC_CUTOFF_EXPR))] = 0
    print(gene_bias)

    scanpy_bias = np.sign(log2fc_scanpy)
    scanpy_bias[~((padj < PADJ_CUTOFF_EXPR) & (abs(log2fc_scanpy) > LFC_CUTOFF_EXPR))] = 0
    print(scanpy_bias)

    marker_info = pd.read_hdf(marker_file, key="info")
    top_markers = pd.read_hdf(marker_file, key="top")
    marker_info.index = marker_info.index.map(
        lambda x: 'nan' if pd.isna(x) else x
    )


    print(top_markers)

    print(count_bias)
    count_bias = count_bias.merge(
            top_markers, left_index=True,
            right_index=True
    )
    print(count_bias)

    gene_order = gene_bias.index
    cluster_order = gene_bias.columns

    print(sum(pd.isna(gene_order)))
    print(sum(pd.isna(marker_info.index)))


    # make sure that the columns in marker_info are in cluster order
    marker_info = marker_info.reindex(cluster_order, axis=1, level='cluster')

    print(gene_order)
    marker_info = marker_info.loc[gene_order, :]
    print(marker_info)

    is_marker = (
        (marker_info['marker_padj'] < PADJ_CUTOFF_MARKER) &
        (marker_info['marker_log2fc'] > LFC_CUTOFF_MARKER)
    )

    gene_meta = (
        adata.var.loc[gene_order, :]
        .assign(female_cls = (gene_bias > 0).sum(axis=1))
        .assign(male_cls = (gene_bias < 0).sum(axis=1))
    )

    cluster_meta = (
        count_bias.loc[cluster_order, :]
        .assign(female_gene = (gene_bias > 0).sum(axis=0))
        .assign(male_gene = (gene_bias < 0).sum(axis=0))
    )

    layers = OrderedDict()
    for stat in expr_bias.columns.get_level_values("stats").unique():
        layers[stat] = expr_bias[stat].values

    layers["bias_scanpy"] = sparse.csr_matrix(scanpy_bias)
    layers["is_marker"] = sparse.csr_matrix(is_marker.astype(int))

    for stat in marker_info.columns.get_level_values("stats").unique():
        layers[stat] = marker_info[stat].values

    uns = OrderedDict()
    uns["info"] = info
    uns["stats"] = tissue_stats

    adata = ad.AnnData(
        sparse.csr_matrix(gene_bias),
        obs = gene_meta,
        var = cluster_meta,
        layers = layers,
        uns = uns,
    )
    print(adata)

    adata.write_h5ad(outfile)

do_all(exprfile, cellsfile, resol, count_bias_file, expr_bias_file, marker_file, outfile)

