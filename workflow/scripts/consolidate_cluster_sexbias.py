#!/usr/bin/env python3
import pandas as pd
import numpy as np
import anndata as ad
from scipy import sparse
from collections import OrderedDict


from define_sexbias_cutoffs import *
from get_bias_info import *

exprfile = snakemake.input["expr"]
count_bias_file = snakemake.input["cntbias"]
expr_bias_file = snakemake.input["exprbias"]
marker_file = snakemake.input["mark"]
resol = snakemake.wildcards['resol']
cell_filter = snakemake.wildcards['cellfilt']
outfile = snakemake.output[0]

print(exprfile)
print(outfile)


def do_all(exprfile, resol, count_bias_file, expr_bias_file, marker_file, outfile):

    adata = ad.read_h5ad(exprfile)
    tissue_stats = adata.uns["stats"]

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

    # make sure that the columns in marker_info are in cluster order
    marker_info = marker_info.reindex(cluster_order, axis=1, level='cluster')

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

do_all(exprfile, resol, count_bias_file, expr_bias_file, marker_file, outfile)

