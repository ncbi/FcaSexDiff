#!/usr/bin/env python3
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad

from define_sexbias_cutoffs import *

exprfile = snakemake.input[0]
cellsfile = snakemake.input[1]
resol = snakemake.wildcards['resol']
outfile = snakemake.output[0]

print(exprfile)
print(cellsfile)
print(resol)
print(outfile)

def get_sexbiased_expression_in_cluster(adata):

    # use normalized expression all through
    both = np.expm1(adata.X)
    male = np.expm1(adata[adata.obs.sex == "male"].X)
    female = np.expm1(adata[adata.obs.sex == "female"].X)

    def all_mean(mat):
        if mat.shape[0] == 0:
            # empty matrix, return a vector of zeros
            return np.zeros(mat.shape[1])
        return np.ravel(mat.mean(axis=0))

    def non_zero_mean(mat):
        if mat.shape[0] == 0:
            # empty matrix, return a vector of zeros
            return np.zeros(mat.shape[1])
        # force number of non-zeros in row be a positive
        # number to avoid divide by error
        num_nz = (mat!=0).sum(axis=0)
        num_nz[num_nz == 0] = 1
        return np.ravel(np.true_divide(mat.sum(axis=0),num_nz))

    def fraction_non_zero(mat):
        if mat.shape[0] == 0:
            # empty matrix, return a vector of zeros
            return np.zeros(mat.shape[1])
        # force number of non-zeros in row be a positive
        # number to avoid divide by error
        num_nz = (mat!=0).sum(axis=0)
        return np.ravel(num_nz/mat.shape[0])

    res = pd.DataFrame({
        'avg_nz_all'   : non_zero_mean(both),
        'avg_nz_female': non_zero_mean(female),
        'avg_nz_male'  : non_zero_mean(male),
        'avg_all'      : all_mean(both),
        'avg_female'   : all_mean(female),
        'avg_male'     : all_mean(male),
        'log2fc'       : np.log2((all_mean(female)+PSEUDO_EXPR)/(all_mean(male)+PSEUDO_EXPR)),
    }, index = adata.var.index)

    ngene = both.shape[1]

    res_de = pd.DataFrame({
        'log2fc_scanpy' : np.zeros(ngene),
        'pval'          : np.zeros(ngene) + 1.0,
        'padj'          : np.zeros(ngene) + 1.0,
        'score'         : np.zeros(ngene),
        'frac_female'    : fraction_non_zero(female),
        'frac_male'      : fraction_non_zero(male),
    }, index = adata.var.index)

    error = False

    # workaround for scanpy bug:
    # https://bytemeta.vip/repo/scverse/scanpy/issues/2239
    if "log1p" in adata.uns.keys():
        adata.uns["log1p"]["base"] = None

    # scanpy assumes expression is already in log scale

    try:
        sc.tl.rank_genes_groups(
            adata, 'sex', groups=['female'], reference='male',
            pts=True, method='wilcoxon', key_added = "wilcoxon"
        )
    except AttributeError:
        print('\tSample doesn\'t contain one group')
        error = True
    except IndexError:
        print('\tNot enough cells in one group')
        error = True
    except ZeroDivisionError:
        print('\tZero division!')
        error = True
    except ValueError:
        print('\tSample doesn\'t contain one group')
        error = True
    if not error:
        res_de = (
            pd.DataFrame({
                col: [t[0] for t in adata.uns['wilcoxon'][rec]]
                for col,rec in [
                    ('symbol','names'),
                    ('log2fc_scanpy', 'logfoldchanges'),
                    ('pval', 'pvals'),
                    ('padj', 'pvals_adj'),
                    ('score', 'scores'),
                ]
            })
            .set_index('symbol')
        )
        # rows in pts are in different order of symbols, handle separately
        pts = pd.DataFrame(adata.uns['wilcoxon']['pts'])
        pts.columns = [f"frac_{col}" for col in pts.columns]
        res_de = res_de.merge(pts, left_index=True, right_index=True)

    return res.merge(res_de, left_index=True, right_index=True)

def do_all(exprfile, cellsfile, resol, outfile):
    adata = ad.read_h5ad(exprfile)
    cells = pd.read_csv(cellsfile)
    adata = adata[cells.CellID]


    assert("cluster" not in adata.obs.columns)
    meta_data = (
        adata.obs
        # anndata saves samples as categorical variable but our code
        # assumes that it is of string type, not categorical
        .assign(cluster = lambda df: df[resol])
        .assign(sample = lambda df: df["sample"].astype(str))
    )
    print(adata)


    clusters = meta_data.cluster.unique()
    expr_bias = pd.concat(
        [
            get_sexbiased_expression_in_cluster(adata[meta_data.cluster == x])
            for x in clusters
        ],
        keys = clusters,
        names = ['cluster', 'stats'],
        axis = 1
    )
    print(expr_bias)

    # there is some issue with saving nan in index to hdf
    expr_bias.index = expr_bias.index.map(lambda x: f"dmel_{x}")

    expr_bias.to_hdf(outfile, key="bias")

do_all(exprfile, cellsfile, resol, outfile)

