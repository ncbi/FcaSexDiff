#!/usr/bin/env python3
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
from ast import literal_eval
from tqdm import tqdm
from scipy import sparse
from scipy.stats import binom_test
from scipy.stats import fisher_exact
from collections import OrderedDict
from statsmodels.stats.multitest import multipletests

PSEUDO_EXPR = 1e-9
PSEUDO_COUNT = 1e-256

PADJ_CUTOFF_EXPR = 0.001
LFC_CUTOFF_EXPR = 1

exprfile = snakemake.input[0]
resol = snakemake.wildcards['resol']
outfile = snakemake.output[0]

print(exprfile)
print(resol)
print(outfile)

info = dict(
  cluster = (
      "Cluster name (cell type if annotation is used)"
  ),
  major_annotation = (
      "Major cell type in the cluster, useful if clusters are"
      " not by annotation"
  ),
  cell_count_female = (
      "Number of female cells in the cluster"
  ),
  cell_count_male = (
      "Number of male cells in the cluster"
  ),
  sample_cell_count_female = (
      "Number of female cells in the tissue"
  ),
  sample_cell_count_male = (
      "Number of male cells in the tissue"
  ),
  log2_count_bias = (
      "log2((NF + eps)/(NM + eps)) where NF (NM) is the normalized (by total"
      " cells in the tissue) female (male) cell count for the cluster and"
      f" eps = {PSEUDO_COUNT}"
  ),
  count_bias_padj = (
      "BH corrected pvalue from Binomial test for cluster-count-bias"
  ),
  count_bias_type = (
      "Cluster is female- or male-biased based on count_bias_padj < 0.05 and "
      "log2_count_bias > 1 or < -1 (2-fold change)"
  ),
  female_gene = (
      "Number of genes female biased in the cluster"
  ),
  male_gene = (
      "Number of genes male biased in the cluster"
  ),
  stats = (
      "Statistics for a gene within cluster"
  ),
  avg_all = (
      "Average expression of gene in all cells in cluster"
  ),
  avg_male = (
      "Average expression of gene in all male cells in cluster"
  ),
  avg_female = (
      "Average expression of gene in all female cells in cluster"
  ),
  avg_nz_all = (
      "Average expression of gene in all cells in cluster where it expressed"
  ),
  avg_nz_male = (
      "Average expression of gene in all male cells in cluster where it"
      " expressed"
  ),
  avg_nz_female = (
      "Average expression of gene in all female cells in cluster where it"
      " expressed"
  ),
  frac_female = (
      "Fraction of female cells where gene is expressed (atleast one UMI)"
  ),
  frac_male = (
      "Fraction of male cells where gene is expressed (atleast one UMI)"
  ),
  log2fc = (
      "Log2 of fold change of average expression in favor of female cells"
      " against male cells"
  ),
  log2fc_scanpy = (
      "Log2 of fold change of average expression in favor of female cells"
      " against male cells as computed by scanpy.tl.rank_genes_groups"
  ),
  padj = (
      "BH corrected pvalue from Wilcoxon test for male cells vs female cells"
      " expression as computed by scanpy.tl.rank_genes_groups"
  ),
  bias = (
      f"Gene is female- or male-biased based on padj < {PADJ_CUTOFF_EXPR}"
      f" and log2fc > {LFC_CUTOFF_EXPR} or < -{LFC_CUTOFF_EXPR}"
      f" ({2**LFC_CUTOFF_EXPR}-fold change)"
  ),
  bias_scanpy = (
      f"gene is female- or male-biased based on padj < {PADJ_CUTOFF_EXPR}"
      f" and log2fc_scanpy > {LFC_CUTOFF_EXPR} or < -{LFC_CUTOFF_EXPR}"
      f" ({2**LFC_CUTOFF_EXPR}-fold change)"
  ),
  chr = (
      "Chromosome where gene is located"
  ),
  symbol = "Gene symbol",
  FBgn = "Flybase ID of gene",
  female_cls = (
      "Number of clusters where gene is female-biased"
  ),
  male_cls = (
      "Number of clusters where gene is male-biased"
  ),
  umi_tissue = (
      "Average number of UMIs for the gene in the whole tissue"
  ),
  nz_umi_tissue = (
      "Average number of UMIs for the gene in cells in the whole tissue where"
      " it expressed"
  ),
  norm_tissue = (
      "Normalized expression of gene (before log1p) in the whole tissue"
  ),
  nz_norm_tissue = (
      "Normalized expression of gene (before log1p) in cells in whole tissue"
      " where it expressed"
  ),
)


print(info)


def merge_dict(male, female):
    male = {} if male == 0 else literal_eval(male)
    female = {} if female == 0 else literal_eval(female)

    # return a tuple (#female, #male) for each key
    return {
        k: (female.get(k, 0), male.get(k, 0))
        for k in set(female) | set(male)
    }


def compute_sexbiased_counts_male_to_female(meta):
    meta = meta[["sex", "cluster", "annotation"]]

    # count totals cells in sample for each sex
    total  = (
        meta[['sex']].assign(sample_cell_count=1)
        .groupby('sex')
        .count()
    )

    # compute cell_fraction (normalized cell count) for each sex
    # also keep count of annotated cell types in each cluster 
    res = (
        meta.assign(cell_count=1)
        .groupby(['cluster', 'sex'])
        .agg({
            'cell_count': 'sum',
            'annotation': lambda x: str(
                {k:v for k,v in x.value_counts().items() if v > 0}
            )
        })
        .unstack('sex', fill_value=0)
        .stack('sex')
        .merge(total, left_on='sex', right_on='sex', left_index=True,
               right_index=True)
        .assign(cell_fraction = lambda df: df.cell_count / df.sample_cell_count)
        .unstack('sex', fill_value=0)
    )

    # flatten columns for both sexes
    res.columns = ['_'.join(x) for x in res.columns]

    # flattening annotations require merging two dictionaries for counts
    res['annotation'] = res.apply(
        lambda x: merge_dict(x['annotation_male'], x['annotation_female']),
        axis = 1
    )

    # major annotation is one which has maximum sum of (#female, #male) tuple
    res['major_annotation'] = res['annotation'].apply(
        lambda x: max(x, key=lambda k: sum(x.get(k)))
    )

    res['cluster_type'] = res.apply(
        lambda x: 'male_only' if x['cell_count_female'] == 0 else
                  'female_only' if x['cell_count_male'] == 0 else
                  'has_both_sex',
        axis = 1
    )

    res['log2_count_bias'] = np.log2(
        (PSEUDO_COUNT + res.cell_fraction_female) /
        (PSEUDO_COUNT + res.cell_fraction_male)
    )

    def test_significance(x):
        nsucc = x['cell_count_male']
        nfail = x['cell_count_female']
        sample_succ = x['sample_cell_count_male']
        sample_fail = x['sample_cell_count_female']
        prob = sample_succ / (sample_succ + sample_fail)
        binom_res = binom_test([nsucc, nfail], p = prob, alternative = 'two-sided')
        fisher_res = fisher_exact([[nsucc, sample_succ - nsucc],
                                   [nfail, sample_fail - nfail]],
                                  alternative = 'two-sided')[1]
        return pd.Series({'pval_binom': binom_res, 'pval_fisher':fisher_res})

    res = res.merge(res.apply(test_significance, axis=1),
                    left_index=True, right_index=True)

    # pvals_corrected 2nd item in returned tuple from multipletests
    res['padj_binom'] = multipletests(res['pval_binom'], method='fdr_bh')[1]
    res['padj_fisher'] = multipletests(res['pval_fisher'], method='fdr_bh')[1]

    return res



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


def do_all(exprfile, reosl, outfile):
    adata = ad.read_h5ad(exprfile)
    adata = adata[adata.obs.sex.isin(["female", "male"])]

    assert("cluster" not in adata.obs.columns)
    adata.obs = adata.obs.assign(cluster = lambda df: df[resol])
    print(adata)

    count_bias = compute_sexbiased_counts_male_to_female(adata.obs)
    print(count_bias)

    clusters = adata.obs.cluster.unique()
    expr_bias = pd.concat(
        [
            get_sexbiased_expression_in_cluster(adata[adata.obs.cluster == x])
            for x in clusters
        ],
        keys = clusters,
        names = ['cluster', 'stats'],
        axis = 1
    )
    expr_bias.columns = expr_bias.columns.swaplevel()
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


    gene_meta = (
        adata.var.loc[gene_bias.index, :]
        .assign(female_cls = (gene_bias > 0).sum(axis=1))
        .assign(male_cls = (gene_bias < 0).sum(axis=1))
    )

    cluster_meta = (
        count_bias.loc[gene_bias.columns, :]
        .assign(female_gene = (gene_bias > 0).sum(axis=0))
        .assign(male_gene = (gene_bias < 0).sum(axis=0))
    )

    layers = OrderedDict()
    for stat in expr_bias.columns.get_level_values("stats").unique():
        layers[stat] = expr_bias[stat].values

    uns = OrderedDict()
    uns["info"] = info

    layers["bias_scanpy"] = sparse.csr_matrix(scanpy_bias)

    adata = ad.AnnData(
        sparse.csr_matrix(gene_bias),
        obs = gene_meta,
        var = cluster_meta,
        layers = layers,
        uns = uns,
    )
    print(adata)

    adata.write_h5ad(outfile)

do_all(exprfile, resol, outfile)

