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

exprfile = snakemake.input[0]
resol = snakemake.wildcards['resol']
newexprfile = snakemake.output[0]
markerfile = snakemake.output[1]

print(exprfile)
print(resol)
print(newexprfile)
print(markerfile)

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



def find_markers(adata):

    # use normalized expression all through
    adata.X = adata.layers['norm']

    error = False

    markers = pd.DataFrame()

    # scanpy assumes expression is already in log scale
    adata.X = np.log1p(adata.X)
    try:
        sc.tl.rank_genes_groups(
            adata, 'cluster', groups="all", reference="rest",
            pts=True, method='wilcoxon', key_added = "wilcoxon"
        )
    except AttributeError as e:
        print('\tSample doesn\'t contain one group')
        print(e)
        error = True
    except IndexError:
        print('\tNot enough cells in one group')
        error = True
    except ZeroDivisionError:
        print('\tZero division!')
        error = True
    except ValueError as e:
        print('\tSample doesn\'t contain one group')
        print(e)
        error = True
    if not error:
        frames = []
        clusters = adata.obs.cluster.unique()
        for cls in clusters:
            frames.append(pd.DataFrame({
                # t is a rec.array which can be indexed by col name
                col : [t[cls] for t in adata.uns['wilcoxon'][rec]]
                for col,rec in [
                    ('symbol','names'),
                    ('log2fc', 'logfoldchanges'),
                    ('pval', 'pvals'),
                    ('padj', 'pvals_adj'),
                    ('score', 'scores'),
                ]
            }).set_index('symbol'))
        markers = pd.concat(frames, axis=1, keys = clusters, names =
                            ['cluster', 'stats'])
        markers.columns = markers.columns.swaplevel()

        print(markers)
        # rows in pts are in different order of symbols, handle separately
        pts = pd.DataFrame(adata.uns['wilcoxon']['pts'])
        pts = pd.concat([pts], axis=1, keys=["pct"], names=["stats", "cluster"])

        markers = (
            pd.concat([markers, pts], axis=1)
            .sort_index(1)
        )
        print(markers)

    return adata, markers


def do_all(exprfile, reosl, newexprfile, markerfile):
    adata = ad.read_h5ad(exprfile)
    adata = adata[adata.obs.sex.isin(["female", "male"])]

    assert("cluster" not in adata.obs.columns)
    adata.obs = adata.obs.assign(cluster = lambda df: df[resol].astype(str))
    adata.obs = adata.obs.assign(cluster = lambda df: df["cluster"].apply(lambda x: f"C{x}"))
    print(adata)
    print(adata.obs[["cluster"]])


    count_bias = compute_sexbiased_counts_male_to_female(adata.obs)
    print(count_bias)

    res_ad, markers = find_markers(adata)
    print(res_ad)
    print(markers)

    is_marker = (markers['padj'] < 0.05) & (markers['log2fc'] > 1)

    gene_meta = (
        adata.var.loc[markers.index, :]
        .assign(marker_cls = is_marker.sum(axis=1))
    )

    cluster_meta = (
        count_bias.loc[is_marker.columns, :]
    )

    layers = OrderedDict()
    for stat in markers.columns.get_level_values("stats").unique():
        layers[stat] = markers[stat].values

    ad.AnnData(
        sparse.csr_matrix(is_marker),
        obs = gene_meta,
        var = cluster_meta,
        layers = layers,
    ).write_h5ad(markerfile)

    res_ad.write_h5ad(newexprfile)

do_all(exprfile, resol, newexprfile, markerfile)

