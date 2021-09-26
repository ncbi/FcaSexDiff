#!/usr/bin/env python3
import pandas as pd
import numpy as np
from tqdm import tqdm
import scanpy as sc
import anndata as ad
from scipy.stats import binom_test
from statsmodels.stats.multitest import multipletests
from scipy.stats import fisher_exact
from ast import literal_eval
import sys, os

PSEUDO_EXPR = 1e-9
PSEUDO_COUNT = 1e-256

exprfile = snakemake.input[0]
tissue = snakemake.wildcards['tissue']
resol = snakemake.wildcards['resol']
outfile = snakemake.output[0]

print(exprfile)
print(tissue)
print(resol)
print(outfile)

adata = ad.read_h5ad(exprfile)
adata = adata[adata.obs.sex.isin(["female", "male"])]
print(adata)

def merge_dict(male, female):
    male = {} if male == 0 else literal_eval(male)
    female = {} if female == 0 else literal_eval(female)
    # return a tuple (#female, #male) for each key
    return {
        k: (female.get(k, 0), male.get(k, 0))
        for k in set(female) | set(male)
    }


def compute_sexbiased_counts_male_to_female(meta):
    use_cols = ["sex", resol]
    if resol != "annotation":
        use_cols += ["annotation"]
    meta = (
        meta[use_cols]
        .assign(cluster = lambda df: df[resol])
    )
    print(meta)
    # count totals cells in sample for each sex
    total  = (
        meta[['sex']].assign(sample_cell_count=1)
        .groupby('sex')
        .count()
    )
    print(total)

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
    print(res)

    res.columns = ['_'.join(x) for x in res.columns]

    res['annotation'] = res.apply(
        lambda x: merge_dict(x['annotation_male'], x['annotation_female']),
        axis = 1
    )
    # major annotation is one which has maximum sum of (#female, #male) tuple
    res['major_annotation'] = res['annotation'].apply(
        lambda x: max(x, key=lambda k: sum(x.get(k)))
    )

    res['cluster_type'] = res.apply(lambda x: 'male_only' if x['cell_count_female'] ==
                               0 else 'female_only' if x['cell_count_male'] ==
                               0 else 'has_both_sex', axis = 1)

    res['log2_male_to_female'] = np.log2((PSEUDO_COUNT + res.cell_fraction_male) /
                                         (PSEUDO_COUNT + res.cell_fraction_female))
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
    return res.reset_index()



def get_sexbiased_expression_in_cluster(x):
    print(x)
    # for computing average consider only non zero values
    nz = x.replace(0.0, np.NaN)
    nz_avg_male = nz.loc[:, nz.columns.get_level_values('sex')=='male'].mean(skipna=True, axis=1)
    nz_avg_female = nz.loc[:, nz.columns.get_level_values('sex')=='female'].mean(skipna=True, axis=1)
    nz_avg_all = nz.mean(skipna=True, axis=1)

    avg_male = x.loc[:, x.columns.get_level_values('sex')=='male'].mean(skipna=True, axis=1)
    avg_female = x.loc[:, x.columns.get_level_values('sex')=='female'].mean(skipna=True, axis=1)
    avg_all = x.mean(skipna=True, axis=1)

    res = pd.DataFrame({
        'avg_nz_male': nz_avg_male,
        'avg_nz_female': nz_avg_female,
        'avg_nz_all': nz_avg_all,
        'avg_male': avg_male,
        'avg_female': avg_female,
        'avg_all': avg_all,
        'myl2fc': np.log2((avg_male+PSEUDO_EXPR)/(avg_female+PSEUDO_EXPR)),
    }).fillna(0).reset_index()

    res_de = pd.DataFrame({
        'log2fc': [0]*x.shape[0],
        'pval': [1.0]*x.shape[0],
        'padj': [1.0]*x.shape[0]
    }, index = x.index).reset_index()

    adata = sc.AnnData(
        np.log1p(x.values.T),
        obs=x.columns.to_frame(index=False),
        var=x.index.to_frame(index=False).set_index('symbol')
    )
    error = False
    try:
        sc.tl.rank_genes_groups(adata, 'sex', groups=['male'], pts=False,
                                reference='female', method='wilcoxon',
                                key_added = "wilcoxon")
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
        res_de = pd.DataFrame({
            col: [t[0] for t in adata.uns['wilcoxon'][rec]]
            for col,rec in [('symbol','names'), ('log2fc', 'logfoldchanges'),
                            ('pval', 'pvals'), ('padj', 'pvals_adj')]
        })
    res = res.merge(res_de)
    print(res['symbol avg_all avg_male avg_female log2fc padj'.split()])

    return res.set_index(['symbol'])




count_bias = compute_sexbiased_counts_male_to_female(adata.obs)
print(count_bias)

quit()

count_bias.to_hdf(outfile, key='count_bias')


clusters = expr.columns.unique('cluster')
print(clusters)

expr_bias = pd.concat(
    [
        get_sexbiased_expression_in_cluster(
            expr.loc[:,(expr.columns.get_level_values('cluster') == x)]
        )
        for x in clusters
    ],
    keys = clusters,
    names = ['cluster', 'stats'],
    axis = 1
).loc[avg.symbol, :]
print(expr_bias)

expr_bias.index = pd.MultiIndex.from_frame(avg)

expr_bias.to_hdf(outfile, key='expr_bias')

