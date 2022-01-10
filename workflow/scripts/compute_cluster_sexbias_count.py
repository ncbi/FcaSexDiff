#!/usr/bin/env python3
import json
import pandas as pd
import numpy as np
import anndata as ad
from scipy.stats import binom_test
from scipy.stats import fisher_exact
from scipy.stats import mannwhitneyu
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from numpyencoder import NumpyEncoder # Needed for numpy -> json conversion

from define_sexbias_cutoffs import *

exprfile = snakemake.input[0]
resol = snakemake.wildcards['resol']
outfile = snakemake.output[0]

print(exprfile)
print(resol)
print(outfile)

def merge_annotations(female, male):
    if female == 0: female = "{}"
    if male == 0: male = "{}"
    return json.dumps(
        pd.concat(
            {
                "f": pd.DataFrame(json.loads(female)),
                "m": pd.DataFrame(json.loads(male)),
            },
            names = ["sex", "sample"],
            axis=1,
            sort=True,
        )
        .T
        .fillna(0)
        .groupby("sex")
        .agg(lambda x: {
            k:v for k,v in dict(x.reset_index("sex", drop=True)).items()
            if v > 0
        })
        .agg(lambda x: dict(x))
        .to_dict(),
        cls=NumpyEncoder
    )

def get_major_annotation(x):
    return (
        pd.concat(
            {k: pd.DataFrame(v).T for k, v in json.loads(x).items()},
            names = ["annotation", "sex"],
            sort=True,
            axis=0
        )
        .fillna(0)
        .reset_index()
        .melt(id_vars=["annotation", "sex"],
              var_name="sample",
              value_name="cell_count")
        .groupby("annotation")
        .agg("sum")
        .idxmax()
        ["cell_count"]
     )

def compute_sexbiased_counts_male_to_female(meta):
    meta = (
        meta[["sex", "cluster", "annotation", "sample"]]
        .rename(columns={"annotation": "annotations"})
        .set_index("sex")
    )
    print(meta)



    def get_list(x):
        return [] if pd.isna(x) or (x == 0) else list(json.loads(x).values())

    def process_single_sex(meta):
        tissue_count = meta.shape[0]

        # count totals cells in each replicate in tissue
        replicate_counts  = (
            meta[["sample"]]
            .assign(tissue_rep_counts=1)
            .groupby(['sample'])
            .agg("sum")
        )
        print(replicate_counts)

        def agg_normal(x):
            x = x.reset_index(level=['cluster'], drop=True)
            return json.dumps(dict(x), cls=NumpyEncoder)

        def agg_nested(x):
            x = x.reset_index(level=['cluster'], drop=True)
            x = x.apply(lambda r: json.loads(r))
            return json.dumps(dict(x), cls=NumpyEncoder)

        # compute cluster_rep_fracs (normalized cell count
        # for each replicate within the tissue separately)
        # also keep count of annotated cell types in each cluster 
        res = (
            meta.assign(cluster_rep_counts=1)
            .groupby(['sample', 'cluster'])
            .agg({
                'cluster_rep_counts': 'sum',
                'annotations': lambda x: json.dumps({
                    k:v for k,v in x.value_counts().items() if v > 0
                })
            })
            .merge(replicate_counts, left_index=True, right_index=True)
            .fillna({
                "annotations" : "{}",
                "cluster_rep_counts" : 0.0,
                "tissue_rep_counts" : 0.0,
            })
            .assign(cluster_rep_fracs = lambda df: (
                df.cluster_rep_counts / df.tissue_rep_counts
            ))
            .groupby(['cluster'])
            # now save each entry as a dictionary keyed by replicate
            .agg({
                "annotations" : agg_nested,
                "cluster_rep_fracs": agg_normal,
                "cluster_rep_counts": agg_normal,
                "tissue_rep_counts": agg_normal,
            })
        )

        print(res)

        # also add information combining all replicates together

        res_reps_together = (
            meta[["cluster"]]
            .assign(cluster_count=1)
            .groupby(['cluster'])
            .agg('sum')
            .assign(tissue_count = tissue_count)
            .assign(cluster_frac = lambda df: df.cluster_count / df.tissue_count)
        )
        print(res_reps_together)

        res = (
            res.merge(res_reps_together, left_index=True, right_index=True)
            .assign(cluster_rep_fracs_mean = lambda df: (
                df.cluster_rep_fracs.apply(lambda x: np.mean(get_list(x)))
            ))
            .assign(cluster_rep_fracs_sd = lambda df: (
                df.cluster_rep_fracs.apply(lambda x: np.std(get_list(x), ddof = 1))
            ))
        )

        print(res)

        return res

    female_meta = meta.xs("female")
    male_meta = meta.xs("male")

    res = pd.concat(
        [
            process_single_sex(female_meta),
            process_single_sex(male_meta),
        ],
        axis = 1,
        keys = ["female", "male"],
        names = ["sex", "stats"],
        sort = True,
    )
    # it may happen that all replicates of a single sex are
    # are not present, then fill up with na values as appropriate
    na_values = pd.Series(dict(
        annotations            = "{}",
        tissue_count           = 0.0,
        cluster_count          = 0.0,
        cluster_frac           = 0.0,
        cluster_rep_fracs      = "{}",
        cluster_rep_counts     = "{}",
        tissue_rep_counts      = "{}",
        cluster_rep_fracs_mean = 0.0,
        cluster_rep_fracs_sd   = 0.0,
    ))
    # fillna using a dataframe whose index matched with the
    # columns (multi-index) of the dataframe being filled
    res = res.fillna(pd.concat(
        [na_values, na_values],
        keys=["female", "male"],
        names=["sex", "stats"]
    ))
    print(res)

    # move stats at the top of column index
    # and rearrange columns so that stats for the two sexes are together
    res = res.swaplevel(axis=1).sort_index(axis=1)

    # flatten columns for both sexes
    res.columns = ['_'.join(x) for x in res.columns]

    # flattening annotations require merging two dictionaries for counts
    res['annotations'] = res.apply(
        lambda x: merge_annotations(x['annotations_female'], x['annotations_male']),
        axis = 1
    )
    print(res[['annotations']])


    # major annotation is one which has maximum sum from male and female cells
    res['major_annotation'] = res['annotations'].apply(get_major_annotation)

    res['cluster_type'] = res.apply(
        lambda x: 'male_only' if x['cluster_count_female'] == 0 else
                  'female_only' if x['cluster_count_male'] == 0 else
                  'has_both_sex',
        axis = 1
    )


    def test_significance(x):
        female_cell_fracs = get_list(x['cluster_rep_fracs_female'])
        male_cell_fracs = get_list(x['cluster_rep_fracs_male'])

        nsucc = x['cluster_count_female']
        nfail = x['cluster_count_male']
        sample_succ = x['tissue_count_female']
        sample_fail = x['tissue_count_male']

        prob = sample_succ / (sample_succ + sample_fail)

        pval_binom = binom_test([nsucc, nfail], p = prob, alternative = 'two-sided')
        pval_fisher = fisher_exact([[nsucc, sample_succ - nsucc],
                                    [nfail, sample_fail - nfail]],
                                   alternative = 'two-sided')[1]
        pval_wilcox = mannwhitneyu(female_cell_fracs, male_cell_fracs)[1]
        pval_ttest = ttest_ind(female_cell_fracs, male_cell_fracs,
                               equal_var=False)[1]
        # t-test returns nan when standard deviation is zero
        if np.isnan(pval_ttest):
            pval_ttest = 1.0

        log2bias = np.log2(
            (PSEUDO_COUNT + np.mean(female_cell_fracs)) /
            (PSEUDO_COUNT + np.mean(male_cell_fracs))
        )

        return pd.Series({
            'log2_count_bias': log2bias,
            'pval_binom'     : pval_binom,
            'pval_fisher'    : pval_fisher,
            'pval_wilcox'    : pval_wilcox,
            'pval_ttest'     : pval_ttest,
        })

    res = res.merge(res.apply(test_significance, axis=1),
                    left_index=True, right_index=True)
    print(res)

    # pvals_corrected 2nd item in returned tuple from multipletests
    res['padj_binom'] = multipletests(res['pval_binom'], method='fdr_bh')[1]
    res['padj_fisher'] = multipletests(res['pval_fisher'], method='fdr_bh')[1]
    res['padj_wilcox'] = multipletests(res['pval_wilcox'], method='fdr_bh')[1]
    res['padj_ttest'] = multipletests(res['pval_ttest'].fillna(1.0), method='fdr_bh')[1]

    print(res)

    n_female_reps = len(female_meta["sample"].unique())
    n_male_reps = len(male_meta["sample"].unique())

    padj_use = (
        "padj_wilcox" if ((n_female_reps > 1) and (n_male_reps > 1)) else
        "padj_binom"
    )

    def encode_count_bias(log2bias, padj, cluster_type):
        return (
            "FemaleOnly" if cluster_type == "female_only" else
            "MaleOnly" if cluster_type == "male_only" else
            "FemaleSignificant" if ((log2bias > LFC_CUTOFF_COUNT) &
                                    (padj < PADJ_CUTOFF_COUNT)) else
            "FemaleNonsignificant" if log2bias > LFC_CUTOFF_COUNT else
            "MaleSignificant" if ((log2bias < -LFC_CUTOFF_COUNT) &
                                  (padj < PADJ_CUTOFF_COUNT)) else
            "MaleNonsignificant" if log2bias < -LFC_CUTOFF_COUNT else
            "Unbiased"
        )

    print(res[["log2_count_bias", padj_use, "cluster_type"]])

    res = (
        res.assign(count_bias_padj = lambda df: df[padj_use]) #.apply(pvalstr))
        .assign(count_bias_type = lambda df: df.apply(lambda x: encode_count_bias(
            x["log2_count_bias"], x["count_bias_padj"], x["cluster_type"]
        ), axis=1))
    )
    print(res[["count_bias_type", "log2_count_bias", "count_bias_padj", "cluster_type"]])

    return res


def do_all(exprfile, resol, outfile):
    adata = ad.read_h5ad(exprfile)

    assert("cluster" not in adata.obs.columns)
    meta_data = (
        adata.obs
        # anndata saves samples as categorical variable but our code
        # assumes that it is of string type, not categorical
        .assign(cluster = lambda df: df[resol].astype(str))
        .assign(sample = lambda df: df["sample"].astype(str))
    )
    print(adata)

    count_bias = compute_sexbiased_counts_male_to_female(meta_data)
    print(count_bias)

    count_bias.to_hdf(outfile, key='bias')

do_all(exprfile, resol, outfile)

