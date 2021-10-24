#!/usr/bin/env python3
import json
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
from ast import literal_eval
from tqdm import tqdm
from scipy import sparse
from scipy.stats import binom_test
from scipy.stats import fisher_exact
from scipy.stats import mannwhitneyu
from scipy.stats import ttest_ind
from collections import OrderedDict
from statsmodels.stats.multitest import multipletests
from numpyencoder import NumpyEncoder # Needed for numpy -> jason conversion

PSEUDO_EXPR = 1e-9
PSEUDO_COUNT = 1e-256

PADJ_CUTOFF_COUNT = 0.05
LFC_CUTOFF_COUNT = 1

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
  annotations_female = (
    "Replicate-wise counts of annotated celltypes among female cells in the cluster"
  ),
  annotations_male = (
    "Replicate-wise counts of annotated celltypes among male cells in the cluster"
  ),
  annotations = (
    "Sex and replicate -wise counts of annotated celltypes among all cells in the cluster"
  ),
  major_annotation = (
      "Major annotated celltype in the cluster, useful if clusters are"
      " not by annotation"
  ),
  tissue_count_female = (
    "Number of female cells in the tissue (combining all replicates)"
  ),
  tissue_count_male = (
    "Number of male cells in the tissue (combining all replicates)"
  ),
  cluster_type = (
    "Classify cluster based on whether there is at least one cell in the cluster from either sex"
  ),
  cluster_count_female = (
      "Number of female cells in the cluster (combining all replicates)"
  ),
  cluster_count_male = (
      "Number of male cells in the cluster (combining all replicates)"
  ),
  tissue_rep_counts_female = (
      "Replicate-wise counts of female cells in the tissue"
  ),
  tissue_rep_counts_male = (
      "Replicate-wise counts of male cells in the tissue"
  ),
  cluster_rep_counts_female = (
      "Replicate-wise counts of female cells in the cluster"
  ),
  cluster_rep_counts_male = (
      "Replicate-wise counts of male cells in the cluster"
  ),
  cluster_frac_female = (
    "Fraction of cells in the cluster that are female (combining all replicates)"
  ),
  cluster_frac_male = (
    "Fraction of cells in the cluster that are male (combining all replicates)"
  ),
  cluster_rep_fracs_female = (
    "Replicate-wise fractions of cells in the cluster that are female"
  ),
  cluster_rep_fracs_male = (
    "Replicate-wise fractions of cells in the cluster that are male"
  ),
  cluster_rep_fracs_mean_female = (
    "Mean of replicate-wise fractions of cells in the cluster that are female "
    "(should differ from fraction of cells in the cluster that are female if the "
    "replicates have different number of cells)"
  ),
  cluster_rep_fracs_mean_male = (
    "Mean of replicate-wise fractions of cells in the cluster that are male "
    "(should differ from fraction of cells in the cluster that are male if the "
    "replicates have different number of cells)"
  ),
  cluster_rep_fracs_sd_female = (
    "Standard deviation of replicate-wise fractions of cells in the cluster that are female"
  ),
  cluster_rep_fracs_sd_male = (
    "Standard deviation of replicate-wise fractions of cells in the cluster that are male"
  ),
  pval_binom = (
    "Pvalue from bionamial test on the counts considering all replicates together"
  ),
  pval_fisher = (
    "Pvalue from Fisher's test on the counts considering all replicates together"
  ),
  pval_wilcox = (
    "Pvalue from Wilcoxon-Mann-Whitney U test on the replicate-wise fractions"
  ),
  pval_ttest = (
    "Pvalue from 2-sample t test on the replicate-wise fractions"
  ),
  padj_binom = (
      "BH corrected pvalue from binomial test on the counts considering all replicates together"
  ),
  padj_fisher = (
    "BH corrected pvalue from Fisher's test on the counts considering all replicates together"
  ),
  padj_wilcox = (
    "BH corrected pvalue from Wilcoxon-Mann-Whitney U test on the replicate-wise fractions"
  ),
  padj_ttest = (
    "BH corrected pvalue from s-sample t test on the replicate-wise fractions"
  ),
  log2_count_bias = (
      "log2((NF + eps)/(NM + eps)) where NF (NM) is the replicate-average "
      " fraction of cells in cluster that are female (male) and"
      f" eps = {PSEUDO_COUNT}"
  ),
  count_bias_padj = (
      "BH corrected pvalue from Binomial test for tissues having single "
      "replicate, otherwise from Wilcoxon-Mann-Whitney U test"
  ),
  count_bias_type = (
      f"Cluster is female- or male-biased based on count_bias_padj < "
      f"{PADJ_CUTOFF_COUNT} and log2_count_bias > {LFC_CUTOFF_COUNT} or "
      f"< -{LFC_CUTOFF_COUNT} ({2**LFC_CUTOFF_COUNT}-fold change)"
  ),
  female_gene = (
      "Number of genes female-biased in the cluster"
  ),
  male_gene = (
      "Number of genes male-biased in the cluster"
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


def merge_annotations(female, male):
    if female == 0: female = "{}"
    if male == 0: male = "{}"
    return json.dumps(
        pd.concat(
            {
                "f": pd.DataFrame(json.loads(female)),
                "m": pd.DataFrame(json.loads(male)),
            },
            names = ["sex", "replicate"],
            axis=1
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
            axis=0
        )
        .fillna(0)
        .reset_index()
        .melt(id_vars=["annotation", "sex"],
              var_name="replicate",
              value_name="cell_count")
        .groupby("annotation")
        .agg("sum")
        .idxmax()
        ["cell_count"]
     )

def compute_sexbiased_counts_male_to_female(meta):
    meta = (
        meta[["sex", "cluster", "annotation", "sample_id"]]
        .rename(columns={"annotation": "annotations", "sample_id": "replicate"})
        #.assign(replicate = lambda df: df.replicate.str.split("__").str[0])
        .assign(replicate = lambda df: df.replicate.str.split("__").str[1].str.split("_").str[0])
        .set_index("sex")
    )
    print(meta)

    def get_list(x):
        return [] if pd.isna(x) or (x == 0) else list(json.loads(x).values())

    def process_single_sex(meta):
        tissue_count = meta.shape[0]

        # count totals cells in each replicate in tissue
        replicate_counts  = (
            meta[["replicate"]]
            .assign(tissue_rep_counts=1)
            .groupby(['replicate'])
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
            .groupby(['replicate', 'cluster'])
            .agg({
                'cluster_rep_counts': 'sum',
                'annotations': lambda x: json.dumps({
                    k:v for k,v in x.value_counts().items() if v > 0
                })
            })
            .merge(replicate_counts, left_index=True, right_index=True)
            .fillna(0)
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
    )
    # it may happen that all replicates of a single sex are
    # are not present
    res = res.fillna(0)

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

    # pvals_corrected 2nd item in returned tuple from multipletests
    res['padj_binom'] = multipletests(res['pval_binom'], method='fdr_bh')[1]
    res['padj_fisher'] = multipletests(res['pval_fisher'], method='fdr_bh')[1]
    res['padj_wilcox'] = multipletests(res['pval_wilcox'], method='fdr_bh')[1]
    res['padj_ttest'] = multipletests(res['pval_ttest'].fillna(1.0), method='fdr_bh')[1]

    n_female_reps = len(female_meta.replicate.unique())
    n_male_reps = len(male_meta.replicate.unique())

    padj_use = (
        "padj_wilcox" if ((n_female_reps > 1) and (n_male_reps > 1)) else
        "padj_binomial"
    )

    res = (
        res.assign(count_bias_padj = lambda df: df[padj_use]) #.apply(pvalstr))
        .assign(count_bias_type = lambda df: df.apply(
            lambda x: "Female" if (((x['count_bias_padj'] < PADJ_CUTOFF_COUNT) &
                                    (x["log2_count_bias"] > LFC_CUTOFF_COUNT)) |
                                   (x["cluster_type"] == "female_only")) else
                      "Male"  if (((x['count_bias_padj'] < PADJ_CUTOFF_COUNT) &
                                    (x["log2_count_bias"] < -LFC_CUTOFF_COUNT)) |
                                   (x["cluster_type"] == "male_only")) else
                      "Unbiased",
            axis=1
        ))
    )
    print(res)

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
    adata.obs = adata.obs.assign(cluster = lambda df: (
        df[resol] if resol == "annotaion" else df[resol].apply(
            lambda x: f"{resol}C{x}"
        )
    ))
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

