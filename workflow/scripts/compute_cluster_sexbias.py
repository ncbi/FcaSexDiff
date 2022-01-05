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
from numpyencoder import NumpyEncoder # Needed for numpy -> json conversion

PSEUDO_EXPR = 1e-9
PSEUDO_COUNT = 1e-256

PADJ_CUTOFF_COUNT = 0.05
LFC_CUTOFF_COUNT = 1

PADJ_CUTOFF_MARKER = 0.001
LFC_CUTOFF_MARKER = 1

PADJ_CUTOFF_EXPR = 0.001
LFC_CUTOFF_EXPR = 1

exprfile = snakemake.input[0]
sex_specific_annotations_file = snakemake.input[1]
resol = snakemake.wildcards['resol']
cell_filter = snakemake.wildcards['cellfilt']
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
  frac_all = (
      "Fraction of cells where gene is expressed (atleast one UMI)"
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
  marker_log2fc = (
      "Log2 of fold change of average expression in favor of cells"
      " in the cluster against rest of the cells"
  ),
  marker_padj = (
      "BH corrected pvalue from Wilcoxon test for cluster cells vs rest cells"
      " expression as computed by scanpy.tl.rank_genes_groups"
  ),
  marker_score = (
      "Marker score computed by scanpy.tl.rank_genes_groups on clusters"
  ),
  is_marker = (
      f"Gene is cluster marker based on marker_padj < {PADJ_CUTOFF_MARKER}"
      f" and marker_log2fc > {LFC_CUTOFF_MARKER}"
      f" ({2**LFC_CUTOFF_MARKER}-fold change)"
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


def find_markers(adata, groupby):

    # use normalized expression all through
    error = False
    marker_info = pd.DataFrame()
    top_markers = {}

    try:
        sc.tl.rank_genes_groups(
            adata, groupby, groups = "all", reference = "rest",
            pts = True, method = 'wilcoxon', key_added = "marker_wilcoxon"
        )
    except AttributeError as e:
        print('\tSample doesn\'t contain one group')
        print(e)
        error = True
    except IndexError:
        print('\tNot enough cells in one group')
        error = True
    except ZeroDivisionError as e:
        print('\tZero division!')
        print(e)
        error = True
    except ValueError as e:
        print('\tSample doesn\'t contain one group')
        print(e)
        error = True
    if not error:
        frames = []
        clusters = adata.obs[groupby].unique()
        for cls in clusters:
            tmp_frame = pd.DataFrame({
                # t is a rec.array which can be indexed by col name
                col : [t[cls] for t in adata.uns['marker_wilcoxon'][rec]]
                for col,rec in [
                    ('symbol','names'),
                    ('marker_log2fc', 'logfoldchanges'),
                    ('marker_pval', 'pvals'),
                    ('marker_padj', 'pvals_adj'),
                    ('marker_score', 'scores'),
                ]
            }).set_index('symbol')
            # scanpy already arrange the genes as per descending marker scores
            top_markers[cls] = ','.join(tmp_frame.head(10).index)
            frames.append(tmp_frame)

        marker_info = pd.concat(frames, axis=1, keys = clusters, names =
                            ['cluster', 'stats'])
        marker_info.columns = marker_info.columns.swaplevel()

        print(marker_info)
        # rows in pts are in different order of symbols, handle separately
        pts = pd.DataFrame(adata.uns['marker_wilcoxon']['pts'])
        pts = pd.concat([pts], axis=1, keys=["frac_all"], names=["stats", "cluster"])

        marker_info = (
            pd.concat([marker_info, pts], axis=1)
            .sort_index(1)
        )
        print(marker_info)

    return top_markers, marker_info

def get_tissue_stats(adata, name):
    return pd.DataFrame(dict(
        n_gene               = [adata.shape[1]],
        n_cluster            = [len(adata.obs[resol].unique())],
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


def do_all(exprfile, sex_specific_annotations_file, resol, outfile):
    adata = ad.read_h5ad(exprfile)

    # ideally we should call find_markers after filtering the cells.
    # however scanpy gives ZeroDivisionError, calling here does not give
    top_markers, marker_info = find_markers(adata, resol)
    print(top_markers)
    print(marker_info)

    # first discard cells that have 'mix' sex
    adata = adata[adata.obs.sex.isin(["female", "male"])]

    tissue_stats = get_tissue_stats(adata, "AllCells")

    if cell_filter in ["NoSexspecArtef", "NoSexspecArtefMuscle"]:
        # next discard cells that have sex-specific annotations
        # we get the information from sex_specific_annotations_file
        # that has a column sex_specific. Ensure it does not conflict 
        assert("sex_specific" not in adata.obs.columns)
        annots = (
            adata.obs[["annotation"]]
            .merge(pd.read_csv(sex_specific_annotations_file), how = "left")
            .drop_duplicates()
        )
        # we also make sure there is no annotation that is not
        # known in sex_specific_annotations_file
        print("\n".join(annots[annots.sex_specific.isna()]["annotation"].unique().tolist()))
        assert(sum(annots.sex_specific.isna()) == 0)
        to_remove = annots.annotation[annots.sex_specific != "no"].tolist()
        # remove cells annotated as artefacts too
        to_remove += ["artefact"]
        if cell_filter == "NoSexspecArtefMuscle":
            to_remove += ["muscle cell"]
        # now do the actual filtering
        adata = adata[~adata.obs.annotation.isin(to_remove)]
        print(adata)
        print(adata.obs.annotation.unique().tolist())

    if cell_filter != "AllCells":
        tissue_stats_filtered = get_tissue_stats(adata, cell_filter)
        tissue_stats = pd.concat([tissue_stats, tissue_stats_filtered])
    print(tissue_stats)


    assert("cluster" not in adata.obs.columns)
    meta_data = (
        adata.obs
        .assign(cluster = lambda df: df[resol])
        # anndata saves samples as categorical variable but our code
        # assumes that it is of string type, not categorical
        .assign(sample = lambda df: df["sample"].astype(str))
    )
    print(adata)


    count_bias = compute_sexbiased_counts_male_to_female(meta_data)
    print(count_bias)

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



    # make sure that the columns in marker_info are in cluster order
    marker_info = marker_info.reindex(clusters, axis=1, level='cluster')

    is_marker = (
        (marker_info['marker_padj'] < PADJ_CUTOFF_MARKER) &
        (marker_info['marker_log2fc'] > LFC_CUTOFF_MARKER)
    )


    gene_meta = (
        adata.var.loc[gene_bias.index, :]
        .assign(female_cls = (gene_bias > 0).sum(axis=1))
        .assign(male_cls = (gene_bias < 0).sum(axis=1))
    )

    cluster_meta = (
        count_bias.loc[gene_bias.columns, :]
        .assign(top_markers = pd.Series(top_markers))
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

do_all(exprfile, sex_specific_annotations_file, resol, outfile)

