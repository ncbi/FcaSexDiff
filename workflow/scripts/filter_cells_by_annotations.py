#!/usr/bin/env python3
import re
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


exprfile = snakemake.input[0]
sex_specific_annotations_file = snakemake.input[1]
cell_filter = snakemake.wildcards['cellfilt']
outfile = snakemake.output[0]

print(exprfile)
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


def do_all(exprfile, sex_specific_annotations_file, outfile):
    adata = ad.read_h5ad(exprfile)

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

    adata.uns["stats"] = tissue_stats

    print(adata)

    adata.write_h5ad(outfile)

do_all(exprfile, sex_specific_annotations_file, outfile)

