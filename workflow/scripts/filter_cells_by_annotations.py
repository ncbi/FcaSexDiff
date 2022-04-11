#!/usr/bin/env python3
import re
import json
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad


exprfile = snakemake.input[0]
sex_specific_annotations_file = snakemake.input[1]
cell_filter = snakemake.wildcards['cellfilt']
outfile = snakemake.output[0]

print(exprfile)
print(outfile)


def do_all(exprfile, sex_specific_annotations_file, outfile):
    adata = ad.read_h5ad(exprfile)

    # first discard cells that have 'mix' sex
    adata = adata[adata.obs.sex.isin(["female", "male"])]

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

    print(adata)
    pd.Series(adata.obs_names).to_csv(outfile, index=False)


do_all(exprfile, sex_specific_annotations_file, outfile)

