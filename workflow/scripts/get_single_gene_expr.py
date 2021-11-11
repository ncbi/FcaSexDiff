#!/usr/bin/env python3
import pandas as pd
import numpy as np
import anndata as ad

exprfile = snakemake.input[0]
sex_specific_annotations_file = snakemake.input[1]
resol = snakemake.wildcards['resol']
gene = snakemake.wildcards['gene']
outfile = snakemake.output[0]

print(exprfile)
print(resol)
print(gene)
print(outfile)

def do_all(exprfile, sex_specific_annotations_file, reosl, outfile):
    adata = ad.read_h5ad(exprfile)

    # first discard cells that have 'mix' sex
    adata = adata[adata.obs.sex.isin(["female", "male"])]

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
    # now do the actual filtering
    adata = adata[~adata.obs.annotation.isin(to_remove)]
    print(adata)

    assert("cluster" not in adata.obs.columns)
    adata.obs = adata.obs.assign(cluster = lambda df: (
        df[resol] if resol == "annotation" else df[resol].apply(
            lambda x: f"{resol}C{x}"
        )
    ))
    print(adata)

    expr = adata[:, [gene, "dsx"]].to_df()

    df = (
        adata.obs[["tissue", "sex", "cluster", "annotation"]]
        .assign(**{gene: expr["AkhR"]})
        .assign(dsx = expr["dsx"])
    )
    print(df)

    df.to_csv(outfile, sep="\t", index=True)


do_all(exprfile, sex_specific_annotations_file, resol, outfile)

