import anndata as ad
import pandas as pd

outfile = snakemake.output[0]
resol = snakemake.wildcards.resol

tissues = [
#    "all",
    "antenna",
    "body",
    "body_wall",
    "fat_body",
    "gut",
    "haltere",
    "head",
    "heart",
    "leg",
    "malpighian_tubule",
    "oenocyte",
    "proboscis_and_maxillary_palps",
    "trachea",
    "wing",
]

def combine_biased_genes(tissues, resolution, outfile):
    bias = pd.concat(
        [
            pd.DataFrame(
                adata.X.todense(),
                index = adata.obs_names,
                columns = adata.var_names
            )
            for tissue in tissues
            for infile in [
                # note this is a single element array
                f"exports/sexdiff/cellfilter~NoSexspecArtef/"
                + f"resolution~{resolution}/{tissue}/"
                + f"sexdiff_{tissue}_stringent_{resolution}_NoSexspecArtef.h5ad"
            ]
            for adata in [ad.read_h5ad(infile)]
        ],
        axis = 1,
        keys = tissues,
        names = ["tissue", "cluster"]
    )

    print(bias)

    bias = bias.fillna(0)
    bias = bias.loc[abs(bias).sum(axis = 1) > 0]

    print(bias)

    bias = pd.concat([bias], keys=["bias"], axis = 1)

    print(bias)

    bias = bias.stack("cluster").stack("tissue")

    print(bias)

    bias = bias.loc[abs(bias.bias) > 0, :]

    print(bias)

    bias.to_csv(outfile, sep = "\t")

if resol == "L6.0":
    tissues = ["body", "head"]

combine_biased_genes(tissues, resol, outfile)

