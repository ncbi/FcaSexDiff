import anndata as ad
import pandas as pd
from os.path import basename

paths = snakemake.input
tissues = ['_'.join(basename(x).split("_")[1:]) for x in snakemake.input]

print(paths)

summary = pd.concat(
    [
        (
            ad.read_h5ad(path)
            .var[["female_gene", "male_gene"]]
            .rename(columns=lambda x: x.split("_")[0])
        )
        for path in paths
    ],
    axis=1, keys=tissues, names=["tissue", "gene_bias"]
)
print(summary)


summary.to_excel("exports/extras/summary_sexbiased_genes_annotation.xlsx")

