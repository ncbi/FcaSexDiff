import pandas as pd
import anndata as ad
import numpy as np
import scanpy as sc

adata = ad.read_h5ad(snakemake.input[0])
print(adata)

adata = adata[~adata.obs.scpopcorn_cluster.isnull()]
print(adata)



adata.obs = (
    adata.obs.assign(subcls = lambda df: df.apply(
        lambda x: f'C{int(x["scpopcorn_cluster"])}_{x["sex"]}', axis=1)
    )
    .assign(subcls = lambda df: pd.Categorical(df["subcls"]))
)

print(adata.obs)

print("Copying normalized exp to X ...")
adata.X = np.log1p(adata.layers["norm"])

#sc.tl.dendrogram(adata, "subcls", n_pcs=30)

print("Going to create correlation matrix ...")
ax = sc.pl.correlation_matrix(adata, "subcls", show=False, figsize=(10,8))

ax[0].figure.savefig(snakemake.output[0])
