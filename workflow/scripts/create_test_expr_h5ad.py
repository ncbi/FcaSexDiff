import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from numpy.random import normal

def format_test_data(infile, outfile, popcornfile):
  cell_meta = (
      pd.read_csv(infile, index_col=0, nrows=5, header=None, prefix="cell_")
      .T
  )
  cell_meta.index = cell_meta.index.rename("CellID")
  cell_meta.columns = cell_meta.columns.rename("")
  print(cell_meta)

  # generate random tsne coordinates for the cells in clusters
  clusters = cell_meta["L6.0"].unique().tolist()
  centroids = normal(loc=[0,0], scale=[10,10], size=(len(clusters),2))
  print(centroids)

  cell_meta = cell_meta.merge(
      cell_meta["L6.0"].apply(
          lambda x: pd.Series(
              np.ravel(normal(centroids[clusters.index(x)], [1, 1], [1, 2])),
              index=["tSNE1", "tSNE2"]
          )
      ),
      left_index=True, right_index=True
  )
  print(cell_meta)

  df = pd.read_csv(infile, index_col=0, skiprows=5, header=None, prefix="cell_")
  df.index = df.index.rename("symbol")
  df.columns = df.columns.rename("CellID")
  print(df)

  adata = ad.AnnData(
      sparse.csr_matrix(df.values.T),
      obs=cell_meta.drop(columns=["popcorn"]),
      var=pd.DataFrame(index=df.index),
  )
  print(adata)

  adata.write_h5ad(outfile)
  cell_meta[["popcorn"]].to_csv(popcornfile, sep="\t", header=False)

format_test_data(snakemake.input[0], snakemake.output[0], snakemake.output[1])
