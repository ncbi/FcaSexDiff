import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

def format_test_data(infile, outfile):
  cell_meta = (
      pd.read_csv(infile, index_col=0, nrows=2, header=None, prefix="cell_")
      .T
  )
  cell_meta.index = cell_meta.index.rename("CellID")
  cell_meta.columns = cell_meta.columns.rename("")
  print(cell_meta)

  df = pd.read_csv(infile, index_col=0, skiprows=2, header=None, prefix="cell_")
  df.index = df.index.rename("symbol")
  df.columns = df.columns.rename("CellID")
  print(df)

  adata = ad.AnnData(
      sparse.csr_matrix(df.values.T),
      obs=cell_meta,
      var=pd.DataFrame(index=df.index),
  )
  print(adata)

  adata.write_h5ad(outfile)

format_test_data(snakemake.input[0], snakemake.output[0])
