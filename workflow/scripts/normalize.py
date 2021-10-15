import numpy as np
import pandas as pd
import scanpy as sc
import loompy
import json
import anndata as ad

chr_file = 'resources/flybase/gene_chr_FB2019_06_r6.31.tsv'

# force pandas to not convert gene 'nan' by setting na_filter to False
g2chr = (
  pd.read_table(chr_file, header=None, names=['FBgn', 'symbol', 'chr'],
                dtype='str', na_filter=False)
  .set_index("symbol")
)
print(g2chr)


def log_normalize(infile, outfile):

  adata = ad.read_h5ad(infile)
  print(adata)

  from collections import OrderedDict
  adata.layers = OrderedDict()
  adata.layers["umi"] = adata.X.copy()
  umi = adata.layers["umi"]

  sc.pp.normalize_total(adata, target_sum=1e4)
  norm = adata.X.copy()

  sc.pp.log1p(adata)

  adata.var = (
    adata.var
    # augment FBgn id and chromosome with gene symbol in the rows
    .merge(g2chr, how='left', left_index=True, right_index=True)
    # keep average umi of all cells in the tissue
    .assign(umi_tissue = umi.mean(axis=0).T)
    .assign(nz_umi_tissue = (umi.sum(axis=0)/(umi!=0).sum(0)).T)
    # keep average normalized expression of all cells in the tissue
    .assign(norm_tissue = norm.mean(axis=0).T)
    .assign(nz_norm_tissue = (norm.sum(axis=0)/(norm!=0).sum(0)).T)
  )

  adata.write_h5ad(outfile)


log_normalize(snakemake.input[0], snakemake.output[0])
