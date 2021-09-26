import numpy as np
import pandas as pd
import loompy
import json
import anndata as ad


chr_file = 'resources/flybase/gene_chr_FB2019_06_r6.31.tsv'

# force pandas to not convert gene 'nan' by setting na_filter to False
g2chr = pd.read_table(chr_file, header=None, names=['FBgn', 'symbol', 'chr'],
                      dtype='str', na_filter=False)
print(g2chr)


def loom_to_h5ad(loom_path, h5ad_path):
  with loompy.connect(loom_path, validate=False) as ds:

    # get resolution level names for 'Clusterings'
    a = json.loads(ds.attrs['MetaData'])
    resol = (
      pd.DataFrame(a['clusterings'])
      [['id', 'name']]
      .assign(label = lambda df: df.name.replace(
        regex = {
            r'Leiden resolution (\d+\.\d+).*' : r'L\1',
        }
      ))
    )
    print(resol)

    # get names for 'Embeddings'
    embed = (
        pd.DataFrame(a['embeddings'])
        [['id', 'name']]
        .assign(label = lambda df: df.name.replace(
            regex = {
                r'HVG t-SNE' : r'tSNE',
                r'HVG UMAP' : r'UMAP',
                r'HVG PC1/PC2' : r'PC',
                r'SCENIC AUC UMAP' : r'SCENIC_AUC_UMAP',
                r'SCENIC AUC t-SNE' : r'SCENIC_AUC_tSNE',
            }
        ))
    )
    print(embed)


    cell_meta = pd.DataFrame.from_dict({x:ds.ca[x] for x in ds.ca.keys()})

    # breakup 'Clusterings' column to individual resolution levels
    cell_meta = (
        cell_meta
        .merge(cell_meta['Clusterings'].apply(
            lambda x: pd.Series(x, index=resol['label'])
        ), left_index=True, right_index=True)
        .drop(columns = ['Clusterings'])
    )

    # breakup 'Embeddings_X, _Y' columns to individual embeddings
    cell_meta = (
        cell_meta
        .merge(cell_meta['Embeddings_X'].apply(
            lambda x: pd.Series(x, index=embed['label']+'1')
        ), left_index=True, right_index=True)
        .merge(cell_meta['Embeddings_Y'].apply(
            lambda x: pd.Series(x, index=embed['label']+'2')
        ), left_index=True, right_index=True)
        .drop(columns=['Embeddings_X', 'Embeddings_Y'])
    )

    cell_meta = cell_meta.set_index('CellID')

    print(cell_meta)

    gene_meta = (
        pd.DataFrame({'symbol': ds.ra['Gene']})
        # augment FBgn id and chromosome with gene symbol in the rows
        .merge(g2chr, how='left')
        .set_index('symbol')
    )

    umi = ds.sparse().T.tocsr()

    # keep in gene_meta the average umi of all cells in the tissue
    gene_meta['umi_tissue'] = (umi.sum(axis=0)/(umi!=0).sum(0)).T
    print(gene_meta)

    # now normalize gene expression data
    from sklearn.preprocessing import normalize
    norm = normalize(umi, norm='l1', axis=1)*1e6

    # keep average normalized expression of all cells in the tissue
    gene_meta['norm_tissue'] = (norm.sum(axis=0)/(norm!=0).sum(0)).T

    print(gene_meta)


    from collections import OrderedDict
    layers = OrderedDict()
    layers["norm"] = norm

    adata = ad.AnnData(
        umi,
        obs=cell_meta,
        var=gene_meta,
        layers=layers,
    )
    print(adata)

    adata.write_h5ad(h5ad_path)



#loom_to_h5ad(
#    loom_path = 'resources/loom/male_reproductive_glands_relaxed.loom',
#    h5ad_path = 'male_reproductive_glands_relaxed.h5ad',
#)

loom_to_h5ad(snakemake.input[0], snakemake.output[0])
