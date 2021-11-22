import numpy as np
import pandas as pd
import loompy
import json
import anndata as ad


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

    # create a new cell_meta column representing short FCA sample id
    assert("sample" not in cell_meta.columns)
    # extract FCA2 from d54b72e4__FCA2_MaleFemale_Body
    cell_meta["sample"] = (
        cell_meta["sample_id"]
        .str.split("__").str[1]
        .str.split("_").str[0]
    )

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
        .set_index('symbol')
    )

    umi = ds.sparse().T.tocsr()

    adata = ad.AnnData(
        umi,
        obs=cell_meta,
        var=gene_meta,
    )
    print(adata)

    adata.write_h5ad(h5ad_path)


loom_to_h5ad(snakemake.input[0], snakemake.output[0])
