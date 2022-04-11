import numpy as np
import pandas as pd
import loompy
import json
import anndata as ad


def loom_to_h5ad(loom_path, h5ad_path):

  with loompy.connect(loom_path, validate=False) as ds:

    # get resolution level names for 'Clusterings'
    a = json.loads(ds.attrs['MetaData'])
    print(a["clusterings"])

    for clustering in a['clusterings']:
        print(f'Clustering: {clustering["name"]} has ID: {clustering["id"]}')

    clustering_id_to_use = 1

    print(ds.ra[f'ClusterMarkers_{clustering_id_to_use}'].dtype.names)
    print(pd.DataFrame(ds.ra[f'ClusterMarkers_{clustering_id_to_use}']))
    print(pd.DataFrame(ds.ra[f'ClusterMarkers_{clustering_id_to_use}_avg_logFC']))

    resol = (
      pd.DataFrame(a['clusterings'])
      [['id', 'name']]
      .assign(label = lambda df: df.name.replace(
        regex = {
            r'Leiden resolution (\d+\.\d+).*' : r'L\1',
        }
      ))
      .assign(id = lambda df: df["id"].map(str))
      .set_index("id")
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
        .assign(id = lambda df: df.id.map(str))
        .set_index("id")
    )
    print(embed)

    exclude_list = [
        "Clusterings", "MotifRegulonsAUC", "TrackRegulonsAUC",
        "Embedding", "Embeddings_X", "Embeddings_Y",
    ]

    # make column meta data from loom column attributes
    # except for Clusterings which will be handled separately
    cell_meta = pd.DataFrame.from_dict({
        x:ds.ca[x]
        for x in ds.ca.keys()
        if x not in exclude_list
    })
    print(cell_meta)
    print(cell_meta.columns)
    print(exclude_list)

    # breakup 'Embeddings_X, _Y' columns to individual embeddings

    embed_x = pd.DataFrame(ds.ca["Embeddings_X"])
    embed_x.columns = (
        embed.loc[embed_x.columns.tolist()]
        ["label"].map(lambda x: f"{x}1")
    )
    print(embed_x)

    embed_y = pd.DataFrame(ds.ca["Embeddings_Y"])
    embed_y.columns = (
        embed.loc[embed_y.columns.tolist()]
        ["label"].map(lambda x: f"{x}2")
    )
    print(embed_y)

    # breakup 'Clusterings' column to individual resolution levels
    clusterings = pd.DataFrame(ds.ca["Clusterings"])
    print(clusterings)

    resol = resol.loc[clusterings.columns.tolist()]
    clusterings.columns = resol["label"].tolist()
    print(clusterings)
    print(clusterings.iloc[:,0:9])

    cell_meta = pd.concat([cell_meta, embed_x, embed_y, clusterings], axis=1)
    print(cell_meta)
    print(cell_meta.columns)

    # change integer cluster number to formatted strings
    for res_lab in resol['label']:
        if res_lab[0] == "L":
            # do it for leiden levels only 
            ndigit = np.ceil(np.log10(cell_meta[res_lab].max()+1)).astype(int)
            cell_meta[res_lab] = cell_meta[res_lab].apply(
                lambda x: f"{res_lab}C{x:0>{ndigit}d}"
            )
    print(cell_meta)


    # create a new cell_meta column representing short FCA sample id
    assert("sample" not in cell_meta.columns)
    # extract FCA2 from d54b72e4__FCA2_MaleFemale_Body
    cell_meta["sample"] = (
        cell_meta["sample_id"]
        .str.split("__").str[1]
        .str.split("_").str[0]
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
