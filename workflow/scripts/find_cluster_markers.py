#!/usr/bin/env python3
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad

exprfile = snakemake.input[0]
resol = snakemake.wildcards['resol']
markerfile = snakemake.output[0]

print(exprfile)
print(resol)
print(markerfile)

def find_markers(adata, groupby):

    # use normalized expression all through
    error = False
    marker_info = pd.DataFrame()
    top_markers = {}

    # workaround for scanpy bug:
    # https://bytemeta.vip/repo/scverse/scanpy/issues/2239
    if "log1p" in adata.uns.keys():
        adata.uns["log1p"]["base"] = None

    try:
        sc.tl.rank_genes_groups(
            adata, groupby, groups = "all", reference = "rest",
            pts = True, method = 'wilcoxon', key_added = "marker_wilcoxon"
        )
    except AttributeError as e:
        print('\tSample doesn\'t contain one group')
        print(e)
        error = True
    except IndexError:
        print('\tNot enough cells in one group')
        error = True
    except ZeroDivisionError as e:
        print('\tZero division!')
        print(e)
        error = True
    except ValueError as e:
        print('\tSample doesn\'t contain one group')
        print(e)
        error = True
    if not error:
        frames = []
        clusters = adata.obs[groupby].unique()
        for cls in clusters:
            tmp_frame = pd.DataFrame({
                # t is a rec.array which can be indexed by col name
                col : [t[cls] for t in adata.uns['marker_wilcoxon'][rec]]
                for col,rec in [
                    ('symbol','names'),
                    ('marker_log2fc', 'logfoldchanges'),
                    ('marker_pval', 'pvals'),
                    ('marker_padj', 'pvals_adj'),
                    ('marker_score', 'scores'),
                ]
            }).set_index('symbol')
            # scanpy already arrange the genes as per descending marker scores
            top_markers[cls] = ','.join(tmp_frame.head(10).index)
            frames.append(tmp_frame)

        marker_info = pd.concat(frames, axis=1, keys = clusters, names =
                            ['cluster', 'stats'])
        marker_info.columns = marker_info.columns.swaplevel()

        print(marker_info)
        # rows in pts are in different order of symbols, handle separately
        pts = pd.DataFrame(adata.uns['marker_wilcoxon']['pts'])
        pts = pd.concat([pts], axis=1, keys=["frac_all"], names=["stats", "cluster"])

        marker_info = (
            pd.concat([marker_info, pts], axis=1)
            .sort_index(1)
        )
        print(marker_info)

    return top_markers, marker_info



def do_all(exprfile, resol, outfile):
    adata = ad.read_h5ad(exprfile)

    # ideally we should call find_markers after filtering the cells.
    # however scanpy gives ZeroDivisionError, calling here does not give
    top_markers, marker_info = find_markers(adata, resol)
    print(top_markers)
    print(marker_info)

    print(top_markers.values())
    print(top_markers.keys())
    top_markers = pd.DataFrame(dict(
        markers = list(top_markers.values()),
        cluster = list(top_markers.keys())
    )).set_index("cluster")
    print(top_markers)
        
    marker_info.to_hdf(outfile, key = "info")
    top_markers.to_hdf(outfile, key = "top")


do_all(exprfile, resol, markerfile)

