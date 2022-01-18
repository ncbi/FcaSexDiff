#!/usr/bin/env python3

import pandas as pd
import anndata as ad
from io import StringIO
from utils import DF2Excel, pvalstr

bias_file = snakemake.input[0]
outfile = snakemake.output["big"]
gene_info_xlsx = snakemake.output["genes"]
cluster_info_xlsx = snakemake.output["clusters"]
chunk_file_prefix = snakemake.params["chunk_prefix"]
tissue = snakemake.wildcards.tissue
resol = snakemake.wildcards.resol

to_hide = (tissue != "test")

print(bias_file)
print(outfile)
print(tissue)
print(resol)
print(to_hide)


ordered_stats = [
    "avg_all", "avg_female", "avg_male",
    "avg_nz_female", "avg_nz_male",
    "frac_all", "frac_female", "frac_male",
    "padj", "log2fc", "bias",
    "log2fc_scanpy", "bias_scanpy",
    #"marker_padj", "marker_log2fc", "is_marker",
    "marker_score",
]

ordered_rows = [
    "chr", "symbol", "FBgn", "umi_tissue", "norm_tissue",
    "nz_umi_tissue", "nz_norm_tissue", "female_cls", "male_cls",
]

ordered_cols = [
    'cluster', 'markers',
#    'annotations_female', 'annotations_male',
    'annotations', 'major_annotation',
#    'tissue_count_female', 'tissue_count_male',
    'cluster_count_female', 'cluster_count_male', 'cluster_type',
#    'cluster_frac_female', 'cluster_frac_male',
#    'tissue_rep_counts_female', 'tissue_rep_counts_male',
    'cluster_rep_counts_female', 'cluster_rep_counts_male',
#    'cluster_rep_fracs_female', 'cluster_rep_fracs_male',
#    'cluster_rep_fracs_mean_female', 'cluster_rep_fracs_mean_male',
#    'cluster_rep_fracs_sd_female', 'cluster_rep_fracs_sd_male',
#    'pval_binom', 'pval_fisher', 'pval_wilcox', 'pval_ttest',
#    'padj_binom', 'padj_fisher', 'padj_wilcox', 'padj_ttest',
    'count_bias_padj', 'log2_count_bias', 'count_bias_type',
    'female_gene', 'male_gene',
    'stats',
]

if to_hide:
    # skip some of the data
    criteria = lambda x: not (("nz_" in x) or ("_scanpy" in x))
    ordered_stats = list(filter(criteria, ordered_stats))
    ordered_rows = list(filter(criteria, ordered_rows))
    ordered_cols = list(filter(criteria, ordered_cols))


def get_stats_frame(adata, layer = None, is_bias = False, is_marker = False):
    tmp = adata.layers[layer] if layer else adata.X
    if is_bias:
        tmp = tmp.todense()
        dat = tmp.astype(str)
        dat[tmp > 0] = "Female"
        dat[tmp < 0] = "Male"
        dat[tmp == 0] = ""
    elif is_marker:
        tmp = tmp.todense()
        dat = tmp.astype(str)
        dat[tmp > 0] = "Marker"
        dat[tmp <= 0] = ""
    else:
        dat = tmp
    return (
        pd.DataFrame(dat, index=adata.obs_names, columns=adata.var_names)
    )


def prepare_bias_table(adata):

    expr_bias = pd.concat(
        [
            get_stats_frame(
                adata,
                layer = stats if stats != "bias" else None,
                is_bias = stats in ["bias", "bias_scanpy"],
                is_marker = stats in ["is_marker"],
            )
            for stats in ordered_stats
        ],
        keys = ordered_stats,
        names = ['stats', 'cluster'],
        axis = 1
    )

    for  col in adata.var.columns:
        print(str(col))
        print(adata.var[col].tolist())
        print("\n")

    expr_bias.columns = pd.MultiIndex.from_frame(
        expr_bias.columns.to_frame(index=False)
        .merge(adata.var.reset_index(), how="left")
        [ordered_cols]
    )

    expr_bias.index = pd.MultiIndex.from_frame(
        expr_bias.index.to_frame(index=False)
        .merge(adata.obs.reset_index(), how="left")
        [ordered_rows]
    )

    return expr_bias


def write_to_excel(detailed, info, tissue_stats, outfile):

    d2x = DF2Excel(outfile)

    readme = d2x.writer.book.add_worksheet("ReadMe")

    cell_format = d2x.writer.book.add_format()
    cell_format.set_align('left')
    cell_format.set_align('vcenter')

    bold_format = d2x.writer.book.add_format()
    bold_format.set_align('left')
    bold_format.set_align('vcenter')
    bold_format.set_bold()

    readme.write(
        0, 0,
        "Cluster level sexbias information both by counts and by gene expression",
        bold_format
    )
    #readme.write(
    #    2, 0,
    #    "Sheet 'Summary' shows  (only for genes that are biased in atleast"
    #    " one cluster)",
    #    cell_format
    #)
    readme.write(
        2, 0,
        "Sheet 'Detailed' shows for each gene G and cluster C if G is sex-biased"
        " (Male/Female) in C and the details of average expression, log fold change,"
        " pvalue etc.",
        cell_format
    )
    readme.write(
        4, 0,
        "The details about rows and columns are as given below.",
        cell_format
    )
    readme.write(6, 0, "Column information:", bold_format)
    row = 7
    for col in ordered_cols:
        readme.write(row, 0, col, cell_format)
        readme.write(row, 1, info.get(col, "unknown"), cell_format)
        row += 1
    for col in ordered_stats:
        readme.write(row, 1, col, cell_format)
        readme.write(row, 2, info.get(col, "unknown"), cell_format)
        row += 1
    row += 1
    readme.write(row, 0, "Row information:", bold_format)
    row += 1
    for col in ordered_rows:
        readme.write(row, 0, col, cell_format)
        readme.write(row, 1, info.get(col, "unknown"), cell_format)
        row += 1

    #d2x.write(summary, 'Summary', write_index=True)

    print(tissue_stats)

    tissue_stats.loc["comment", "tissue_rep_counts_female"] = "replicates"
    tissue_stats.loc["comment", "tissue_rep_counts_male"] = "replicates"
    tissue_stats.T.to_excel(d2x.writer, sheet_name='Detailed')
    d2x.write(detailed, 'Detailed')

    d2x.close()


def export_excel(expr_bias, info, tissue_stats, outfile):

    detailed = (
        expr_bias.sort_index(0)
        .sort_index(1)
        .reindex(ordered_stats, axis=1, level='stats')
    )

    #if to_hide:
    #    # trim non-interesting genes
    #    detailed = detailed.loc[
    #        detailed.index.get_level_values('male_cls') +
    #        detailed.index.get_level_values('female_cls') > 0,
    #        :
    #    ]

    print(detailed)

    #summary = detailed.loc[
    #    :,
    #    detailed.columns.get_level_values('stats') == 'bias'
    #]
    #summary = pd.DataFrame(
    #    summary.values,
    #    index = summary.index.get_level_values("symbol"),
    #    columns = summary.columns.get_level_values("cluster"),
    #)
    #print(summary)

    #d2x.write(summary, 'Summary', write_index=True)

    write_to_excel(detailed, info, tissue_stats, outfile)


    cluster_info = detailed.columns.to_frame(index=False)
    gene_info = detailed.index.to_frame(index=False)

    old_columns = cluster_info.columns.tolist()

    ncluster_in_chunk = 20
    clusters = cluster_info.cluster.unique().tolist()
    chunks = [
        clusters[i:i+ncluster_in_chunk]
        for i in range(0, len(clusters), ncluster_in_chunk)
    ]
    chunks = pd.DataFrame(
        [[cl, ch+1] for ch, cls in enumerate(chunks) for cl in cls],
        columns = ['cluster', 'chunk']
    )
    cluster_info = cluster_info.merge(chunks, how="left")
    print(cluster_info)

    cluster_info.set_index('cluster').to_excel(cluster_info_xlsx)
    gene_info.set_index('symbol').to_excel(gene_info_xlsx)

    detailed.columns = pd.MultiIndex.from_frame(
        cluster_info[["chunk"] + old_columns]
    )
    # there is some bug in snakemake, following line gives error
    # for ch, grp in detailed.groupby("cluster"):
    for ch in cluster_info.chunk.unique():
        grp = detailed.loc[
            :,
            detailed.columns.get_level_values("chunk") == ch
        ]
        chunk_file = f"{chunk_file_prefix}{ch}.xlsx"
        write_to_excel(grp, info, tissue_stats, chunk_file)



# main function

adata = ad.read_h5ad(bias_file)
bias = prepare_bias_table(adata)
export_excel(bias, adata.uns["info"], adata.uns["stats"], outfile)

