#!/usr/bin/env python3

import pandas as pd
import anndata as ad
from io import StringIO
from utils import DF2Excel, pvalstr

bias_file = snakemake.input[0]
outfile = snakemake.output[0]

print(bias_file)
print(outfile)

readme_text = """
Cluster level sexbias information both by counts and by gene expression

Sheet 'Summary' shows for each gene G and cluster C if G is sex-biased (Male/Female) in C
Sheet 'Detailed' shows details of average expression, log fold change, pvalue in addition
In both sheets the details about rows and columns are as given below.

Column information:
cluster|Cluster name (cell type if annotation is used)
major_annotation|Major cell type in the cluster, useful if clusters are not by annotation
cluster_type|Whether cluster contains cells of a single sex or both
cell_count_male|Number of male cells in the cluster
cell_count_female|Number of female cells in the cluster
log2_count_bias|log2((NM + eps)/(FM + eps)) where NM (NF) is the normalized (by total cells in the tissue sample) male (female) cell count for the cluster
count_bias_sig|Significance of the count bias based on Fisher's test
stats|Statistics for a gene within cluster
|avg_all|Average expression of gene in all cells in cluster
|avg_male|Average expression of gene in all male cells in cluster
|avg_female|Average expression of gene in all female cells in cluster
|log2fc|Log2 of fold change of average expression in favor of male cells against female cells
|padj|BH corrected pvalue from Wilcoxon test for male cells vs female cells expression
|bias|gene is male- or female-biased based on padj<0.05 and log2fc >1 or <-1 (2-fold change)


Row information:
chr|Chromosome where gene is located
symbol|Gene symbol
FBgn|Flybase ID of gene
male_cls|Number of clusters where gene is male-biased
female_cls|Number of clusters where gene is female-biased
biased_cls|Number of clusters where gene is either male- or female-biased
"""

readme = pd.read_csv(StringIO(readme_text), names=['short', 'long', 'more'],
                     sep='|', engine='python', skip_blank_lines=False)

print(readme)

adata = ad.read_h5ad(bias_file)

ordered_stats = 'avg_all avg_female avg_male log2fc padj'.split()

tmp = adata.X.todense()
bias = tmp.astype(str)
bias[tmp > 0] = "Male"
bias[tmp < 0] = "Female"
bias[tmp == 0] = ""
bias = (
    pd.DataFrame(bias, index=adata.obs_names, columns=adata.var_names)

)

frames = [
    pd.DataFrame(adata.layers[stat], index=adata.obs_names, columns=adata.var_names)
    for stat in ordered_stats
]

expr_bias = pd.concat(
    [bias] + frames,
    keys = ['bias'] + ordered_stats,
    names = ['stats', 'cluster'],
    axis = 1
)
print(expr_bias)


ordered_rows = [
    "chr", "symbol", "FBgn", "umi_tissue", "norm_tissue", "female_cls",
    "male_cls",
]

ordered_cols = [
    'cluster', 'major_annotation', 'cluster_type',
    'cell_count_female', 'cell_count_male',
    'log2_count_bias', 'count_bias_sig',
    'female_gene', 'male_gene', 'stats',
]

expr_bias.columns = pd.MultiIndex.from_frame(
    expr_bias.columns.to_frame(index=False)
    .merge(adata.var.reset_index(), how="left")
    .assign(count_bias_sig = lambda df: df['padj_binom'].apply(pvalstr))
    [ordered_cols]
)
print(expr_bias)

expr_bias.index = pd.MultiIndex.from_frame(
    expr_bias.index.to_frame(index=False)
    .merge(adata.obs.reset_index(), how="left")
    [ordered_rows]
)
print(expr_bias)

detailed = (
    expr_bias.sort_index(0)
    .sort_index(1)
    .reindex(ordered_stats + ['bias'], axis=1, level='stats')
)
print(detailed)

# trim non-interesting genes
if (detailed.shape[0] > 100):
    detailed = detailed.loc[
        detailed.index.get_level_values('male_cls') +
        detailed.index.get_level_values('female_cls') > 0,
        :
    ]
print(detailed)

summary = detailed.loc[:, detailed.columns.get_level_values('stats') == 'bias']
print(summary)


d2x = DF2Excel(outfile)
d2x.write(readme, 'ReadMe', write_header=False)
d2x.write(summary, 'Summary')
d2x.write(detailed, 'Detailed')
d2x.close()

