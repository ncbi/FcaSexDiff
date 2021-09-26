#!/usr/bin/env python3

import pandas as pd
import sys, os
sys.path.append(os.path.abspath("../utils"))
from utils import *
from df2excel import *
from io import StringIO

bias_file = snakemake.input[0]
gene_meta_file = snakemake.input[1]
outfile = snakemake.output[0]

print(bias_file)
print(gene_meta_file)
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


# get cluster information
count_bias = (
    pd.read_hdf(bias_file, key='count_bias')
    [['cluster', 'major_annotation', 'cluster_type',
      'cell_count_female', 'cell_count_male',
      #'sample_cell_count_female', 'sample_cell_count_male',
      'log2_male_to_female', 'padj_binom']]
    .rename(columns={
        'log2_male_to_female': 'log2_count_bias',
        'padj_binom': 'count_bias_sig'
    })
    .assign(count_bias_sig = lambda df: df['count_bias_sig'].apply(pvalstr))
)
print(count_bias)

use_stats = 'avg_all avg_male avg_female log2fc padj bias'.split()

expr_bias = pd.read_hdf(bias_file, key='expr_bias')
expr_bias = expr_bias.loc[:, expr_bias.columns.get_level_values('stats').isin(use_stats)]
print(expr_bias)


padj = expr_bias.loc[:, expr_bias.columns.get_level_values('stats') == 'padj']
log2fc = expr_bias.loc[:, expr_bias.columns.get_level_values('stats') == 'log2fc']

padj.columns = padj.columns.droplevel('stats')
log2fc.columns = log2fc.columns.droplevel('stats')

male_bias = (padj < 0.05) & (log2fc > 1)
female_bias = (padj < 0.05) & (log2fc < -1)

bias = pd.DataFrame(index=male_bias.index, columns=male_bias.columns)
bias[male_bias == True] = 'Male'
bias[female_bias == True] = 'Female'
bias = pd.concat([bias], keys=['bias'], names=['stats'], axis=1)
bias.columns = bias.columns.reorder_levels(list(range(1, bias.columns.nlevels))+[0])

print(bias)

expr_bias = (
    pd.concat([expr_bias, bias], axis=1)
    .stack('stats')
    [count_bias['cluster']]
)
print(expr_bias)

expr_bias.columns = pd.MultiIndex.from_frame(count_bias)
print(expr_bias)

expr_bias = (
    expr_bias.unstack('stats')
    .reindex(use_stats, axis=1, level='stats')
)
print(expr_bias)


# compute number of biased clusters for each gene
male_bias_count = male_bias.sum(1)
female_bias_count = female_bias.sum(1)

# augment biased cluster number for each gene on the rows
expr_bias.index = pd.MultiIndex.from_frame(
    expr_bias.index.to_frame()
    .merge(pd.read_table(gene_meta_file, index_col='symbol'),
           how='left', left_index=True, right_index=True)
    .assign(male_cls = male_bias_count)
    .assign(female_cls = female_bias_count)
    .assign(biased_cls = male_bias_count + female_bias_count)
)

# trim non-interesting genes
expr_bias = expr_bias.loc[expr_bias.index.get_level_values('biased_cls')>0,:]

expr_bias.index = expr_bias.index.reorder_levels('chr symbol FBgn avg_base male_cls female_cls biased_cls'.split())
print(expr_bias)

expr_bias = expr_bias.sort_index()
print(expr_bias)

summary = expr_bias.loc[:, expr_bias.columns.get_level_values('stats') == 'bias']
print(summary)


d2x = DF2Excel(outfile)
d2x.write(readme, 'ReadMe', write_header=False)
d2x.write(summary, 'Summary')
d2x.write(expr_bias, 'Detailed')
d2x.close()
