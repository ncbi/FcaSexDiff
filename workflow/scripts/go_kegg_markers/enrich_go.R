library(anndata)
library(tidyverse)
library(clusterProfiler)


infile = snakemake@input[[1]]
outfile = snakemake@output[[1]]
tissue = snakemake@wildcards[["tissue"]]
cluster = snakemake@wildcards[["cluster"]]

print(infile)
print(tissue)
print(cluster)

# SET THE DESIRED ORGANISM HERE
organism = "org.Dm.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

symbol2fbgn <- (
  read.delim("resources/flybase/gene_chr_FB2019_06_r6.31.tsv", header=FALSE)
  %>% pull(var = V1, name = V2)
)
print(head(symbol2fbgn))

adata <- read_h5ad(infile)

score <- adata$layers["marker_score"]
padj <- adata$layers["padj"]

fbgns <- symbol2fbgn[row.names(score)]

# we want the log2 fold change
original_gene_list <- score[, cluster]

# name the vector
names(original_gene_list) <- fbgns

# omit any NA values
gene_list <- na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

# Exctract significant results (padj < 0.05)
sig_genes = padj[, cluster] < 0.05

# From significant results, we want to filter on log2fold change
genes <- score[, cluster]

# Name the vector
names(genes) <- fbgns

# omit NA values
genes <- na.omit(genes)

# filter on min log2fold change (log2FoldChange > 1)
genes <- names(genes)[abs(genes) > 1]

gse <- tryCatch(
  expr = {
    # Your code...
    # goes here...
    # ...
    enrichGO(
      gene = genes
      , universe = names(gene_list)
      , OrgDb = organism
      , keyType = 'ENSEMBL'
      , readable = T
      , ont = "ALL"
      , pvalueCutoff = 0.25
      , qvalueCutoff = 0.50
    )
  }
  , error = function(e){
    # (Optional)
    # Do this if an error is caught...
    return(NA)
  }
  # , warning = function(w){
  #   # (Optional)
  #   # Do this if an warning is caught...
  # }
  # , finally = {
  #   # (Optional)
  #   # Do this at the end before quitting the tryCatch structure...
  # }
)


saveRDS(gse, outfile)
