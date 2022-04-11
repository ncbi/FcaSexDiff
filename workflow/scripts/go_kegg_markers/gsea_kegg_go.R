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

log2fc <- adata$layers["log2fc"]
fbgns <- symbol2fbgn[row.names(log2fc)]

# we want the log2 fold change
original_gene_list <- log2fc[, cluster]

# name the vector
names(original_gene_list) <- fbgns

# omit any NA values
gene_list <- na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

gse <- tryCatch(
  expr = {
    # Your code...
    # goes here...
    # ...
    gseGO(
      geneList = gene_list
      , ont ="ALL"
      , keyType = "ENSEMBL"
      , nPermSimple = 10000
      #, minGSSize = 3
      #, maxGSSize = 800
      , pvalueCutoff = 0.5
      , verbose = TRUE
      , OrgDb = organism
      , pAdjustMethod = "none"
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
