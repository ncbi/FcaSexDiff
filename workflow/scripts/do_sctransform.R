library(anndata)
library(Seurat)
library(Matrix)
library(tidyverse)

infile <- "scraps/lognorm_h5ad/lognorm_body_stringent_female.h5ad"
seuratfile <- "scraps/lognorm_h5ad/seurat_sct_body_stringent_female.rds"
outfile <- "scraps/lognorm_h5ad/sct_lognorm_body_stringent_female.h5ad"

infile <- snakemake@input[[1]]
outfile <- snakemake@output[[1]]

print(infile)
print(outfile)



transform <- function(ad) {

  meta <- ad$obs
  print(head(meta))

  counts <- as(t(as.matrix(ad$X)), "sparseMatrix")
  row_names <- rownames(ad$var)
  rownames(counts) <- rownames(ad$var)
  colnames(counts) <- rownames(ad$obs)
  print(counts[1:5,1:5])

  so <- CreateSeuratObject(
    counts = counts
    , project = 'FCATissueBody'
    , meta.data = meta
  )

  so <- SCTransform(
    so
    , vars.to.regress = c('batch')
    , return.only.var.genes = FALSE
    # change it to zero to match number of genes in RNA & SCT
    # but that gives some NA error
    , min_cells = 5
    , verbose = TRUE
  )

  print(dim(so))
  print(dim(so@assays$RNA))
  dim(so@assays$SCT)

  sct_counts <- Matrix::t(so@assays$SCT@data)
  print(sct_counts[1:5, 1:5])

  # SCTransform has removed some genes
  # so truncate scanpy object
  ad <- ad[,colnames(sct_counts)]

  ad$layers[["umi"]] <- ad$X
  ad$X <- sct_counts

  ad
}

ad <- read_h5ad(infile)
ad

male <- ad[ad$obs$sex == "male",]
female <- ad[ad$obs$sex == "female",]

ad <- concat(list(transform(female), transform(male)), join = "outer")
ad

g2chr <- read.delim(
  "resources/flybase/gene_chr_FB2019_06_r6.31.tsv"
  , header = FALSE
  , row.names = 2
)
names(g2chr) <- c("FBgn", "chr")
g2chr = g2chr[ad$var_names,]

head(g2chr)

umi <- ad$layers[["umi"]]
norm <- ad$X

var_name <- attr(ad$var_names,"name")

ad$var <- (
  cbind(ad$var, g2chr)
  %>% mutate(umi_tissue = colMeans(umi))
  %>% mutate(nz_umi_tissue = colSums(umi)/colSums(umi != 0))
  %>% mutate(norm_tissue = colMeans(norm))
  %>% mutate(nz_norm_tissue = colSums(norm)/colSums(norm != 0))
)

attr(ad$var_names,"name") <- var_name

write_h5ad(ad, outfile)
