library(ComplexHeatmap)
library(tidyverse)
library(anndata)

infile <- snakemake@input[[1]]
outfile <- snakemake@output[[1]]
tissue <- snakemake@wildcards[["tissue"]]

print(infile)
print(outfile)
print(tissue)


process <- function(use_tissue, infile, outfile) {
  outbase <- tools::file_path_sans_ext(basename(outfile))
  print(outbase)

  outfile1 <- tempfile(outbase)
  outfile2 <- tempfile(outbase)

  ad = read_h5ad(infile)
  print(ad)

  mat  = as.matrix(ad$X)
  zscore = ad$layers["score"]

  biased_genes = rowSums(abs(mat)) > 0 

  mat  = mat[biased_genes, ]
  zscore = zscore[biased_genes, ]


  pdf(outfile1, height=10, width=20)

  h1 <- Heatmap(
    mat,
    show_row_names = FALSE
  )
  h <- draw(h1)

  mat = mat[row_order(h), column_order(h)]
  zscore = zscore[row_order(h), column_order(h)]

  h1 <- Heatmap(
    mat,
    show_row_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE
  )

  h2 <- Heatmap(
    zscore,
    show_row_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE
  )

  draw(h1 + h2)
  dev.off()

  pdf(outfile2, height=250, width=30)

  h3 <- Heatmap(
    mat
    , cluster_rows = FALSE
    , cluster_columns = FALSE
  )

  h4 <- Heatmap(
    zscore
    , cluster_rows = FALSE
    , cluster_columns = FALSE
  )

  draw(h3 + h4)
  dev.off()

  system(paste("mergepdf", outfile1, outfile2, outfile))
}



process(tissue, infile, outfile)
