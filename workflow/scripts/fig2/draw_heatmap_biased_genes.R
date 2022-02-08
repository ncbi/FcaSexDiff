library(ComplexHeatmap)
library(tidyverse)

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
  df <- (
    read.delim(infile)
    %>% filter(tissue == use_tissue)
  )

  head(df)
  unique(df$cluster)

  mat  <- (
    select(df, symbol, cluster, bias)
    %>% spread("cluster", "bias", fill = 0)
    %>% column_to_rownames("symbol")
  )


  pdf(outfile1, height=250, width=15)
  h1 <- Heatmap(mat)
  draw(h1)
  dev.off()

  pdf(outfile2, height=10, width=10)
  h2 <- Heatmap(
    mat,
    show_row_names = FALSE
  )
  draw(h2)
  dev.off()

  system(paste("mergepdf", outfile2, outfile1, outfile))
}

process(tissue, infile, outfile)
