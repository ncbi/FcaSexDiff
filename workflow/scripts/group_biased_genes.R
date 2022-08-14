library(tidyverse)

process <- function(infile, outfile) {
bias <- read.csv(infile, row.names=1)

ncluster <- ncol(bias)

group_gene <- function(x) {
  x <- (
      stack(x)
      %>% separate(
        ind
        , c("annotation", "female_cells", "male_cells", "label")
        , remove = FALSE
        , sep = "[_]"
      )
      %>% separate(label, c("resol", "cluster"), remove = FALSE, sep = "C")
  )
  males <- x %>% filter(values == -1)
  females <- x %>% filter(values == 1)
  nmale <- nrow(males)
  nfemale <- nrow(females)
  if (min(nmale, nfemale) > 0) {
    group <- "MixedBiased"
  } else if (nmale > nfemale) {
    group <- "MaleOnlyBiased"
  } else if (nmale < nfemale) {
    group <- "FemaleOnlyBiased"
  } else {
    group <- "Undecided"
  }
  data.frame(
    n_female_bias = nfemale
    , n_male_bias = nmale
    , bias_group = group
    , n_bias = nmale + nfemale
    , short_female_clusters = paste(females$cluster, collapse=",")
    , short_male_clusters = paste(males$cluster, collapse=",")
    , long_female_clusters = paste(females$ind, collapse=",")
    , long_male_clusters = paste(males$ind, collapse=",")
  )
}

res <- (
  bias
  %>% rownames_to_column("symbol")
  %>% group_by(symbol)
  %>% summarise(info = list(group_gene(cur_data())))
  %>% unnest(info)
  %>% arrange(bias_group, desc(n_bias))
)

print(res)

write.table(res, outfile, sep="\t", row.names = FALSE)
}

process(snakemake@input[[1]], snakemake@output[[1]])
