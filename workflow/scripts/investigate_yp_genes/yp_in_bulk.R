library(tidyverse)

tissues <- c(
  "AC" = "Abdomen carcass",
  "GO" = "Gonad",
  "HD" = "Head",
  "TX" = "Thorax",
  "VS" = "Viscera"
)

alldf <- data.frame()

for (ts in names(tissues)) {
male_file = paste0("/panfs/pan1/scanalysis/flydata/filtered-normalized-counts/count~norm,replicate~all,sex~M,species~all,strain~rep,tissue~", ts, "/filtered.tsv")

yp_fbgns <- c(
  "Yp1" = "FBgn0004045"
  , "Yp2" = "FBgn0005391"
  , "Yp3" = "FBgn0004047"
)


male <- (
  read.delim(male_file)
  %>% filter(ref_dmel %in% yp_fbgns)
  %>% gather("sample", "expr", -ref_dmel)
  %>% separate(sample, c("species", "strain", "sex", "tissue", "replicate"))
)

head(male)
dim(male)

female_file = paste0("/panfs/pan1/scanalysis/flydata/filtered-normalized-counts/count~norm,replicate~all,sex~F,species~all,strain~rep,tissue~", ts, "/filtered.tsv")

female <- (
  read.delim(female_file)
  %>% filter(ref_dmel %in% yp_fbgns)
  %>% gather("sample", "expr", -ref_dmel)
  %>% separate(sample, c("species", "strain", "sex", "tissue", "replicate"))
)

head(female)
dim(female)

df <- (
  bind_rows(male, female)
  %>% mutate(gene = factor(
    ref_dmel
    , levels = yp_fbgns
    , labels = names(yp_fbgns)
    , ordered = TRUE
  ))
)

alldf <- rbind(alldf, df %>% mutate(tissue = tissues[ts]))
}


p <- (
  ggplot(alldf, aes(x = gene, y = log(1+expr), color = sex))
  + geom_boxplot(outlier.shape=NA)
  #+ facet_wrap(~species, ncol=3)
  + facet_grid(tissue~species)
  #+ coord_flip()
)

pdf("Yp_genes_in_bulk.pdf", width = 10, height = 10)
p
dev.off()

