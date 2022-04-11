library(tidyverse)

yp_fbgns <- c(
  "Yp1" = "FBgn0004045"
  , "Yp2" = "FBgn0005391"
  , "Yp3" = "FBgn0004047"
)

tissues <- c(
  "AC" = "Abdomen carcass",
  "GO" = "Gonad",
  "HD" = "Head",
  "IR" = "Internal reproductive tract",
  "TE" = "Terminalia",
  "TX" = "Thorax",
  "VS" = "Viscera",
  "WO" = "Whole body"
)

df <- (
  read.delim("dmel_yp_genes.bsv", sep = " ", header=FALSE)
)
names(df) <- c("species", "strain", "sex", "ts", "replicate", "Yp1", "Yp3", "Yp2")

print(df)


df <- (
  gather(df, "gene", "expr", -species, -strain, -sex, -ts, -replicate)
  %>% mutate(tissue = tissues[ts])
#  %>% mutate(gene = factor(
#    ref_dmel
#    , levels = yp_fbgns
#    , labels = names(yp_fbgns)
#    , ordered = TRUE
#  ))
)

p <- (
  ggplot(df, aes(x = gene, y = log(1+expr), color = sex))
  + geom_boxplot(outlier.shape=NA)
  + facet_wrap(~tissue, ncol=4)
  #+ coord_flip()
)

pdf("Yp_genes_in_bulk_dmel.pdf", width = 7, height = 5)
p
dev.off()

quit()

