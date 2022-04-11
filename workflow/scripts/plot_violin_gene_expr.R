library(tidyverse)
library(ggbeeswarm)

resol = snakemake@wildcards[['resol']]
symbol = snakemake@wildcards[['gene']]

infile = snakemake@input[[1]]
outfile = snakemake@output[[1]]

df = (
  read.delim(infile, sep='\t')
  %>% mutate(gene = .[, !!symbol])
)

head(df)

annot = (
  df
  %>% count(cluster, annotation)
  %>% group_by(cluster)
  %>% summarize(major_annotation = annotation[which.max(n)])
)

print(annot)


info = (
  df
  %>% mutate(count=1)
  %>% group_by(cluster, sex)
  %>% summarize(count = sum(count), nz = sum(gene > 0))
  %>% mutate(label = paste0(count, "(", nz, ")"))
  %>% select(cluster, sex, label)
  %>% spread(sex, label)
  %>% left_join(annot)
)
print(info)

if (resol == "annotation") {
  info = info  %>% mutate(label = paste0(cluster, " ", female, ",", male))
} else {
  info = info  %>% mutate(label = paste0(cluster, "|", major_annotation, " ", female, ",", male))
}

df = left_join(df, info)
head(df)

ncluster = length(unique(df$label))

pdf(outfile, height=2*ncluster, width=15)
(
 ggplot(df, aes(x=sex, y=gene, fill=sex, color=sex, group=batch))
 + geom_violin(alpha = 0.25)
 + geom_boxplot(outlier.colour=NA, position=position_dodge(width=0.8), alpha=0.25)
 + geom_point(position=position_jitterdodge(dodge.width=0.9), color="black", alpha=0.25)
 + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
 + labs(x="cluster", y=paste(symbol, "expression (log normalized)"))
 + facet_wrap(~label, ncol=2)
 + coord_flip()
)
