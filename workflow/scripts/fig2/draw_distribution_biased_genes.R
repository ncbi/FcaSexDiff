library(tidyverse)

bias = read.delim("all_biased_gene.tsv", sep = "\t")

head(bias)

female_bias = bias[bias$bias > 0,]
head(female_bias)

df = female_bias %>% group_by(cluster) %>% summarize(ngene = sum(bias))

print(df)

ggplot(df, aes(x=ngene)) +
  geom_histogram() +
  labs(x="Number of female biased genes (n)", y="Number of clusters with n female biased genes")


male_bias = bias[bias$bias < 0,]
head(female_bias)

df = male_bias %>% group_by(cluster) %>% summarize(ngene = -sum(bias))

print(df)

ggplot(df, aes(x=ngene)) +
  geom_histogram() +
  labs(x="Number of male biased genes (n)", y="Number of clusters with n male biased genes")




