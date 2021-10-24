#!/usr/bin/env Rscript
library(tidyverse)
library(ggthemes)
library(ggrepel)
library(patchwork)
library(ggrastr)
library(anndata)

print(getwd())
source("workflow/scripts/utils.R")

RED = "#e8402b"
GREEN = "#3c811b"
BLUE = "#1f45f7"

#bias_file = "exports/sexdiff_h5ad/resol~annotation/sexdiff_body_stringent_annotation.h5ad"
#expr_file = "scraps/lognorm_h5ad/lognorm_body_stringent.h5ad"


expr_file = snakemake@input[[1]]
bias_file = snakemake@input[[2]]
out_file = snakemake@output[[1]]

resol = snakemake@wildcards[['resol']]

print(expr_file)
print(bias_file)
print(out_file)
print(resol)
 

bias <- (
  get_var_anndata(read_h5ad(bias_file))
  %>% select(cluster, major_annotation, cluster_rep_fracs_mean_female,
             cluster_rep_fracs_mean_male, count_bias_padj,
             log2_count_bias, count_bias_type)
  %>% mutate(pct_female = 100*cluster_rep_fracs_mean_female)
  %>% mutate(pct_male = 100*cluster_rep_fracs_mean_male)
  %>% mutate(log2_count_bias = ifelse(log2_count_bias > 25, 25,
                                      ifelse(log2_count_bias < -25, -25,
                                             log2_count_bias)))
)

head(bias)

df <- (
  read_h5ad(expr_file)$obs
  %>% select(tSNE1, tSNE2, !!resol, sex)
  %>% filter(sex %in% c("female", "male"))
  %>% mutate(cluster = .[,resol]) # rename column named by resol as cluster
  %>% mutate(cluster = paste0(!!resol, "C", cluster))
  %>% left_join(bias)
)

head(df)


ff = df %>% filter(sex=='female')
mf = df %>% filter(sex=='male')


clusters <-
(
    ggplot(df, aes(tSNE1, tSNE2, color=cluster))
    + ggrastr::rasterise(geom_point(size=1, alpha=1))
    + theme_void()
    + labs(title='clusters')
    + theme(plot.title = element_text(color = "black", hjust=0.5))
    + coord_fixed()
)

count_bias <-
(
    ggplot(df, aes(tSNE1, tSNE2, color=log2_count_bias))
    + ggrastr::rasterise(geom_point(size=1, alpha=1))
    + theme_void()
    + scale_colour_gradient2(
      low = "blue",
      mid = "gray",
      high = "red",
    )
    #+ sc_blue
    + labs(title='count bias in clusters (normalized, replicate average, log2 scale)')
    + theme(plot.title = element_text(color = "black", hjust=0.5))
    + coord_fixed()
)


count_bias_female <-
(
    ggplot(ff, aes(tSNE1, tSNE2, color=pct_female))
    + ggrastr::rasterise(geom_point(size=1, alpha=1))
    + theme_void()
    + scale_colour_gradient(
      low = "gray",
      high = "red",
    )
    + labs(title='% female cells (replicate average) in clusters (shown on female cells only)')
    + theme(plot.title = element_text(color = "red", hjust=0.5))
    + coord_fixed()
)

count_bias_male <-
(
    ggplot(mf, aes(tSNE1, tSNE2, color=pct_male))
    + ggrastr::rasterise(geom_point(size=1, alpha=1))
    + theme_void()
    + scale_colour_gradient(
      low = "gray",
      high = "blue",
    )
    + labs(title='% male cells (replicate average) in clusters (shown on male cells only)')
    + theme(plot.title = element_text(color = "blue", hjust=0.5))
    + coord_fixed()
)

p = (
     clusters
     + count_bias
     + count_bias_female
     + count_bias_male
     + plot_layout(ncol = 2, widths = 1, guides = 'collect')
)

ggsave(out_file, p, width=20, height=15)

