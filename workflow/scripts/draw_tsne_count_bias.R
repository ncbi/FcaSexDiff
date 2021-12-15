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
rds_file = snakemake@output[[2]]

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
  %>% mutate(bias_label = as.character(count_bias_type))
  %>% mutate(bias_label = factor(bias_label,
                                 levels = c("MaleOnly", "MaleSignificant", "MaleNonsignificant", 
                                            "Unbiased",
                                            "FemaleNonsignificant", "FemaleSignificant",  "FemaleOnly"),
                                 ordered=TRUE))
  %>% mutate(cluster_label = ifelse(bias_label != "Unbiased", as.character(cluster), ''))
)

if (resol != "annotation") {
  # label clusters without the prefix showing cluster resolution
  bias <- mutate(bias, cluster_label = substring(cluster_label, nchar(!!resol)+2))
}

head(bias)

df <- (
  read_h5ad(expr_file)$obs
  %>% filter(sex %in% c("female", "male"))
  %>% mutate(cluster = .[,resol]) # copy column named by resol as cluster
  %>% select(tSNE1, tSNE2, cluster, sex)
)

head(df)

df <- right_join(df, bias)

head(df)


centroids <- (
  select(df, cluster, tSNE1, tSNE2, cluster_label, bias_label)
  %>% group_by(cluster)
  %>% summarize(tSNE1 = mean(tSNE1), tSNE2 = mean(tSNE2),
                cluster_label = dplyr::first(cluster_label),
                bias_label = dplyr::first(bias_label))
)
head(centroids)


ff = df %>% filter(sex=='female')
mf = df %>% filter(sex=='male')


clusters <-
(
    ggplot(df, aes(tSNE1, tSNE2, color=cluster))
    + ggrastr::rasterise(geom_point(size=.05, alpha=1))
    + theme_void()
    + labs(title='clusters')
    + theme(plot.title = element_text(color = "black", hjust=0.5))
    + coord_fixed()
    + theme(legend.position = "none")
)

count_bias <-
(
    ggplot(df, aes(tSNE1, tSNE2, color=bias_label))
    + ggrastr::rasterise(geom_point(size=.05, alpha=1))
    + theme_void()
    #+ scale_colour_gradient2(
    #  low = "blue",
    #  mid = "gray",
    #  high = "red",
    #)
    #+ sc_blue
    + geom_text_repel(data=centroids, aes(label=cluster_label), color="black",
                      min.segment.length = 0, max.overlaps=Inf, size=4)
    + scale_color_manual(values = c(
        "FemaleOnly" = "red",
        "FemaleSignificant" = "red",
        "FemaleNonsignificant" = "#BC544B", #"pink",
        "MaleOnly" = "blue",
        "MaleSignificant" = "blue",
        "MaleNonsignificant" = "#73C2FB", #maya
        "Unbiased" = "gray"
    ), drop = FALSE)
    + labs(title='count bias in clusters')
    + theme(plot.title = element_text(color = "black", hjust=0.5))
    + theme(legend.position = "bottom")
    + guides(color = guide_legend(override.aes = list(size=8)))
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
     #+ count_bias_female
     #+ count_bias_male
     + plot_layout(ncol = 2, widths = 1)
)

ggsave(out_file, p, width=20, height=10)
saveRDS(count_bias, rds_file)

