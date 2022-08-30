library(tidyverse)
library(patchwork)

infiles <- snakemake@input
pdf_file <- snakemake@output[[1]]

input <- (
    stack(infiles)
    %>% filter(ind != "")
    %>% rename(tissue := ind, path := values)
)

tissues <- c(
    "head",
    "body",
    "antenna",
    "proboscis_and_maxillary_palps",
    "haltere",
    "wing",
    "body_wall",
    "leg",
    "heart",
    "gut",
    "malpighian_tubule",
    "fat_body",
    "oenocyte",
    "trachea"
)

keep_groups <- c(
    "AllExpressed",
    "FemaleOnlyBiased",
    "MixedBiased",
    "MaleOnlyBiased",
    "AllBiased"
)

enrichment_types <- c(
    "dsx_occupancy",
    "dsx_one_hop",
    "dsx_two_hop",
    "dsx_ff_y",
    "dsx_ff_x",
    "dsx_ff_x_high",
#    "putative_tf",
    "trusted_tf"
)

group_colors <- c(
    "AllExpressed" = "gray50",
    "FemaleOnlyBiased" = "#ED746D", # "red",
    "MixedBiased" = "#C77BF8",  #"purple",
    "MaleOnlyBiased" = "#4D71BE", #"blue"
    "AllBiased" = "#868049" #"blue"
)
print(group_colors)

df <- (
    input
    %>% group_by(tissue)
    %>% summarize(res = read.delim(path))
    %>% unnest(cols = c(res))
    %>% gather("enrichment_type", "enrichment_value", -tissue, -bias_group)
    %>% filter(enrichment_type != "expressed")
    %>% filter(bias_group %in% keep_groups)
)

head(df)

pvals <- (
    input
    %>% mutate(path = paste0(tools::file_path_sans_ext(path), "_pvals.tsv"))
    %>% group_by(tissue)
    %>% summarize(res = read.delim(path))
    %>% unnest(cols = c(res))
    %>% gather("enrichment_type", "enrichment_pvalue", -tissue, -bias_group)
    %>% filter(enrichment_type != "expressed")
    %>% mutate(sig = ifelse(enrichment_pvalue < 0.001, "***",
                            ifelse(enrichment_pvalue < 0.01, "**",
                                   ifelse(enrichment_pvalue < 0.05, "*", ""))))
    %>% filter(bias_group %in% keep_groups)
)

head(pvals)

df <- (
    merge(df, pvals)
    %>% mutate(tissue = factor(tissue, levels = tissues, ordered = TRUE))
    %>% mutate(enrichment_type = factor(
        enrichment_type, levels = enrichment_types, ordered = TRUE
    ))
    %>% mutate(GeneGroup = factor(bias_group, levels = keep_groups))
)

p1 <- (
    df
    %>% filter(enrichment_type != "trusted_tf")
    %>% ggplot(aes(x = GeneGroup, y = enrichment_value, fill=GeneGroup))
    + geom_bar(stat = "identity", position = "dodge")
    + geom_text(aes(label = sig), vjust = 1, size = 5)
    + facet_grid(tissue ~ enrichment_type)
    #+ ggh4x::facet_grid2(
    #    tissue ~ enrichment_type, scales = "free_y", independent = "y"
    #)
    + scale_fill_manual(values = group_colors)
    + labs(
      x = "Biased gene group",
      y = "Fraction of genes within group"
    )
    + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    + theme(legend.position = "top")
)

expected <- df %>% filter(GeneGroup == "AllExpressed")
femaleonly <- df %>% filter(GeneGroup == "FemaleOnlyBiased")
maleonly <- df %>% filter(GeneGroup == "MaleOnlyBiased")
mixed <- df %>% filter(GeneGroup == "MixedBiased")
biased <- df %>% filter(GeneGroup == "AllBiased")

femaleonly <- (
    left_join(femaleonly, expected,
              by = c("tissue", "enrichment_type"), suffix=c("", ".exp"))
    %>% mutate(enrichment_value = enrichment_value / enrichment_value.exp)
)
maleonly <- (
    left_join(maleonly, expected,
              by = c("tissue", "enrichment_type"), suffix=c("", ".exp"))
    %>% mutate(enrichment_value = enrichment_value / enrichment_value.exp)
)
mixed <- (
    left_join(mixed, expected,
              by = c("tissue", "enrichment_type"), suffix=c("", ".exp"))
    %>% mutate(enrichment_value = enrichment_value / enrichment_value.exp)
)
biased <- (
    left_join(biased, expected,
              by = c("tissue", "enrichment_type"), suffix=c("", ".exp"))
    %>% mutate(enrichment_value = enrichment_value / enrichment_value.exp)
)

p3 <- (
    rbind(femaleonly, maleonly, mixed) #, biased)
    %>% mutate(enrichment_value = ifelse(
        enrichment_value > 10, 10,  enrichment_value
    ))
    %>% filter(enrichment_type != "trusted_tf")
    %>% ggplot(aes(x = GeneGroup, y = enrichment_value, fill = GeneGroup))
    + geom_bar(stat = "identity", position = "dodge")
    + geom_text(aes(label = sig), vjust = 1, size = 5)
    + facet_grid(tissue ~ enrichment_type)
    + scale_fill_manual(values = group_colors)
    + labs(
      x = "Biased gene group",
      y = "Fold change in comparison to AllExpressed"
    )
    + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    + theme(legend.position = "top")
)

p2 <- (
    df %>% filter(enrichment_type == "trusted_tf")
    %>% ggplot(aes(x = GeneGroup, y = enrichment_value, fill = GeneGroup))
    + geom_bar(stat = "identity", position = "dodge")
    + geom_text(aes(label = sig), vjust = 1, size = 5)
    + facet_grid(tissue ~ enrichment_type)
    + scale_fill_manual(values = group_colors)
    + labs(
        x = "Biased gene group",
        y = "Fraction of genes within group"
    )
    + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
)

pdf("dsx_targets_tfs_enrichments.pdf", height=12, width=18)
(
    p1 + p3 + p2 + plot_layout(nrow = 1, widths = c(4,4,1), guides = "collect")
    & theme(legend.position = "top")
)
dev.off()

