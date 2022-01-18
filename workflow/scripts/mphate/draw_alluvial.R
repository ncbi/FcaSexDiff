library(tidyverse)
library(ggalluvial)
library(ggrepel)

infile <- snakemake@input[[1]]
outfile <- snakemake@output[[1]]

clusters <- (
  read.delim(infile)
  %>% mutate(Sex=sex)
)
head(clusters)
dim(clusters)

lodes <- (
  to_lodes_form(
    clusters,
    axes = setdiff(names(clusters), c('Sex', 'freq')),
    key = "resolution"
  )
  %>% mutate(stratum = stringr::str_wrap(stratum, 25))
)
head(lodes)
dim(lodes)

pdf(outfile, width=15, height=10)

ggplot(lodes, aes(x = resolution, y=freq, stratum = stratum, alluvium = alluvium)) +
  geom_label_repel(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    hjust = 0,
    nudge_x = 0.1,
    direction = "y",
    # Repel away from the left edge, not from the right.
    xlim = c(NA, Inf)
  ) +
  geom_alluvium(aes(fill = Sex), width = 1/12, alpha = 0.5) +
  geom_stratum(width = 1/12, fill = "gray25", color = "grey") +
  coord_cartesian(clip = "off") +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("Cell composition of multi-scale PHATE clusters") +
  theme_minimal(base_size=10)
