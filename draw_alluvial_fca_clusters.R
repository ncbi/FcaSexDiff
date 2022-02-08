library(tidyverse)
library(ggalluvial)
library(ggrepel)
library(anndata)
source("workflow/scripts/utils.R")

#infile <- snakemake@input[[1]]
#outfile <- snakemake@output[[1]]

process <- function(tissue) {

cluster_info_file <- paste0("exports/fig1/fig1_data_", tissue, ".csv")
expr_file <- paste0("scraps/filtered_lognorm_h5ad/cellfilter~NoSexspecArtef",
                    "/filtered_lognorm_",
                    tissue,
                    "_stringent_NoSexspecArtef.h5ad")
cell_meta <- (
  get_obs_anndata(read_h5ad(expr_file))
  %>% select(CellID, sex, L6.0, annotation)
  %>% count(sex, L6.0, annotation)
  %>% rename(freq = n)
)

head(cell_meta)
names(cell_meta)

cluster_info <- (
  read.csv(cluster_info_file)
  %>% select(cluster, count_bias_type, expr_bias_type, rp_bias_type, nonrp_bias_type)
  %>% add_row(cluster = c("female", "male"), count_bias_type = c("female", "male"), expr_bias_type = c("female", "male"), rp_bias_type = c("female", "male"), nonrp_bias_type = c("female", "male"))
  %>% rename(stratum = cluster)
)

head(cluster_info)
names(cluster_info)

clusters <- (
  cell_meta
  %>% mutate(Sex=sex)
)
head(clusters)
dim(clusters)

lodes <- (
  to_lodes_form(
    clusters,
    axes = c("sex", "L6.0", "annotation"), #setdiff(names(clusters), c('Sex', 'freq')),
    key = "resolution"
  )
  %>% left_join(cluster_info)
  %>% replace_na(list(count_bias_type="unknown2", expr_bias_type="unknown1"))
  %>% mutate(stratum = stringr::str_wrap(stratum, 25))
)
head(lodes)
dim(lodes)


bias_colors <- c(
    "female" = "red",
    "FemaleOnly" = "#A60000", #"red",
    "FemaleSignificant" = "#FF5233", #"red",
    "FemaleNonsignificant" = "#FFbbbc", # "#BC544B", #"pink",
    "Unbiased" = "gray",
    "MaleNonsignificant" = "#B2D8FF", #"#73C2FB", #maya
    "MaleSignificant" = "#023672", #"blue",
    "MaleOnly" = "#021732", #"blue"
    "male" = "blue",
    "unknown1" = "black",
    "unknown2" = "gray95"
)
outfile <- paste0("alluvial_", tissue, ".pdf")

pdf(outfile, width=20, height=50)

p <- ggplot(lodes, aes(x = resolution, y=freq, stratum = stratum, alluvium = alluvium)) +
  geom_alluvium(aes(fill = Sex), width = 1/12, alpha = 0.5) +
  geom_stratum(aes(fill=count_bias_type), width = 1/12, color = "grey25") +
  geom_label_repel(
    stat = "stratum",
    aes(label = after_stat(stratum), color=expr_bias_type),
    hjust = 0,
    nudge_x = 0.2,
    direction = "y",
    alpha = 0.8,
    # Repel away from the left edge, not from the right.
    xlim = c(NA, Inf)
  ) +
  coord_cartesian(clip = "off") +
  #scale_fill_brewer(type = "qual", palette = "Set1") +
  scale_fill_manual(values = bias_colors, name="") +
  scale_color_manual(values = bias_colors, name=NULL) +
  ggtitle("Cell composition FCA clusters, box colored by count bias, label colored by expression bias") +
  guides(colour = "none") +
  theme_minimal(base_size=11)

print(p)

dev.off()

}

#process("body")
process("head")
