library(plotgardener)
library(patchwork)
library(tidyverse)
library(ggthemes)
library(ggrepel)
library(anndata)
library(glue)
source("workflow/scripts/utils.R")

expr_h5ad_fmt <- "scraps/lognorm_h5ad/lognorm_{tissue}_stringent.h5ad"
sexdiff_h5ad_fmt <- (
  "exports/sexdiff/cellfilter~NoSexspecArtef/resolution~L6.0/{tissue}\\
  /sexdiff_{tissue}_stringent_L6.0_NoSexspecArtef.h5ad"
)
rp_expr_fmt <- (
  "exports/extra/cellfilter~NoSexspecArtef/resolution~L6.0/{tissue}\\
  /ribosomal_clusters_{tissue}_stringent_L6.0_NoSexspecArtef.csv"
)
out_csv_body <- "exports/fig1/fig1_data_body.csv"
out_csv_head <- "exports/fig1/fig1_data_head.csv"
out_pdf_file <- "exports/fig1/fig1.pdf"

get_bias <- function(tissue) {
  sexdiff_h5ad_file = glue(sexdiff_h5ad_fmt)
  ad <- read_h5ad(sexdiff_h5ad_file)
  ngene <- ad$uns[["stats"]]["NoSexspecArtef", "n_gene"]
  print(ngene)
  bias <- (
    get_var_anndata(ad)
    %>% select(cluster, major_annotation, cluster_rep_fracs_mean_female,
               cluster_rep_fracs_mean_male, cluster_rep_fracs_sd_female,
               cluster_rep_fracs_sd_male, count_bias_padj,
               log2_count_bias, count_bias_type, female_gene, male_gene)
    %>% mutate(cluster_label = substring(cluster, nchar("L6.0")+2))
    %>% mutate(pct_cells_female = 100*cluster_rep_fracs_mean_female)
    %>% mutate(pct_cells_male = 100*cluster_rep_fracs_mean_male)
    %>% mutate(sd_cells_female = 100*cluster_rep_fracs_sd_female)
    %>% mutate(sd_cells_male = 100*cluster_rep_fracs_sd_male)
    %>% mutate(xmin = pct_cells_male - sd_cells_male)
    %>% mutate(xmax = pct_cells_male + sd_cells_male)
    %>% mutate(ymin = pct_cells_female - sd_cells_female)
    %>% mutate(ymax = pct_cells_female + sd_cells_female)
    %>% mutate(count_bias_label = ifelse(count_bias_type != "Unbiased",
                              as.character(cluster_label), ''))
    %>% mutate(pct_genes_female = 100*female_gene/ngene)
    %>% mutate(pct_genes_male = 100*male_gene/ngene)
    %>% mutate(expr_bias_type = ifelse(count_bias_type == "FemaleOnly", "FemaleOnly",
                                       ifelse(count_bias_type == "MaleOnly", "MaleOnly",
                                              ifelse(pct_genes_female > 2*pct_genes_male, "FemaleBiased",
                                                     ifelse(pct_genes_male > 2*pct_genes_female, "MaleBiased", "Unbiased")))))
    %>% mutate(expr_bias_label = ifelse((pct_genes_female > 2*pct_genes_male) | (pct_genes_male > 2*pct_genes_female),
                                       as.character(cluster_label), ''))
    %>% left_join(
      read.csv(glue(rp_expr_fmt))
      %>% replace_na(list(rp_avg_female=0, rp_avg_male=0))
    )
    %>% mutate(rp_bias_type = ifelse(count_bias_type == "FemaleOnly", "FemaleOnly",
                                     ifelse(count_bias_type == "MaleOnly", "MaleOnly",
                                            ifelse(rp_avg_female > 2*rp_avg_male, "FemaleBiased",
                                                   ifelse(rp_avg_male > 2*rp_avg_female, "MaleBiased", "Unbiased")))))
    %>% mutate(rp_bias_label = ifelse((rp_avg_female > 2*rp_avg_male) | (rp_avg_male > 2*rp_avg_female),
                                      as.character(cluster_label), ''))
  )
  print(bias)
}

get_cell_data <- function(tissue) {
  expr_file = glue(expr_h5ad_fmt)
  (
    read_h5ad(expr_file)$obs
    %>% filter(sex %in% c("female", "male"))
    %>% mutate(cluster = .[,"L6.0"])
    %>% select(tSNE1, tSNE2, cluster, sex)
  )
}

body_cells <- get_cell_data("body")
body_bias <- get_bias("body")
# since some cells were filtered, use only the cells
# in the clusters from the bias table
body <- left_join(body_bias, body_cells)
head(body)
write.csv(body_bias, out_csv_body)

head_cells <- get_cell_data("head")
head_bias <- get_bias("head")
# since some cells were filtered, use only the cells
# in the clusters from the bias table
head <- left_join(head_bias, head_cells)
head(head)
write.csv(head_bias, out_csv_head)

count_bias_colors <- c(
    "FemaleOnly" = "#A60000", #"red",
    "FemaleSignificant" = "#FF5233", #"red",
    "FemaleNonsignificant" = "#FFbbbc", # "#BC544B", #"pink",
    "Unbiased" = "gray",
    "MaleNonsignificant" = "#B2D8FF", #"#73C2FB", #maya
    "MaleSignificant" = "#023672", #"blue",
    "MaleOnly" = "#021732" #"blue"
)

expr_bias_colors <- c(
    "FemaleOnly" = "#A60000", #"red",
    "FemaleBiased" = "#FF5233",
    "Unbiased" = "gray",
    "MaleBiased" = "#023672",
    "MaleOnly" = "#021732" #"blue"
)

rp_bias_colors <- c(
    "FemaleOnly" = "#A60000", #"red",
    "FemaleBiased" = "#FF5233",
    "Unbiased" = "gray",
    "MaleBiased" = "#023672",
    "MaleOnly" = "#021732" #"blue"
)


get_tsne <- function(df, color_col, label_col, color_values, title, legend_name) {
  (
    ggplot(data=df, mapping=aes_string(x='tSNE1', y='tSNE2', color=color_col))
    + ggrastr::rasterise(geom_point(size=.05, alpha=1))
    + theme_void(base_size=10)
    #+ geom_text_repel(data=centroids, aes_string(label=label_col), color="black",
    #                  min.segment.length = 0, max.overlaps=Inf, size=4)
    + scale_color_manual(values=color_values)
    + labs(title=title)
    + theme(plot.title = element_text(color = "black", hjust=0.5))
    + theme(legend.position = "bottom")
    + guides(color = guide_legend(legend_name, override.aes = list(size=3), ncol=1))
    + coord_fixed()
  )
}

get_scatter_count <- function(df, title) {
  (
    ggplot(df, aes(y=pct_cells_female, x=pct_cells_male, color=count_bias_type))
    + geom_abline(intercept = 0, slope = 2, linetype='dotted')
    + geom_abline(intercept = 0, slope = 0.5, linetype='dotted')
    + geom_abline(intercept = 0, slope = 1, linetype='dashed', color="gray75")
    + geom_errorbar(aes(ymin=ymin, ymax=ymax), width=0.5, size=0.5, color="gray")
    + geom_errorbarh(aes(xmin=xmin, xmax=xmax), height=0.5, size=0.5, color="gray")
    #+ geom_text_repel(aes(label=count_bias_label), min.segment.length = 0, max.overlaps=Inf)
    + geom_point(size=1)
    + labs(
        title=title,
        y='% cells female (replicate average)',
        x='% cells male (replicate average)'
    )
    #+ coord_fixed()
    + theme_few(base_size=10)
    + scale_color_manual(values = count_bias_colors, drop=FALSE)
    + theme(legend.position='none')
    + theme(plot.title = element_text(color = "black", hjust=0.5))
  )
}

get_scatter_expr <- function(df, title) {
  (
    ggplot(df, aes(y=pct_genes_female, x=pct_genes_male, color=expr_bias_type))
    + geom_abline(intercept = 0, slope = 2, linetype='dotted')
    + geom_abline(intercept = 0, slope = 0.5, linetype='dotted')
    + geom_abline(intercept = 0, slope = 1, linetype='dashed', color="gray75")
    #+ geom_text_repel(aes(label=expr_bias_label), min.segment.length = 0, max.overlaps=Inf)
    + geom_point(size=1)
    + labs(
        title=title,
        y='% genes female biased',
        x='% genes male biased'
    )
    #+ coord_fixed()
    + theme_few(base_size=10)
    + scale_color_manual(values = expr_bias_colors, drop=FALSE)
    + theme(legend.position='none')
    + theme(plot.title = element_text(color = "black", hjust=0.5))
  )
}

get_scatter_rp <- function(df, title) {
  (
    ggplot(df, aes(y=rp_avg_female, x=rp_avg_male, color=rp_bias_type))
    + geom_abline(intercept = 0, slope = 2, linetype='dotted')
    + geom_abline(intercept = 0, slope = 0.5, linetype='dotted')
    + geom_abline(intercept = 0, slope = 1, linetype='dashed', color="gray75")
    #+ geom_text_repel(aes(label=rp_bias_label), min.segment.length = 0, max.overlaps=Inf)
    + geom_point(size=1)
    + labs(
        title=title,
        y='Avg RP gene expr in female',
        x='Avg RP gene expr in male'
    )
    #+ coord_fixed()
    + theme_few(base_size=10)
    + scale_color_manual(values = rp_bias_colors, drop=FALSE)
    + guides(color=NULL) #legend.position='none')
    + theme(plot.title = element_text(color = "black", hjust=0.5))
  )
}

body_tsne_count_bias <- get_tsne(
  body,
  color_col='count_bias_type',
  label_col='count_bias_label',
  color_values = count_bias_colors,
  title='count bias in body',
  legend_name="Count bias"
)

body_tsne_expr_bias <- get_tsne(
  body,
  color_col='expr_bias_type',
  label_col='expr_bias_label',
  color_values = expr_bias_colors,
  title='expression bias in body',
  legend_name="Expression bias"
)

body_tsne_rp_bias <- get_tsne(
  body,
  color_col='rp_bias_type',
  label_col='rp_bias_label',
  color_values = expr_bias_colors,
  title='RP gene expression bias in body',
  legend_name="Expression bias"
)

head_tsne_count_bias <- get_tsne(
  head,
  color_col='count_bias_type',
  label_col='count_bias_label',
  color_values = count_bias_colors,
  title='count bias in head',
  legend_name="Count bias"
)

head_tsne_expr_bias <- get_tsne(
  head,
  color_col='expr_bias_type',
  label_col='expr_bias_label',
  color_values = expr_bias_colors,
  title='expression bias in head',
  legend_name="Expression bias"
)

head_tsne_rp_bias <- get_tsne(
  head,
  color_col='rp_bias_type',
  label_col='rp_bias_label',
  color_values = rp_bias_colors,
  title='RP gene expression bias in head',
  legend_name="Expression bias"
)

body_scatter_count_bias <- get_scatter_count(body_bias, title='count bias in body')
body_scatter_expr_bias <- get_scatter_expr(body_bias, title='expression bias in body')
body_scatter_rp_bias <- get_scatter_rp(body_bias, title='RP expression bias in body')

head_scatter_count_bias <- get_scatter_count(head_bias, title='count bias in head')
head_scatter_expr_bias <- get_scatter_expr(head_bias, title='expression bias in head')
head_scatter_rp_bias <- get_scatter_rp(head_bias, title='RP expression bias in head')


tsnes <- (
  body_tsne_count_bias + body_tsne_rp_bias + body_tsne_expr_bias 
  + head_tsne_count_bias + head_tsne_rp_bias + head_tsne_expr_bias 
  + plot_layout(ncol=3, guides="collect")
  + plot_annotation(tag_levels = list(c('A', 'B', 'C', 'D', 'E', 'F')))
  & theme(plot.tag = element_text(face = 'bold'))
  & theme(legend.position="right")
)

scatters <- (
  body_scatter_count_bias + body_scatter_rp_bias + body_scatter_expr_bias 
 + head_scatter_count_bias + head_scatter_rp_bias + head_scatter_expr_bias 
  + plot_layout(ncol=3, guides="collect")
  + plot_annotation(tag_levels = list(c('G', 'H', 'I', 'J', 'K', 'L')))
  & theme(plot.tag = element_text(face = 'bold'))
  & theme(legend.position="none")
)

pdf(out_pdf_file, width=10, height=13)

fig_width = 8.5
fig_height = 11

nrow = 4
ncol = 3

panel_width = fig_width / ncol
panel_height = fig_height / nrow

pageCreate(
  width = fig_width, height = fig_height, default.units = "inches"
)

plotGG(
  plot = tsnes,
  x = 0*panel_width, y = 0*panel_height,
  width = 3*panel_width, height = 2*panel_height - 0.5,
  just = c("left", "top")
)

plotGG(
  plot = scatters,
  x = 0*panel_width, y = 2*panel_height - 0.5,
  width = 3*panel_width, height = 2*panel_height + 0.5,
  just = c("left", "top")
)

dev.off()

