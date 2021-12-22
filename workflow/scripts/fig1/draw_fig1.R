library(plotgardener)
library(patchwork)
library(tidyverse)
library(ggthemes)
library(ggrepel)
library(anndata)
library(glue)
source("workflow/scripts/utils.R")

expr_h5ad_fmt <- "scraps/lognorm_h5ad/lognorm_{tissue}_stringent.h5ad"
out_csv_body <- "exports/fig1/fig1_data_body.csv"
out_csv_head <- "exports/fig1/fig1_data_head.csv"
out_pdf_file <- "exports/fig1/fig1.pdf"


body_bias <- read.csv(out_csv_body, row.names=1)
head(body_bias)

head_bias <- read.csv(out_csv_head, row.names=1)
head(head_bias)

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
# since some cells were filtered, use only the cells
# in the clusters from the bias table
body <- left_join(body_bias, body_cells)
head(body)

head_cells <- get_cell_data("head")
# since some cells were filtered, use only the cells
# in the clusters from the bias table
head <- left_join(head_bias, head_cells)
head(head)

theme_base_size = 9

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

nonrp_bias_colors <- c(
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
    + theme_void(base_size=theme_base_size)
    #+ geom_text_repel(data=centroids, aes_string(label=label_col), color="black",
    #                  min.segment.length = 0, max.overlaps=Inf, size=4)
    + scale_color_manual(values=color_values)
    + labs(title=title)
    + theme(plot.title = element_text(color = "black", hjust=0.5))
    + theme(legend.position = "bottom")
    + theme(legend.box.margin = margin(t = 0, r = 5, b = 0, l = 5, unit = "cm"))
    + guides(color = guide_legend(
        legend_name,
        title.position = "top",
        title.hjust = 0.5,
        override.aes = list(size=3), nrow=2, byrow=TRUE)
    )
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
    + theme_few(base_size=theme_base_size)
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
    + theme_few(base_size=theme_base_size)
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
    + theme_few(base_size=theme_base_size)
    + scale_color_manual(values = rp_bias_colors, drop=FALSE)
    + guides(color=NULL) #legend.position='none')
    + theme(plot.title = element_text(color = "black", hjust=0.5))
  )
}

get_scatter_nonrp <- function(df, title) {
  (
    ggplot(df, aes(y=pct_nonrp_female, x=pct_nonrp_male, color=nonrp_bias_type))
    + geom_abline(intercept = 0, slope = 2, linetype='dotted')
    + geom_abline(intercept = 0, slope = 0.5, linetype='dotted')
    + geom_abline(intercept = 0, slope = 1, linetype='dashed', color="gray75")
    #+ geom_text_repel(aes(label=nonrp_bias_label), min.segment.length = 0, max.overlaps=Inf)
    + geom_point(size=1)
    + labs(
        title=title,
        y='% non-RP genes female biased',
        x='% non-RP genes male biased'
    )
    #+ coord_fixed()
    + theme_few(base_size=theme_base_size)
    + scale_color_manual(values = nonrp_bias_colors, drop=FALSE)
    + theme(legend.position='none')
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
  title='RP gene expr bias in body',
  legend_name="Expression bias"
)

body_tsne_nonrp_bias <- get_tsne(
  body,
  color_col='nonrp_bias_type',
  label_col='nonrp_bias_label',
  color_values = nonrp_bias_colors,
  title='NonRP expr bias in body',
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
  title='RP gene expr bias in head',
  legend_name="Expression bias"
)

head_tsne_nonrp_bias <- get_tsne(
  head,
  color_col='nonrp_bias_type',
  label_col='nonrp_bias_label',
  color_values = nonrp_bias_colors,
  title='NonRP expr bias in head',
  legend_name="Expression bias"
)

body_scatter_count_bias <- get_scatter_count(body_bias, title='count bias in body')
body_scatter_expr_bias <- get_scatter_expr(body_bias, title='expression bias in body')
body_scatter_rp_bias <- get_scatter_rp(body_bias, title='RP expr bias in body')
body_scatter_nonrp_bias <- get_scatter_nonrp(body_bias, title='NonRP expr bias in body')

head_scatter_count_bias <- get_scatter_count(head_bias, title='count bias in head')
head_scatter_expr_bias <- get_scatter_expr(head_bias, title='expression bias in head')
head_scatter_rp_bias <- get_scatter_rp(head_bias, title='RP expr bias in head')
head_scatter_nonrp_bias <- get_scatter_nonrp(head_bias, title='NonRP expr bias in head')


tsnes <- (
  body_tsne_count_bias + body_tsne_expr_bias + body_tsne_rp_bias + body_tsne_nonrp_bias
  + head_tsne_count_bias + head_tsne_expr_bias + head_tsne_rp_bias + head_tsne_nonrp_bias
  + plot_layout(ncol=4, guides="collect")
  + plot_annotation(tag_levels = list(c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H')))
  & theme(plot.tag = element_text(face = 'bold'))
  & theme(legend.position="bottom")
)

scatters <- (
  body_scatter_count_bias + body_scatter_expr_bias + body_scatter_rp_bias + body_scatter_nonrp_bias 
 + head_scatter_count_bias + head_scatter_expr_bias + head_scatter_rp_bias + head_scatter_nonrp_bias 
  + plot_layout(ncol=4, guides="collect")
  + plot_annotation(tag_levels = list(c('I', 'J', 'K', 'L', 'M', 'N', 'O', 'P')))
  & theme(plot.tag = element_text(face = 'bold'))
  & theme(legend.position="none")
)

pdf(out_pdf_file, width=11, height=13)

fig_width = 8.5
fig_height = 11

nrow = 4
ncol = 4

panel_width = fig_width / ncol
panel_height = fig_height / nrow

pageCreate(
  width = fig_width, height = fig_height, default.units = "inches"
)

plotGG(
  plot = tsnes,
  x = 0*panel_width, y = 0*panel_height,
  width = 4*panel_width, height = 2*panel_height,
  just = c("left", "top")
)

plotGG(
  plot = scatters,
  x = 0*panel_width, y = 2*panel_height,
  width = 4*panel_width, height = 2*panel_height,
  just = c("left", "top")
)

dev.off()

