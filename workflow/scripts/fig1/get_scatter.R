library(tidyverse)
library(ggthemes)
library(ggrepel)
library(anndata)
library(glue)

scatter_data_file <- snakemake@input[["scatter"]]
outfile <- snakemake@output[[1]]
resol <- snakemake@wildcards[["resol"]]
bias_type <- snakemake@wildcards[["bias"]]

print(bias_type)



bias_data <- read.csv(scatter_data_file)
head(bias_data)

theme_base_size = 9

get_scatter <- function(df, x_col, y_col, color_col, label_col, x_label, y_label, title, legend_name) {
    df <- (
        df
        %>% rename(x := !!x_col)
        %>% rename(y := !!y_col)
        %>% rename(color := !!color_col)
        %>% rename(label := !!label_col)
        %>% select(cluster, x, y, color, label, major_annotation)
    )
    color_values <- unique(df %>% pull(color))
    color_values <- setNames(color_values, color_values)
    print(color_values)

    axis_limit <- max(max(df$x), max(df$y))
    p1 <- (
        ggplot(df, aes(x, y, color=color))
        + geom_abline(intercept = 0, slope = 2, linetype='dotted')
        + geom_abline(intercept = 0, slope = 0.5, linetype='dotted')
        + geom_abline(intercept = 0, slope = 1, linetype='dashed', color="gray75")
        #+ geom_errorbar(aes(ymin=ymin, ymax=ymax), width=0.5, size=0.5, color="gray")
        #+ geom_errorbarh(aes(xmin=xmin, xmax=xmax), height=0.5, size=0.5, color="gray")
        + geom_point(size=1)
        + labs(title = title, y = y_label, x = x_label)
        + coord_fixed()
        + xlim(0, axis_limit)
        + ylim(0, axis_limit)
        + theme_few(base_size=theme_base_size)
        + scale_color_manual(values=color_values)
        #+ scale_color_manual(values = count_bias_colors, drop=FALSE)
        + theme(legend.position='none')
        + theme(plot.title = element_text(color = "black", hjust=0.5))
   )
    p2 <- (
        p1
        + geom_text_repel(aes(label=label), min.segment.length = 0, max.overlaps = Inf)
    )
   list(p1, p2)
}

x_col <- list(
    count = "pct_cells_male",
    expr = "pct_genes_male",
    rp = "pct_rp_male",
    rpavg = "rp_avg_male",
    nonrp = "pct_nonrp_male"
)

x_label <- list(
    count = "% cells male (replicate average)",
    expr = "% genes male biased",
    rp = "% RP genes male biased",
    rpavg = "Avg RP expression in male cells",
    nonrp = "% non-RP genes male biased"
)

y_col <- list(
    count = "pct_cells_female",
    expr = "pct_genes_female",
    rp = "pct_rp_female",
    rpavg = "rp_avg_female",
    nonrp = "pct_nonrp_female"
)

y_label <- list(
    count = "% cells female (replicate average)",
    expr = "% genes female biased",
    rp = "% RP genes female biased",
    rpavg = "Avg RP expression in female cells",
    nonrp = "% non-RP genes female biased"
)

color_col <- list(
    count = "count_color",
    expr = "expr_color",
    rp = "rp_color",
    rpavg = "rp_color",
    nonrp = "nonrp_color"
)

label_col <- list(
    count = "count_label",
    expr = "expr_label",
    rp = "rp_label",
    rpavg = "rp_label",
    nonrp = "nonrp_label"
)

title <- list(
    count = "Count bias",
    expr = "Expression bias",
    rp = "RP expr bias",
    rpavg = "Avg RP expr",
    nonrp = "NonRP expr bias"
)

legend <- list(
    count = "Count bias",
    expr = "Expression bias",
    rp = "RP expr bias",
    rpavg = "Avg RP expr",
    nonrp = "NonRP expr bias"
)

plots <- get_scatter(
  bias_data,
  x_col[[bias_type]],
  y_col[[bias_type]],
  color_col[[bias_type]],
  label_col[[bias_type]],
  x_label[[bias_type]],
  y_label[[bias_type]],
  title[[bias_type]],
  legend[[bias_type]]
)

pdf(outfile, width = 8, height = 8)
print(plots[[1]])
print(plots[[2]])
dev.off()
