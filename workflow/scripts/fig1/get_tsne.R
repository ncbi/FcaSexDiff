library(tidyverse)
library(ggthemes)
library(ggrepel)
library(anndata)
library(glue)

tsne_data_file <- snakemake@input[["tsne"]]
outfile <- snakemake@output[[1]]
resol <- snakemake@wildcards[["resol"]]
bias_type <- snakemake@wildcards[["bias"]]

print(bias_type)



bias_data <- readRDS(tsne_data_file)
head(bias_data)


theme_base_size = 9

get_tsne <- function(df, color_col, label_col, title, legend_name) {
    df <- (
        df
        %>% rename(color := !!color_col)
        %>% rename(label := !!label_col)
        %>% select(cluster, tSNE1, tSNE2, color, label, major_annotation)
    )

    color_values <- unique(df %>% pull(color))
    color_values <- setNames(color_values, color_values)
    print(color_values)

    centroids <- (
        df %>% group_by(cluster)
        %>% summarize(
            tSNE1 = mean(tSNE1)
            , tSNE2 = mean(tSNE2)
            , color = dplyr::first(color)
            , label = dplyr::first(label)
            , major_annotation = dplyr::first(major_annotation)
        )
    )
    print(head(centroids))


    p1 <- (
        ggplot(data = df, mapping = aes(x = tSNE1, y = tSNE2, color = color))
        #+ ggrastr::rasterise(geom_point(size=.0001, alpha=1))
        + geom_point(size=.00001, alpha=1)
        + theme_void(base_size = theme_base_size)
        + scale_color_manual(values=color_values)
        + labs(title=title)
        + theme(plot.title = element_text(color = "black", hjust=0.5))
        + theme(legend.position = "bottom")
        + theme(legend.box.margin = margin(t = 0, r = 5, b = 0, l = 5, unit = "cm"))
        #+ guides(color = guide_legend(
        #    legend_name
        #    , title.position = "top"
        #    , title.hjust = 0.5
        #    , override.aes = list(size = 3)
        #    #, nrow = 1
        #    #, byrow = TRUE
        #))
        #+ scale_x_continuous(limits = c(xmin, xmax), expand = expansion(mult = 0.75))
        #+ scale_y_continuous(limits = c(ymin, ymax), expand = expansion(mult = 0.75))
        + coord_cartesian(clip = "off")
        + theme(legend.position = "none")
        + coord_fixed()
    )

    cent1 <- centroids %>% filter(tSNE1 >= 0, tSNE2 >= 0)
    cent2 <- centroids %>% filter(tSNE1 < 0, tSNE2 >= 0)
    cent3 <- centroids %>% filter(tSNE1 < 0, tSNE2 < 0)
    cent4 <- centroids %>% filter(tSNE1 >= 0, tSNE2 < 0)

    xmax <- max(df$tSNE1)
    xmin <- min(df$tSNE1)
    ymax <- max(df$tSNE2)
    ymin <- min(df$tSNE2)

    #xmax <- xmax + 0.1 * (xmax - xmin)
    #xmin <- xmin - 0.1 * (xmax - xmin)
    #ymax <- ymax + 0.1 * (ymax - ymin)
    #ymin <- ymin - 0.1 * (ymax - ymin)
    
    p2 <- (
        p1
        + geom_label_repel(
            data = cent1
            , mapping = aes(label = label)
            #, color="black"
            , xlim = c(xmax, NA)
            , ylim = c(0, NA)
            , hjust = 0
            , min.segment.length = 0
            , max.overlaps = 3
            , size = 4
            , alpha = 0.7
            , show.legend = FALSE
        )
        + geom_label_repel(
            data = cent2
            , mapping = aes(label = label)
            #, color="black"
            , xlim = c(NA, xmin)
            , ylim = c(0, NA)
            , hjust = 1
            , min.segment.length = 0
            , max.overlaps = 3
            , size = 4
            , alpha = 0.7
            , show.legend = FALSE
        )
        + geom_label_repel(
            data = cent4
            , mapping = aes(label = label)
            #, color="black"
            , xlim = c(xmax, NA)
            , ylim = c(NA, 0)
            , hjust = 0
            , min.segment.length = 0
            , max.overlaps = 3
            , size = 4
            , alpha = 0.7
            , show.legend = FALSE
        )
        + geom_label_repel(
            data = cent3
            , mapping = aes(label = label)
            #, color="black"
            , xlim = c(NA, xmin)
            , ylim = c(NA, 0)
            , hjust = 1
            , min.segment.length = 0
            , max.overlaps = 3
            , size = 4
            , alpha = 0.7
            , show.legend = FALSE
        )
    )
    list(p1, p2)
}

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

plots <- get_tsne(
  bias_data,
  color_col[[bias_type]],
  label_col[[bias_type]],
  title[[bias_type]],
  legend[[bias_type]]
)

pdf(outfile, width = 8, height = 8)
print(plots[[1]])
print(plots[[2]])
dev.off()
