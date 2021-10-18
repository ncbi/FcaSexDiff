library(anndata)
library(tidyverse)
library(Matrix)
library(ggrepel)

infile <- snakemake@input[[1]]
outfile <- snakemake@output[[1]]

print(infile)
print(outfile)

adata <- read_h5ad(infile)

adata

bias = (
    as.data.frame(as.matrix(adata$X))
    %>% rownames_to_column("symbol")
    %>% gather("cluster", "bias", -symbol)
)
head(bias)


avg_all = (
    as.data.frame(adata$layers[["avg_all"]])
    %>% rownames_to_column("symbol")
    %>% gather("cluster", "avg_all", -symbol)
)
head(avg_all)

avg_male = (
    as.data.frame(adata$layers[["avg_male"]])
    %>% rownames_to_column("symbol")
    %>% gather("cluster", "avg_male", -symbol)
)
head(avg_male)


avg_female = (
    as.data.frame(adata$layers[["avg_female"]])
    %>% rownames_to_column("symbol")
    %>% gather("cluster", "avg_female", -symbol)
)
head(avg_female)

log2fc = (
    as.data.frame(adata$layers[["log2fc"]])
    %>% rownames_to_column("symbol")
    %>% gather("cluster", "log2fc", -symbol)
)
head(log2fc)


padj = (
    as.data.frame(adata$layers[["padj"]])
    %>% rownames_to_column("symbol")
    %>% gather("cluster", "padj", -symbol)
)
head(padj)

df = (
    inner_join(bias, avg_all)
    %>% inner_join(avg_male)
    %>% inner_join(avg_female)
    %>% inner_join(log2fc)
    %>% inner_join(padj)
    %>% filter(!(cluster %in% c("unannotated", "artefact")))
    %>% dplyr::mutate(label = ifelse(grepl(pattern="RpL", symbol), symbol, ""))
    %>% mutate(cluster = gsub("adult", "", cluster))
    %>% mutate(cluster = gsub("Malpighian tubule", "MT", cluster))
    %>% mutate(bias = ifelse(bias > 0, "female_biased",
                             ifelse(bias < 0, "male_biased", "unbiased")))
    %>% mutate(bias = factor(bias, levels=c("unbiased", "male_biased", "female_biased"), ordered=TRUE))
    %>% arrange(bias)
)

head(df)

# Because of pseudo count 1e-9, there are blobs around +/- 25 for log2fc
# Generate same sets of plots removing those blobs
df2 = filter(df, abs(log2fc) < 10)


# identify specific genes on these plot I particular ribosomal proteins will be informative (RpL…. and RpS)

ncluster = length(unique(df$cluster))
n_plot_cols = 6
n_plot_rows = ceiling(ncluster/n_plot_cols)

pdf(outfile, height=2*n_plot_rows, width=2*n_plot_cols)


#1. scatter plots for the expression in the in the major Malpighian clusters. genes (use both log and no log expression since I am not sure what would be better for the visualization here).

(
    ggplot(df, aes(avg_male, avg_female))
    + ggrastr::rasterise(geom_point(aes(color=bias), size=0.2, alpha=0.2))
    + facet_wrap(~cluster, ncol=n_plot_cols)
    + ggtitle("Female expression vs male expression")
    + scale_color_manual(values=c("female_biased"="red", "male_biased"="blue", "unbiased"="black"))
    + theme_minimal()
    + theme(legend.position="bottom")
)
(
    ggplot(df2, aes(avg_male, avg_female))
    + ggrastr::rasterise(geom_point(aes(color=bias), size=0.2, alpha=0.2))
    + facet_wrap(~cluster, ncol=n_plot_cols)
    + ggtitle("Female expression vs male expression (genes with fold change 2^10 or less)")
    + scale_color_manual(values=c("female_biased"="red", "male_biased"="blue", "unbiased"="black"))
    + theme_minimal()
    + theme(legend.position="bottom")
)

# 1 but on log scale

(
    ggplot(df, aes(log1p(avg_male), log1p(avg_female)))
    + ggrastr::rasterise(geom_point(aes(color=bias), size=0.2, alpha=0.2))
    + facet_wrap(~cluster, ncol=n_plot_cols)
    + ggtitle("Female expression vs male expression both in log scale")
    + scale_color_manual(values=c("female_biased"="red", "male_biased"="blue", "unbiased"="black"))
    + theme_minimal()
    + theme(legend.position="bottom")
)
(
    ggplot(df2, aes(log1p(avg_male), log1p(avg_female)))
    + ggrastr::rasterise(geom_point(aes(color=bias), size=0.2, alpha=0.2))
    + facet_wrap(~cluster, ncol=n_plot_cols)
    + ggtitle("Female expression vs male expression both in log scale (genes with fold change 2^10 or less)")
    + scale_color_manual(values=c("female_biased"="red", "male_biased"="blue", "unbiased"="black"))
    + theme_minimal()
    + theme(legend.position="bottom")
)

#2. suggested by Brian Sex-biased female:male expression ratio (log2) versus average expression intensity (log2) 

(
    ggplot(df, aes(log1p(avg_all), log2fc))
    + ggrastr::rasterise(geom_point(aes(color=bias), size=0.2, alpha=0.2))
    + facet_wrap(~cluster, ncol=n_plot_cols)
    + ggtitle("Fold change vs average expression")
    + scale_color_manual(values=c("female_biased"="red", "male_biased"="blue", "unbiased"="black"))
    + theme_minimal()
    + theme(legend.position="bottom")
)
(
    ggplot(df2, aes(log1p(avg_all), log2fc))
    + ggrastr::rasterise(geom_point(aes(color=bias), size=0.2, alpha=0.2))
    + facet_wrap(~cluster, ncol=n_plot_cols)
    + ggtitle("Fold change vs average expression (genes with fold change 2^10 or less)")
    + scale_color_manual(values=c("female_biased"="red", "male_biased"="blue", "unbiased"="black"))
    + theme_minimal()
    + theme(legend.position="bottom")
)

#3. Same as 2 but rather than expression ratio use log p-value for difference with positive sign for female and negative sign for males

ndf = df %>% mutate(log10padj = ifelse(bias == "female_biased", -log10(padj), log10(padj)))
ndf2 = df2 %>% mutate(log10padj = ifelse(bias == "female_biased", -log10(padj), log10(padj)))

(
    ggplot(ndf, aes(log1p(avg_all), log10padj))
    + ggrastr::rasterise(geom_point(aes(color=bias), size=0.2, alpha=0.2))
    + facet_wrap(~cluster, ncol=n_plot_cols)
    + ggtitle("Pvalue vs average expression")
    + scale_color_manual(values=c("female_biased"="red", "male_biased"="blue", "unbiased"="black"))
    + theme_minimal()
    + theme(legend.position="bottom")
)
(
    ggplot(ndf2, aes(log1p(avg_all), log10padj))
    + ggrastr::rasterise(geom_point(aes(color=bias), size=0.2, alpha=0.2))
    + facet_wrap(~cluster, ncol=n_plot_cols)
    + ggtitle("Pvalue vs average expression (genes with fold change 2^10 or less)")
    + scale_color_manual(values=c("female_biased"="red", "male_biased"="blue", "unbiased"="black"))
    + theme_minimal()
    + theme(legend.position="bottom")
)

#4. Volcano plots

(
    ggplot(df, aes(log2fc, -log10(padj)))
    + ggrastr::rasterise(geom_point(aes(color=bias), size=0.2, alpha=0.2))
    + facet_wrap(~cluster, ncol=n_plot_cols)
    + ggtitle("Volcano plot")
    + geom_text_repel(aes(label=label))
    + scale_color_manual(values=c("female_biased"="red", "male_biased"="blue", "unbiased"="black"))
    + theme_minimal()
    + theme(legend.position="bottom")
)
(
    ggplot(df2, aes(log2fc, -log10(padj)))
    + ggrastr::rasterise(geom_point(aes(color=bias), size=0.2, alpha=0.2))
    + facet_wrap(~cluster, ncol=n_plot_cols)
    + ggtitle("Volcano plot (genes with fold change 2^10 or less)")
    + geom_text_repel(aes(label=label))
    + scale_color_manual(values=c("female_biased"="red", "male_biased"="blue", "unbiased"="black"))
    + theme_minimal()
    + theme(legend.position="bottom")
)

