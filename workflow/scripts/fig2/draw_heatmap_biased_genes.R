library(ComplexHeatmap)
library(tidyverse)
library(anndata)
library(circlize)
library(RColorBrewer)
library(svglite)

source("workflow/scripts/utils.R")

infile <- snakemake@input[[1]]
outfile <- snakemake@output[[1]]
tissue <- snakemake@wildcards[["tissue"]]

print(infile)
print(outfile)
print(tissue)

rp_genes <- scan("resources/rp_genes_flybase_FBgg0000141.txt", what="", sep="\n")
head(rp_genes)

special_genes <- c(
    "Myc",
    "NELF-B",
    "Nelf-A",
    "Nelf-E",
    "tor",
    "Tor",
    "Sfp23F",
    "Sfp24Ba",
    "Sfp24Bb",
    "Sfp24Bc",
    "Sfp24Bd",
    "Sfp24C1",
    "Sfp24F",
    "Sfp26Ac",
    "Sfp26Ad",
    "Sfp33A1",
    "Sfp33A2",
    "Sfp33A3",
    "Sfp33A4",
    "Sfp35C",
    "Sfp36F",
    "Sfp38D",
    "Sfp51E",
    "Sfp53D",
    "Sfp60F",
    "Sfp65A",
    "Sfp70A4",
    "Sfp77F",
    "Sfp78E",
    "Sfp79B",
    "Sfp84E",
    "Sfp87B",
    "Sfp93F",
    "Sfp96F"
  )

trans_genes <- scan("resources/translation_genes.txt", what="", sep="\n")

special_genes <- c(trans_genes, special_genes)



biased_genes_small_pdf <- gsub(".pdf$", "_small.pdf", outfile)
biased_genes_small_svg <- gsub(".pdf$", ".svg", biased_genes_small_pdf)

biased_genes_tall_pdf <- gsub(".pdf$", "_tall.pdf", outfile)
biased_genes_tall_svg <- gsub(".pdf$", ".svg", biased_genes_tall_pdf)

translation_genes_pdf <- gsub(".pdf$", "_translation_genes.pdf", outfile)
translation_genes_svg <- gsub(".pdf$", ".svg", translation_genes_pdf)


ad = read_h5ad(infile)
print(ad)

cluster_info <- (
  get_var_anndata(ad)
  %>% filter(cluster_count_female > 10, cluster_count_male > 10)
  %>% select(cluster, major_annotation, cluster_count_female, cluster_count_male)
  %>% mutate(long_name = paste0(major_annotation, "_f", cluster_count_female, "_m", cluster_count_male, "_", cluster))
)

print(head(cluster_info))

mat  = as.matrix(ad$X)
zscore = ad$layers["score"]

keep_cluster = colnames(mat) %in% cluster_info$cluster

mat  = mat[, keep_cluster]
zscore = zscore[, keep_cluster]

colnames(mat)  = cluster_info$long_name
colnames(zscore)  = cluster_info$long_name

biased_genes = rowSums(abs(mat)) > 0 

mat  = mat[biased_genes, ]
zscore_biased = zscore[biased_genes, ]

color_zscore = colorRamp2(seq(-20, 20, length = 3), c("blue", "#EEEEEE", "red"))




h1 <- Heatmap(
    mat,
    show_row_names = FALSE,
    column_names_max_height = unit(10, "cm")
)
h <- draw(h1)

mat = mat[row_order(h), column_order(h)]
zscore_biased = zscore_biased[row_order(h), column_order(h)]

rp_zscore <- rbind(zscore_biased[rownames(zscore_biased) %in% rp_genes, ],
                 zscore[rownames(zscore) %in% special_genes,])

h1 <- Heatmap(
    mat,
    show_row_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_names_max_height = unit(10, "cm")
)

h2 <- Heatmap(
    zscore_biased,
    col = color_zscore,
    show_row_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_names_max_height = unit(10, "cm")
)

pdf(biased_genes_small_pdf, height=10, width=20)
draw(h1 + h2)
dev.off()

svglite(biased_genes_small_svg, height=10, width=20)
draw(h1 + h2)
dev.off()

h3 <- Heatmap(
    mat
    , cluster_rows = FALSE
    , cluster_columns = FALSE
    , column_names_max_height = unit(10, "cm")
)

h4 <- Heatmap(
    zscore_biased
    , col = color_zscore
    , cluster_rows = FALSE
    , cluster_columns = FALSE
    , column_names_max_height = unit(10, "cm")
)

pdf(biased_genes_tall_pdf, height=250, width=30)
draw(h3 + h4)
dev.off()

svglite(biased_genes_tall_svg, height=250, width=30)
draw(h3 + h4)
dev.off()

h5 <- Heatmap(
    rp_zscore
    , col = color_zscore
    , cluster_rows = FALSE
    , cluster_columns = FALSE
    , column_names_max_height = unit(10, "cm")
)

pdf(translation_genes_pdf, height=30, width=20)
draw(h5)
dev.off()

svglite(translation_genes_svg, height=30, width=20)
draw(h5)
dev.off()

system(paste("mergepdf", biased_genes_small_pdf, biased_genes_tall_pdf, translation_genes_pdf, outfile))

