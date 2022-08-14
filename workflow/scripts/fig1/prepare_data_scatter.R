library(tidyverse)
library(anndata)
library(glue)
library(colorspace)
source("workflow/scripts/utils.R")

sexdiff_file <- snakemake@input[["sexdiff"]]
rp_expr_file <- snakemake@input[["rp_expr"]]
scatter_data_file <- snakemake@output[[1]]
tissue <- snakemake@wildcards[["tissue"]]
resol <- snakemake@wildcards[["resol"]]

print(rp_expr_file)
print(sexdiff_file)
print(tissue)
print(scatter_data_file)

rpexpr <- read_h5ad(rp_expr_file)
rp_genes <- rpexpr$var_names
head(rp_genes)

rp_avg <- (
    get_obs_anndata(rpexpr)
    %>% mutate(rp_avg = Matrix::rowMeans(rpexpr$X))
    %>% rename(cluster := !!resol)
    %>% select(sex, cluster, rp_avg)
    %>% group_by(cluster, sex)
    %>% summarize(rp_avg = mean(rp_avg))
    %>% mutate(sex = paste0("rp_avg_", sex))
    %>% spread(sex, rp_avg)
)
head(rp_avg)

get_pct <- function(mat) {
    100 * Matrix::colSums(mat) / nrow(mat)
}

sexdiff <- read_h5ad(sexdiff_file)
ngene <- sexdiff$uns[["stats"]]["NoSexspecArtef", "n_gene"]
print(ngene)

expr_bias <- sexdiff$X

expr_bias_rp <- expr_bias[(row.names(expr_bias) %in% rp_genes), ]
rp_bias <- data.frame(
    pct_rp_female = get_pct(expr_bias_rp > 0),
    pct_rp_male = get_pct(expr_bias_rp < 0),
    cluster = colnames(expr_bias)
)
head(rp_bias)

expr_bias_nonrp <- expr_bias[!(row.names(expr_bias) %in% rp_genes), ]
non_rp_bias <- data.frame(
    pct_nonrp_female = get_pct(expr_bias_nonrp > 0),
    pct_nonrp_male = get_pct(expr_bias_nonrp < 0),
    cluster = colnames(expr_bias)
)
head(non_rp_bias)

col_fun_female <- circlize::colorRamp2(c(0, 1), c("gray75", "red"))
col_fun_male <- circlize::colorRamp2(c(0, 1), c("gray75", "blue"))

blend_color <- function(pct_female, pct_male) {
    col_female <- col_fun_female(pct_female)
    col_male <- col_fun_male(pct_male)
    rgb_female <- t(col2rgb(col_female))
    rgb_male <- t(col2rgb(col_male))
    col <- mixcolor(0.5, sRGB(rgb_male/255), sRGB(rgb_female/255))
    rgb(col@coords)
}

fix_rank <- function(pct_female, pct_male) {
    pct_max <- pmax(pct_female, pct_male)
    order(pct_max)
}

show_label <- function(pct_female, pct_male, label) {
    female_lb <- quantile(pct_female, 0.10)
    male_lb <- quantile(pct_male, 0.10)
    female_ub <- quantile(pct_female, 0.75)
    male_ub <- quantile(pct_male, 0.75)
    interesting <- (
        (pct_female > 1.25 * pct_male)
        | (pct_male > 1.25 * pct_female)
        | (pct_female > female_ub)
        | (pct_male > male_ub)
    )
    keep <- (
        interesting
        & (pct_female > female_lb)
        & (pct_male > male_lb)
    )
    ifelse(keep, label, "")
}

shorten_label <- function(label) {
    if (resol == "annotation") {
        label
    } else {
        substring(label, nchar(resol) + 2)
    }
}

bias <- (
    get_var_anndata(sexdiff)
    %>% select(cluster, major_annotation, cluster_rep_fracs_mean_female,
               cluster_count_female, cluster_count_male,
               cluster_rep_fracs_mean_male, cluster_rep_fracs_sd_female,
               cluster_rep_fracs_sd_male, count_bias_padj,
               log2_count_bias, count_bias_type, female_gene, male_gene)
    %>% mutate(cluster_label = shorten_label(cluster))
    %>% mutate(pct_cells_female = 100*cluster_rep_fracs_mean_female)
    %>% mutate(pct_cells_male = 100*cluster_rep_fracs_mean_male)
    %>% mutate(sd_cells_female = 100*cluster_rep_fracs_sd_female)
    %>% mutate(sd_cells_male = 100*cluster_rep_fracs_sd_male)
    %>% mutate(xmin = pct_cells_male - sd_cells_male)
    %>% mutate(xmax = pct_cells_male + sd_cells_male)
    %>% mutate(ymin = pct_cells_female - sd_cells_female)
    %>% mutate(ymax = pct_cells_female + sd_cells_female)
    %>% mutate(pct_genes_female = 100*female_gene/ngene)
    %>% mutate(pct_genes_male = 100*male_gene/ngene)
    %>% left_join(rp_avg)
    %>% replace_na(list(rp_avg_female=0, rp_avg_male=0))
    %>% left_join(rp_bias)
    %>% left_join(non_rp_bias)
    %>% mutate(count_color = blend_color(pct_cells_female, pct_cells_male))
    %>% mutate(expr_color = blend_color(pct_genes_female, pct_genes_male))
    %>% mutate(rp_color = blend_color(pct_rp_female, pct_rp_male))
    %>% mutate(nonrp_color = blend_color(pct_nonrp_female, pct_nonrp_male))
    %>% mutate(count_label = show_label(pct_cells_female, pct_cells_male, cluster_label))
    %>% mutate(expr_label = show_label(pct_genes_female, pct_genes_male, cluster_label))
    %>% mutate(rp_label = show_label(pct_rp_female, pct_rp_male, cluster_label))
    %>% mutate(nonrp_label = show_label(pct_nonrp_female, pct_nonrp_male, cluster_label))
    %>% mutate(count_rank = fix_rank(pct_cells_female, pct_cells_male))
    %>% mutate(expr_rank = fix_rank(pct_genes_female, pct_genes_male))
    %>% mutate(rp_rank = fix_rank(pct_rp_female, pct_rp_male))
    %>% mutate(nonrp_rank = fix_rank(pct_nonrp_female, pct_nonrp_male))
)
head(bias)

write.csv(bias, scatter_data_file)

