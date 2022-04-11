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

rp_genes <- scan('resources/rp_genes_flybase_FBgg0000141.txt', what="", sep='\n')

get_bias <- function(tissue) {
  sexdiff_h5ad_file = glue(sexdiff_h5ad_fmt)
  ad <- read_h5ad(sexdiff_h5ad_file)
  ngene <- ad$uns[["stats"]]["NoSexspecArtef", "n_gene"]
  print(ngene)
  expr_bias = ad$X
  expr_bias_nonrp = expr_bias[!(row.names(expr_bias) %in% rp_genes),]
  non_rp <- data.frame(
    pct_nonrp_female = 100*Matrix::colSums(expr_bias_nonrp > 0)/nrow(expr_bias_nonrp),
    pct_nonrp_male = 100*Matrix::colSums(expr_bias_nonrp < 0)/nrow(expr_bias_nonrp),
    cluster = colnames(expr_bias)
  )
  bias <- (
    get_var_anndata(ad)
    %>% select(cluster, major_annotation, cluster_rep_fracs_mean_female,
               cluster_count_female, cluster_count_male,
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
    %>% left_join(non_rp)
    %>% mutate(nonrp_bias_type = ifelse(count_bias_type == "FemaleOnly", "FemaleOnly",
                                       ifelse(count_bias_type == "MaleOnly", "MaleOnly",
                                              ifelse(pct_nonrp_female > 2*pct_nonrp_male, "FemaleBiased",
                                                     ifelse(pct_nonrp_male > 2*pct_nonrp_female, "MaleBiased", "Unbiased")))))
    %>% mutate(nonrp_bias_label = ifelse((pct_nonrp_female > 2*pct_nonrp_male) | (pct_nonrp_male > 2*pct_nonrp_female),
                                       as.character(cluster_label), ''))
  )
  print(bias)
}

body_bias <- get_bias("body")
write.csv(body_bias, out_csv_body)

head_bias <- get_bias("head")
write.csv(head_bias, out_csv_head)

