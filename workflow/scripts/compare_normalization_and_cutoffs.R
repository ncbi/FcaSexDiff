library(anndata)
library(tidyverse)


infile <- "exports/sexdiff_h5ad/cellfilt~NoSexspecArtef/expr~LogNorm/fcaver~stringent/resol~L4.0/tissue~body/sexdiff.h5ad"

infile <- snakemake@input[[1]]
outfile <- snakemake@output[[1]]

print(infile)
print(outfile)

ad <- read_h5ad(infile)

padj <- ad$layers[["padj"]]
log2fc <- ad$layers[["log2fc"]]
zscore <- ad$layers[["score"]]

get_bias_score <- function(pval_cutoff, log2fc_cutoff, mostly_frac) {
  male_biased <- (padj < pval_cutoff) & (log2fc < -log2fc_cutoff)
  male_biased <- male_biased[, colSums(male_biased) > 0]
  male_biased_genes <- rowSums(male_biased)
  mostly_male_biased_genes <- male_biased_genes[
    male_biased_genes > mostly_frac*ncol(male_biased)
  ]
  male_pearson <- cor(t(zscore[names(mostly_male_biased_genes), colnames(male_biased)]), method = "pearson")
  male_pearson <- male_pearson[lower.tri(male_pearson, diag=FALSE)]
  male_spearman <- cor(t(zscore[names(mostly_male_biased_genes), colnames(male_biased)]), method = "spearman")
  male_spearman <- male_spearman[lower.tri(male_spearman, diag=FALSE)]
  
  female_biased <- (padj < pval_cutoff) & (log2fc > log2fc_cutoff)
  female_biased <- female_biased[, colSums(female_biased) > 0]
  female_biased_genes <- rowSums(female_biased)
  mostly_female_biased_genes <- female_biased_genes[
    female_biased_genes > mostly_frac*ncol(female_biased)
  ]
  female_pearson <- cor(t(zscore[names(mostly_female_biased_genes), colnames(female_biased)]), method = "pearson")
  female_pearson <- female_pearson[lower.tri(female_pearson, diag=FALSE)]
  female_spearman <- cor(t(zscore[names(mostly_female_biased_genes), colnames(female_biased)]), method = "spearman")
  female_spearman <- female_spearman[lower.tri(female_spearman, diag=FALSE)]
  
  df <- data.frame(
    bias_mean_spearman_sexwise = mean(c(mean(female_spearman), mean(male_spearman)))
    , bias_mean_spearman_all = mean(c(female_spearman, male_spearman))
    , bias_mean_pearson_sexwise = mean(c(mean(female_pearson), mean(male_pearson)))
    , bias_mean_pearson_all = mean(c(female_pearson, male_pearson))
    , n_mostly_female = length(mostly_female_biased_genes)
    , n_mostly_male = length(mostly_male_biased_genes)
    , genes_mostly_female = paste(names(mostly_female_biased_genes), collapse = ",")
    , genes_mostly_male = paste(names(mostly_male_biased_genes), collapse = ",")
  )
}

df <- expand.grid(
  padj_cutoff = c(0.05, 0.01, 0.001)
  , fc_cutoff = c(1.25, 1.5, 2)
  , mostly_frac = c(0.67, 0.75, 0.8, 0.9)
)

df <- (
  df
  %>% rowwise()
  %>% mutate(newdf = list(get_bias_score(padj_cutoff, log2(fc_cutoff), mostly_frac)))
  %>% unnest(newdf)
)
df

write.table(df, outfile, row.names = FALSE, sep = "\t")

