library(tidyverse)
library(ggthemes)
library(ggrepel)
library(anndata)

source("workflow/scripts/utils.R")


infile <- snakemake@input[[1]]
outfile <- snakemake@output[[1]]
csvfile <- snakemake@output[[2]]
resol <- snakemake@wildcards[['resol']]

print(infile)
print(outfile)

#infile = "exports/sexdiff_h5ad/resol~annotation/sexdiff_body_stringent_annotation.h5ad"

bias <- (
  get_var_anndata(read_h5ad(infile))
  %>% select(cluster, major_annotation, cluster_rep_fracs_mean_female,
             cluster_rep_fracs_mean_male, cluster_rep_fracs_sd_female,
             cluster_rep_fracs_sd_male, count_bias_padj,
             log2_count_bias, count_bias_type)
  %>% mutate(pct_female = 100*cluster_rep_fracs_mean_female)
  %>% mutate(pct_male = 100*cluster_rep_fracs_mean_male)
  %>% mutate(sd_female = 100*cluster_rep_fracs_sd_female)
  %>% mutate(sd_male = 100*cluster_rep_fracs_sd_male)
  %>% mutate(xmin = pct_male - sd_male)
  %>% mutate(xmax = pct_male + sd_male)
  %>% mutate(ymin = pct_female - sd_female)
  %>% mutate(ymax = pct_female + sd_female)
  %>% mutate(label = ifelse(((pct_female > 2) | (count_bias_type == "Female")) | 
                            ((pct_male > 2) | (count_bias_type == "Male")),
                            as.character(cluster), ''))
)

if (resol != "annotation") {
  bias <- mutate(bias, label = substring(label, nchar(!!resol)+2))
}
head(bias)

(
  select(bias, cluster, major_annotation, count_bias_type, count_bias_padj,
         pct_male, pct_female, sd_male, sd_female)
  %>% write.csv(csvfile, row.names=F)
)

base <- (
    ggplot(bias, aes(y=pct_female, x=pct_male, color=count_bias_type))
    + geom_abline(intercept = 0, slope = 2, linetype='dotted')
    + geom_abline(intercept = 0, slope = 0.5, linetype='dotted')
    + geom_errorbar(aes(ymin=ymin, ymax=ymax), width=0.5, size=0.5, color="gray")
    + geom_errorbarh(aes(xmin=xmin, xmax=xmax), height=0.5, size=0.5, color="gray")
    + geom_text_repel(aes(label=label), min.segment.length = 0, max.overlaps=Inf)
    + geom_point(size=1)
    + theme_few()
    + scale_color_manual(values=c("Female"="red", "Male"="blue", "Unbiased"="black"))
)

main <- (
    base
    + theme(legend.position = "bottom")
    + labs(
        y='% female cells (replicate average)',
        x='% male cells (replicate average)'
    )
)


inset <- (
  base
  + lims(x=c(0,2), y=c(0,2))
  + labs(x="", y="")
  + theme(
    legend.position = "none",
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )
)

vp <- grid::viewport(width = 0.4, height = 0.4, x = 0.05, y = 1, just=c("left", "top"))

pdf(outfile, height=11, width=10)
print(main)
print(inset, vp = vp)

