library(tidyverse)
library(grid)
library(gridExtra)

process <- function(tissue) {
    df = read.csv(paste0("ribosomal_proteins_", tissue, ".csv"))
    
    cor_vals = (
      df %>% group_by(batch)
      %>% summarize(pearson = round(cor(avg_rp, n_counts), 4))
      %>% mutate(label=paste("R:", pearson))
    )
    
    (
      ggplot(df, aes(x=log(avg_rp+1), y=log(n_counts+1), color=sex)) 
      + geom_point(size=0.1)
      + geom_label(data=cor_vals, aes(label=label), x=Inf, y=Inf, color="black", vjust=1, hjust=1)
      + facet_wrap(~batch, ncol=4)
    )
}
pdf("rp_body_head.pdf", width=12, height=9)

print(process("body"))
print(process("head"))
