library(tidyverse)
library(anndata)

expr_file <- snakemake@input[["expr"]]
scatter_data_file <- snakemake@input[["scatter"]]
tsne_data_file <- snakemake@output[[1]]
resol <- snakemake@wildcards[["resol"]]

print(expr_file)
print(scatter_data_file)

bias <- (
  read.csv(scatter_data_file, row.names=1)
  %>% rowwise()
  #%>% filter(max(cluster_count_male, cluster_count_female) >= 10)
  %>% filter(min(cluster_count_male, cluster_count_female) > 10)
)

get_cell_data <- function(tissue) {
  (
    read_h5ad(expr_file)$obs
    %>% filter(sex %in% c("female", "male"))
    %>% rename(cluster := !!resol)
    %>% select(tSNE1, tSNE2, cluster, sex)
  )
}

cells <- get_cell_data("body")

# since some cells were filtered, use only the cells
# in the clusters from the bias table
df <- left_join(bias, cells)
head(df)

saveRDS(df, tsne_data_file)
