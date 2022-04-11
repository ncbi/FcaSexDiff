library(tidyverse)
library(anndata)
source("workflow/scripts/utils.R")

get_data <- function(infile) {
  if (file.exists(infile)) {
    return(read.csv(infile))
  }

  fat_body <- (
    read_h5ad("imports/expr_h5ad/expr_fat_body_stringent.h5ad")
    %>% get_obs_anndata()
    %>% select(CellID, annotation, n_counts, n_genes)
    %>% rename(nUMI = n_counts, nGene = n_genes)
    %>% mutate(tissue = "fat_body")
  )

  heart <- (
    read_h5ad("imports/expr_h5ad/expr_heart_stringent.h5ad")
    %>% get_obs_anndata()
    %>% select(CellID, annotation, n_counts, n_genes)
    %>% rename(nUMI = n_counts, nGene = n_genes)
    %>% mutate(tissue = "heart")
  )

  library(Seurat)
  testis <- readRDS("/net/intdev/flycellatlas/FcaTestisFollowUp/data/scratch/seurat_adult_sn_all.rds")
  testis <- (
    testis@meta.data
    %>% rownames_to_column("CellID")
    %>% select(CellID, annotation, nCount_RNA, nFeature_RNA)
    %>% rename(nUMI = nCount_RNA, nGene = nFeature_RNA)
    %>% mutate(tissue = "testis")
  )
  df <- (
    bind_rows(fat_body, heart, testis)
  )
  write.csv(df, infile, row.names = FALSE)
  df
}

df <- get_data("data_for_umi_box_plot.csv")

cluster_order <- (
  group_by(df, annotation, tissue)
  %>% summarize(avg_umi = median(nUMI))
  %>% ungroup()
  %>% group_by(annotation)
  %>% summarize(avg_umi = max(avg_umi))
  %>% ungroup()
  %>% arrange(avg_umi)
  %>% pull(annotation)
)
head(cluster_order)

df <- (
  df
  %>% gather("metric", "count", -CellID, -annotation, -tissue)
  %>% mutate(metric = factor(metric, levels = c("nUMI", "nGene"), ordered = TRUE))
  %>% mutate(annotation = factor(annotation, levels = cluster_order, ordered = TRUE))
)

p <- (
  ggplot(df, aes(x=annotation, y=count, color=tissue))
  + geom_boxplot(outlier.size = 0.5)
  + coord_flip()
  + facet_wrap(~metric, ncol=2, scales="free_x")
  + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
)

pdf("umi_box_plot.pdf", height=10, width=7)

p

dev.off()
