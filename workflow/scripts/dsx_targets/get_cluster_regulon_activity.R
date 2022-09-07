library(SCENIC)
library(AUCell)
library(tidyverse)
library(anndata)

source("workflow/scripts/utils.R")

expr_path <- snakemake@input[["adata"]]
auc_path <- snakemake@input[["auc"]]
cells_path <- snakemake@input[["cells"]]
resol <- snakemake@wildcards[["resol"]]

outfile <- snakemake@output[[1]]


cells <- read.csv(cells_path) %>% pull(CellID)

adata <- read_h5ad(expr_path)
adata <- adata[cells]

cell_meta <- (
    get_obs_anndata(adata)
    %>% rename(cluster := !!resol)
    %>% select(CellID, cluster)
)

head(cell_meta)

# This function is supposed to be included in the next version of AUCell
binarizeAUC <- function(auc, thresholds)
{
  thresholds <- thresholds[intersect(names(thresholds), rownames(auc))]
  regulonsCells <- setNames(lapply(names(thresholds),
                                   function(x) {
                                     trh <- thresholds[x]
                                     names(which(getAUC(auc)[x,]>trh))
                                   }),names(thresholds))
  regulonActivity <- reshape2::melt(regulonsCells)
  binaryRegulonActivity <- t(table(regulonActivity[,1], regulonActivity[,2]))
  class(binaryRegulonActivity) <- "matrix"

  return(binaryRegulonActivity)
}

obj <- readRDS(auc_path)

thresholds_df <- data.frame(
    regulon = names(obj$thresholds),
    threshold = unname(obj$thresholds)
)

head(thresholds_df)



df <- (
    as.data.frame(getAUC(obj$auc))
    %>% rownames_to_column("regulon")
    %>% gather("CellID", "activity", -regulon)
    %>% inner_join(cell_meta)
    %>% group_by(regulon, cluster)
    %>% summarise(avg_activity = mean(activity))
    %>% ungroup()
    %>% left_join(thresholds_df)
    %>% mutate(binarized_avg_activity = ifelse(avg_activity > threshold, 1, 0))
)

head(df)

binaryRegulonActivity <- binarizeAUC(obj$auc, obj$thresholds)
dim(binaryRegulonActivity)


binAct_subset <- binaryRegulonActivity[, which(colnames(binaryRegulonActivity) %in% cells)]
dim(binAct_subset)

binAct_subset[1:5,1:5]

df2 <- (
    as.data.frame(binAct_subset)
    %>% rownames_to_column("regulon")
    %>% gather("CellID", "activity", -regulon)
    %>% inner_join(cell_meta)
    %>% group_by(regulon, cluster)
    %>% summarise(avg_binarized_activity = mean(activity))
)

head(df2)

df <- left_join(df, df2)

head(df)

write.table(df, outfile, sep = "\t", row.names = FALSE)

