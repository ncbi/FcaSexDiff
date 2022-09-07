library(SCopeLoomR)
library(SCENIC)
library(tidyverse)

scenicLoomPath <- snakemake@input[["loom"]]
scenicLoomPath

targets_path <- snakemake@output[["targets"]]
auc_path <- snakemake@output[["auc"]]

loom <- open_loom(scenicLoomPath)

# Read information from loom file:
#exprMat <- get_dgem(loom)
#exprMat_log <- log2(exprMat+1) # Better if it is logged/normalized
regulons_incidMat <- get_regulons(loom, column.attr.name="MotifRegulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom)
regulonAucThresholds <- get_regulon_thresholds(loom)
#embeddings <- get_embeddings(loom)
#cellClusters <- get_clusterings(loom)

close_loom(loom)



# 2022-09-06: It seems there is a bug in the get_regulon_thresholds
# function for which the names and values of the thresholds are swapped
regulonAucThresholds <- setNames(
    as.numeric(names(regulonAucThresholds)),
    unname(regulonAucThresholds)
)
saveRDS(list(auc = regulonAUC, thresholds = regulonAucThresholds), auc_path)


regulons_df <- bind_rows(lapply(names(regulons), function(tf) {
    data.frame(x = tf, y = regulons[[tf]])
}))
write.table(regulons_df, targets_path, sep = "\t")

