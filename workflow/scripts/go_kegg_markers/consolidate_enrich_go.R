library(anndata)
library(tidyverse)
library(clusterProfiler)

process <- function(tissue) {

pdf(
  paste0("enrich_go_dotplots_marker_score_", tissue, ".pdf"),
  height = 7, width = 10
)

go <- data.frame()

clusters <- (
  read.delim(paste0("clusters_", tissue, ".csv"))
  %>% arrange(cluster)
  %>% pull(cluster)
)

print(clusters)

for (cls in clusters) {
  gse <-readRDS(paste0("GO_enrich_", tissue, "/", tissue, "_", cls, ".rds"))
  if (!is.na(gse)) {
  #go <- rbind(go, as.data.frame(gse) %>% mutate(cluster = cls))
  print(go)

  p <- (
    dotplot(gse) #, showCategory=20, split=".sign")
    #+ facet_grid(.~.sign)
    + ggtitle(paste("GO Enrich for", tissue, cls))
    #+ theme(plot.margin = unit(c(1,1,1,10), "cm"))
    + scale_y_discrete(labels = function(x) formatC(x, width = 30))
  )
  print(p)
} else {
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.5, y = 0.5, paste(
    "Either there was some error in running Enrich GO on this cluster,\n",
    cls,
    "or there were no significant terms."
    )
  , cex = 1.6, col = "black"
  )
}
}
dev.off()

write.csv(go, paste0("GO_", tissue, ".csv"))


}

#process("head")
process("body")

#import anndata as ad
#
#adata = (
#    ad.read_h5ad("../data/imports/expr_h5ad/expr_body_stringent.h5ad")
#    .obs[["sex", "L6.0", "annotation", "Annotation"]]
#)
#adata = adata.loc[adata["L6.0"] == "L6.0C60"]
#print(adata)
#adata = adata.query("sex != 'mix'")
#print(adata)
#quit()
#
#adata = (
#    ad.read_h5ad("scraps/lognorm_h5ad/lognorm_body_stringent.h5ad")
#    .obs[["sex", "L6.0", "annotation"]]
#)
#adata = adata.query("sex != 'mix'")
#adata = adata.loc[adata["L6.0"] == "L6.0C60"]
#print(adata)
#print(adata.groupby("annotation").count())
#
#adata = (
#    ad.read_h5ad("scraps/filtered_lognorm_h5ad/cellfilter~NoSexspecArtef/filtered_lognorm_body_stringent_NoSexspecArtef.h5ad")
#    .obs[["sex", "L6.0", "annotation", "S_annotation"]]
#)
#adata = adata.query("sex != 'mix'")
#adata = adata.loc[adata["L6.0"] == "L6.0C60"]
#print(adata)
#print(adata.columns)
