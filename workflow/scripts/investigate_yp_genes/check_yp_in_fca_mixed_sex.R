library(anndata)
library(tidyverse)

process <- function(tissue) {
adata <- read_h5ad(paste0(
  "scraps/lognorm_h5ad/lognorm_",
  tissue,
  "_stringent.h5ad"
))
adata <- adata[adata$obs$sex == 'mix']

df <- (
  adata$obs
  %>% select(sex, L4.0, annotation)
  %>% mutate(
      rox1 = adata$X[, "lncRNA:roX1"],
      rox2 = adata$X[, "lncRNA:roX2"],
      yp1 = adata$X[, "Yp1"],
      yp2 = adata$X[, "Yp2"],
      yp3 = adata$X[, "Yp3"]
  )
  %>% mutate(
      predicted_sex = ifelse(rox1 + rox2 > 1, "male", "female")
  )
  %>% select(-rox1, -rox2)
)

print(head(df))

print(dim(df))

df2 <- df %>% filter(!grepl("fat body", annotation, fixed = TRUE))
print(dim(df2))

print(dim(df %>% filter(grepl("fat body", annotation, fixed=TRUE))))

print(df %>% count(annotation))

df2 <- df2 %>% mutate(L4.0 = "NON FATBODY")

df = bind_rows(df, df2)

annot = (
  df
  %>% count(L4.0, annotation)
  %>% group_by(L4.0)
  %>% summarize(major_annotation = annotation[which.max(n)])
)
print(annot)


cells <- (
  df
  %>% count(L4.0, predicted_sex)
  %>% mutate(label = paste0(substring(predicted_sex, 1, 1), n))
  %>% group_by(L4.0)
  %>% summarize(label = paste0(label, collapse=","), n = sum(n))
  %>% left_join(annot)
  %>% mutate(label = paste(L4.0, label, major_annotation, sep=","))
  %>% mutate(label = fct_reorder(label, n, .desc=TRUE))
)
print(cells)

df <- (
  left_join(df, cells)
  %>% gather("gene", "lognorm", -label, -sex, -predicted_sex, -L4.0, -annotation, -major_annotation, -n)
)
head(df)

pdf(
  paste0("boxplot_mixed_sex_rox_yp_", tissue, ".pdf"),
  height = 50, width = 10
)

p <- (
  ggplot(df, aes(x = gene, y = lognorm, color = predicted_sex))
  + geom_boxplot(outlier.shape=NA)
  #+ facet_wrap(~label, ncol=4)
  + coord_flip()
  + facet_grid(label~gene, scale="free_x", switch = "y")
  + theme(strip.text.y.left = element_text(angle = 0))
  #+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
)
print(p)

dev.off()
}

process("head")
process("body")
