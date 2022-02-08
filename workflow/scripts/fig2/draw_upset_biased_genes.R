library(tidyverse)
library(ComplexUpset)
library(svglite)
library(ggplot2)
#source('../colors.R')

infile <- snakemake@input[[1]]
outfile <- snakemake@output[[1]]
tissue <- snakemake@wildcards[["tissue"]]

print(infile)
print(outfile)
print(tissue)



process <- function(use_tissue, infile, outfile) {

biased = (
  read.delim(infile)
  %>% filter(tissue == use_tissue)
  %>% mutate(bias = ifelse(bias > 0, "female",
                           ifelse(bias < 0, "male", "unbiased")))
  %>% filter(bias != "unbiased")
  %>% mutate(set = paste(cluster, bias, sep="|"))
  %>% mutate(val = TRUE)
  %>% select(symbol, set, val)
  %>% spread("set", "val", fill = FALSE)
)
print(head(biased))
print(dim(biased))

keep_chromosomes = c("2L", "2R", "3L", "3R", "4", "X", "Y")

chrnames <- (
    read.table("resources/flybase/gene_chr_FB2019_06_r6.31.tsv")
    %>% select(V2,V3)
    %>% mutate(symbol=V2, chr=V3)
    %>% select(symbol, chr)
    %>% mutate(chr = ifelse(chr %in% keep_chromosomes, chr, "other"))
)
print(head(chrnames))

chr_wt <- (
  chrnames 
  %>% count(chr)
  %>% mutate(wt = 100.0/n)
)
print(chr_wt)

clusters = setdiff(names(biased), c("symbol"))
print(clusters)

biased = (
    biased 
    %>% left_join(chrnames)
    %>% left_join(chr_wt)
    %>% filter(chr != "other")
    %>% mutate(chr = paste0(chr, " (", n, ")"))
)
print(head(biased))
print(dim(biased))


#outfile = paste0("upset_biased_genes_L6.0_", use_tissue, ".pdf")
#outfile2 = paste0("upset_biased_genes_L6.0_", use_tissue, ".tsv")

minsize = ifelse(use_tissue == "body", 2, 1)

pdf(outfile, width=70, height=70)

query_by_bias_type = function(data, groups, params_by_bias_type,
                              shared, ...) {
  intersections = unique(
    upset_data(data, groups, ...)$plot_intersections_subset
  )
  lapply(
    intersections,
    FUN = function(x) {
      members = ComplexUpset:::get_intersection_members(x)[[1]]
      tdf = (
        data.frame(member = members)
        %>% separate(member, c("name", "bias"), sep="[|]", remove=FALSE)
      )
      biases = unique(tdf$bias)
      bias_type = ifelse(length(biases) > 1, "mixed", biases[1])
      if (!(bias_type %in% names(params_by_bias_type))) {
        param = list(fill='gray')
      } else {
        param = params_by_bias_type[[bias_type]]
      }
      args = c(
        list(intersect=members),
        shared,
        param
      )
      do.call(upset_query, args)
    }
  )
}

cluster_meta <- (
    data.frame(set = clusters)
    %>% separate(set, c("name", "bias"), sep = "[|]", remove = FALSE)
)
print(cluster_meta)

show_genes <- function(intersections) {
  res <- data.frame(isect = intersections)
  print(head(res))

  lapply(
    1:nrow(res),
    FUN = function(i) {
      print("------")
      x = res$isect[i]
      members = ComplexUpset:::get_intersection_members(x)[[1]]
      others = setdiff(clusters, members)
      df = biased[rowSums(biased[,members, drop=FALSE]) == length(members) & 
              rowSums(biased[,others,drop=FALSE]) == 0,]
      print(nrow(df))
      if (nrow(df) > 5) {
        return ("")
      }
      paste(df$symbol, collapse=",")
    }
  )
}

mystat <- ggproto("mystat", Stat,
  compute_group = function(data, scales) {
    data$label = nrow(data)
    print(data)
    ret = unique(data)
    print(ret)
    ret
  },
  required_aes = c("x") #c("x", "y")
)

mystat2 <- ggproto("mystat2", Stat,
  compute_group = function(data, scales) {
    print(data)
    ret = unique(data)
    print(ret)
    ret
  },
  required_aes = c("x") #c("x", "y")
)

p <- upset(
  biased,
  clusters,
  name = 'Biased genes in clusters',
  #height_ratio = 9,
  width_ratio = 0.1,
  themes = upset_default_themes(text=element_text(size=30)),
  min_size = minsize,
  mode = "exclusive_intersection",
  encode_set = FALSE,
  base_annotations=list(
    'Intersection size'=intersection_size(
      text_colors=c(
        on_background='brown', on_bar='yellow'
      )
    )
    + geom_text(
                stat=mystat2,
                mapping=aes(label=show_genes(intersection)),
      #mapping = aes(label = ifelse(length(symbol) > 100, "",
      #                             paste(symbol, collapse=","))),
      #mapping = aes(label = show_genes(intersection)),
      angle = 90,
      hjust = 0,
      size = 5 
    )
    + ylab('Intersection size')
  ),
  annotations = list(
    'chromosome'=(
      ggplot(mapping=aes(fill=chr, weight=wt))
      + geom_bar(stat='count', position='stack')
      + geom_text(stat=mystat, y=-Inf, vjust=0)
      #+ geom_text(aes(label=inclusive_intersection_size))
      + scale_fill_brewer(type = "qual", palette = "Dark2") 
      + facet_wrap(~chr, ncol=1, , strip.position="left")
      + ylab('Percent of genes in chromosome')
      + theme(
        strip.placement = "outside"
      )
    )
  ),
  stripes=upset_stripes(
    mapping=aes(color=bias),
    colors=c(
      'female'='lightpink',
      'male'='lightblue'
    ),
    data=cluster_meta
  ),
  matrix=intersection_matrix(
    geom=geom_point(shape='circle filled', size=8)
  ),
  queries = query_by_bias_type(
    biased,
    clusters,
    min_size = minsize,
    params_by_bias_type=list(
      'female'=list(fill='lightpink'),
      'male'=list(fill='lightblue'),
      'mixed'=list(fill='purple')
    ),
    shared=list(
      only_components=c("intersections_matrix", "Intersection size"),
      color='black'
    )
  )
) + patchwork::plot_layout(heights=c(3,1,8))
print(p)

dev.off()

}

process(tissue, infile, outfile)

