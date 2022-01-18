# Load package
library(networkD3)
library(tidyverse)

infile = "links.tsv"
outfile = "sankey.html"

plot_sankey <- function(infile, outfile) {
  links = read.delim(infile, header=T, sep='\t')

  print(head(links))


  # From these flows we need to create a node data frame: it lists every entities involved in the flow
  nodes <- data.frame(
                      name=c(as.character(links$src), 
                             as.character(links$tgt)) %>% unique()
  )
  print(nodes)

  # With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
  links$IDsource <- match(links$src, nodes$name)-1 
  links$IDtarget <- match(links$tgt, nodes$name)-1



  # Make the Network
  p <- sankeyNetwork(Links = links, Nodes = nodes,
                     Source = "IDsource", Target = "IDtarget",
                     Value = "count", NodeID = "name", 
                     units = "cells", fontSize = 20, nodeWidth = 30,
                     iterations=500, height=3500, width=4000)
  p


  # save the widget
  library(htmlwidgets)
  saveWidget(p, file=outfile)

  #library(webshot)
  #webshot("sankey.html", "sankey.pdf")
}

plot_sankey(infile, outfile)
plot_sankey("links_filtered.tsv", "sankey_filtered.html")
