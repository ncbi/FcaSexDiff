
# get the var table in add anndata, including the varnames
# as a column in the table
get_var_anndata <- function(adata) {
  df = adata$var
  df[,attr(adata$var_names,"name")] = adata$var_names
  df
}


## convert integer clusters to strings
#format_cluster <- function(cluster_ids, prefix="", ndigit=-1) {
#  if (ndigit < 0) ndigit = ceiling(log10(max(cluster_ids)+1))
#  cluster_ids = formatC(cluster_ids, width = ndigit, format = "d", flag = "0")
#  paste0(prefix, "C", cluster_ids)
#}
