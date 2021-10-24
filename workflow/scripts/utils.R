
# get the var table in add anndata, including the varnames
# as a column in the table
get_var_anndata <- function(adata) {
  df = adata$var
  df[,attr(adata$var_names,"name")] = adata$var_names
  df
}
