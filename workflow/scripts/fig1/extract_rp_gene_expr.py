import anndata as ad
import pandas as pd
import sys

#infile = f"imports/expr_h5ad/expr_{tissue}_stringent.h5ad"


infile = snakemake.input["expr"]
rp_genes_file = snakemake.input["rpgenes"]
filtered_cells_csv = snakemake.input["cellfilt"]

outfile = snakemake.output[0]

print(infile)
print(outfile)
print(rp_genes_file)

rp_genes = (
    pd.read_table(rp_genes_file, header = None)
    [0].tolist()
)

keep_cells = (
    pd.read_table(filtered_cells_csv)
    .CellID.tolist()
)

adata = ad.read_h5ad(infile)
adata = adata[adata.obs_names.isin(keep_cells), adata.var_names.isin(rp_genes)]

adata.write_h5ad(outfile)

