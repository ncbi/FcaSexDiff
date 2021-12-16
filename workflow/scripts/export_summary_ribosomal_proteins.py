import anndata as ad
import pandas as pd
import sys

infiles = snakemake.input
outfile = snakemake.output[0]

print(infiles)
print(outfile)

writer = pd.ExcelWriter(outfile, engine='xlsxwriter')
for tissue,infile in infiles.items():
    print(infile)
    print(tissue)
    df = pd.read_hdf(infile)
    df.to_excel(writer, sheet_name=tissue)

writer.save()
