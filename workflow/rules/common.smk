import glob

import pandas as pd
#from snakemake.utils import validate
#validate(config, schema="../schemas/config.schema.yaml")

samples = (
    pd.read_csv(config["samples"], sep="\t", comment = "#",
                dtype={"tissue": str, "fcaver":str})
    .assign(resols = lambda df: df.resols.apply(
        lambda x: ["annotation"] + ([] if pd.isna(x) else x.split(','))
    ))
    .assign(cellfilts = lambda df: df.cellfilts.apply(
        lambda x: ["NoSexspecArtef"] + ([] if pd.isna(x) else x.split(','))
    ))
    .set_index(["tissue", "fcaver"], drop=False)
    .sort_index()
)
#validate(samples, schema="../schemas/samples.schema.yaml")


final_output = []

def get_final_output():
  global final_output
  return final_output

def append_final_output(new_output):
  global final_output
  final_output += new_output

