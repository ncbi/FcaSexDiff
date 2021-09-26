import glob

import pandas as pd
#from snakemake.utils import validate
#validate(config, schema="../schemas/config.schema.yaml")

samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"tissue": str, "fcaver":str})
    .set_index(["tissue", "fcaver"], drop=False)
    .sort_index()
)

sex_specific = ['testis', 'ovary', 'male_reproductive_glands', 'all']

samples = samples.query("tissue not in @sex_specific")
samples = samples.query("tissue == 'leg'")

#validate(samples, schema="../schemas/samples.schema.yaml")

final_output = []

def get_final_output():
  global final_output
  return final_output

def append_final_output(new_output):
  global final_output
  final_output += new_output

