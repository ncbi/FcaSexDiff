import glob

import pandas as pd
#from snakemake.utils import validate
#validate(config, schema="../schemas/config.schema.yaml")

samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"tissue": str, "fcaver":str})
    .assign(resols = lambda df: df.resols.apply(
        lambda x: ["annotation"] + ([] if pd.isna(x) else x.split())
    ))
    .set_index(["tissue", "fcaver"], drop=False)
    .sort_index()
)
#validate(samples, schema="../schemas/samples.schema.yaml")

exclude = ['testis', 'ovary', 'male_reproductive_glands',
           'gonad',
           'all',
#           'test',
#           'antenna', 'head',
#           'body',
           ]

samples = samples.query("tissue not in @exclude")

##samples = samples.query("tissue in ['leg', 'malpighian_tubule']")
##samples = samples.query("tissue in ['leg']")
##samples = samples.query("tissue in ['testis', 'ovary']")
##samples = samples.query("tissue in ['gonad']")
#samples = samples.query("tissue in ['malpighian_tubule']")
#samples = samples.query("tissue in ['test']")
#samples = samples.query("tissue in ['body', 'malpighian_tubule']")
#samples = samples.query("tissue in ['body', 'head']")
#samples = samples.query("tissue in ['body']")


final_output = []

def get_final_output():
  global final_output
  return final_output

def append_final_output(new_output):
  global final_output
  final_output += new_output

