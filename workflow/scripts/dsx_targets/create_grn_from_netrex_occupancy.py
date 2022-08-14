import pandas as pd
import networkx as nx

grn_female = (
    pd.read_table("../data/DSX_targets/NetRex-Female.txt")
    [["TF", "Gene"]]
)
print(grn_female)

grn_male = (
    pd.read_table("../data/DSX_targets/NetRex-Male.txt")
    [["TF", "Gene"]]
)
print(grn_male)

dsx_targets = (
    pd.read_table("DSX_targets.tsv")[["Gene_Name"]]
    .rename(columns={"Gene_Name": "Gene"})
    .assign(TF = "dsx")
)
print(dsx_targets)
print(dsx_targets.query("Gene == 'dsx' & TF == 'dsx'"))

grn = pd.concat([grn_female, grn_male, dsx_targets], axis =
        0).drop_duplicates()
print(grn)

grn.to_csv("resources/GRN.tsv", sep="\t", index=False)

