import pandas as pd
import networkx as nx

grn  = pd.read_table("../data/DSX_targets/NetRex-Female.txt")
print(grn)

G = nx.from_pandas_edgelist(grn,
        source = "TF",
        target = "Gene",
        edge_attr = "Rank",
        create_using = nx.DiGraph)

biased_genes = pd.read_table("head_biased_gene_groups.tsv")
print(biased_genes)

def shortest_path(symbol):
    if symbol in G.nodes():
      d = nx.shortest_path_length(G, source = "dsx", target = symbol)
    else:
      d = float('Inf')
    return d

biased_genes["shortest"] = biased_genes.symbol.apply(shortest_path)

columns = biased_genes.columns.tolist()

biased_genes = (
        biased_genes.sort_values("shortest")
        [columns[0:4] + ["shortest"] + columns[5:-1]]
)
biased_genes.to_csv("head_biased_gene_groups_dsx_dist.tsv", index=False, sep="\t")

