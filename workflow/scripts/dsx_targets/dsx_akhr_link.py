import pandas as pd
import networkx as nx

grn  = pd.read_table("../data/DSX_targets/NetRex-Female.txt")
print(grn)

G = nx.from_pandas_edgelist(grn,
        source = "TF",
        target = "Gene",
        edge_attr = "Rank",
        create_using = nx.DiGraph)

dsx_akhr_paths = {}
for path in nx.all_simple_paths(G, source='dsx', target='AkhR', cutoff = 4):
  print(path)
  path_str = " -> ".join([f"{node}" for node in path])
  max_rank = max([
      grn.query(f"(TF == '{h}') & (Gene == '{t}')").iloc[0,:]["Rank"]
      for h,t in zip(path, path[1:])
  ])
  dsx_akhr_paths[path_str] = max_rank

dsx_akhr_paths = pd.DataFrame.from_dict(dsx_akhr_paths.items())
dsx_akhr_paths.columns = ["path", "max_rank"]
dsx_akhr_paths = dsx_akhr_paths.sort_values("max_rank")
print(dsx_akhr_paths)
dsx_akhr_paths.to_csv("dsx_akhr_paths_in_female.csv", index = False)

grn  = pd.read_table("../data/DSX_targets/NetRex-Male.txt")
print(grn)

G = nx.from_pandas_edgelist(grn,
        source = "TF",
        target = "Gene",
        edge_attr = "Rank",
        create_using = nx.DiGraph)

dsx_akhr_paths = {}
for path in nx.all_simple_paths(G, source='dsx', target='AkhR', cutoff = 4):
  print(path)
  path_str = " -> ".join([f"{node}" for node in path])
  max_rank = max([
      grn.query(f"(TF == '{h}') & (Gene == '{t}')").iloc[0,:]["Rank"]
      for h,t in zip(path, path[1:])
  ])
  dsx_akhr_paths[path_str] = max_rank

dsx_akhr_paths = pd.DataFrame.from_dict(dsx_akhr_paths.items())
dsx_akhr_paths.columns = ["path", "max_rank"]
dsx_akhr_paths = dsx_akhr_paths.sort_values("max_rank")
print(dsx_akhr_paths)
dsx_akhr_paths.to_csv("dsx_akhr_paths_in_male.csv", index = False)

