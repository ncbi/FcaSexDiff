import pandas as pd
body_clusters = (
  pd.read_csv("clusters_body.csv")
  ["cluster"]
  .tolist()
)
print(body_clusters)
print(expand("body_{cluster}.rds", cluster = body_clusters))

rule all:
  #input: expand("GO_gsea_body/body_{cluster}.rds", cluster = body_clusters)
  input: expand("GO_enrich_body/body_{cluster}.rds", cluster = body_clusters)

rule do_go_gsea:
  input:
    "exports/sexdiff/cellfilter~NoSexspecArtef/resolution~L6.0/"
    "{tissue}/sexdiff_{tissue}"
    "_stringent_L6.0_NoSexspecArtef.h5ad"
  output:
    "GO_gsea_{tissue}/{tissue}_{cluster}.rds"
  script:
    "gsea_go.R"

rule do_go_enrich:
  input:
    "exports/sexdiff/cellfilter~NoSexspecArtef/resolution~L6.0/"
    "{tissue}/sexdiff_{tissue}"
    "_stringent_L6.0_NoSexspecArtef.h5ad"
  output:
    "GO_enrich_{tissue}/{tissue}_{cluster}.rds"
  script:
    "enrich_go.R"
