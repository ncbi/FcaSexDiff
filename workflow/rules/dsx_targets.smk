rule get_biased_genes:
  input:
    sexdiff_h5ad.path
  output:
    biased_genes.path,
  script:
    "../scripts/get_biased_genes_filtered_clusters.R"


rule group_biased_genes:
  input:
    biased_genes.path,
  output:
    biased_gene_groups.path,
  script:
    "../scripts/group_biased_genes.R"


rule compute_enrichment_dsx_targets_tfs:
  input:
    expr = normalized_h5ad.path,
    bias_grps = biased_gene_groups.path,
  output:
    dsx_targets_tfs_enrich.path,
  script:
    "../scripts/dsx_targets/compute_enrichment_dsxtargets_tfs.py"

tasks = samples.explode("resols").explode("cellfilts")

append_final_output(
    expand(
    expand(
      [
        dsx_targets_tfs_enrich.path
      ],
      zip,
      allow_missing=True,
      tissue=tasks["tissue"],
      fcaver=tasks["fcaver"],
      resol=tasks["resols"],
      cellfilt=tasks["cellfilts"],
    ),
      expr = "LogNorm",
    )
)

