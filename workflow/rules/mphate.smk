rule run_mphate:
  input:
    "imports/expr_h5ad/expr_{tissue}_{fcaver}.h5ad",
  output:
    "scraps/mphate_pkl/mphate_{tissue}_{fcaver}.pkl",
    "scraps/mphate_clusters/mphate_clusters_{tissue}_{fcaver}.tsv",
  script:
    "../scripts/mphate/run_mphate.py"


rule prepare_sankey_data:
  input:
    "imports/expr_h5ad/expr_{tissue}_{fcaver}.h5ad",
    "scraps/mphate_clusters/mphate_clusters_{tissue}_{fcaver}.tsv",
  output:
    "scraps/mphate_sankey/mphate_sankey_{tissue}_{fcaver}.h5",
    "scraps/mphate_sankey/mphate_alluvial_{tissue}_{fcaver}.tsv",
  script:
    "../scripts/mphate/prepare_sankey_data.py"


rule draw_sankey:
  input:
    "scraps/mphate_sankey/mphate_sankey_{tissue}_{fcaver}.h5",
  output:
    "sankey_{tissue}_{fcaver}.html",
  script:
    "../scripts/mphate/draw_sankey.py"


rule draw_alluvial:
  input:
    "scraps/mphate_sankey/mphate_alluvial_{tissue}_{fcaver}.tsv",
  output:
    "alluvial_{tissue}_{fcaver}.pdf",
  script:
    "../scripts/mphate/draw_alluvial.R"


append_final_output(
    expand(
      [
          "scraps/mphate_clusters/mphate_clusters_{tissue}_{fcaver}.tsv",
          "scraps/mphate_sankey/mphate_sankey_{tissue}_{fcaver}.h5",
          "sankey_{tissue}_{fcaver}.html",
          "sankey_{tissue}_{fcaver}.pdf",
      ],
      zip,
      allow_missing = True,
      tissue=samples["tissue"],
      fcaver=samples["fcaver"],
    )
)

