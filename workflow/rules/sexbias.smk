
rule compute_sexbias:
  input:
    "hoards/int_h5ad/int_{tissue}_{fcaver}.h5ad",
    "resources/sex_specific_annotations.csv",
  output:
    "exports/sexdiff_h5ad_{cellfilt}/resol~{resol}/sexdiff_{tissue}_{fcaver}_{resol}_{cellfilt}.h5ad",
  script:
    "../scripts/compute_cluster_sexbias.py"


rule export_excel:
  input:
    "exports/sexdiff_h5ad_{cellfilt}/resol~{resol}/sexdiff_{tissue}_{fcaver}_{resol}_{cellfilt}.h5ad",
  output:
    "exports/sexdiff_xlsx_{cellfilt}/resol~{resol}/sexdiff_{tissue}_{fcaver}_{resol}_{cellfilt}.xlsx",
  script:
    "../scripts/export_sexbias_to_excel.py"


rule summarize_sexbias:
  input:
    expand(
      expand(
        "exports/sexdiff_h5ad_{cellfilt}/resol~{resol}/sexdiff_{tissue}_{fcaver}_{resol}_{cellfilt}.h5ad",
        zip,
        allow_missing=True,
        tissue=samples["tissue"],
        fcaver=samples["fcaver"],
      ),
      resol=["annotation"],
      cellfilt=["AllCells","NoSexspecArtef"]
    )
  output:
    "exports/extras/summary_sexbiased_genes_annotation.xlsx",
  script:
    "../scripts/summarize_sexbias.py"


rule plot_ma_plots:
  input:
    "exports/sexdiff_h5ad_{cellfilt}/resol~{resol}/sexdiff_{tissue}_{fcaver}_{resol}_{cellfilt}.h5ad",
  output:
    "exports/extras/maplots_{cellfilt}/resol~{resol}/maplots_{tissue}_{fcaver}_{resol}_{cellfilt}.pdf",
  script:
    "../scripts/plot_ma_plots.R"


rule plot_scatter_count_bias:
  input:
    "exports/sexdiff_h5ad_{cellfilt}/resol~{resol}/sexdiff_{tissue}_{fcaver}_{resol}_{cellfilt}.h5ad",
  output:
    "exports/extras/scatter_plots_{cellfilt}/resol~{resol}/scatter_count_bias_{tissue}_{fcaver}_{resol}_{cellfilt}.pdf",
    "exports/extras/scatter_plots_{cellfilt}/resol~{resol}/scatter_count_bias_{tissue}_{fcaver}_{resol}_{cellfilt}.csv",
  script:
    "../scripts/plot_scatter_count_bias.R"


rule draw_tsne_count_bias:
  input:
    "scraps/lognorm_h5ad/lognorm_{tissue}_{fcaver}.h5ad",
    "exports/sexdiff_h5ad_{cellfilt}/resol~{resol}/sexdiff_{tissue}_{fcaver}_{resol}_{cellfilt}.h5ad",
  output:
    "exports/extras/scatter_plots_{cellfilt}/resol~{resol}/tsne_count_bias_{tissue}_{fcaver}_{resol}_{cellfilt}.pdf",
  script:
    "../scripts/draw_tsne_count_bias.R"


rule get_single_gene_expr:
  input:
    "hoards/int_h5ad/int_{tissue}_{fcaver}.h5ad",
    "resources/sex_specific_annotations.csv",
  output:
    "exports/extras/violin_plots_{cellfilt}/resol~{resol}/{gene}_expr_{tissue}_{fcaver}_{resol}_{cellfilt}.tsv",
  script:
    "../scripts/get_single_gene_expr.py"



rule get_violin_gene_expr:
  input:
    "exports/extras/violin_plots_{cellfilt}/resol~{resol}/{gene}_expr_{tissue}_{fcaver}_{resol}_{cellfilt}.tsv",
  output:
    "exports/extras/violin_plots_{cellfilt}/resol~{resol}/{gene}_violin_{tissue}_{fcaver}_{resol}_{cellfilt}.pdf",
  script:
    "../scripts/plot_violin_gene_expr.R"



tasks = samples.explode("resols").explode("cellfilts")
specials = tasks.query("tissue in ['body', 'head', 'test']")

append_final_output(
    expand(
      [
        "exports/sexdiff_h5ad_{cellfilt}/resol~{resol}/sexdiff_{tissue}_{fcaver}_{resol}_{cellfilt}.h5ad",
        "exports/sexdiff_xlsx_{cellfilt}/resol~{resol}/sexdiff_{tissue}_{fcaver}_{resol}_{cellfilt}.xlsx",
        "exports/extras/maplots_{cellfilt}/resol~{resol}/maplots_{tissue}_{fcaver}_{resol}_{cellfilt}.pdf",
        "exports/extras/violin_plots_{cellfilt}/resol~{resol}/AkhR_violin_{tissue}_{fcaver}_{resol}_{cellfilt}.pdf",
      ],
      zip,
      allow_missing=True,
      tissue=tasks["tissue"],
      fcaver=tasks["fcaver"],
      resol=tasks["resols"],
      cellfilt=tasks["cellfilts"],
    ) + expand(
      [
        "exports/extras/scatter_plots_{cellfilt}/resol~{resol}/tsne_count_bias_{tissue}_{fcaver}_{resol}_{cellfilt}.pdf",
        "exports/extras/scatter_plots_{cellfilt}/resol~{resol}/scatter_count_bias_{tissue}_{fcaver}_{resol}_{cellfilt}.pdf",
      ],
      zip,
      allow_missing=True,
      tissue=specials["tissue"],
      fcaver=specials["fcaver"],
      resol=specials["resols"],
      cellfilt=specials["cellfilts"],
    ) + [
      "exports/extras/summary_sexbiased_genes_annotation.xlsx",
    ]
)

