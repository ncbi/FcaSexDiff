
rule compute_sexbias:
  input:
    "hoards/int_h5ad/int_{tissue}_{fcaver}.h5ad",
  output:
    "exports/sexdiff_h5ad/resol~{resol}/sexdiff_{tissue}_{fcaver}_{resol}.h5ad",
  script:
    "../scripts/compute_cluster_sexbias.py"


rule export_excel:
  input:
    "exports/sexdiff_h5ad/resol~{resol}/sexdiff_{tissue}_{fcaver}_{resol}.h5ad",
  output:
    "exports/sexdiff_xlsx/resol~{resol}/sexdiff_{tissue}_{fcaver}_{resol}.xlsx",
  script:
    "../scripts/export_sexbias_to_excel.py"


rule summarize_sexbias:
  input:
    expand(
      expand(
        "exports/sexdiff_h5ad/resol~{resol}/sexdiff_{tissue}_{fcaver}_{resol}.h5ad",
        zip,
        allow_missing=True,
        tissue=samples["tissue"],
        fcaver=samples["fcaver"],
      ),
      resol=["annotation"]
    )
  output:
    "exports/extras/summary_sexbiased_genes_annotation.xlsx",
  script:
    "../scripts/summarize_sexbias.py"


rule plot_ma_plots:
  input:
    "exports/sexdiff_h5ad/resol~{resol}/sexdiff_{tissue}_{fcaver}_{resol}.h5ad",
  output:
    "exports/extras/maplots/resol~{resol}/maplots_{tissue}_{fcaver}_{resol}.pdf",
  script:
    "../scripts/plot_ma_plots.R"


rule plot_scatter_count_bias:
  input:
    "exports/sexdiff_h5ad/resol~{resol}/sexdiff_{tissue}_{fcaver}_{resol}.h5ad",
  output:
    "exports/extras/scatter_plots/resol~{resol}/scatter_count_bias_{tissue}_{fcaver}_{resol}.pdf",
  script:
    "../scripts/plot_scatter_count_bias.R"


rule draw_tsne_count_bias:
  input:
    "scraps/lognorm_h5ad/lognorm_{tissue}_{fcaver}.h5ad",
    "exports/sexdiff_h5ad/resol~{resol}/sexdiff_{tissue}_{fcaver}_{resol}.h5ad",
  output:
    "exports/extras/scatter_plots/resol~{resol}/tsne_count_bias_{tissue}_{fcaver}_{resol}.pdf",
  script:
    "../scripts/draw_tsne_count_bias.R"


tasks = samples.explode("resols")
specials = tasks.query("tissue in ['body', 'head']")

append_final_output(
    expand(
      [
        "exports/sexdiff_h5ad/resol~{resol}/sexdiff_{tissue}_{fcaver}_{resol}.h5ad",
        "exports/sexdiff_xlsx/resol~{resol}/sexdiff_{tissue}_{fcaver}_{resol}.xlsx",
        "exports/extras/maplots/resol~{resol}/maplots_{tissue}_{fcaver}_{resol}.pdf",
      ],
      zip,
      allow_missing=True,
      tissue=tasks["tissue"],
      fcaver=tasks["fcaver"],
      resol=tasks["resols"],
    ) + expand(
      [
        "exports/extras/scatter_plots/resol~{resol}/tsne_count_bias_{tissue}_{fcaver}_{resol}.pdf",
        "exports/extras/scatter_plots/resol~{resol}/scatter_count_bias_{tissue}_{fcaver}_{resol}.pdf",
      ],
      zip,
      allow_missing=True,
      tissue=specials["tissue"],
      fcaver=specials["fcaver"],
      resol=specials["resols"],
    ) + [
      "exports/extras/summary_sexbiased_genes_annotation.xlsx",
    ]
)

