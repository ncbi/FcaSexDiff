
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

rule plot_ma_plots:
  input:
    "exports/sexdiff_h5ad/resol~{resol}/sexdiff_{tissue}_{fcaver}_{resol}.h5ad",
  output:
    "exports/extras/maplots/resol~{resol}/maplots_{tissue}_{fcaver}_{resol}.pdf",
  script:
    "../scripts/plot_ma_plots.R"

from itertools import cycle

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


append_final_output(
  expand(
    expand(
      [
        "exports/sexdiff_h5ad/resol~{resol}/sexdiff_{tissue}_{fcaver}_{resol}.h5ad",
        "exports/sexdiff_xlsx/resol~{resol}/sexdiff_{tissue}_{fcaver}_{resol}.xlsx",
        "exports/extras/maplots/resol~{resol}/maplots_{tissue}_{fcaver}_{resol}.pdf",
      ],
      zip,
      allow_missing=True,
      tissue=samples["tissue"],
      fcaver=samples["fcaver"],
    ),
    resol=['annotation']
  ) + [
    "exports/extras/summary_sexbiased_genes_annotation.xlsx",
  ]
)

