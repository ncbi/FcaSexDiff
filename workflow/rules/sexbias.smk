
rule compute_sexbias:
  input:
    "resources/int_h5ad/int_{tissue}_{fcaver}.h5ad",
  output:
    "resources/sexbias/resol~{resol}/sexbias_{tissue}_{fcaver}_{resol}.h5ad",
  script:
    "../scripts/compute_cluster_sexbias.py"


rule export_excel:
  input:
    "resources/sexbias/resol~{resol}/sexbias_{tissue}_{fcaver}_{resol}.h5ad",
  output:
    "results/sexbias/resol~{resol}/sexbias_{tissue}_{fcaver}_{resol}.xlsx",
  script:
    "../scripts/export_sexbias_to_excel.py"


append_final_output(
  expand(
    expand(
      "results/sexbias/resol~{resol}/sexbias_{tissue}_{fcaver}_{resol}.xlsx",
      zip,
      allow_missing=True,
      tissue=samples["tissue"],
      fcaver=samples["fcaver"],
    ),
    resol=['annotation']
  )
)
get_final_output()

