
rule compute_sexbias:
  input:
    "scraps/lognorm_h5ad/lognorm_{tissue}_{fcaver}.h5ad",
    "resources/sex_specific_annotations.csv",
  output:
    "exports/sexdiff/cellfilter~{cellfilt}/resolution~{resol}/{tissue}/sexdiff_{tissue}_{fcaver}_{resol}_{cellfilt}.h5ad",
  script:
    "../scripts/compute_cluster_sexbias.py"


rule export_excel:
  input:
    "exports/sexdiff/cellfilter~{cellfilt}/resolution~{resol}/{tissue}/sexdiff_{tissue}_{fcaver}_{resol}_{cellfilt}.h5ad",
  output:
    "exports/sexdiff/cellfilter~{cellfilt}/resolution~{resol}/{tissue}/sexdiff_{tissue}_{fcaver}_{resol}_{cellfilt}.xlsx",
  script:
    "../scripts/export_sexbias_to_excel.py"


tasks = samples.explode("resols").explode("cellfilts")
specials = tasks.query("tissue in ['body', 'head', 'test']")

append_final_output(
    expand(
      [
        "exports/sexdiff/cellfilter~{cellfilt}/resolution~{resol}/{tissue}/sexdiff_{tissue}_{fcaver}_{resol}_{cellfilt}.h5ad",
        "exports/sexdiff/cellfilter~{cellfilt}/resolution~{resol}/{tissue}/sexdiff_{tissue}_{fcaver}_{resol}_{cellfilt}.xlsx",
      ],
      zip,
      allow_missing=True,
      tissue=tasks["tissue"],
      fcaver=tasks["fcaver"],
      resol=tasks["resols"],
      cellfilt=tasks["cellfilts"],
    )
)
