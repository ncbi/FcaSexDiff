rule find_markers:
  input:
    "scraps/lognorm_h5ad/lognorm_{tissue}_{fcaver}.h5ad",
  output:
    "scraps/markers_h5/resolution~{resol}/markers_{tissue}_{fcaver}_{resol}.h5",
  script:
    "../scripts/find_cluster_markers.py"

append_final_output(
    [
        "scraps/markers_h5/resolution~annotation/markers_body_stringent_annotation.h5",
    ]
)

rule compute_sexbias_count:
  input:
    "scraps/filtered_lognorm_h5ad/cellfilter~{cellfilt}/filtered_lognorm_{tissue}_{fcaver}_{cellfilt}.h5ad"
  output:
    "scraps/sexdiff/cellfilter~{cellfilt}/resolution~{resol}/{tissue}/sexdiff_count_{tissue}_{fcaver}_{resol}_{cellfilt}.h5",
  script:
    "../scripts/compute_cluster_sexbias_count.py"


rule compute_sexbias_expr:
  input:
    "scraps/filtered_lognorm_h5ad/cellfilter~{cellfilt}/filtered_lognorm_{tissue}_{fcaver}_{cellfilt}.h5ad"
  output:
    "scraps/sexdiff/cellfilter~{cellfilt}/resolution~{resol}/{tissue}/sexdiff_expr_{tissue}_{fcaver}_{resol}_{cellfilt}.h5",
  script:
    "../scripts/compute_cluster_sexbias_expression.py"

rule consolidate_sexbias:
  input:
    expr = "scraps/filtered_lognorm_h5ad/cellfilter~{cellfilt}/filtered_lognorm_{tissue}_{fcaver}_{cellfilt}.h5ad",
    mark = "scraps/markers_h5/resolution~{resol}/markers_{tissue}_{fcaver}_{resol}.h5",
    cntbias = "scraps/sexdiff/cellfilter~{cellfilt}/resolution~{resol}/{tissue}/sexdiff_count_{tissue}_{fcaver}_{resol}_{cellfilt}.h5",
    exprbias = "scraps/sexdiff/cellfilter~{cellfilt}/resolution~{resol}/{tissue}/sexdiff_expr_{tissue}_{fcaver}_{resol}_{cellfilt}.h5",
  output:
    "exports/sexdiff/cellfilter~{cellfilt}/resolution~{resol}/{tissue}/sexdiff_{tissue}_{fcaver}_{resol}_{cellfilt}.h5ad",
  script:
    "../scripts/consolidate_cluster_sexbias.py"


rule export_excel:
  input:
    "exports/sexdiff/cellfilter~{cellfilt}/resolution~{resol}/{tissue}/sexdiff_{tissue}_{fcaver}_{resol}_{cellfilt}.h5ad",
  output:
    "exports/sexdiff/cellfilter~{cellfilt}/resolution~{resol}/{tissue}/sexdiff_{tissue}_{fcaver}_{resol}_{cellfilt}.xlsx",
  script:
    "../scripts/export_sexbias_to_excel.py"


tasks = samples.explode("resols").explode("cellfilts")

append_final_output(
    expand(
      [
        "scraps/sexdiff/cellfilter~{cellfilt}/resolution~{resol}/{tissue}/sexdiff_count_{tissue}_{fcaver}_{resol}_{cellfilt}.h5",
    "scraps/sexdiff/cellfilter~{cellfilt}/resolution~{resol}/{tissue}/sexdiff_expr_{tissue}_{fcaver}_{resol}_{cellfilt}.h5",
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
