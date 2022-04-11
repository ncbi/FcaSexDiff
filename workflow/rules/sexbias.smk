from files.sexbias_files import *

rule find_markers:
  input:
    normalized_h5ad.path
  output:
    marker_info.path
  script:
    "../scripts/find_cluster_markers.py"

rule compute_sexbias_count:
  input:
    h5ad_from_loom.path,
    filtered_cells_csv.path,
  output:
    sexdiff_count.path
  script:
    "../scripts/compute_cluster_sexbias_count.py"


rule compute_sexbias_expr:
  input:
    normalized_h5ad.path,
    filtered_cells_csv.path,
  output:
    sexdiff_expr.path
  script:
    "../scripts/compute_cluster_sexbias_expression.py"

rule consolidate_sexbias:
  input:
    expr = normalized_h5ad.path,
    cells = filtered_cells_csv.path,
    mark = marker_info.path,
    cntbias = sexdiff_count.path,
    exprbias = sexdiff_expr.path,
  output:
    sexdiff_h5ad.path
  script:
    "../scripts/consolidate_cluster_sexbias.py"


rule export_excel:
  input:
    sexdiff_h5ad.path
  output:
    sexdiff_excel.path
  script:
    "../scripts/export_sexbias_to_excel.py"

#rule combine_biased_genes:
#  output:
#    "scraps/biased_genes/all_biased_genes_{resol}.tsv",
#  script:
#    "../scripts/combine_biased_genes.py"
#    
#
#
#tasks = samples.explode("resols").explode("cellfilts")
#
append_final_output(
#    expand(
#      [
#        "scraps/sexdiff/cellfilter~{cellfilt}/resolution~{resol}/{tissue}/sexdiff_count_{tissue}_{fcaver}_{resol}_{cellfilt}.h5",
#    "scraps/sexdiff/cellfilter~{cellfilt}/resolution~{resol}/{tissue}/sexdiff_expr_{tissue}_{fcaver}_{resol}_{cellfilt}.h5",
#        "exports/sexdiff/cellfilter~{cellfilt}/resolution~{resol}/{tissue}/sexdiff_{tissue}_{fcaver}_{resol}_{cellfilt}.h5ad",
#        "exports/sexdiff/cellfilter~{cellfilt}/resolution~{resol}/{tissue}/sexdiff_{tissue}_{fcaver}_{resol}_{cellfilt}.xlsx",
#      ],
#      zip,
#      allow_missing=True,
#      tissue=tasks["tissue"],
#      fcaver=tasks["fcaver"],
#      resol=tasks["resols"],
#      cellfilt=tasks["cellfilts"],
#    ) +
     [
       sexdiff_excel.path.format(
         tissue = "body", fcaver = "stringent",
         cellfilt = "NoSexspecArtef",  resol = "L4.0",
         expr = "LogNorm",
       ),
       sexdiff_excel.path.format(
         tissue = "body", fcaver = "stringent",
         cellfilt = "NoSexspecArtef",  resol = "L4.0",
         expr = "SCTrans",
       ),
       sexdiff_excel.path.format(
         tissue = "head", fcaver = "stringent",
         cellfilt = "NoSexspecArtef",  resol = "L4.0",
         expr = "LogNorm",
       ),
       sexdiff_excel.path.format(
         tissue = "head", fcaver = "stringent",
         cellfilt = "NoSexspecArtef",  resol = "L4.0",
         expr = "SCTrans",
       ),
#      "scraps/biased_genes/all_biased_genes_annotation.tsv",
#      "scraps/biased_genes/all_biased_genes_L4.0.tsv",
    ]
)

print(expand(
  sexdiff_excel.path,
  tissue = ["head", "body"],
  fcaver = "stringent",
  cellfilt = "NoSexspecArtef",
  resol = "L4.0",
  expr = ["LogNorm", "SCTrans"]
))

quit()

