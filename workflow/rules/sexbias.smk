
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

rule compare_normalization_cutoffs:
  input:
    sexdiff_h5ad.path
  output:
    compare_params_sexdiff.path
  script:
    "../scripts/compare_normalization_and_cutoffs.R"


#rule combine_biased_genes:
#  output:
#    "scraps/biased_genes/all_biased_genes_{resol}.tsv",
#  script:
#    "../scripts/combine_biased_genes.py"
#    
#
#
tasks = samples.explode("resols").explode("cellfilts")
print(tasks)

append_final_output(
    expand(
        expand(
            [sexdiff_excel.path],
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

