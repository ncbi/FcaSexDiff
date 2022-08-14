
#ruleorder: create_test_h5ad > combine_gonads > loom_to_h5ad

rule loom_to_h5ad:
  input:
    "imports/loom/{tissue}_{fcaver}.loom",
  output:
    h5ad_from_loom.path,
  script:
    "../scripts/loom_to_h5ad.py"


#rule create_test_h5ad:
#  input:
#    "resources/expr_test_relaxed.csv",
#  output:
#    "imports/expr_h5ad/expr_test_relaxed.h5ad",
#    "hoards/scpopcorn/scpopcorn_test_relaxed.txt"
#  script:
#    "../scripts/create_test_expr_h5ad.py"
#
#
#rule combine_gonads:
#  input:
#    "imports/expr_h5ad/expr_testis_relaxed.h5ad",
#    "imports/expr_h5ad/expr_ovary_relaxed.h5ad",
#  output:
#    "imports/expr_h5ad/expr_gonad_relaxed.h5ad",
#  run:
#    adata = ad.AnnData.concatenate(
#        ad.read_h5ad(input[0]),
#        ad.read_h5ad(input[1]),
#        index_unique = None,
#        join = "outer",
#    )
#    adata.write_h5ad(output[0])
#
#
rule filter_cells_by_annotations:
  input:
    h5ad_from_loom.path,
    "resources/sex_specific_annotations.csv",
  output:
    filtered_cells_csv.path,
  script:
    "../scripts/filter_cells_by_annotations.py"

ruleorder: normalize_expr > sctransform_expr

rule normalize_expr:
  input:
    h5ad_from_loom.path,
  output:
    normalized_h5ad.path
  script:
    "../scripts/normalize.py"

rule sctransform_expr:
  input:
    h5ad_from_loom.path,
  output:
    normalized_h5ad.path
  script:
    "../scripts/do_sctransform.R"


append_final_output(
  expand(
    normalized_h5ad.path,
    zip,
    tissue = samples["tissue"],
    fcaver = samples["fcaver"],
    expr = ["LogNorm"],
  )
)
#
#tasks = samples.explode("cellfilts")
#
#append_final_output(expand(
#  "scraps/filtered_lognorm_h5ad/cellfilter~{cellfilt}/filtered_lognorm_{tissue}_{fcaver}_{cellfilt}.h5ad",
#  zip,
#  tissue = tasks["tissue"],
#  fcaver = tasks["fcaver"],
#  cellfilt=tasks["cellfilts"],
#))
#
