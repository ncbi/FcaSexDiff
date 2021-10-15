ruleorder: create_test_h5ad > combine_gonads > loom_to_h5ad

rule loom_to_h5ad:
  input:
    "imports/loom/{tissue}_{fcaver}.loom",
  output:
    "imports/expr_h5ad/expr_{tissue}_{fcaver}.h5ad",
  log:
    "logs/loom_to_h5ad_{tissue}_{fcaver}.log",
  script:
    "../scripts/loom_to_h5ad.py"


rule create_test_h5ad:
  input:
    "resources/expr_test_relaxed.csv",
  output:
    "imports/expr_h5ad/expr_test_relaxed.h5ad",
  script:
    "../scripts/create_test_expr_h5ad.py"


rule combine_gonads:
  input:
    "imports/expr_h5ad/expr_testis_relaxed.h5ad",
    "imports/expr_h5ad/expr_ovary_relaxed.h5ad",
  output:
    "imports/expr_h5ad/expr_gonad_relaxed.h5ad",
  run:
    adata = ad.AnnData.concatenate(
        ad.read_h5ad(input[0]),
        ad.read_h5ad(input[1]),
        index_unique = None,
        join = "outer",
    )
    adata.write_h5ad(output[0])


rule normalize_expr:
  input:
    "imports/expr_h5ad/expr_{tissue}_{fcaver}.h5ad",
  output:
    "scraps/lognorm_h5ad/lognorm_{tissue}_{fcaver}.h5ad",
  script:
    "../scripts/normalize.py"


append_final_output(expand(
  "scraps/lognorm_h5ad/lognorm_{tissue}_{fcaver}.h5ad",
  zip,
  tissue = samples["tissue"],
  fcaver = samples["fcaver"],
))


