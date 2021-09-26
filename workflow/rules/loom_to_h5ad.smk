rule loom_to_h5ad:
  input:
    "resources/loom/{tissue}_{fcaver}.loom",
  output:
    "resources/exp_h5ad/exp_{tissue}_{fcaver}.h5ad",
  log:
    "logs/loom_to_h5ad_{tissue}_{fcaver}.log",
  script:
    "../scripts/loom_to_h5ad.py"

downloads = samples.query("loom != ''")

append_final_output(expand(
  "resources/exp_h5ad/exp_{tissue}_{fcaver}.h5ad",
  zip,
  tissue=downloads["tissue"],
  fcaver=downloads["fcaver"],
))


