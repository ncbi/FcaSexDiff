rule get_loom:
  output:
    "resources/loom/{tissue}_{fcaver}.loom",
  log:
    "logs/get_loom_{tissue}_{fcaver}.log",
  run:
    loom = samples.loom[(wildcards.tissue, wildcards.fcaver)]
    loom = loom + "/download"
    shell("wget -O {output} {loom}")

downloads = samples.query("loom != ''")


append_final_output(expand(
  "resources/loom/{tissue}_{fcaver}.loom",
  zip,
  tissue=downloads["tissue"],
  fcaver=downloads["fcaver"],
))


