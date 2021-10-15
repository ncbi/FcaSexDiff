rule get_loom:
  output:
    "imports/loom/{tissue}_{fcaver}.loom",
  log:
    "logs/get_loom_{tissue}_{fcaver}.log",
  run:
    loom = samples.loom[(wildcards.tissue, wildcards.fcaver)]
    loom = loom + "/download"
    shell("wget -O {output} {loom}")

# download looms for samples that have non-empty links
downloads = samples.query("~(loom.isna() | (loom == ''))", engine="python")

append_final_output(expand(
  "imports/loom/{tissue}_{fcaver}.loom",
  zip,
  tissue=downloads["tissue"],
  fcaver=downloads["fcaver"],
))


