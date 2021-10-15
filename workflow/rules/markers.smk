
rule find_markers:
  input:
    "resources/int_h5ad/int_{tissue}_{fcaver}.h5ad",
  output:
    "resources/mrk_h5ad/mrk_{tissue}_{fcaver}_{resol}.h5ad",
    "resources/markers/resol~{resol}/markers_{tissue}_{fcaver}_{resol}.h5ad",
  script:
    "../scripts/find_cluster_markers.py"


rule export_markers_to_excel:
  input:
    "resources/markers/resol~{resol}/markers_{tissue}_{fcaver}_{resol}.h5ad",
  output:
    "results/markers/resol~{resol}/markers_{tissue}_{fcaver}_{resol}.xlsx",
  script:
    "../scripts/export_markers_to_excel.py"

rule visualize_markers:
  input:
    "resources/markers/resol~{resol}/markers_{tissue}_{fcaver}_{resol}.h5ad",
  output:
    "results/markers/resol~{resol}/visualize_markers_{tissue}_{fcaver}_{resol}.ipynb",
  shell:
    "papermill workflow/notebooks/visualizing-marker-genes.ipynb "
    "{output} -p tissue {wildcards.tissue} -p infile {input} "


append_final_output(
  expand(
    expand(
      #"resources/markers/resol~{resol}/markers_{tissue}_{fcaver}_{resol}.h5ad",
      "results/markers/resol~{resol}/visualize_markers_{tissue}_{fcaver}_{resol}.ipynb",
      #"results/markers/resol~{resol}/markers_{tissue}_{fcaver}_{resol}.xlsx",
      zip,
      allow_missing=True,
      tissue=samples["tissue"],
      fcaver=samples["fcaver"],
    ),
    resol=['scpopcorn_cluster']
  )
)
get_final_output()

