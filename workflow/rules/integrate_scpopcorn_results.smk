import sys
sys.path.append('workflow/scripts')

from scpopcorn import MergeSingleCell
from scpopcorn import SingleCellData
import anndata as ad
import pickle

rule augment_h5ad_obs:
  input:
    "resources/exp_h5ad/exp_{tissue}_{fcaver}.h5ad",
    "resources/scpopcorn/scpopcorn_{tissue}_{fcaver}.txt",
  output:
    "resources/int_h5ad/int_{tissue}_{fcaver}.h5ad",
  run:
    adata = ad.read_h5ad(input[0])
    pcres = (
      pd.read_table(input[1], names=["CellID", "scpopcorn_cluster"])
      .set_index("CellID")
    )
    # keep cluster number -1 for the cells (artefacts) not used in
    # scpopcorn integration
    adata.obs = (
      adata.obs.merge(pcres, how="left", left_index=True, right_index=True)
      .fillna({"scpopcorn_cluster": -1})
    )
    adata.write_h5ad(output[0])


rule plot_sankey:
  input:
    "resources/int_h5ad/int_{tissue}_{fcaver}.h5ad",
  output:
    "results/sankey_scpopcorn/{tissue}_{fcaver}_sankey_scpopcorn.html",
  script:
    "../scripts/plot_sankey_integration.py"


append_final_output(expand(
  "results/sankey_scpopcorn/{tissue}_{fcaver}_sankey_scpopcorn.html",
  zip,
  tissue=samples["tissue"],
  fcaver=samples["fcaver"],
))


