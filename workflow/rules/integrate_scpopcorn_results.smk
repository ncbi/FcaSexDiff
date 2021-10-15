import sys
sys.path.append('workflow/scripts')

from scpopcorn import MergeSingleCell
from scpopcorn import SingleCellData
import anndata as ad
import pickle

rule augment_h5ad:
  input:
    "scraps/lognorm_h5ad/lognorm_{tissue}_{fcaver}.h5ad",
    "hoards/scpopcorn/scpopcorn_{tissue}_{fcaver}.txt",
  output:
    "hoards/int_h5ad/int_{tissue}_{fcaver}.h5ad",
  run:
    adata = ad.read_h5ad(input[0])
    print(adata)
    pcres = (
      pd.read_table(input[1], names=["CellID", "scpopcorn_cluster"])
      .set_index("CellID")
    )
    print(pcres)
    adata.obs = (
      adata.obs.merge(pcres, how="left", left_index=True, right_index=True)
    )
    print(adata)
    adata.write_h5ad(output[0])

rule plot_sankey:
  input:
    "hoards/int_h5ad/int_{tissue}_{fcaver}.h5ad",
  output:
    "exports/scpopcorn_integration/{tissue}/{tissue}_{fcaver}_sankey_scpopcorn.html",
  script:
    "../scripts/plot_sankey_integration.py"

rule plot_cluster_similarity:
  input:
    "hoards/int_h5ad/int_{tissue}_{fcaver}.h5ad",
  output:
    "exports/scpopcorn_integration/{tissue}/{tissue}_{fcaver}_cluster_similarity_scpopcorn.pdf",
  script:
    "../scripts/plot_cluster_similarity.py"


append_final_output(expand(
  [
    "hoards/int_h5ad/int_{tissue}_{fcaver}.h5ad",
    #"exports/scpopcorn_integration/{tissue}/{tissue}_{fcaver}_sankey_scpopcorn.html",
    #"exports/scpopcorn_integration/{tissue}/{tissue}_{fcaver}_cluster_similarity_scpopcorn.pdf",
  ],
  zip,
  tissue=samples["tissue"],
  fcaver=samples["fcaver"],
))


