import sys
sys.path.append('workflow/scripts')

from scpopcorn import MergeSingleCell
from scpopcorn import SingleCellData
import anndata as ad
import pickle

rule make_single_cell_obj:
  input:
    "resources/exp_h5ad/exp_{tissue}_{fcaver}.h5ad",
  output:
    "resources/scpopcorn/tmp_{tissue}_{fcaver}/sc_merged.pkl",
    "resources/scpopcorn/tmp_{tissue}_{fcaver}/pcp.pkl",
  run:
    adata = ad.read_h5ad(input[0])

    # skip the cells annotated as artefacts
    adata = adata[adata.obs["annotation"] != "artefact"]

    male = adata[adata.obs["sex"] == "male"]
    female = adata[adata.obs["sex"] == "female"]

    n_male_celltype = len(male.obs['annotation'].unique())
    n_female_celltype = len(female.obs['annotation'].unique())
    n_celltype = len(adata.obs['annotation'].unique())

    pcparams = {}

    pcparams['nsuper1'] = min(1000, int(male.obs.shape[0]/20))
    pcparams['nsuper2'] = min(1000, int(female.obs.shape[0]/20))

    #pcparams['kmin'] = min(n_male_celltype, n_female_celltype, n_celltype-5)
    #pcparams['kmax'] = max(n_male_celltype, n_female_celltype, n_celltype+5)
    #pcparams['k'] = min(n_celltype,(pcparams['kmin'] + pcparams['kmax'])//2)

    pcparams['kmin'] = n_celltype
    pcparams['kmax'] = n_celltype+10
    pcparams['k'] = n_celltype

    Test1 = SingleCellData()
    Test1.SetData_H5AD(male)

    Test2 = SingleCellData()
    Test2.SetData_H5AD(female)

    Test1.Normalized_per_Cell()
    Test1.FindHVG()
    Test1.Log1P()

    Test2.Normalized_per_Cell()
    Test2.FindHVG()
    Test2.Log1P()

    MSingle = MergeSingleCell(Test1, Test2)

    pickle.dump(MSingle, open(output[0], 'wb'))
    pickle.dump(pcparams, open(output[1], 'wb'))

rule build_withion_between:
  input:
    "resources/scpopcorn/tmp_{tissue}_{fcaver}/sc_merged.pkl",
    "resources/scpopcorn/tmp_{tissue}_{fcaver}/pcp.pkl",
  output:
    "resources/scpopcorn/tmp_{tissue}_{fcaver}/sc_within_between.pkl",
  run:
    MSingle = pickle.load(open(input[0], 'rb'))
    pcparams = pickle.load(open(input[1], 'rb'))

    MSingle.MultiDefineSuperCell(pcparams["nsuper1"], pcparams["nsuper1"])
    MSingle.ConstructWithinSimiarlityMat_SuperCellLevel()
    MSingle.ConstructBetweenSimiarlityMat_SuperCellLevel()

    pickle.dump(MSingle, open(output[0], 'wb'))

rule joint_partition_NKcut:
  input:
    "resources/scpopcorn/tmp_{tissue}_{fcaver}/sc_within_between.pkl",
    "resources/scpopcorn/tmp_{tissue}_{fcaver}/pcp.pkl",
  output:
    "resources/scpopcorn/tmp_{tissue}_{fcaver}/sc_nkcut.pkl",
  run:
    MSingle = pickle.load(open(input[0], 'rb'))
    pcparams = pickle.load(open(input[1], 'rb'))

    MSingle.SDP_NKcut(pcparams["k"])

    pickle.dump(MSingle, open(output[0], 'wb'))


rule NKcut_rounding:
  input:
    "resources/scpopcorn/tmp_{tissue}_{fcaver}/sc_nkcut.pkl",
    "resources/scpopcorn/tmp_{tissue}_{fcaver}/pcp.pkl",
  output:
    "resources/scpopcorn/scpopcorn_{tissue}_{fcaver}.txt",
    "resources/scpopcorn/scpopcorn_{tissue}_{fcaver}.pkl",
  run:
    MSingle = pickle.load(open(input[0], 'rb'))
    pcparams = pickle.load(open(input[1], 'rb'))

    CResult = MSingle.NKcut_Rounding(pcparams["kmin"], pcparams["kmax"])

    MSingle.OutputResult(output[0])
    pickle.dump(MSingle, open(output[1], 'wb'))


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
    # increase cluster number by 1, keep cluster 0 as artefacts
    pcres["scpopcorn_cluster"] += 1
    adata.obs = (
      adata.obs.merge(pcres, how="left", left_index=True, right_index=True)
      .fillna({"scpopcorn_cluster": 0})
    )
    adata.write_h5ad(output[0])


rule plot_sankey:
  input:
    "resources/int_h5ad/int_{tissue}_{fcaver}.h5ad",
  output:
    "results/sankey_scpopcorn/{tissue}_{fcaver}_sankey_scpopcorn.html",
  script:
    "../scripts/plot_sankey_integration.py"


rule get_stats:
  input:
    "resources/scpopcorn/scpopcorn_{tissue}_{fcaver}.pkl",
  output:
    "resources/scpopcorn/scpopcoen_stats_{tissue}_{fcaver}.txt",
  run:
    MSingle = pickle.load(open(input[0], 'rb'))
    MSingle.StatResult()

print(samples)

append_final_output(expand(
  "results/sankey_scpopcorn/{tissue}_{fcaver}_sankey_scpopcorn.html",
  #"resources/scpopcorn/scpopcoen_stats_{tissue}_{fcaver}.txt",
  #"resources/scpopcorn/tmp_{tissue}_{fcaver}/sc_within_between.pkl",
  #"resources/scpopcorn/tmp_{tissue}_{fcaver}/sc_nkcut.pkl",
  #"resources/int_h5ad/int_{tissue}_{fcaver}.h5ad",
  zip,
  tissue=samples["tissue"],
  fcaver=samples["fcaver"],
))


