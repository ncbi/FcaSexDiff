import anndata as ad
import pandas as pd
import sys

writer = pd.ExcelWriter('ribosomal_proteins_L6.0_clusterwise.xlsx', engine='xlsxwriter')

sex_specific_annotations_file = "resources/sex_specific_annotations.csv"

def process(tissue):
    infile = f"imports/expr_h5ad/expr_{tissue}_stringent.h5ad"

    adata = ad.read_h5ad(infile)

    adata = adata[:, adata.var_names.str.startswith('RpS') |
                     adata.var_names.str.startswith('RpL')]

    print(", ".join(adata.var_names.tolist()))

    adata = adata[adata.obs.sex.isin(["male", "female"])]

    cell_filter = "NoSexspecArtef"
    if cell_filter in ["NoSexspecArtef", "NoSexspecArtefMuscle"]:
        # next discard cells that have sex-specific annotations
        # we get the information from sex_specific_annotations_file
        # that has a column sex_specific. Ensure it does not conflict 
        assert("sex_specific" not in adata.obs.columns)
        annots = (
            adata.obs[["annotation"]]
            .merge(pd.read_csv(sex_specific_annotations_file), how = "left")
            .drop_duplicates()
        )
        # we also make sure there is no annotation that is not
        # known in sex_specific_annotations_file
        print("\n".join(annots[annots.sex_specific.isna()]["annotation"].unique().tolist()))
        assert(sum(annots.sex_specific.isna()) == 0)
        to_remove = annots.annotation[annots.sex_specific != "no"].tolist()
        # remove cells annotated as artefacts too
        to_remove += ["artefact"]
        if cell_filter == "NoSexspecArtefMuscle":
            to_remove += ["muscle cell"]
        # now do the actual filtering
        adata = adata[~adata.obs.annotation.isin(to_remove)]
        print(adata)
        print(adata.obs.annotation.unique().tolist())


    meta_data = (
        adata.obs[["sex", "annotation", "batch", "L6.0", 'n_counts', 'n_genes']]
        .rename(columns={"n_counts":"n_umi"})
        .assign(cluster = lambda df: df["L6.0"].apply(lambda x: f"L6.0C{int(x):0>3d}"))
        .assign(batch = lambda df: df.batch.apply(lambda x: f"R{int(x):0>2d}"))
    )
    print(meta_data)
    meta_data["batch"] = meta_data["batch"].astype(str)

    annot = (
        meta_data[["cluster", "annotation"]]
        .assign(n=1)
        .groupby(["cluster", "annotation"])
        .agg("sum")
        .dropna()
        .reset_index()
        .groupby("cluster")
        .apply(lambda df: pd.Series({"major_annotation":
                                     df.annotation[df.n.idxmax()]}))
        .reset_index()
    )
    print(annot)

    expr = pd.DataFrame(
        adata.X.todense(),
        index = adata.obs_names,
        columns = adata.var_names
    )

    meta_data["rp_avg"] = expr.mean(axis=1)

    meta_data = meta_data.merge(expr, left_index=True, right_index=True)

    print(meta_data)

    meta_data.to_csv(f"ribosomal_proteins_{tissue}.csv")

    male = (
        meta_data.query("sex == 'male'")
        [["batch", "cluster", "n_umi", "rp_avg"]]
        .assign(n_cell = 1)
        .groupby(["cluster", "batch"])
        .agg({"n_umi": "mean", "rp_avg": "mean", "n_cell":"sum"})
        .dropna()
    )
    print(male)
    male_batch_total = (
        male
        .groupby(["cluster"])
        .agg({"n_umi": "mean", "rp_avg": "mean", "n_cell":"sum"})
        .dropna()
    )
    print(male_batch_total)
    male_batch_total = pd.concat([male_batch_total], keys=["RAll"],
                                      names=["batch"]).reorder_levels([1,0])
    print(male_batch_total)
    male = (
        pd.concat([male, male_batch_total])
        .unstack("batch")
    )
    print(male)
    
    
    female = (
        meta_data.query("sex == 'female'")
        [["batch", "cluster", "n_umi", "rp_avg"]]
        .assign(n_cell = 1)
        .groupby(["cluster", "batch"])
        .agg({"n_umi": "mean", "rp_avg": "mean", "n_cell":"sum"})
        .dropna()
    )
    print(female)
    female_batch_total = (
        female
        .groupby(["cluster"])
        .agg({"n_umi": "mean", "rp_avg": "mean", "n_cell":"sum"})
        .dropna()
    )
    print(female_batch_total)
    female_batch_total = pd.concat([female_batch_total], keys=["RAll"],
                                      names=["batch"]).reorder_levels([1,0])
    print(female_batch_total)
    female = (
        pd.concat([female, female_batch_total])
        .unstack("batch")
    )
    print(female)

    df = pd.concat(
        [female, male],
        keys=["female", "male"],
        names=["sex", "stats", "batch"],
        axis=1
    ).rename_axis("cluster", axis=0)
    print(df)
    df = df.reorder_levels([0,2,1], axis=1).sort_index(axis=1)
    print(df)
    df.index = pd.MultiIndex.from_frame(
        df.index.to_frame(index=False)
        .merge(annot)
    )
    print(df)

    df.to_excel(writer, sheet_name=tissue)


process("body")
process("head")

writer.save()
