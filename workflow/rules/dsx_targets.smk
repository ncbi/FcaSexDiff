rule create_grn_from_netrex_occupancy:
    input:
        netrex_male =  "imports/DSX_targets/NetRex-Male.txt",
        netrex_female = "imports/DSX_targets/NetRex-Female.txt",
        dsx_occupancy = "imports/DSX_targets/DSX_occupancy_targets.tsv",
    output:
        "resources/GRN.tsv"
    script:
        "../scripts/dsx_targets/create_grn_from_netrex_occupancy.py"


rule extract_feed_forward_motifs:
    input:
        "resources/GRN.tsv"
    output:
        ff_motif = "exports/dsx_targets/dsx_ff_motifs.tsv",
        one_hop = "exports/dsx_targets/dsx_one_hop.txt",
        two_hop = "exports/dsx_targets/dsx_two_hop.txt",
    script:
        "../scripts/dsx_targets/extract_feed_forward_motifs_from_grn.py"



rule get_biased_genes:
    input:
        sexdiff_h5ad.path
    output:
        biased_genes.path,
    script:
        "../scripts/get_biased_genes_filtered_clusters.R"


rule group_biased_genes:
    input:
        biased_genes.path,
    output:
        biased_gene_groups.path,
    script:
        "../scripts/group_biased_genes.R"


rule compute_enrichment_dsx_targets_tfs:
    input:
        expr = normalized_h5ad.path,
        bias_grps = biased_gene_groups.path,
        dsx_occupancy = "imports/DSX_targets/DSX_occupancy_targets.tsv",
        putative_tfs = "imports/FlyMine/FlyMine_FlyTF_putative.tsv",
        trusted_tfs = "imports/FlyMine/FlyMine_FlyTF_trusted.tsv",
        ff_motif = "exports/dsx_targets/dsx_ff_motifs.tsv",
        one_hop = "exports/dsx_targets/dsx_one_hop.txt",
        two_hop = "exports/dsx_targets/dsx_two_hop.txt",
    output:
        dsx_targets_tfs_enrich.path
    script:
        "../scripts/dsx_targets/compute_enrichment_dsxtargets_tfs.py"

def prepare_input_paths_enrichment(wc):
    return {
        task["tissue"]: expand(
            dsx_targets_tfs_enrich.path,
            tissue = task["tissue"],
            fcaver = task["fcaver"],
            resol = task["resols"],
            cellfilt = task["cellfilts"],
            expr = "LogNorm",
        )
        for i, task in tasks.iterrows()
    }

rule plot_dsx_targets_tfs_enrichments:
    input:
        unpack(prepare_input_paths_enrichment)
    output:
        "exports/dsx_targets/dsx_targets_tfs_enrichments.pdf"
    script:
        "../scripts/dsx_targets/plot_dsx_targets_enrichment.R"

rule get_regulon_info:
    input:
        loom = "imports/loom/{tissue}_{fcaver}.loom",
    output:
        targets = regulon_targets.path,
        auc = regulon_auc.path,
    script:
        "../scripts/dsx_targets/get_regulon_info.R"

rule get_cluster_regulon_activity:
    input:
        auc = regulon_auc.path,
        cells = filtered_cells_csv.path,
        adata = h5ad_from_loom.path
    output:
        cluster_regulon_activity.path
    script:
        "../scripts/dsx_targets/get_cluster_regulon_activity.R"


rule check_dsx_ffl_status:
    input:
        expr = normalized_h5ad.path,
        bias = sexdiff_h5ad.path,
        biased_groups = biased_gene_groups.path,
        ff_motif = "exports/dsx_targets/dsx_ff_motifs.tsv",
        regulon_targets = regulon_targets.path,
        regulon_activity = cluster_regulon_activity.path,
    output:
        summary = dsx_ffl_status.path,
        detailed = dsx_ffl_detailed_status.path,
    script:
        "../scripts/dsx_targets/check_ffl_status_in_tissue.py"

rule export_dsx_ffl_status:
    input:
        ffl = dsx_ffl_detailed_status.path,
        bias = sexdiff_h5ad.path,
    output:
        dsx_ffl_detailed_status_excel.path,
    script:
        "../scripts/dsx_targets/export_ffl_status_to_excel.py"


append_final_output([
    "resources/GRN.tsv",
    "exports/dsx_targets/dsx_ff_motifs.tsv",
    "exports/dsx_targets/dsx_targets_tfs_enrichments.pdf",
])

append_final_output(
    expand(
        expand(
          [
            dsx_targets_tfs_enrich.path,
            dsx_ffl_detailed_status_excel.path,
          ],
          zip,
          allow_missing=True,
          tissue=tasks["tissue"],
          fcaver=tasks["fcaver"],
          resol=tasks["resols"],
          cellfilt=tasks["cellfilts"],
        ),
        expr = "LogNorm",
    )
)

append_final_output(
    expand(
          [
            regulon_auc.path,
          ],
          zip,
          allow_missing=True,
          tissue=tasks["tissue"],
          fcaver=tasks["fcaver"],
    )
)


append_final_output(
    expand(
          [
            cluster_regulon_activity.path,
          ],
          zip,
          allow_missing=True,
          tissue=tasks["tissue"],
          fcaver=tasks["fcaver"],
          resol=tasks["resols"],
          cellfilt=tasks["cellfilts"],
    )
)


