def expanded_path(wc, inp):
    todo = tasks.query("tissue == @wc.tissue")
    return expand(
        [inp],
        zip,
        allow_missing = True,
        fcaver = todo["fcaver"],
        resol = todo["resols"],
        cellfilt = todo["cellfilts"],
        expr = "LogNorm",
    )

rule copy_bias_scatter_plots:
    input:
        lambda wc: expanded_path(wc, bias_tsne_pdf.path)
    output:
        touch("{tissue}_scatter_{bias}_bias")
    params:
        box = "BoxNIH:FcaSexDiff/FollowupPaper/Figures/BiasScatterPlots/"
    shell:
        "rclone copyto {input} {params.box}/{output}.pdf"

rule copy_bias_tsne_plots:
    input:
        lambda wc: expanded_path(wc, bias_tsne_pdf.path)
    output:
        touch("{tissue}_tsne_{bias}_bias")
    params:
        box = "BoxNIH:FcaSexDiff/FollowupPaper/Figures/BiasTsnePlots/"
    shell:
        "rclone copyto {input} {params.box}/{output}.pdf"

rule copy_dsx_targets_enrichment:
    input:
        lambda wc: expanded_path(wc, dsx_targets_tfs_enrich.path)
    output:
        touch("{tissue}_dsx_targets_tfs_enrichment")
    params:
        box = "BoxNIH:FcaSexDiff/FollowupPaper/Data/DsxTargetsTFsEnrichments"
    shell:
        "rclone copyto {input} {params.box}/{output}.tsv"

rule copy_biased_gene_groups:
    input:
        lambda wc: expanded_path(wc, biased_gene_groups.path)
    output:
        touch("{tissue}_biased_gene_groups")
    params:
        box = "BoxNIH:FcaSexDiff/FollowupPaper/Data/BiasedGeneGroups"
    shell:
        "rclone copyto {input} {params.box}/{output}.tsv"

append_final_output(
    expand(
        [
            "{tissue}_scatter_{bias}_bias",
            "{tissue}_tsne_{bias}_bias",
        ],
         tissue = tasks["tissue"],
         bias = ["count", "expr", "rp", "nonrp"],
    )
)
append_final_output(
    expand(
        [
            "{tissue}_dsx_targets_tfs_enrichment",
            "{tissue}_biased_gene_groups",
        ],
         tissue = tasks["tissue"],
    )
)
