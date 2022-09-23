BOX_ROOT = "BoxNIH:FcaSexDiff"

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
        lambda wc: expanded_path(wc, bias_scatter_pdf.path)
    output:
        touch("FollowupPaper/Figures/BiasScatterPlots/{tissue}_scatter_{bias}_bias.pdf")
    shell:
        "echo rclone copyto {input} {BOX_ROOT}/{output}"

rule copy_bias_tsne_plots:
    input:
        lambda wc: expanded_path(wc, bias_tsne_pdf.path)
    output:
        touch("FollowupPaper/Figures/BiasTsnePlots/{tissue}_tsne_{bias}_bias.pdf")
    shell:
        "rclone copyto {input} {BOX_ROOT}/{output}"

rule copy_dsx_targets_enrichment:
    input:
        lambda wc: expanded_path(wc, dsx_targets_tfs_enrich.path)
    output:
        touch("FollowupPaper/Data/DsxTargetsTFsEnrichments/{tissue}_dsx_targets_tfs_enrichment.tsv")
    shell:
        "rclone copyto {input} {BOX_ROOT}/{output}"

rule copy_dsx_ffl_status:
    input:
        lambda wc: expanded_path(wc, dsx_ffl_status.path)
    output:
        touch("FollowupPaper/Data/DsxFFL/Summary/{tissue}_dsx_ffl_motif_status.tsv")
    shell:
        "rclone copyto {input} {BOX_ROOT}/{output}"

rule copy_dsx_ffl_detailed_status:
    input:
        lambda wc: expanded_path(wc, dsx_ffl_detailed_status_excel.path)
    output:
        touch("FollowupPaper/Data/DsxFFL/{tissue}_dsx_ffl_motif_detailed_status.xlsx")
    shell:
        "rclone copyto {input} {BOX_ROOT}/{output}"

rule copy_biased_gene_groups:
    input:
        lambda wc: expanded_path(wc, biased_gene_groups.path)
    output:
        touch("FollowupPaper/Data/BiasedGeneGroups/{tissue}_biased_gene_groups.tsv")
    shell:
        "rclone copyto {input} {BOX_ROOT}/{output}"

rule copy_sexdiff_excel:
    input:
        lambda wc: expanded_path(wc, sexdiff_excel.path)
    output:
        touch("FollowupPaper/Data/SexDiffExcel/{tissue}_sexdiff.xlsx")
    shell:
        "rclone copyto {input} {BOX_ROOT}/{output}"

rule copy_sexdiff_h5ad:
    input:
        lambda wc: expanded_path(wc, sexdiff_h5ad.path)
    output:
        touch("FollowupPaper/Data/SexDiffH5ad/{tissue}_sexdiff.h5ad")
    shell:
        "rclone copyto {input} {BOX_ROOT}/{output}"

append_final_output(
    expand(
        [
            "FollowupPaper/Figures/BiasScatterPlots/{tissue}_scatter_{bias}_bias.pdf",
            "FollowupPaper/Figures/BiasTsnePlots/{tissue}_tsne_{bias}_bias.pdf",
        ],
         tissue = tasks["tissue"],
         bias = ["count", "expr", "rp", "rpavg", "nonrp"],
    )
)
append_final_output(
    expand(
        [
            "FollowupPaper/Data/SexDiffExcel/{tissue}_sexdiff.xlsx",
            "FollowupPaper/Data/SexDiffH5ad/{tissue}_sexdiff.h5ad",
            "FollowupPaper/Data/DsxTargetsTFsEnrichments/{tissue}_dsx_targets_tfs_enrichment.tsv",
            "FollowupPaper/Data/DsxFFL/Summary/{tissue}_dsx_ffl_motif_status.tsv",
            "FollowupPaper/Data/DsxFFL/{tissue}_dsx_ffl_motif_detailed_status.xlsx",
            "FollowupPaper/Data/BiasedGeneGroups/{tissue}_biased_gene_groups.tsv",
        ],
         tissue = tasks["tissue"],
    )
)
