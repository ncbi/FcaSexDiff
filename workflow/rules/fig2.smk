
rule draw_heatmap:
    input:
        sexdiff_h5ad.path
    output:
        heatmap_biased_genes.path
    script:
        "../scripts/fig2/draw_heatmap_biased_genes.R"

rule draw_upset:
    input:
        "scraps/biased_genes/all_biased_genes_{resol}.tsv",
    output:
        "scraps/fig2/upset_biased_genes_{tissue}_{resol}.pdf",
    script:
        "../scripts/fig2/draw_upset_biased_genes.R"


append_final_output(
    expand(
        expand(
            [heatmap_biased_genes.path],
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

