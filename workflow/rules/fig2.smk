rule draw_heatmap:
  input:
    "scraps/biased_genes/all_biased_genes_{resol}.tsv",
  output:
    "scraps/fig2/heatmap_biased_genes_{tissue}_{resol}.pdf",
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
  [
    "scraps/fig2/heatmap_biased_genes_body_L6.0.pdf",
    "scraps/fig2/heatmap_biased_genes_head_L6.0.pdf",
    "scraps/fig2/upset_biased_genes_body_L6.0.pdf",
    "scraps/fig2/upset_biased_genes_head_L6.0.pdf",
  ]
)

