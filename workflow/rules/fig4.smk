# find genes that vary across pseudotime using graph_test function in monocle3
# which in turn runs Moran's I test. this step does not filter any gene, keeps
# results of Moran's I test for all genes.

rule get_degs_across_pseudotime:
  input:
    "imports/monocle_cds/monocle_cds_adult_{rnaseq}_germline.rds",
  output:
    "scratch/warped/morans_test_adult_{rnaseq}.rds",
  log:
    "logs/get_degs_across_pseudotime_{rnaseq}.log",
  script:
    "../scripts/fig4/get_degs_across_pseudotime.R"


# generate a matrix of gene expression smoothened (using spline)
# and scaled (Z-score) across pseudotime

rule get_smooth_expr:
  input:
    "imports/monocle_cds/monocle_cds_adult_{rnaseq}_germline.rds",
  output:
    "scratch/warped/smooth_expr_adult_{rnaseq}.rds"
  log:
    "logs/get_smooth_expr_adult_{rnaseq}.log",
  script:
    "../scripts/fig4/get_smooth_expr.R"


# use dynamic time warping to align sc and sn pseudotime

rule warp:
  input:
    sc_cds = "imports/monocle_cds/monocle_cds_adult_sc_germline.rds",
    sn_cds = "imports/monocle_cds/monocle_cds_adult_sn_germline.rds",
    sc_morans = "scratch/warped/morans_test_adult_sc.rds",
    sn_morans = "scratch/warped/morans_test_adult_sn.rds",
    sc_smooth = "scratch/warped/smooth_expr_adult_sc.rds",
    sn_smooth = "scratch/warped/smooth_expr_adult_sn.rds",
  output:
    "scratch/warped/warped_adult_sc_sn.rds",
  log:
    "logs/warp_adult_sc_sn.log",
  script:
    "../scripts/fig4/warp_adult_sc_sn.R"

# Generate fig4 panels separately

# panel showing dimplots on nuclei data

rule get_dim_plots_nuclei:
  input:
    cds = "imports/monocle_cds/monocle_cds_adult_sn_germline.rds",
    seurat = "imports/seurat/seurat_adult_sn_all.rds",
  output:
    "scratch/fig4/dim_plots_nuclei.rds",
  log:
    "logs/get_dim_plots_nuclei.log",
  script:
    "../scripts/fig4/get_dim_plots_nuclei.R"


# panel showing dimplots on cell data

rule get_dim_plots_cell:
  input:
    "imports/monocle_cds/monocle_cds_adult_sc_germline.rds",
  output:
    "scratch/fig4/dim_plots_cell.rds",
  log:
    "logs/get_dim_plots_cell.log",
  script:
    "../scripts/fig4/get_dim_plots_cell.R"


# panel showing feature plots on cell data

rule get_feature_plots_cell:
  input:
    "imports/monocle_cds/monocle_cds_adult_sc_germline.rds",
  output:
    "scratch/fig4/feature_plots_cell.rds",
  log:
    "logs/get_feature_plots_cell.log",
  script:
    "../scripts/fig4/get_feature_plots_cell.R"


# panel showing warped gene expression

rule get_warped_plots:
  input:
    "scratch/warped/warped_adult_sc_sn.rds",
  output:
    "scratch/fig4/warped_plots.rds",
  log:
    "logs/get_warped_plots.log",
  script:
    "../scripts/fig4/get_warped_plots.R"


# panel showing heatmap of smoothened gene expression

rule get_heatmaps:
  input:
    sc_cds = "imports/monocle_cds/monocle_cds_adult_sc_germline.rds",
    sn_cds = "imports/monocle_cds/monocle_cds_adult_sn_germline.rds",
    sc_smooth = "scratch/warped/smooth_expr_adult_sc.rds",
    sn_smooth = "scratch/warped/smooth_expr_adult_sn.rds",
  output:
    "scratch/fig4/heatmaps.rds",
  log:
    "logs/get_heatmaps.log",
  script:
    "../scripts/fig4/get_heatmaps.R"


# stitch together fig4 panels using plotgardener

rule draw_fig4:
  input:
    dim_plots_cell = "scratch/fig4/dim_plots_cell.rds",
    dim_plots_nuclei = "scratch/fig4/dim_plots_nuclei.rds",
    feature_plots_cell = "scratch/fig4/feature_plots_cell.rds",
    heatmaps = "scratch/fig4/heatmaps.rds",
    warped_plots = "scratch/fig4/warped_plots.rds",
  output:
    "exports/figures/fig4.pdf",
  log:
    "logs/draw_fig4.log",
  script:
    "../scripts/fig4/draw_fig4.R"


append_final_output(
  [
    "scratch/warped/warped_adult_sc_sn.rds",
    "exports/figures/fig4.pdf",
  ]
)

