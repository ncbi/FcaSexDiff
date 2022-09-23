from define_sexbias_cutoffs import *

info = dict(
  cluster = (
    "Cluster name (cell type if annotation is used)"
  ),
  markers = (
    "Top 10 marker genes that are differentially highly expressed in"
    " the cluster than the rest of the clusters"
  ),
  annotations_female = (
    "Replicate-wise counts of annotated celltypes among female cells in the cluster"
  ),
  annotations_male = (
    "Replicate-wise counts of annotated celltypes among male cells in the cluster"
  ),
  annotations = (
    "Sex and replicate -wise counts of annotated celltypes among all cells in the cluster"
  ),
  major_annotation = (
      "Major annotated celltype in the cluster, useful if clusters are"
      " not by annotation"
  ),
  tissue_count_female = (
    "Number of female cells in the tissue (combining all replicates)"
  ),
  tissue_count_male = (
    "Number of male cells in the tissue (combining all replicates)"
  ),
  cluster_type = (
    "Classify cluster based on whether there is at least one cell in the cluster from either sex"
  ),
  cluster_count_female = (
      "Number of female cells in the cluster (combining all replicates)"
  ),
  cluster_count_male = (
      "Number of male cells in the cluster (combining all replicates)"
  ),
  tissue_rep_counts_female = (
      "Replicate-wise counts of female cells in the tissue"
  ),
  tissue_rep_counts_male = (
      "Replicate-wise counts of male cells in the tissue"
  ),
  cluster_rep_counts_female = (
      "Replicate-wise counts of female cells in the cluster"
  ),
  cluster_rep_counts_male = (
      "Replicate-wise counts of male cells in the cluster"
  ),
  cluster_frac_female = (
    "Fraction of cells in the cluster that are female (combining all replicates)"
  ),
  cluster_frac_male = (
    "Fraction of cells in the cluster that are male (combining all replicates)"
  ),
  cluster_rep_fracs_female = (
    "Replicate-wise fractions of cells in the cluster that are female"
  ),
  cluster_rep_fracs_male = (
    "Replicate-wise fractions of cells in the cluster that are male"
  ),
  cluster_rep_fracs_mean_female = (
    "Mean of replicate-wise fractions of cells in the cluster that are female "
    "(should differ from fraction of cells in the cluster that are female if the "
    "replicates have different number of cells)"
  ),
  cluster_rep_fracs_mean_male = (
    "Mean of replicate-wise fractions of cells in the cluster that are male "
    "(should differ from fraction of cells in the cluster that are male if the "
    "replicates have different number of cells)"
  ),
  cluster_rep_fracs_sd_female = (
    "Standard deviation of replicate-wise fractions of cells in the cluster that are female"
  ),
  cluster_rep_fracs_sd_male = (
    "Standard deviation of replicate-wise fractions of cells in the cluster that are male"
  ),
  pval_binom = (
    "Pvalue from bionamial test on the counts considering all replicates together"
  ),
  pval_fisher = (
    "Pvalue from Fisher's test on the counts considering all replicates together"
  ),
  pval_wilcox = (
    "Pvalue from Wilcoxon-Mann-Whitney U test on the replicate-wise fractions"
  ),
  pval_ttest = (
    "Pvalue from 2-sample t test on the replicate-wise fractions"
  ),
  padj_binom = (
      "BH corrected pvalue from binomial test on the counts considering all replicates together"
  ),
  padj_fisher = (
    "BH corrected pvalue from Fisher's test on the counts considering all replicates together"
  ),
  padj_wilcox = (
    "BH corrected pvalue from Wilcoxon-Mann-Whitney U test on the replicate-wise fractions"
  ),
  padj_ttest = (
    "BH corrected pvalue from s-sample t test on the replicate-wise fractions"
  ),
  log2_count_bias = (
      "log2((NF + eps)/(NM + eps)) where NF (NM) is the replicate-average "
      " fraction of cells in cluster that are female (male) and"
      f" eps = {PSEUDO_COUNT}"
  ),
  count_bias_padj = (
      "BH corrected pvalue from Binomial test for tissues having single "
      "replicate, otherwise from Wilcoxon-Mann-Whitney U test"
  ),
  count_bias_type = (
      f"Cluster is female- or male-biased based on count_bias_padj < "
      f"{PADJ_CUTOFF_COUNT} and log2_count_bias > {LFC_CUTOFF_COUNT} or "
      f"< -{LFC_CUTOFF_COUNT} ({2**LFC_CUTOFF_COUNT}-fold change)"
  ),
  female_gene = (
      "Number of genes female-biased in the cluster"
  ),
  male_gene = (
      "Number of genes male-biased in the cluster"
  ),
  stats = (
      "Statistics for a gene within cluster"
  ),
  avg_all = (
      "Average expression of gene in all cells in cluster"
  ),
  avg_male = (
      "Average expression of gene in all male cells in cluster"
  ),
  avg_female = (
      "Average expression of gene in all female cells in cluster"
  ),
  avg_nz_all = (
      "Average expression of gene in all cells in cluster where it expressed"
  ),
  avg_nz_male = (
      "Average expression of gene in all male cells in cluster where it"
      " expressed"
  ),
  avg_nz_female = (
      "Average expression of gene in all female cells in cluster where it"
      " expressed"
  ),
  frac_all = (
      "Fraction of cells where gene is expressed (atleast one UMI)"
  ),
  frac_female = (
      "Fraction of female cells where gene is expressed (atleast one UMI)"
  ),
  frac_male = (
      "Fraction of male cells where gene is expressed (atleast one UMI)"
  ),
  log2fc = (
      "Log2 of fold change of average expression in favor of female cells"
      " against male cells"
  ),
  log2fc_scanpy = (
      "Log2 of fold change of average expression in favor of female cells"
      " against male cells as computed by scanpy.tl.rank_genes_groups"
  ),
  padj = (
      "BH corrected pvalue from Wilcoxon test for male cells vs female cells"
      " expression as computed by scanpy.tl.rank_genes_groups"
  ),
  bias = (
      f"Gene is female- or male-biased based on padj < {PADJ_CUTOFF_EXPR}"
      f" and log2fc > {LFC_CUTOFF_EXPR} or < -{LFC_CUTOFF_EXPR}"
      f" ({2**LFC_CUTOFF_EXPR}-fold change)"
  ),
  bias_scanpy = (
      f"gene is female- or male-biased based on padj < {PADJ_CUTOFF_EXPR}"
      f" and log2fc_scanpy > {LFC_CUTOFF_EXPR} or < -{LFC_CUTOFF_EXPR}"
      f" ({2**LFC_CUTOFF_EXPR}-fold change)"
  ),
  marker_log2fc = (
      "Log2 of fold change of average expression in favor of cells"
      " in the cluster against rest of the cells"
  ),
  marker_padj = (
      "BH corrected pvalue from Wilcoxon test for cluster cells vs rest cells"
      " expression as computed by scanpy.tl.rank_genes_groups"
  ),
  marker_score = (
      "Marker score computed by scanpy.tl.rank_genes_groups on clusters"
  ),
  is_marker = (
      f"Gene is cluster marker based on marker_padj < {PADJ_CUTOFF_MARKER}"
      f" and marker_log2fc > {LFC_CUTOFF_MARKER}"
      f" ({2**LFC_CUTOFF_MARKER}-fold change)"
  ),
  chr = (
      "Chromosome where gene is located"
  ),
  symbol = "Gene symbol",
  FBgn = "Flybase ID of gene",
  female_cls = (
      "Number of clusters where gene is female-biased"
  ),
  male_cls = (
      "Number of clusters where gene is male-biased"
  ),
  umi_tissue = (
      "Average number of UMIs for the gene in the whole tissue"
  ),
  nz_umi_tissue = (
      "Average number of UMIs for the gene in cells in the whole tissue where"
      " it expressed"
  ),
  norm_tissue = (
      "Normalized expression of gene (before log1p) in the whole tissue"
  ),
  nz_norm_tissue = (
      "Normalized expression of gene (before log1p) in cells in whole tissue"
      " where it expressed"
  ),
  x = (
    "dsx target that coregulates another dsx target y"
  ),
  y = (
    "dsx target that is coregulated by another dsx target x"
  ),
  y_in_x_regulon = (
      "Whether y is among the regulon named x in FCA SCENIC data"
  ),
  n_y_biased = (
      "Number of clusters where y is female or male biased"
  ),
  n_y_female_biased = (
      "Number of clusters where y is female biased"
  ),
  ffelement = (
      "Detailed information about the feed forward loop involving dsx and its"
      " two targets x and y where x also coragulates y. "
  ),
  ffl_dsx= (
      "dsx expression in the cluster: 'F' if femaled biased, else 'M' if male"
      " biased, else E if atleast 10 percent of cells have dsx UMI atleast 1,"
      " else N implying not expressed"
  ),
  ffl_x= (
      "x expression in the cluster: 'F' if femaled biased, else 'M' if male"
      " biased, else E if atleast 10 percent of cells have x UMI atleast 1,"
      " else N implying not expressed"
  ),
  ffl_y= (
      "y expression in the cluster: 'F' if femaled biased, else 'M' if male"
      " biased, else E if atleast 10 percent of cells have y UMI atleast 1,"
      " else N implying not expressed"
  ),
  scenic_baa = (
      "binarized average activity score of the FCA SCENIC regulon x in the"
      " cluster computed by applying AUC threshold on the average AUC scores"
      " of all cells in the cluster"
  ),
  scenic_aba = (
      "average binarized activity score of the FCA SCENIC regulon x in the"
      " cluster computed by taking average of binarized AUC score of each"
      " cell in the cluster"
  ),
)

