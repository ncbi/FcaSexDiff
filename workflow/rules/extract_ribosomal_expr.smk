
rule extract_rp_expr:
  input:
    "imports/expr_h5ad/expr_{tissue}_{fcaver}.h5ad",
    "resources/sex_specific_annotations.csv",
    "resources/all_rp_genes.txt",
  output:
    "scraps/extra/cellfilter~{cellfilt}/resolution~{resol}/{tissue}/ribosomal_details_{tissue}_{fcaver}_{resol}_{cellfilt}.h5",
    "scraps/extra/cellfilter~{cellfilt}/resolution~{resol}/{tissue}/ribosomal_cells_{tissue}_{fcaver}_{resol}_{cellfilt}.csv",
    "exports/extra/cellfilter~{cellfilt}/resolution~{resol}/{tissue}/ribosomal_clusters_{tissue}_{fcaver}_{resol}_{cellfilt}.csv",
  script:
    "../scripts/extract_expression_ribosomal_proteins.py"

rule:
  input:
    body="scraps/extra/cellfilter~NoSexspecArtef/resolution~L6.0/body/ribosomal_details_body_stringent_L6.0_NoSexspecArtef.h5",
    head="scraps/extra/cellfilter~NoSexspecArtef/resolution~L6.0/head/ribosomal_details_head_stringent_L6.0_NoSexspecArtef.h5",
  output:
    "exports/extra/cellfilter~{cellfilt}/resolution~{resol}/ribosomal_proteins_L6.0_clusterwise.xlsx",
  script:
    "../scripts/export_summary_ribosomal_proteins.py"

append_final_output(
  [
    "scraps/extra/cellfilter~NoSexspecArtef/resolution~L6.0/body/ribosomal_details_body_stringent_L6.0_NoSexspecArtef.h5",
    "scraps/extra/cellfilter~NoSexspecArtef/resolution~L6.0/head/ribosomal_details_head_stringent_L6.0_NoSexspecArtef.h5",
    "exports/extra/cellfilter~NoSexspecArtef/resolution~L6.0/ribosomal_proteins_L6.0_clusterwise.xlsx",
  ]
)

