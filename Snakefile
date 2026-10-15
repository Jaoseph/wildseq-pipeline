# Snakefile
import pandas as pd

configfile: "/scratche/slow/jblab/chong02/wildseq/config/main.yaml"

SLX_IDs = config["SLX_IDs"]
metadata = pd.read_csv(config["sample_metadata"], dtype=str).set_index("sample_index").to_dict(orient="index")
sample_indices = list(metadata.keys())


rule all:
    input:
        expand("results/{SLX_ID}/{sample_index}/prefilter_qc_plots/{sample_index}_featureScatter_percent.mt_prefilter.png",
              SLX_ID=SLX_IDs,
              sample_index=sample_indices),
        expand("results/{SLX_ID}/{sample_index}/clone_assignment/{SLX_ID}_{sample_index}_all_assigned_seurat.rds",
               SLX_ID=SLX_IDs,
               sample_index=sample_indices),
        expand("results/{SLX_ID}/{sample_index}/clone_assignment/{SLX_ID}_{sample_index}_cellNumbers.tsv",
               SLX_ID=SLX_IDs,
               sample_index=sample_indices)

# Include the individual rule files for each step of the workflow
include: "workflow/rules/crukci_to_illumina.smk"
include: "workflow/rules/cellranger_count.smk"
include: "workflow/rules/barcode_extraction_10X_tdTom.smk"
include: "workflow/rules/barcode_umi_match.smk"
include: "workflow/rules/enrich_pcr.smk"
include: "workflow/rules/barcode_assignment.smk"
include: "workflow/rules/prefilter_qc.smk"
include: "workflow/rules/clone_assignment.smk"
include: "workflow/rules/normalization_qc_analysis.smk"
# include: "workflow/rules/differential_expression_analysis.smk"
# include: "workflow/rules/demo_integration.smk"
# include: "workflow/rules/control_analysis.smk"
# include: "workflow/rules/scTransform.smk"
# include: "workflow/rules/scTransformSLX26635.smk"