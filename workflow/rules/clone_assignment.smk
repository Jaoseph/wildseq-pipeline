# clone_assignment.smk


rule clone_assignment:
    input:
        filtered_feature_bc_matrix="results/{SLX_ID}/{sample_index}/cellranger_count/outs/filtered_feature_bc_matrix",
        barcode_assigned_table="results/{SLX_ID}/{sample_index}/barcode_assignment/barcode_assignment_{SLX_ID}_{sample_index}.txt"
    output:
        full_seurat="results/{SLX_ID}/{sample_index}/clone_assignment/{SLX_ID}_{sample_index}_prefiltered_seurat.rds",
        all_assigned_seurat="results/{SLX_ID}/{sample_index}/clone_assignment/{SLX_ID}_{sample_index}_all_assigned_seurat.rds",
        cell_numbers_table="results/{SLX_ID}/{sample_index}/clone_assignment/{SLX_ID}_{sample_index}_cellNumbers.tsv"
    params:
        filter_nFeatureRNA_lower=lambda wildcards: metadata[wildcards.sample_index]["lower"],
        filter_nFeatureRNA_upper=lambda wildcards: metadata[wildcards.sample_index]["upper"],
        filter_percent_mt=lambda wildcards: metadata[wildcards.sample_index]["mt"]
    shell:
        """
        echo "Running seurat clone assignment on {wildcards.SLX_ID} {wildcards.sample_index}"
        Rscript workflow/scripts/clone_assignment.R \
            --slx_id {wildcards.SLX_ID} \
            --sample {wildcards.sample_index} \
            --filtered_feature_bc_matrix {input.filtered_feature_bc_matrix} \
            --barcode_assigned_table {input.barcode_assigned_table} \
            --filter_nFeatureRNA_lower {params.filter_nFeatureRNA_lower} \
            --filter_nFeatureRNA_upper {params.filter_nFeatureRNA_upper} \
            --filter_percent_mt {params.filter_percent_mt} \
            --full_seurat {output.full_seurat} \
            --all_assigned_seurat {output.all_assigned_seurat} \
            --cell_numbers_table {output.cell_numbers_table}
        echo "Completed seurat clone assignment on {wildcards.SLX_ID} {wildcards.sample_index}"
        """
