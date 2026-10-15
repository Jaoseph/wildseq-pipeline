# prefilter_qc.smk

rule prefilter_qc:
    input:
        filtered_feature_bc_matrix="results/{SLX_ID}/{sample_index}/cellranger_count/outs/filtered_feature_bc_matrix"
        #filtered_feature_bc_matrix=rules.cellranger_count.output.filtered_feature_bc_matrix
    output:
        prefilter_mt="results/{SLX_ID}/{sample_index}/prefilter_qc_plots/{sample_index}_featureScatter_percent.mt_prefilter.png",
        prefilter_nFeature="results/{SLX_ID}/{sample_index}/prefilter_qc_plots/{sample_index}_featureScatter_nFeature_RNA_prefilter.png"
    shell:
        """
        echo "Plotting prefilter QC plots for {wildcards.SLX_ID} {wildcards.sample_index}"
        Rscript workflow/scripts/prefilter_qc.R \
            --slx_id {wildcards.SLX_ID} \
            --sample_index {wildcards.sample_index} \
            --filtered_feature_bc_matrix {input.filtered_feature_bc_matrix}
        echo "Finished plotting for {wildcards.SLX_ID} {wildcards.sample_index}"
        """