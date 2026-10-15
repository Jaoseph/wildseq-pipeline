# normalization_qc_analysis.smk

rule normalization_qc_analysis:
    output:
        libsize_detectedgenes_beforenorm="results/{SLX_ID}/seurat_analyses/QC_plots/{SLX_ID}_combined_lib_size_detectedgenes_relationship_before_norm.pdf",
        libsize_detectedgenes_afternorm="results/{SLX_ID}/seurat_analyses/QC_plots/{SLX_ID}_combined_lib_size_detectedgenes_relationship_after_norm.pdf",
        mean_variance="results/{SLX_ID}/seurat_analyses/QC_plots/{SLX_ID}_combined_mean_variance_relationship_before_norm.pdf",
        mean_detection="results/{SLX_ID}/seurat_analyses/QC_plots/{SLX_ID}_combined_mean_detection_vs_poisson_before_norm.pdf",
        residual_var="results/{SLX_ID}/seurat_analyses/QC_plots/{SLX_ID}_combined_residual_variance_vs_geomean_after_norm.pdf"
    shell:
        """
        echo "Running Seurat normalization visualization on {wildcards.SLX_ID}"
        Rscript workflow/scripts/seurat/normalization_qc_analysis.R \
            --slx_id {wildcards.SLX_ID} \
            --output_dir results/{wildcards.SLX_ID}/seurat_analyses/QC_plots
        """
