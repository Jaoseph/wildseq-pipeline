# barcode_umi_match.smk

rule barcode_umi_match:
    input:
        bowtie_output=rules.barcode_extraction_10X_tdTom.output.bowtie_output,
        cb_ub_txt=rules.barcode_extraction_10X_tdTom.output.cb_ub_txt
    output:
        cb_table="results/{SLX_ID}/{sample_index}/barcode_umi_match/{sample_index}_CBC_table_full.txt"
    threads: 8
    shell:
        """
        echo "Running barcode UMI match on {wildcards.SLX_ID} {wildcards.sample_index}"
        Rscript workflow/scripts/barcode_10X.R \
            --sample {wildcards.sample_index} \
            --theoretical_bowtie {input.bowtie_output} \
            --cb_ub_txt {input.cb_ub_txt} \
            --output_file {output.cb_table}
        """
