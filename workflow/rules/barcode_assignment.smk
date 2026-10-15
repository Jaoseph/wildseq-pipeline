# barcode_assignment.smk


rule barcode_assignment:
    input:
        cb_table=rules.barcode_umi_match.output.cb_table,
        table_full_filtered=rules.enrich_pcr.output.table_full_filtered
    output:
        barcode_assigned_table="results/{SLX_ID}/{sample_index}/barcode_assignment/barcode_assignment_{SLX_ID}_{sample_index}.txt"
    threads: 8
    shell:
        """
        echo "Running barcode assignment on {wildcards.SLX_ID} {wildcards.sample_index}"
        Rscript workflow/scripts/barcode_assignment.R \
            --slx_id {wildcards.SLX_ID} \
            --sample {wildcards.sample_index} \
            --cb_table {input.cb_table} \
            --table_full_filtered {input.table_full_filtered} \
            --output_file {output.barcode_assigned_table}
        """
