# enrich_pcr.smk

rule enrich_pcr:
    input:
        unmapped_txt=rules.barcode_extraction_10X_tdTom.output.unmapped_txt,
        unmapped_count=rules.barcode_extraction_10X_tdTom.output.unmapped_count,
    params:
        enrich_index=lambda wildcards: metadata[wildcards.sample_index]["enrich_index"]
    output:
        table_full="results/{SLX_ID}/{sample_index}/enrich_pcr/{sample_index}_EnrichPCR_table_full.txt",
        table_full_filtered="results/{SLX_ID}/{sample_index}/enrich_pcr/{sample_index}_EnrichPCR_table_full_filtered.txt"
    threads: 8
    shell:
        """
        echo "Processing PCR enrichment library for {wildcards.SLX_ID} {wildcards.sample_index}, {params.enrich_index}"
        bash workflow/scripts/enrich_pcr.sh \
        {input.unmapped_txt} {input.unmapped_count} \
        {wildcards.SLX_ID} {wildcards.sample_index} {params.enrich_index} \
        {output.table_full} {output.table_full_filtered}
        """