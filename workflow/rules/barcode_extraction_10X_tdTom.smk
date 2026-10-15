# barcode_extraction_10X_tdTom.smk

rule barcode_extraction_10X_tdTom:
    input:
        #bam=rules.cellranger_count.output.bam
        bam="results/{SLX_ID}/{sample_index}/cellranger_count/outs/possorted_genome_bam.bam"
    output:
        unmapped_txt="results/{SLX_ID}/{sample_index}/barcode_extraction/unmapped_theoretical.txt",
        unmapped_count="results/{SLX_ID}/{sample_index}/barcode_extraction/unmapped_theoretical_count.txt",
        cb_ub_txt="results/{SLX_ID}/{sample_index}/barcode_extraction/{sample_index}_WILDseq.CB.UB.txt",
        bowtie_output="results/{SLX_ID}/{sample_index}/barcode_extraction/bamtofastq_R2_theoretical.bowtie"
    params:
        bowtie_index="data/theoretical_barcodes2/WILDseq_index2"
    threads: 8
    shell:
        """
        echo "Running barcode extraction on {wildcards.SLX_ID} {wildcards.sample_index}"
        bash workflow/scripts/barcode_extraction_10X_tdTom.sh \
        {input.bam} {wildcards.SLX_ID} {wildcards.sample_index} {params.bowtie_index} \
        {output.unmapped_txt} {output.unmapped_count} {output.cb_ub_txt} {output.bowtie_output}
        """

