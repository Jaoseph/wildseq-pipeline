# Cellranger count Rule
# cellranger_count.smk

rule cellranger_count:
    input:
        fastq_dir="data/raw_data/{SLX_ID}/",
        conversion_log="logs/crukci_to_illumina_conversion_logs/{SLX_ID}/converted_files.log"
    params:
        transcriptome="data/Homo_sapiens_GRCh38_WILDSeq_tdTom/",
        memory=200,
        output_dir=lambda wildcards: f"results/{wildcards.SLX_ID}/{wildcards.sample_index}/cellranger_count"  # parent directory
    output:
        filtered_feature_bc_matrix=directory("results/{SLX_ID}/{sample_index}/cellranger_count/outs/filtered_feature_bc_matrix"),        
        bam="results/{SLX_ID}/{sample_index}/cellranger_count/outs/possorted_genome_bam.bam"
    threads: 20
    shell:
        """
        echo "Running cellranger on {wildcards.SLX_ID} {wildcards.sample_index}"
        rm -rf results/{wildcards.SLX_ID}/{wildcards.sample_index}
        mkdir -p {params.output_dir}
        export PATH=$PATH:/mnt/scratche/slow/jblab/truong01/wild_seq/software/cellranger-9.0.1 && \
        cellranger count \
            --id={wildcards.sample_index} \
            --transcriptome={params.transcriptome} \
            --fastqs={input.fastq_dir} \
            --sample={wildcards.sample_index} \
            --localcores={threads} \
            --localmem={params.memory} \
            --create-bam=true \
            --disable-ui \
            --output-dir={params.output_dir} \
            --tenx-cloud-token-path=""
        """





