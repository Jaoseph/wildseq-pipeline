rule crukci_to_illumina:
    input:
        fastq_dir="data/raw_data/{SLX_ID}/"
    output:
        conversion_log="logs/crukci_to_illumina_conversion_logs/{SLX_ID}/converted_files.log"
    shell:
        """
        echo "Converting file names for {wildcards.SLX_ID}..." > {output.conversion_log}
        python3 workflow/scripts/crukci_to_illumina.py {input.fastq_dir}
        echo "File conversion completed for {wildcards.SLX_ID}." >> {output.conversion_log}
        """
