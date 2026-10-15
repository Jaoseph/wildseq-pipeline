#!/bin/bash

#SBATCH -J wildseq_step1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=200:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jao.chong@cruk.cam.ac.uk
#SBATCH -p epyc
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

snakemake --workflow-profile workflow/profiles/ \
    --use-conda -p -k --snakefile Snakefile \
    --until barcode_assignment prefilter_qc
    # --until cellranger_count
    
