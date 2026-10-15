#!/bin/bash

#SBATCH -J prefilter_slx-26209
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=200:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=an.truong@cruk.cam.ac.uk
#SBATCH -p epyc
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

snakemake --workflow-profile workflow/profiles/ \
    --use-conda -p -k --snakefile Snakefile \
    # --until prefilter_qc 
