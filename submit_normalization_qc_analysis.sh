#!/bin/bash

#SBATCH -J normalization_qc_SLX-26635
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=200:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jao.chong@cruk.cam.ac.uk
#SBATCH -p epyc
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

snakemake \
  --workflow-profile workflow/profiles/ \
  --use-conda -p -k \
  --snakefile Snakefile \
  results/SLX-26635/seurat_analyses/QC_plots/SLX-26635_combined_lib_size_detectedgenes_relationship_before_norm.pdf