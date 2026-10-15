#!/bin/bash

#SBATCH -J name
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --time=200:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jao.chong@cruk.cam.ac.uk
#SBATCH -p epyc
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# snakemake \
  # --workflow-profile workflow/profiles/ \
  # --use-conda -p -k \
  # --snakefile Snakefile \
  # results/SLX-26635/seurat_analyses/nineteenday_analyses/nineteenday_DE_19_617uM_vs_19d_10dWO_wilcox_DEG_list.csv
  # results/SLX-26635/seurat_analyses/ciov6_SLX-26635_merged_SCTransformed_seurat.rds 
  # results/SLX-26209/ciov6_48H_merged_SCTransformed_seurat.rds


# Load Conda environment
source ~/.bashrc
conda activate wildseq

# Move to working directory
cd /scratchc/jblab/chong02/wildseq

# Run your R script
Rscript workflow/scripts/update_thresholds.R