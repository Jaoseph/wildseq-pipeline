#!/bin/bash
#SBATCH --job-name=snake-{rule}-{wildcards}
#SBATCH --cpus-per-task={resources.cpus_per_task}
#SBATCH --mem={resources.mem_mb}
#SBATCH --time={resources.runtime}
#SBATCH --partition={resources.slurm_partition}
#SBATCH --account={resources.slurm_account}
#SBATCH --output={logdir}/{rule}/{rule}-{wildcards}-%j.out
#SBATCH --error={logdir}/{rule}/{rule}-{wildcards}-%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=an.truong@cruk.cam.ac.uk

set -euo pipefail
{exec_job}
