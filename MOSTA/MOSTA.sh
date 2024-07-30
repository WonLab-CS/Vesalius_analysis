#!/bin/bash
#SBATCH --job-name=MOSTA
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --output=array_%A_%a.out
#SBATCH --error=array_%A_%a.err
#SBATCH --array=1-7%2
#SBATCH --mem-per-cpu=200GB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=<patrick.martin@cshs.org>
#SBATCH -p defq


module load R/4.2.1
module load rlibs/4.2.1
module load hdf5
Rscript /common/wonklab/Stereo_seq/pipeline/MOSTA.R $SLURM_ARRAY_TASK_ID
