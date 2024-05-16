#!/bin/bash
#SBATCH --job-name=map
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --output=array_%A_%a.out
#SBATCH --error=array_%A_%a.err
#SBATCH --array=1-9%1
#SBATCH --mem-per-cpu=50GB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=<patrick.martin@cshs.org>
#SBATCH -p defq
module load R/4.2.1
module load rlibs/4.2.1
module load hdf5
Rscript /common/martinp4/benchmarking_out/Vesalius/pipeline/seq_FISH.r $SLURM_ARRAY_TASK_ID
Rscript /common/martinp4/benchmarking_out/Vesalius/pipeline/seq_FISH_plot.r $SLURM_ARRAY_TASK_ID