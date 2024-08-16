#!/bin/bash
#SBATCH --job-name=precast
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=array_%A_%a.out
#SBATCH --error=array_%A_%a.err
#SBATCH --array=1-144%50
#SBATCH --mem-per-cpu=12GB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=<patrick.martin@cshs.org>
#SBATCH -p defq

module load R/4.2.1
module load rlibs/4.2.1
Rscript /common/martinp4/benchmarking_out/precast/pipeline/precast_bench.r $SLURM_ARRAY_TASK_ID "circle"
Rscript /common/martinp4/benchmarking_out/precast/pipeline/precast_bench.r $SLURM_ARRAY_TASK_ID "layered"