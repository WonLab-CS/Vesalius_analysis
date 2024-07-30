#!/bin/bash
#SBATCH --job-name=map
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=array_%A_%a.out
#SBATCH --error=array_%A_%a.err
#SBATCH --array=1-144%100
#SBATCH --mem-per-cpu=12GB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=<patrick.martin@cshs.org>
#SBATCH -p defq

module load R/4.2.1
module load rlibs/4.2.1
Rscript /common/martinp4/benchmarking_out/Vesalius/pipeline/vesalius_combinations.r $SLURM_ARRAY_TASK_ID "circle" 1
Rscript /common/martinp4/benchmarking_out/Vesalius/pipeline/vesalius_combinations.r $SLURM_ARRAY_TASK_ID "layered" 1
Rscript /common/martinp4/benchmarking_out/Vesalius/pipeline/vesalius_combinations.r $SLURM_ARRAY_TASK_ID "circle" 2
Rscript /common/martinp4/benchmarking_out/Vesalius/pipeline/vesalius_combinations.r $SLURM_ARRAY_TASK_ID "layered" 2
Rscript /common/martinp4/benchmarking_out/Vesalius/pipeline/vesalius_combinations.r $SLURM_ARRAY_TASK_ID "circle" 3
Rscript /common/martinp4/benchmarking_out/Vesalius/pipeline/vesalius_combinations.r $SLURM_ARRAY_TASK_ID "layered" 3
Rscript /common/martinp4/benchmarking_out/Vesalius/pipeline/vesalius_combinations.r $SLURM_ARRAY_TASK_ID "circle" 4
Rscript /common/martinp4/benchmarking_out/Vesalius/pipeline/vesalius_combinations.r $SLURM_ARRAY_TASK_ID "layered" 4
Rscript /common/martinp4/benchmarking_out/Vesalius/pipeline/vesalius_combinations.r $SLURM_ARRAY_TASK_ID "circle" 5
Rscript /common/martinp4/benchmarking_out/Vesalius/pipeline/vesalius_combinations.r $SLURM_ARRAY_TASK_ID "layered" 5
Rscript /common/martinp4/benchmarking_out/Vesalius/pipeline/vesalius_combinations.r $SLURM_ARRAY_TASK_ID "circle" 6
Rscript /common/martinp4/benchmarking_out/Vesalius/pipeline/vesalius_combinations.r $SLURM_ARRAY_TASK_ID "layered" 6
Rscript /common/martinp4/benchmarking_out/Vesalius/pipeline/vesalius_combinations.r $SLURM_ARRAY_TASK_ID "circle" 7
Rscript /common/martinp4/benchmarking_out/Vesalius/pipeline/vesalius_combinations.r $SLURM_ARRAY_TASK_ID "layered" 7

Rscript /common/martinp4/benchmarking_out/Vesalius/pipeline/vesalius_combinations_plot.r