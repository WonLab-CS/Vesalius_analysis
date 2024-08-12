#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --job-name=UNIFY
#SBATCH --mem-per-cpu=10GB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=<patrick.martin@cshs.org>
#SBATCH -p defq

module load R/4.2.1
module load rlibs/4.2.1
module load hdf5
Rscript /common/martinp4/benchmarking_out/pipeline/unify_and_score.r
Rscript /common/martinp4/benchmarking_out/pipeline/unify_plot.r

