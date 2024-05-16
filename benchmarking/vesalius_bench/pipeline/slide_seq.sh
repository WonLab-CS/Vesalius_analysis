#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --job-name=SSv2
#SBATCH --mem-per-cpu=80GB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=<patrick.martin@cshs.org>
#SBATCH -p defq

module load R/4.2.1
module load rlibs/4.2.1
module load hdf5
Rscript /common/martinp4/benchmarking_out/Vesalius/pipeline/slide_seq.r
Rscript /common/martinp4/benchmarking_out/Vesalius/pipeline/slide_seq_plot.r
