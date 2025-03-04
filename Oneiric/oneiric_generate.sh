#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=Oneiric
#SBATCH --mem-per-cpu=200GB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=<patrick.martin@cshs.org>
#SBATCH -p defq

module load R/4.2.1 rlibs/4.2.1
script_loc="/home/martinp4/common/Vesalius_analysis/Oneiric/"
output="/common/wonklab/synthetic_spatial/"
seed=1453
Rscript ${script_loc}oneiric_generate.r $output $seed