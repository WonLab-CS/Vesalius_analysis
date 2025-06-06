#!/bin/bash
#SBATCH --job-name=RAZA
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --output=array_%A_%a.out
#SBATCH --error=array_%A_%a.err
#SBATCH --array=1-10000%50
#SBATCH --mem-per-cpu=10GB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=<patrick.martin@cshs.org>
#SBATCH -p defq

module load R/4.2.1 rlibs/4.2.1
input="/common/wonklab/RAZA/split_data/"
output_plots="/common/wonklab/RAZA/output_plots/"
output_data="/common/wonklab/RAZA/output_plots/"
script_loc="/home/martinp4/common/Vesalius_analysis/IMC/"
Rscript ${script_loc}RAZA_balence.R $SLURM_ARRAY_TASK_ID
