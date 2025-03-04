#!/bin/bash
#SBATCH --job-name=StoS
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --output=array_%A_%a.out
#SBATCH --error=array_%A_%a.err
#SBATCH --array=1-3%1
#SBATCH --mem-per-cpu=80GB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=<patrick.martin@cshs.org>
#SBATCH -p defq

module load R/4.2.1 rlibs/4.2.1
input="/common/wonklab/seqFISH/"
output="/common/martinp4/stos/report/stos/"
script_loc="/home/martinp4/common/Vesalius_analysis/tech_to_tech/stos/"
Rscript ${script_loc}stos.r $SLURM_ARRAY_TASK_ID $input $output
Rscript ${script_loc}stos_plot.r $SLURM_ARRAY_TASK_ID $input $output