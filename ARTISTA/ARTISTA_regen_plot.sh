#!/bin/bash
#SBATCH --job-name=AR_RE
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=array_%A_%a.out
#SBATCH --error=array_%A_%a.err
#SBATCH --array=1-1%1
#SBATCH --mem-per-cpu=50GB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=<patrick.martin@cshs.org>
#SBATCH -p defq


module load R/4.2.1 rlibs/4.2.1
input="/common/wonklab/Stereo_seq_arista/report/"
outout="/common/wonklab/Stereo_seq_arista/report/"
script_loc="/home/martinp4/common/Vesalius_analysis/ARTISTA/"
Rscript ${script_loc}ARISTA_regen_plot.r $SLURM_ARRAY_TASK_ID $input $output
