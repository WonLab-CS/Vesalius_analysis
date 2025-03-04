#!/bin/bash
#SBATCH --job-name=seqFISH_format
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --output=array_%A_%a.out
#SBATCH --error=array_%A_%a.err
#SBATCH --array=1-1%1
#SBATCH --mem-per-cpu=15GB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=<patrick.martin@cshs.org>
#SBATCH -p defq
module load R/4.2.1 rlibs/4.2.1
input="/common/wonklab/seqFISH/"
output="/common/wonklab/seqFISH/"
script_loc="/home/martinp4/common/Vesalius_analysis/Spatial_data_sets/seqFISH/"
file_tag="seqFISH_formatted"
Rscript ${script_loc}seqFISH_format.r $SLURM_ARRAY_TASK_ID $input $output $file_tag
