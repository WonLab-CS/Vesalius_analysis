#!/bin/bash
#SBATCH --job-name=SSV2_format
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --output=array_%A_%a.out
#SBATCH --error=array_%A_%a.err
#SBATCH --array=1-1%1
#SBATCH --mem-per-cpu=80GB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=<patrick.martin@cshs.org>
#SBATCH -p defq
module load R/4.2.1 rlibs/4.2.1
input="/common/wonklab/SSv2/"
input_ref="/common/wonklab/mouse_hippocampus_scRNAseq/SCRef_hippocampus.RDS"
output="/common/wonklab/SSv2/"
script_loc="/home/martinp4/common/Vesalius_analysis/Spatial_data_sets/SSv2/"
file_tag="SSv2_formatted"
seed="Puck_200115_08"
query="Puck_190921_21"
Rscript ${script_loc}ssv2_format.r $SLURM_ARRAY_TASK_ID $input $input_ref $output $file_tag $seed $query
