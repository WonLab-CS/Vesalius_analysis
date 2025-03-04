#!/bin/bash
#SBATCH --job-name=Stereo_format
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --output=array_%A_%a.out
#SBATCH --error=array_%A_%a.err
#SBATCH --array=1-1%1
#SBATCH --mem-per-cpu=200GB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=<patrick.martin@cshs.org>
#SBATCH -p defq


# cd /common/wonklab/Stereo_seq/
# wget https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000058/stomics/E16.5_E2S4.MOSTA.h5ad
# wget https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000058/stomics/E16.5_E2S5.MOSTA.h5ad
# wget https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000058/stomics/E16.5_E2S6.MOSTA.h5ad

module load R/4.2.1 rlibs/4.2.1
seed_input="/common/wonklab/Stereo_seq/E16.5_E2S6.MOSTA.h5ad"
query_input="/common/wonklab/Stereo_seq/E16.5_E2S5.MOSTA.h5ad"
output="/common/wonklab/Stereo_seq/"
script_loc="/home/martinp4/common/Vesalius_analysis/Spatial_data_sets/Stereo/"
file_tag="Stereo_formatted"
seed_tag="E2S6"
query_tag="E2S5"
Rscript ${script_loc}Stereo_format.r $SLURM_ARRAY_TASK_ID $seed_input $query_input $output $file_tag $seed_tag $query_tag
