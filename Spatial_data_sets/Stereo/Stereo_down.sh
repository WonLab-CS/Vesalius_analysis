#!/bin/bash
#SBATCH --job-name=Stereo_dwld
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --output=array_%A_%a.out
#SBATCH --error=array_%A_%a.err
#SBATCH --array=1-1%1
#SBATCH --mem-per-cpu=1GB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=<patrick.martin@cshs.org>
#SBATCH -p defq

cd /common/wonklab/Stereo_seq/
wget https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000058/stomics/E16.5_E2S4.MOSTA.h5ad
wget https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000058/stomics/E16.5_E2S5.MOSTA.h5ad
wget https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000058/stomics/E16.5_E2S6.MOSTA.h5ad