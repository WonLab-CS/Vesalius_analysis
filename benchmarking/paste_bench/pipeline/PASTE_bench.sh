#!/bin/bash
#SBATCH --job-name=paste
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=array_%A_%a.out
#SBATCH --error=array_%A_%a.err
#SBATCH --array=1-144%10
#SBATCH --mem-per-cpu=15GB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=<patrick.martin@cshs.org>
#SBATCH -p defq

source /home/martinp4/anaconda3/etc/profile.d/conda.sh
conda activate paste
cd common/
python /common/martinp4/submits/PASTE_bench.py $SLURM_ARRAY_TASK_ID "circle"
python /common/martinp4/submits/PASTE_bench.py $SLURM_ARRAY_TASK_ID "layered"
conda deactivate