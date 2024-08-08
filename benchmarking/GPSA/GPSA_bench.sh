#!/bin/bash
#SBATCH --job-name=gpsa
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=array_%A_%a.out
#SBATCH --error=array_%A_%a.err
#SBATCH --array=1-144%5
#SBATCH --mem-per-cpu=100GB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=<patrick.martin@cshs.org>
#SBATCH -p defq

source /home/martinp4/anaconda3/etc/profile.d/conda.sh
conda activate gpsa
python /common/martinp4/benchmarking_out/GPSA/pipeline/GPSA_bench.py $SLURM_ARRAY_TASK_ID "circle"
python /common/martinp4/benchmarking_out/GPSA/pipeline/GPSA_bench.py $SLURM_ARRAY_TASK_ID "layered"
conda deactivate