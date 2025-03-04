#!/bin/bash
#SBATCH --job-name=SLAT_bench
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=array_%A_%a.out
#SBATCH --error=array_%A_%a.err
#SBATCH --array=1-144%50
#SBATCH --mem-per-cpu=15GB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=<patrick.martin@cshs.org>
#SBATCH -p defq

source /home/martinp4/anaconda3/etc/profile.d/conda.sh
conda activate slat
input="/common/wonklab/synthetic_spatial"
output="/common/martinp4/benchmarking_out/SLAT/report/"
script_loc="/home/martinp4/common/Vesalius_analysis/benchmarking/SLAT_bench/"

python ${script_loc}SLAT_bench.py --slurm_id $SLURM_ARRAY_TASK_ID --type "one_cell" --input $input --output $output
python ${script_loc}SLAT_bench.py --slurm_id $SLURM_ARRAY_TASK_ID --type "two_cell" --input $input --output $output
python ${script_loc}SLAT_bench.py --slurm_id $SLURM_ARRAY_TASK_ID --type "contact_one" --input $input --output $output
python ${script_loc}SLAT_bench.py --slurm_id $SLURM_ARRAY_TASK_ID --type "contact_two" --input $input --output $output
python ${script_loc}SLAT_bench.py --slurm_id $SLURM_ARRAY_TASK_ID --type "circle" --input $input --output $output
python ${script_loc}SLAT_bench.py --slurm_id $SLURM_ARRAY_TASK_ID --type "layered" --input $input --output $output
python ${script_loc}SLAT_bench.py --slurm_id $SLURM_ARRAY_TASK_ID --type "dropped" --input $input --output $output
conda deactivate