#!/bin/bash
#SBATCH --job-name=Tangram_bench
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
conda activate tangram-env

input="/common/wonklab/synthetic_spatial/"
output="/common/martinp4/benchmarking_out/Tangram/report/"
script_loc="/home/martinp4/common/Vesalius_analysis/benchmarking/Tangram_bench/"

python ${script_loc}Tangram_bench.py $SLURM_ARRAY_TASK_ID "one_cell" $input $output
python ${script_loc}Tangram_bench.py $SLURM_ARRAY_TASK_ID "two_cell" $input $output
python ${script_loc}Tangram_bench.py $SLURM_ARRAY_TASK_ID "contact_one" $input $output
python ${script_loc}Tangram_bench.py $SLURM_ARRAY_TASK_ID "contact_two" $input $output
python ${script_loc}Tangram_bench.py $SLURM_ARRAY_TASK_ID "circle" $input $output
python ${script_loc}Tangram_bench.py $SLURM_ARRAY_TASK_ID "layered" $input $output
python ${script_loc}Tangram_bench.py $SLURM_ARRAY_TASK_ID "dropped" $input $output
conda deactivate