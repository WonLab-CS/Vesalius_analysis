#!/bin/bash
#SBATCH --job-name=Vesalius_bench
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=array_%A_%a.out
#SBATCH --error=array_%A_%a.err
#SBATCH --array=1-144%50
#SBATCH --mem-per-cpu=12GB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=<patrick.martin@cshs.org>
#SBATCH -p defq

module load R/4.2.1 rlibs/4.2.1
input="/common/wonklab/synthetic_spatial/"
output="/common/martinp4/benchmarking_out/Vesalius/report/"
script_loc="/home/martinp4/common/Vesalius_analysis/benchmarking/Vesalius_bench/"
for i in {1..14}; do
    Rscript ${script_loc}Vesalius_bench.r $SLURM_ARRAY_TASK_ID "one_cell" $input $output $i
    Rscript ${script_loc}Vesalius_bench.r $SLURM_ARRAY_TASK_ID "two_cell" $input $output $i
    Rscript ${script_loc}Vesalius_bench.r $SLURM_ARRAY_TASK_ID "contact_one" $input $output $i
    Rscript ${script_loc}Vesalius_bench.r $SLURM_ARRAY_TASK_ID "contact_two" $input $output $i
    Rscript ${script_loc}Vesalius_bench.r $SLURM_ARRAY_TASK_ID "circle" $input $output $i
    Rscript ${script_loc}Vesalius_bench.r $SLURM_ARRAY_TASK_ID "layered" $input $output $i
    Rscript ${script_loc}Vesalius_bench.r $SLURM_ARRAY_TASK_ID "dropped" $input $output $i
done