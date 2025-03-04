#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=UNIFY_bench
#SBATCH --mem-per-cpu=10GB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=<patrick.martin@cshs.org>
#SBATCH -p defq

module load R/4.2.1 rlibs/4.2.1
script_loc="/home/martinp4/common/Vesalius_analysis/benchmarking/plotting_and_scoring/"
input_matched="/common/martinp4/benchmarking_out/"
data_type="synthetic"
input_ref="/common/wonklab/synthetic_spatial/"
output="/common/martinp4/benchmarking_out/report/"
ref="sample"
query="sample"

Rscript ${scritp_loc}unify_mapping_score_bench.r $input_matched $input_ref $data_type $output $ref $query
Rscript ${scritp_loc}plot_mapping_score_bench.r $output
Rscript ${script_loc}plot_mapping_event_bench.r $input_matched $input_ref $data_type $output $ref $query
Rscript ${scritp_loc}plot_performance_bench.r $input_matched $input_ref $data_type $output $ref $query

