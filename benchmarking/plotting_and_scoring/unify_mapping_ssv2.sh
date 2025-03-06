#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=UNIFY_ssv2
#SBATCH --mem-per-cpu=10GB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=<patrick.martin@cshs.org>
#SBATCH -p defq

module load R/4.2.1 rlibs/4.2.1
script_loc="/home/martinp4/common/Vesalius_analysis/benchmarking/plotting_and_scoring/"
input_matched="/common/martinp4/benchmarking_out/"
input_ref="/common/wonklab/SSv2/"
data_type="SSv2"
output="/common/martinp4/benchmarking_out/report/"
ref="Puck_200115_08"
query="Puck_190921_21"
Rscript ${script_loc}unify_mapping_score.r $input_matched $input_ref $data_type $output $ref $query
Rscript ${script_loc}plot_mapping_score_bio.r $output $data_type
Rscript ${script_loc}plot_mapping_event_bio.r $input_matched $input_ref $data_type $output $ref $query
Rscript ${script_loc}plot_mapping_contribution_bio.r $input_matched $input_ref $data_type $output $ref $query