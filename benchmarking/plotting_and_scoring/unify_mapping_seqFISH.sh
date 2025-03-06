#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=UNIFY_seqFISH
#SBATCH --mem-per-cpu=10GB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=<patrick.martin@cshs.org>
#SBATCH -p defq

module load R/4.2.1 rlibs/4.2.1
script_loc="/home/martinp4/common/Vesalius_analysis/benchmarking/plotting_and_scoring/"
input_matched="/common/martinp4/benchmarking_out/"
data_type="seqFISH"
input_ref="/common/wonklab/seqFISH/"
output="/common/martinp4/benchmarking_out/report/"
ref="embryo1"
query="embryo3"
Rscript ${scritp_loc}unify_mapping_score.r $input_matched $input_ref $data_type $output $ref $query
Rscript ${scritp_loc}plot_mapping_score_bio.r $output $data_type 
Rscript ${script_loc}plot_mapping_event_bio.r $input_matched $input_ref $data_type $output $ref $query
Rscript ${script_loc}plot_mapping_contribution_bio.r $input_matched $input_ref $data_type $output $ref $query