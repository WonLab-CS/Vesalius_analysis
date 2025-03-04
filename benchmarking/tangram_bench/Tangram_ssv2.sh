#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=15-00:00:00
#SBATCH --cpus-per-task=5
#SBATCH --job-name=Tangram_ssv2
#SBATCH --mem-per-cpu=100GB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=<patrick.martin@cshs.org>
#SBATCH -p defq

source /home/martinp4/anaconda3/etc/profile.d/conda.sh
conda activate tangram-env
data_type="SSv2"
input="/common/wonklab/SSv2/"
output="/common/martinp4/benchmarking_out/Tangram/report/"
script_loc="/home/martinp4/common/Vesalius_analysis/benchmarking/Tangram_bench/"
seed_tag="Puck_200115_08"
query_tag="Puck_190921_21"
python ${script_loc}Tangram_bio_spa.py $data_type $input $output $seed_tag $query_tag

conda deactivate