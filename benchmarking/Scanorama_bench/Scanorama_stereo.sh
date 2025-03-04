#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=Scanorama_stereo
#SBATCH --mem-per-cpu=500GB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=<patrick.martin@cshs.org>
#SBATCH -p defq

source /home/martinp4/anaconda3/etc/profile.d/conda.sh
conda activate scano
data_type="Stereo"
input="/common/wonklab/Stereo_seq/"
output="/common/martinp4/benchmarking_out/Scanorama/report/"
script_loc="/home/martinp4/common/Vesalius_analysis/benchmarking/Scanorama_bench/"
seed_tag="E2S6"
query_tag="E2S5"
python ${script_loc}Scanorama_bio_spa.py $data_type $input $output $seed_tag $query_tag
conda deactivate