#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --job-name=SLAT_stereo
#SBATCH --mem-per-cpu=80GB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=<patrick.martin@cshs.org>
#SBATCH -p defq

source /home/martinp4/anaconda3/etc/profile.d/conda.sh
conda activate slat
data_type="Stereo"
input="/common/wonklab/Stereo_seq/"
output="/common/martinp4/benchmarking_out/SLAT/report/"
script_loc="/home/martinp4/common/Vesalius_analysis/benchmarking/SLAT_bench/"
seed_tag="E2S6"
query_tag="E2S5"
python ${script_loc}SLAT_bio_spa.py $data_type $input $output $seed_tag $query_tag

conda deactivate