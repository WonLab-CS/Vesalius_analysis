#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --job-name=GPSA_seqFISH
#SBATCH --mem-per-cpu=400GB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=<patrick.martin@cshs.org>
#SBATCH -p defq

source /home/martinp4/anaconda3/etc/profile.d/conda.sh
conda activate gpsa
data_type="seqFISH"
input="/common/wonklab/seqFISH/"
output="/common/martinp4/benchmarking_out/GPSA/report/"
script_loc="/home/martinp4/common/Vesalius_analysis/benchmarking/GPSA_bench/"
seed_tag="embryo1"
query_tag="embryo3"
python ${script_loc}GPSA_bio_spa.py $data_type $input $output $seed_tag $query_tag
conda deactivate