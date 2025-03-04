#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --job-name=CytoSpace_seqFISH
#SBATCH --mem-per-cpu=80GB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=<patrick.martin@cshs.org>
#SBATCH -p defq

source /home/martinp4/anaconda3/etc/profile.d/conda.sh
conda activate cytospace

data_type="seqFISH"
input="/common/wonklab/seqFISH/"
output="/common/martinp4/benchmarking_out/CytoSpace/report/"
script_loc="/home/martinp4/common/Vesalius_analysis/benchmarking/CytoSpace_bench/"
seed_tag="embryo1"
query_tag="embryo3"
lab="noLab"
tmp="${output}${data_type}_${lab}/"
mkdir -p $tmp
python ${script_loc}CytoSpace_bio_spa.py $data_type $input $tmp $seed_tag $query_tag

cytospace --single-cell \
   --scRNA-path "${tmp}scRNA_${lab}_${data_type}.txt"\
   --cell-type-path "${tmp}scLabels_${lab}_${data_type}.txt" \
   --st-path "${tmp}stRNA_${lab}_${data_type}.txt" \
   --coordinates-path "${tmp}stCoord_${lab}_${data_type}.txt" \
   --st-cell-type-path "${tmp}stLabels_${lab}_${data_type}.txt" \
   -o $tmp \
   --number-of-selected-spots 10000 \
   --mean-cell-numbers 1 \
   --number-of-processors 1
python ${script_loc}CytoSpace_clean_bio_spa.py $tmp $output $data_type $seed_tag $query_tag
rm -rf "${tmp}"


data_type="seqFISH"
input="/common/wonklab/seqFISH/"
output="/common/martinp4/benchmarking_out/CytoSpace/report/"
script_loc="/home/martinp4/common/Vesalius_analysis/benchmarking/CytoSpace_bench/"
seed_tag="embryo1"
query_tag="embryo3"
lab="trueLab"
tmp="${output}${data_type}_${lab}/"
mkdir -p $tmp
python ${script_loc}CytoSpace_bio_spa.py $data_type $input $tmp $seed_tag $query_tag

cytospace --single-cell \
   --scRNA-path "${tmp}scRNA_${lab}_${data_type}.txt"\
   --cell-type-path "${tmp}scLabels_${lab}_${data_type}.txt" \
   --st-path "${tmp}stRNA_${lab}_${data_type}.txt" \
   --coordinates-path "${tmp}stCoord_${lab}_${data_type}.txt" \
   --st-cell-type-path "${tmp}stLabels_${lab}_${data_type}.txt" \
   -o $tmp \
   --number-of-selected-spots 10000 \
   --mean-cell-numbers 1 \
   --number-of-processors 1
python ${script_loc}CytoSpace_clean_bio_spa.py $tmp $output $data_type $seed_tag $query_tag
rm -rf "${tmp}"