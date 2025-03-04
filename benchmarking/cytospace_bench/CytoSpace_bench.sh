#!/bin/bash
#SBATCH --job-name=CytoSpace_bench
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
conda activate cytospace

cyto_input="/common/martinp4/benchmarking_out/CytoSpace/report/" # loc of temp files for cytospace
input="/common/wonklab/synthetic_spatial/"
output="/common/martinp4/benchmarking_out/CytoSpace/report/"
script_loc="/home/martinp4/common/Vesalius_analysis/benchmarking/CytoSpace_bench/"

type="one_cell"
lab="noLab"
tmp="${cyto_input}${SLURM_ARRAY_TASK_ID}_${lab}_${type}/"
mkdir -p $tmp
python ${script_loc}CytoSpace_bench.py $SLURM_ARRAY_TASK_ID $type $tmp $input 

cytospace --single-cell \
   --scRNA-path "${tmp}scRNA_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt"\
   --cell-type-path "${tmp}scLabels_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --st-path "${tmp}stRNA_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --coordinates-path "${tmp}stCoord_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --st-cell-type-path "${tmp}stLabels_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   -o $tmp \
   --number-of-selected-spots 5000 \
   --mean-cell-numbers 1 \
   --number-of-processors 1
python ${script_loc}CytoSpace_clean.py $tmp $output $type
rm -rf "${tmp}"


type="two_cell"
lab="noLab"
tmp="${cyto_input}${SLURM_ARRAY_TASK_ID}_${lab}_${type}/"
mkdir -p $tmp
python ${script_loc}CytoSpace_bench.py $SLURM_ARRAY_TASK_ID $type $tmp $input 

cytospace --single-cell \
   --scRNA-path "${tmp}scRNA_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt"\
   --cell-type-path "${tmp}scLabels_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --st-path "${tmp}stRNA_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --coordinates-path "${tmp}stCoord_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --st-cell-type-path "${tmp}stLabels_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   -o $tmp \
   --number-of-selected-spots 5000 \
   --mean-cell-numbers 1 \
   --number-of-processors 1
python ${script_loc}CytoSpace_clean.py $tmp $output $type
rm -rf "${tmp}"


type="contact_one"
lab="noLab"
tmp="${cyto_input}${SLURM_ARRAY_TASK_ID}_${lab}_${type}/"
mkdir -p $tmp
python ${script_loc}CytoSpace_bench.py $SLURM_ARRAY_TASK_ID $type $tmp $input 

cytospace --single-cell \
   --scRNA-path "${tmp}scRNA_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt"\
   --cell-type-path "${tmp}scLabels_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --st-path "${tmp}stRNA_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --coordinates-path "${tmp}stCoord_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --st-cell-type-path "${tmp}stLabels_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   -o $tmp \
   --number-of-selected-spots 5000 \
   --mean-cell-numbers 1 \
   --number-of-processors 1
python ${script_loc}CytoSpace_clean.py $tmp $output $type
rm -rf "${tmp}"



type="contact_two"
lab="noLab"
tmp="${cyto_input}${SLURM_ARRAY_TASK_ID}_${lab}_${type}/"
mkdir -p $tmp
python ${script_loc}CytoSpace_bench.py $SLURM_ARRAY_TASK_ID $type $tmp $input 

cytospace --single-cell \
   --scRNA-path "${tmp}scRNA_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt"\
   --cell-type-path "${tmp}scLabels_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --st-path "${tmp}stRNA_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --coordinates-path "${tmp}stCoord_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --st-cell-type-path "${tmp}stLabels_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   -o $tmp \
   --number-of-selected-spots 5000 \
   --mean-cell-numbers 1 \
   --number-of-processors 1
python ${script_loc}CytoSpace_clean.py $tmp $output $type
rm -rf "${tmp}"


type="circle"
lab="noLab"
tmp="${cyto_input}${SLURM_ARRAY_TASK_ID}_${lab}_${type}/"
mkdir -p $tmp
python ${script_loc}CytoSpace_bench.py $SLURM_ARRAY_TASK_ID $type $tmp $input 

cytospace --single-cell \
   --scRNA-path "${tmp}scRNA_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt"\
   --cell-type-path "${tmp}scLabels_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --st-path "${tmp}stRNA_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --coordinates-path "${tmp}stCoord_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --st-cell-type-path "${tmp}stLabels_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   -o $tmp \
   --number-of-selected-spots 5000 \
   --mean-cell-numbers 1 \
   --number-of-processors 1
python ${script_loc}CytoSpace_clean.py $tmp $output $type
rm -rf "${tmp}"


type="layered"
lab="noLab"
tmp="${cyto_input}${SLURM_ARRAY_TASK_ID}_${lab}_${type}/"
mkdir -p $tmp
python ${script_loc}CytoSpace_bench.py $SLURM_ARRAY_TASK_ID $type $tmp $input 

cytospace --single-cell \
   --scRNA-path "${tmp}scRNA_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt"\
   --cell-type-path "${tmp}scLabels_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --st-path "${tmp}stRNA_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --coordinates-path "${tmp}stCoord_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --st-cell-type-path "${tmp}stLabels_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   -o $tmp \
   --number-of-selected-spots 5000 \
   --mean-cell-numbers 1 \
   --number-of-processors 1
python ${script_loc}CytoSpace_clean.py $tmp $output $type
rm -rf "${tmp}"



type="dropped"
lab="noLab"
tmp="${cyto_input}${SLURM_ARRAY_TASK_ID}_${lab}_${type}/"
mkdir -p $tmp
python ${script_loc}CytoSpace_bench.py $SLURM_ARRAY_TASK_ID $type $tmp $input 

cytospace --single-cell \
   --scRNA-path "${tmp}scRNA_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt"\
   --cell-type-path "${tmp}scLabels_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --st-path "${tmp}stRNA_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --coordinates-path "${tmp}stCoord_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --st-cell-type-path "${tmp}stLabels_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   -o $tmp \
   --number-of-selected-spots 5000 \
   --mean-cell-numbers 1 \
   --number-of-processors 1
python ${script_loc}CytoSpace_clean.py $tmp $output $type
rm -rf "${tmp}"



#-----------------------------------------------------------------------------#
# Demonstrating performance when cell labels are given
#-----------------------------------------------------------------------------#

type="one_cell"
lab='trueLab'
tmp="${cyto_input}${SLURM_ARRAY_TASK_ID}_${lab}_${type}/"
mkdir -p $tmp
python ${script_loc}CytoSpace_bench.py $SLURM_ARRAY_TASK_ID $type $tmp $input

cytospace --single-cell \
   --scRNA-path "${tmp}scRNA_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt"\
   --cell-type-path "${tmp}scLabels_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --st-path "${tmp}stRNA_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --coordinates-path "${tmp}stCoord_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --st-cell-type-path "${tmp}stLabels_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   -o $tmp \
   --number-of-selected-spots 5000 \
   --mean-cell-numbers 1 \
   --number-of-processors 1
python ${script_loc}CytoSpace_clean.py $tmp $output $type
rm -rf "${tmp}"


type="two_cell"
lab="trueLab"
tmp="${cyto_input}${SLURM_ARRAY_TASK_ID}_${lab}_${type}/"
mkdir -p $tmp
python ${script_loc}CytoSpace_bench.py $SLURM_ARRAY_TASK_ID $type $tmp $input

cytospace --single-cell \
   --scRNA-path "${tmp}scRNA_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt"\
   --cell-type-path "${tmp}scLabels_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --st-path "${tmp}stRNA_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --coordinates-path "${tmp}stCoord_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --st-cell-type-path "${tmp}stLabels_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   -o $tmp \
   --number-of-selected-spots 5000 \
   --mean-cell-numbers 1 \
   --number-of-processors 1
python ${script_loc}CytoSpace_clean.py $tmp $output $type
rm -rf "${tmp}"


type="contact_one"
lab="trueLab"
tmp="${cyto_input}${SLURM_ARRAY_TASK_ID}_${lab}_${type}/"
mkdir -p $tmp
python ${script_loc}CytoSpace_bench.py $SLURM_ARRAY_TASK_ID $type $tmp $input 

cytospace --single-cell \
   --scRNA-path "${tmp}scRNA_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt"\
   --cell-type-path "${tmp}scLabels_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --st-path "${tmp}stRNA_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --coordinates-path "${tmp}stCoord_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --st-cell-type-path "${tmp}stLabels_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   -o $tmp \
   --number-of-selected-spots 5000 \
   --mean-cell-numbers 1 \
   --number-of-processors 1
python ${script_loc}CytoSpace_clean.py $tmp $output $type
rm -rf "${tmp}"



type="contact_two"
lab="trueLab"
tmp="${cyto_input}${SLURM_ARRAY_TASK_ID}_${lab}_${type}/"
mkdir -p $tmp
python ${script_loc}CytoSpace_bench.py $SLURM_ARRAY_TASK_ID $type $tmp $input 

cytospace --single-cell \
   --scRNA-path "${tmp}scRNA_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt"\
   --cell-type-path "${tmp}scLabels_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --st-path "${tmp}stRNA_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --coordinates-path "${tmp}stCoord_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --st-cell-type-path "${tmp}stLabels_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   -o $tmp \
   --number-of-selected-spots 5000 \
   --mean-cell-numbers 1 \
   --number-of-processors 1
python ${script_loc}CytoSpace_clean.py $tmp $output $type
rm -rf "${tmp}"


type="circle"
lab='trueLab'
tmp="${cyto_input}${SLURM_ARRAY_TASK_ID}_${lab}_${type}/"
mkdir -p $tmp
python ${script_loc}CytoSpace_bench.py $SLURM_ARRAY_TASK_ID $type $tmp $input

cytospace --single-cell \
   --scRNA-path "${tmp}scRNA_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt"\
   --cell-type-path "${tmp}scLabels_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --st-path "${tmp}stRNA_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --coordinates-path "${tmp}stCoord_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --st-cell-type-path "${tmp}stLabels_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   -o $tmp \
   --number-of-selected-spots 5000 \
   --mean-cell-numbers 1 \
   --number-of-processors 1
python ${script_loc}CytoSpace_clean.py $tmp $output $type
rm -rf "${tmp}"


type="layered"
lab="trueLab"
tmp="${cyto_input}${SLURM_ARRAY_TASK_ID}_${lab}_${type}/"
mkdir -p $tmp
python ${script_loc}CytoSpace_bench.py $SLURM_ARRAY_TASK_ID $type $tmp $input

cytospace --single-cell \
   --scRNA-path "${tmp}scRNA_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt"\
   --cell-type-path "${tmp}scLabels_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --st-path "${tmp}stRNA_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --coordinates-path "${tmp}stCoord_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --st-cell-type-path "${tmp}stLabels_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   -o $tmp \
   --number-of-selected-spots 5000 \
   --mean-cell-numbers 1 \
   --number-of-processors 1
python ${script_loc}CytoSpace_clean.py $tmp $output $type
rm -rf "${tmp}"

type="dropped"
lab="trueLab"
tmp="${cyto_input}${SLURM_ARRAY_TASK_ID}_${lab}_${type}/"
mkdir -p $tmp
python ${script_loc}CytoSpace_bench.py $SLURM_ARRAY_TASK_ID $type $tmp $input

cytospace --single-cell \
   --scRNA-path "${tmp}scRNA_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt"\
   --cell-type-path "${tmp}scLabels_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --st-path "${tmp}stRNA_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --coordinates-path "${tmp}stCoord_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   --st-cell-type-path "${tmp}stLabels_${SLURM_ARRAY_TASK_ID}_${lab}_${type}.txt" \
   -o $tmp \
   --number-of-selected-spots 5000 \
   --mean-cell-numbers 1 \
   --number-of-processors 1
python ${script_loc}CytoSpace_clean.py $tmp $output $type
rm -rf "${tmp}"