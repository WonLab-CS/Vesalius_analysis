#!/bin/bash
#SBATCH --job-name=cyto
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=array_%A_%a.out
#SBATCH --error=array_%A_%a.err
#SBATCH --array=1-144%20
#SBATCH --mem-per-cpu=15GB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=<patrick.martin@cshs.org>
#SBATCH -p defq

source /home/martinp4/anaconda3/etc/profile.d/conda.sh
conda activate cytospace
cd /common/martinp4/cytospace/
input='/common/martinp4/benchmarking_out/CytoSpace/report/'
type="circle"
tmp="$input""$SLURM_ARRAY_TASK_ID""_""$type""/"
mkdir -p $tmp
python /common/martinp4/benchmarking_out/CytoSpace/pipeline/CytoSpace_bench.py $SLURM_ARRAY_TASK_ID $type $tmp

cytospace --single-cell \
   --scRNA-path "$tmp""scRNA_""$SLURM_ARRAY_TASK_ID""_""$type"".txt"\
   --cell-type-path "$tmp""scLabels_""$SLURM_ARRAY_TASK_ID""_""$type"".txt" \
   --st-path "$tmp""stRNA_""$SLURM_ARRAY_TASK_ID""_""$type"".txt" \
   --coordinates-path "$tmp""stCoord_""$SLURM_ARRAY_TASK_ID""_""$type"".txt" \
   --st-cell-type-path "$tmp""stLabels_""$SLURM_ARRAY_TASK_ID""_""$type"".txt" \
   -o $tmp \
   --number-of-selected-spots 5000 \
   --mean-cell-numbers 1 \
   --number-of-processors 1
python /common/martinp4/benchmarking_out/CytoSpace/pipeline/CytoSpace_clean.py $tmp $input
rm -rf "$tmp"
# rm -f "$tmp""scRNA_""$SLURM_ARRAY_TASK_ID""_""$type"".txt"
# rm -f "$tmp""scLabels_""$SLURM_ARRAY_TASK_ID""_""$type"".txt"
# rm -f "$tmp""stRNA_""$SLURM_ARRAY_TASK_ID""_""$type"".txt"
# rm -f "$tmp""stCoord_""$SLURM_ARRAY_TASK_ID""_""$type"".txt"
# rm -f "$tmp""stLabels_""$SLURM_ARRAY_TASK_ID""_""$type"".txt"
# rm -f "$tmp""CytoSpace_tmp_query.csv"

type="layered"
tmp="$input""$SLURM_ARRAY_TASK_ID""_""$type""/"
mkdir -p $tmp
python /common/martinp4/benchmarking_out/CytoSpace/pipeline/CytoSpace_bench.py $SLURM_ARRAY_TASK_ID $type $tmp

cytospace --single-cell \
   --scRNA-path "$tmp""scRNA_""$SLURM_ARRAY_TASK_ID""_""$type"".txt"\
   --cell-type-path "$tmp""scLabels_""$SLURM_ARRAY_TASK_ID""_""$type"".txt" \
   --st-path "$tmp""stRNA_""$SLURM_ARRAY_TASK_ID""_""$type"".txt" \
   --coordinates-path "$tmp""stCoord_""$SLURM_ARRAY_TASK_ID""_""$type"".txt" \
   --st-cell-type-path "$tmp""stLabels_""$SLURM_ARRAY_TASK_ID""_""$type"".txt" \
   -o $tmp \
   --number-of-selected-spots 5000 \
   --mean-cell-numbers 1 \
   --number-of-processors 1
python /common/martinp4/benchmarking_out/CytoSpace/pipeline/CytoSpace_clean.py $tmp $input
rm -rf "$tmp"
# rm -f "$tmp""scRNA_""$SLURM_ARRAY_TASK_ID""_""$type"".txt"
# rm -f "$tmp""scLabels_""$SLURM_ARRAY_TASK_ID""_""$type"".txt"
# rm -f "$tmp""stRNA_""$SLURM_ARRAY_TASK_ID""_""$type"".txt"
# rm -f "$tmp""stCoord_""$SLURM_ARRAY_TASK_ID""_""$type"".txt"
# rm -f "$tmp""stLabels_""$SLURM_ARRAY_TASK_ID""_""$type"".txt"
# rm -f "$tmp""CytoSpace_tmp_query.csv"
conda deactivate