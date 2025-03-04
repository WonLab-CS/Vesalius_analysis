#!/bin/bash
#SBATCH --job-name=SLAT_perf
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=array_%A_%a.out
#SBATCH --error=array_%A_%a.err
#SBATCH --array=1-5%2
#SBATCH --mem-per-cpu=500GB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=<patrick.martin@cshs.org>
#SBATCH -p defq



data_input="/common/wonklab/synthetic_spatial/"
script_loc="/home/martinp4/common/Vesalius_analysis/benchmarking/SLAT_bench/"
output="/common/martinp4/benchmarking_out/SLAT/report/"
tool="SLAT"
numeric_values=()


output_csv=${output}${tool}"_comp_performance_metrics.csv"
echo "Creating Output file"
if [[ ! -f "$output_csv" ]]; then
  echo "Counter,TaskID,Method,NumericValue,MaxMemoryGB,RuntimeSeconds" > "$output_csv"
fi

echo "Get Comp Performance files"
# Loop over the files directly
for file in $(find "$data_input" -type f -name "computational_performance_*_cells_*"); do
  numeric_value=$(echo "$file" | sed -E 's/.*_([0-9eE\+\.-]+)_cells.*/\1/')
  
  # Check if numeric_value is not empty before adding it to the array
  if [[ -n "$numeric_value" ]]; then
    numeric_values+=("$numeric_value")
  fi
done

unique_values=($(printf '%s\n' "${numeric_values[@]}" | sort -g | uniq))

counter=0
source /home/martinp4/anaconda3/etc/profile.d/conda.sh
conda activate slat
for value in "${unique_values[@]}"; do
  echo "Processing files with numeric value: $value"

  start_time=$(date +%s)
  MAX_MEM=0
  TIMEOUT=43200  # 12 hours in seconds
  MEM_LIMIT=$((480 * 1024 * 1024))  # 480GB in KB
  runtime=""
  MAX_MEM_GB=""

  python ${script_loc}SLAT_bench.py --slurm_id ${SLURM_ARRAY_TASK_ID} --type "computational_performance_${value}" --input $data_input --output $output &
  PID=$!

  while kill -0 "$PID" 2>/dev/null; do
    current_time=$(date +%s)
    elapsed_time=$((current_time - start_time))

    if [[ $elapsed_time -ge $TIMEOUT ]]; then
      echo "Job exceeded time limit of 12 hours. Killing process $PID."
      kill -9 "$PID"
      runtime="Exceed_limit"
      break
    fi

    MEM=$(ps -o rss= -p "$PID" | awk '{print $1}')
    
    if [[ -n "$MEM" && "$MEM" -gt "$MAX_MEM" ]]; then
      MAX_MEM=$MEM
    fi

    if [[ "$MAX_MEM" -ge "$MEM_LIMIT" ]]; then
      echo "Job exceeded memory limit of 480GB. Killing process $PID."
      kill -9 "$PID"
      runtime="Exceed_limit"
      MAX_MEM_GB="Exceed_limit"
      break
    fi

    sleep 1
  done

  if [[ -z "$runtime" ]]; then
    end_time=$(date +%s)
    runtime=$((end_time - start_time))
  fi

  if [[ -z "$MAX_MEM_GB" ]]; then
    MAX_MEM_GB=$((MAX_MEM / 1048576))
  fi

  ((counter++))

  echo "$counter,$SLURM_JOB_ID,$tool,$value,$MAX_MEM_GB,$runtime" >> "$output_csv"
done

conda deactivate
