#!/bin/bash
#SBATCH --job-name=Vesalius_perf_combi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=array_%A_%a.out
#SBATCH --error=array_%A_%a.err
#SBATCH --array=1-5%5
#SBATCH --mem-per-cpu=50GB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=<patrick.martin@cshs.org>
#SBATCH -p defq



data_input="/common/wonklab/synthetic_spatial/"
script_loc="/home/martinp4/common/Vesalius_analysis/benchmarking/Vesalius_bench/"
output="/common/martinp4/benchmarking_out/Vesalius/report/"
tool="Vesalius"
numeric_values=()


output_csv=${output}${tool}"_comp_performance_combinations.csv"
echo "Creating Output file"
if [[ ! -f "$output_csv" ]]; then
  echo "TaskID,Combination,Method,NumericValue,MaxMemoryGB,RuntimeSeconds" > "$output_csv"
fi
cell_number=(5000 10000)
combi=("f" "n" "t" "c" "fn" "fny" "fc" "ft" "nt" "nc" "fnt" "fnc" "fnct" "fncty")
counter=1
module load R/4.2.1 rlibs/4.2.1
for value in "${cell_number[@]}"; do
  echo "Processing files with numeric value: $value"
    for co in "${combi[@]}" ; do
        start_time=$(date +%s)
        MAX_MEM=0
        TIMEOUT=43200  # 12 hours in seconds
        MEM_LIMIT=$((480 * 1024 * 1024))  # 480GB in KB
        runtime=""
        MAX_MEM_GB=""

        Rscript ${script_loc}Vesalius_bench.r $SLURM_ARRAY_TASK_ID "computational_performance_${value}" $data_input $output $counter &
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

    echo "$SLURM_JOB_ID,$co,$tool,$value,$MAX_MEM_GB,$runtime" >> "$output_csv"
    done
    counter=1
done
