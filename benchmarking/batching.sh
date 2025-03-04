#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=15-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --job-name=BENCHMAIN
#SBATCH --mem-per-cpu=1GB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=<patrick.martin@cshs.org>
#SBATCH -p defq

script_loc="/home/martinp4/common/Vesalius_analysis/benchmarking/"

check_job_running() {
    local job_name="$1"
    
    # Ensure a job name is provided
    if [ -z "$job_name" ]; then
        echo "Error: Please provide a job name."
        return 2
    fi

    while true; do
        # Check if the job with the given name is running under username martinp4
        local running_jobs
        running_jobs=$(squeue -u martinp4 -o "%.25j" | grep -w "$job_name")
        
        if [ -n "$running_jobs" ]; then
            sleep 120
        else
            # Job is no longer running
            echo "Job '$job_name' is not running."
            return 0
        fi
    done
}

# Full list
# tool_list=("Vesalius" "SLAT" "Tangram" "PASTE" "CytoSpace" "Scanorama" "GPSA")
# task_list=("ssv2" "seqFISH" "stereo" "bench" "perf")
# Buffer List - in case of errors and needing to rerun
tool_list=("Vesalius" "SLAT" "Tangram" "PASTE" "CytoSpace" "Scanorama" "GPSA")
task_list=("perf")

#Task exclusion - just in case we want to skip some runs
exclusion_list=("Tangram_bench" "Tangram_ssv2" "Tangram_seqFISH" "Tangram_stereo" \ 
                "PASTE_bench" "PASTE_ssv2" "PASTE_stereo" "PASTE_seqFISH" \ 
                "GPSA_ssv2" "GPSA_stereo" "GPSA_seqFISH" "GPSA_bench" \ 
                "SLAT_bench" "SLAT_ssv2" "SLAT_seqFISH" "SLAT_stereo" \ 
                "CytoSpace_bench" "CytoSpace_ssv2" "CytoSpace_seqFISH" "CytoSpace_stereo" \ 
                "Vesalius_bench" "Vesalius_ssv2" "Vesalius_seqFISH" "Vesalius_stereo" \ 
                "Scanorama_bench" "Scanorama_ssv2" "Scanorama_seqFISH" "Scanorama_stereo")

# Loop over tasks and tools
for task in "${task_list[@]}"; do
    for tool in "${tool_list[@]}"; do
        # Change directory for the current tool
        tool_dir="${script_loc}${tool}_bench"
        if [ -d "$tool_dir" ]; then
            echo "Changing to directory: $tool_dir"
            cd "$tool_dir" || { echo "Failed to change directory to $tool_dir"; exit 1; }
        else
            echo "Directory $tool_dir does not exist. Skipping tool $tool."
            continue
        fi

        job_name="${tool}_${task}"
        echo "Checking exclusion list"
        if [[ " ${exclusion_list[@]} " =~ " $job_name " ]]; then
            echo "Skipping $job_name"
            continue
        fi
        echo "Submitting ${job_name}"

        # Submit the job
        if sbatch "${job_name}.sh"; then
            echo "Job '${job_name}' submitted successfully."
        else
            echo "Failed to submit job '${job_name}'. Skipping."
            continue
        fi

        # Check if the job is running
        check_job_running "${job_name}"
        if [ $? -eq 0 ]; then
            echo "Job '$job_name' completed."
            continue
        fi
    done
    cd "$script_loc" || exit 1
done
