#!/bin/bash
#SBATCH --job-name=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=array_%A_%a.out
#SBATCH --error=array_%A_%a.err
#SBATCH --array=1-1%1
#SBATCH --mem-per-cpu=1GB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=<patrick.martin@cshs.org>
#SBATCH -p defq

job_id="paste"  
tool_path="/common/martinp4/benchmarking_out/"
new_tool=("SLAT","Tangram","Vesalius","precast")
# Function to submit a new job
submit_new_job() {
    new_job=$(ls "$tool_path""$new_tool""/pipeline" | grep ".sh")
    sbatch "$tool_path""$new_tool""/pipeline/""$new_job"
}

for i in "${new_tool}"; do
    while true; do
        jobs=$(squeue -u martinp4 | grep "$job_id" | wc -l)
        if [ "$jobs" -eq 0 ]; then
            echo "Job $job_id no longer exists. Submitting a new job..."
            submit_new_job "$tool_path" $i
            job_id="$i"
            echo "New tool $job_id"
            break  
        else
            echo "Job $job_id is still running. Checking again in 120 seconds..."
        fi
        sleep 120 
    done
done
