#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --job-name=Vesalius_ssv2
#SBATCH --mem-per-cpu=80GB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=<patrick.martin@cshs.org>
#SBATCH -p defq

module load R/4.2.1 rlibs/4.2.1
data_type="SSv2"
input="/common/wonklab/SSv2/"
output="/common/martinp4/benchmarking_out/Vesalius/report/"
script_loc="/home/martinp4/common/Vesalius_analysis/benchmarking/Vesalius_bench/"
seed_tag="Puck_200115_08"
query_tag="Puck_190921_21"
for i in {1..14}; do
    Rscript ${script_loc}Vesalius_bio_spa.r $data_type $input $output $seed_tag $query_tag $i
done

