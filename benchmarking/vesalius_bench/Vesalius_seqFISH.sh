#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --job-name=Vesalius_seqFISH
#SBATCH --mem-per-cpu=80GB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=<patrick.martin@cshs.org>
#SBATCH -p defq

module load R/4.2.1 rlibs/4.2.1
script_loc="/home/martinp4/common/Vesalius_analysis/benchmarking/Vesalius_bench/"
data_type="seqFISH"
input="/common/wonklab/seqFISH/"
output="/common/martinp4/benchmarking_out/Vesalius/report/"
seed_tag="embryo1"
query_tag="embryo3"
for i in {1..14}; do
    Rscript ${script_loc}Vesalius_bio_spa.r $data_type $input $output $seed_tag $query_tag $i
done

