#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --job-name=HTAN_simple_map
#SBATCH --mem-per-cpu=50GB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=<patrick.martin@cshs.org>
#SBATCH -p defq

module load R/4.2.1 rlibs/4.2.1
counts="/common/wonklab/SSv2/slide.seq.raw.counts.rds"
coordinates="/common/wonklab/SSv2/slide.seq.BeadLocations.rds"
annotations="/common/wonklab/SSv2/slide.seq.ano.rds"
output="/common/martinp4/benchmarking_out/Vesalius/report/"
ref_tag="Tumor01"
query_tag="Tumor02"
cores="10"
script_loc="/common/martinp4/Vesalius_analysis/HTAN/"
Rscript ${script_loc}htan_simple_map.r $counts $coordinates $annotations $output $ref_tag $query_tag $cores


