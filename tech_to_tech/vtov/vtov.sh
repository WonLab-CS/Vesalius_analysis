#!/bin/bash
#SBATCH --job-name=VtoV
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=400GB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=<patrick.martin@cshs.org>
#SBATCH -p defq

module load R/4.2.1 rlibs/4.2.1
script_loc="/home/martinp4/common/Vesalius_analysis/tech_to_tech/vtov/"
output="/common/martinp4/stos/report/vtov/"
input_ref="/common/wonklab/mouse_hippocampus_scRNAseq/"
input_viz="/common/wonklab/visium_brain/"
input_hd="/common/wonklab/VisiumHD/Mouse_brain/square_008um/"

Rscript ${script_loc}vtov_annot_visium.r $input_viz $input_ref $output
Rscript ${script_loc}vtov_annot_visiumHD.r $input_hd $input_ref $output
Rscript ${script_loc}vtov.r $output
Rscript ${script_loc}vtov_plot.r $output