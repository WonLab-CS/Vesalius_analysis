# Vesalius analysis 
This repository contains the code used to generate all the plots in the _Vesalius 2.0.0_ manuscript. 

[ADD Link to pre-print/paper]

# Disclaimer
All the analysis was carried on a HPC unit running the (Slurm)[https://slurm.schedmd.com/documentation.html] workload manager. 

We carried out the analyis by submitting jobs as batches. The shell scripts show the batch submission process. 

In the case where many data sets where used (for benchmarking for instance), we use the `$SLURM_ARRAY_TASK_ID` to define which data set
will be selected for analysis. The selection was carried out in R or python. 
The exported file names will always reflect the data sets that were used. 

The process applied to each data set is shown in the R scripts (or python). The same code is applied to all data sets in that batch. 

For individual runs, we suggest modifying the data set up section of the analysis code to account for your specific needs. 
However, in somce case, you will noticed that we load `.rds` files. We save intermediate files during our workflow and use these files in downsteams tasks. This allows to not need to run everything from scratch if it is not required.

# Analysis - Overview

## ARTISTA

This directory contains the analysis related to (Axolotl Brain Regeneration)[https://www.science.org/doi/10.1126/science.abp9444] 

The data can be downloaded from (the STOmics data collection)[https://db.cngb.org/stomics/artista/]

* ARISTA_regen => Code used for processing ARTISTA data and mapping between samples.
* ARTISTA_plot => Code used to generate mapping plots.
* ARTISTA_regen_15_to_20_DPI => Analysis specific to the mapping of 15DPI to 20DPI data sets.

## Benchmarking

This directory contains the analysis related to the benchmarking of Vesalius on Synthetic data sets. 

To generate synthetic data sets, please check out the dedicated repository (onieric)[https://github.com/WonLab-CS/oneiric]

This directory is organized by tool used for benchmarking. 

* [TOOL_NAME]_bench => code used to load and run mapping or spatial alignement task.

The directory called `plotting and scoring` takes the data produced by each tool and unifies all the results into a single data frame that is used for plottin . 

**NOTEs**

* CytoSpace requires and extra data cleaning step. We decided to run it from the command line directly instead of running it through python as it seems that this approach is the preferred approach according to the authors. 

* Vesalius contains extra `R` scripts related to testing different combinations of cost matrices and plotting the outcomes.

## IMC

This directory contains the analysis related to the (_in situ_ Mass Cytometry data)[https://www.nature.com/articles/s41588-022-01041-y]

* build_data => data selection, filtering, and pre-processing.
* RAZA_balence => mapping of samples between each other and export results. 
* RAZA_heat => plotting and clustering of mapping results. 

## MOSTA

This directory contains the analysis related to (Mouse embryonic development)[https://www.sciencedirect.com/science/article/pii/S0092867422003993?via%3Dihub]

The data can be downloaded from (the STOmics data collection)[https://db.cngb.org/stomics/mosta/]

* MOSTA => pre-processing and mapping of embryo data forward in time.
* MOSTA_plot => plotting mapping results .
* MOSTA_cluster => clustering of mapped cells and DEG analysis .

## Spatial Data sets

This directory contains the analysis related to the (Slide-seq V2)[https://www.nature.com/articles/s41587-020-0739-1] and (seqFISH)[https://www.nature.com/articles/s41587-021-01006-2] data.

The data sets can be downloaed from:

(SSV2)[https://singlecell.broadinstitute.org/single_cell/study/SCP815/sensitive-spatial-genome-wide-expression-profiling-at-cellular-resolution#study-summary]

(seqFISH)[https://content.cruk.cam.ac.uk/jmlab/SpatialMouseAtlas2020/]

The Slide-seq sub-directory contains:

* slide_seq => pre-processing and mapping.
* slide_seq_plot => plotting of results. 

The seqFISH sub-directory contains:

* seq_FISH => pre-processing and mapping (using multiple cost matrix combinations - exported).
* seq_FISH_plot => plotting of mapping results. 
* seq_FISH_combination_plot => quantitative scoring of matrix combinations and plotting.

## Tech to Tech
This directory contains the analysis related to cross technology and cross resolution mapping. For this analysis we used data from the following sources:

* (seqFISH)[https://content.cruk.cam.ac.uk/jmlab/SpatialMouseAtlas2020/]
* (Stereo-seq)[https://db.cngb.org/stomics/mosta/]
* (Visium)[https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Mouse_Brain_Rep1/CytAssist_FFPE_Mouse_Brain_Rep1_web_summary.html]
* (VisiumHD)[https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-mouse-brain-he]
* (scRNA brain reference)[https://www.sciencedirect.com/science/article/pii/S0092867418309553?via%3Dihub]

The "stos" (**s**eqFISH to **S**tereo-seq) sub-directory contains:

* stos => pre-processing and mapping
* stos_plot => plotting of mapping results

The "vtov" (**V**isiumHD to **V**isium) sub-directory contains:

* vtov_annot => RCTD annotation of data sets and pre-processing 
* vtov => cross resolution mapping
* vtov_plot => plotting results
