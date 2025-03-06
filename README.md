# Vesalius analysis 
This repository contains the code used to generate all the plots in the _Vesalius 2.0.0_ manuscript. 

[Pre-print](https://www.biorxiv.org/content/10.1101/2024.08.31.610638v2)

# Disclaimer
All the analysis was carried on a HPC unit running the [Slurm](https://slurm.schedmd.com/documentation.html) workload manager. 

We carried out the analyis by submitting jobs as batches. The shell scripts show the batch submission process and the analysis pipeline used in each case. 

Note that in some cases, we create intermediate files which are used later. For instance, during benchmarking, the output of each tool is standardized and then aggregated for scoring and plotting. 

Input directories contain (with exceptions - explictely mentioned) the data in the same form as downloaded from online repositories. 


# Analysis - Overview

## ARTISTA

This directory contains the analysis related to [Axolotl Brain Regeneration](https://www.science.org/doi/10.1126/science.abp9444)

The data can be downloaded from [the STOmics data collection](https://db.cngb.org/stomics/artista/)

* ARISTA_regen => Code used for processing ARTISTA data and mapping between samples.
* ARTISTA_plot => Code used to generate mapping plots.
* ARTISTA_regen_15_to_20_DPI => Analysis specific to the mapping of 15DPI to 20DPI data sets.

Bash files are files used during submission.

## Benchmarking

For this manuscript, we benchmarked our methods against 6 other tools (GPSA, PASTE, SLAT, CytoSpace, Tangram, Scanorama) in synthetic and real data sets.

Below, we present how the synthetic data was generated and how real data was formatted.

### Generating Synthetic Data

For all synthetic data sets, we used the [oneiric](https://github.com/WonLab-CS/oneiric) package. The package contains a dedicated function to generate all synthetic data used in our analysis. The oneric diretory in this repository shows how the data sets were created and plots the output. 

NOTE: To generate the same data sets, please DO NO change the seed that is provided. 

### Spatial Data formatting 

For the real data, we refromatted real spatial transcriptomics data to ensure a consitent input to our benchmarking code. Specifically, we reformatted:

* [seqFISH Mouse embryo](https://content.cruk.cam.ac.uk/jmlab/SpatialMouseAtlas2020/)
* [Slide-seq V2 mouse hippocampus](https://singlecell.broadinstitute.org/single_cell/study/SCP815/sensitive-spatial-genome-wide-expression-profiling-at-cellular-resolution#study-summary)
* [Stereo-seq mouse embryo](https://db.cngb.org/stomics/mosta/)

Each sub-directory shows the reformatting procedure for all data sets. The reformated data is saved in the same directory as the original data and is the expected input of the benchmarking. The formatting consists of splitting data when needed, adjusting coordinates (fix to origin), adding cell type labels (includes deconvolution when required), and cell context (added useing oneiric)

The original publication links are the following:

* [seqFISH)](https://www.nature.com/articles/s41587-021-01006-2)
* [Slide-seq V2](https://www.nature.com/articles/s41587-020-0739-1)
* [Mouse embryonic development](https://www.sciencedirect.com/science/article/pii/S0092867422003993?via%3Dihub)




### Benchmarking

This directory contains all the code related to the benchmarking across all tools in both synthetic and real data. We also include a sub-directory containing the code related to result aggregation, scoring, and plotting.

In brief, each tool use 2 main analysis scripts (except CytoSpace which requires further reformatting) with the appropriate extensions (`.r` or .`py`):

1. `{tool}_{bench}`
2. `{tool}_{bio_spa}`

The first script is use during the mapping of synthetic data sets while the second is used during the mapping of real biological data. 

We used bash scripts to call these analysis files with the appropriate arguments. If you wish to re-run the analysis, please update the bash scripts with the appropriate path to directories, make sure synthetic data is available, and real data has been properly formatted (see above). As noted in the disclaimer, we used a SLURM engine nomenclature. This can removed or replaced with which ever heading suits your needs. 

The bash scripts allow:

1. Synthetic data benchmarking
2. Computational performance benchmarking
3. seqFISH benchmarking
4. Slide-seq (ssv2) benchmarking
5. Stereo-seq (stereo) benchmarking

### Aggreating, Scoring, and Plotting

We aggregate, score and plot all benchmarking results using a set of bash submission files which will call the appropriate R scripts to perform each task. The general pipeline is following:

1. Unify mapping scores
2. Plot Mapping scores
3. Plot Mapping Event

For synthetic data, there is an additional step:

4. Plot computational performance

For real data, there is an addition step:

4. Plot Contribution scores (Vesalius Only)

### Big Batch

We also provide a **big batch** script which will sequentially submit benchmarking tasks across all tools and tasks. 

NOTE: this does require all other scripts to have been updated accordingly.


## Cancer

This directory contains the analysis related to mapping cells across tumor samples in protate cancer (Slide-seq v2). For simplicity, the `r` script called by the bash script performs an end-to-end analysis of samples.

Data can be downloaded on the [GitHub](https://github.com/shenglinmei/ProstateCancerAnalysis) provided by the authors. The original publications is available [here](https://pmc.ncbi.nlm.nih.gov/articles/PMC9905093/#Abs1)


## IMC

This directory contains the analysis related to the [_in situ_ Mass Cytometry data](https://www.nature.com/articles/s41588-022-01041-y)

* build_data => data selection, filtering, and pre-processing.
* RAZA_balence => mapping of samples between each other and export results. 
* RAZA_heat => plotting and clustering of mapping results. 

## MOSTA

This directory contains the analysis related to [Mouse embryonic development](https://www.sciencedirect.com/science/article/pii/S0092867422003993?via%3Dihub)

The data can be downloaded from [the STOmics data collection](https://db.cngb.org/stomics/mosta/)

* MOSTA => pre-processing and mapping of embryo data forward in time.
* MOSTA_plot => plotting mapping results .
* MOSTA_cluster => clustering of mapped cells and DEG analysis .

## Tech to Tech
This directory contains the analysis related to cross technology and cross resolution mapping. For this analysis we used data from the following sources:

* [seqFISH](https://content.cruk.cam.ac.uk/jmlab/SpatialMouseAtlas2020/)
* [Stereo-seq](https://db.cngb.org/stomics/mosta/)
* [Visium](https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Mouse_Brain_Rep1/CytAssist_FFPE_Mouse_Brain_Rep1_web_summary.html)
* [VisiumHD](https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-mouse-brain-he)
* [scRNA brain reference](https://www.sciencedirect.com/science/article/pii/S0092867418309553?via%3Dihub)

The "stos" (**s**eqFISH to **S**tereo-seq) sub-directory contains:

* stos => pre-processing and mapping
* stos_plot => plotting of mapping results

The "vtov" (**V**isiumHD to **V**isium) sub-directory contains:

* vtov_annot => RCTD annotation of data sets and pre-processing 
* vtov => cross resolution mapping
* vtov_plot => plotting results
