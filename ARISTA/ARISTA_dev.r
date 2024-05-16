#-----------------------------------------------------------------------------#
# Loading Libraries and setting seed
#-----------------------------------------------------------------------------#
library(dplyr)
library(Seurat)
library(future)
library(future.apply)
library(ggplot2)
library(dplyr)
library(patchwork)
library(Matrix)
library(deldir)
library(imager)
library(imagerExtra)
library(Morpho)
library(lpSolve, lib = "/common/martinp4/R")
library(TreeDist, lib.loc = "/common/martinp4/R")
library(mclust, lib.loc = "/common/martinp4/R")
library(pwr, lib.loc = "/common/martinp4/R")
library(kohonen, lib.loc = "/common/martinp4/R")
library(registry, lib.loc = "/common/martinp4/R")
library(rngtools, lib.loc = "/common/martinp4/R")
library(NMF, lib.loc = "/common/martinp4/R")
library(vesalius, lib.loc = "/common/martinp4/R")
library(RColorBrewer)


set.seed(1547)
plan(multicore, workers = 2)
max_size <- 100000 * 1024^2
options(future.globals.maxSize = max_size)

#-----------------------------------------------------------------------------#
# directories
#-----------------------------------------------------------------------------#
args <- commandArgs(TRUE)
idx <- as.numeric(args[1])

if (!dir.exists("/common/wonklab/Stereo_seq_arista/report/")) {
    dir.create("/common/wonklab/Stereo_seq_arista/report/")
}
output_plots <- "/common/wonklab/Stereo_seq_arista/report/"
output_data <- "/common/wonklab/Stereo_seq_arista/report/"


#-----------------------------------------------------------------------------#
# load and prepare data
#-----------------------------------------------------------------------------#
file_loc <- list.files("/common/wonklab/Stereo_seq_arista/",
    pattern = ".rds",
    full.names = TRUE)
file_loc <- grep("_arista_dev|_arista_regen",file_loc, value = TRUE, invert = TRUE)
file_vec <- c(grep("Stage44",file_loc, value = TRUE),
    grep("Stage54",file_loc, value = TRUE),
    grep("Stage57",file_loc, value = TRUE),
    grep("Control_Juv",file_loc, value = TRUE),
    grep("Adult",file_loc, value = TRUE))

seed <- file_vec[idx + 1]
stage_seed <- gsub("/common/wonklab/Stereo_seq_arista//","",seed)
stage_seed <- gsub(".rds","",stage_seed)
query <- file_vec[idx]
stage_query <- gsub("/common/wonklab/Stereo_seq_arista//","",query)
stage_query <- gsub(".rds","",stage_query)

seed <- readRDS(seed)
seed_counts <- GetAssayData(seed, layer = "counts")
seed_coord <- GetTissueCoordinates(seed)
seed_coord <- data.frame("barcodes" = rownames(seed_coord),
    "x" = seed_coord$imagerow,
    "y" = seed_coord$imagecol)
seed_cells <- seed@meta.data$Annotation
names(seed_cells) <- rownames(seed@meta.data)
seed <- build_vesalius_assay(coordinates = seed_coord,
    counts = seed_counts)
seed <- add_cells(seed, cells = seed_cells)

query <- readRDS(query)
query_counts <- GetAssayData(query, layer = "counts")
query_coord <- GetTissueCoordinates(query)
query_coord <- data.frame("barcodes" = rownames(query_coord),
    "x" = query_coord$imagerow,
    "y" = query_coord$imagecol)
query_cells <- query@meta.data$Annotation
names(query_cells) <- rownames(query@meta.data)
query <- build_vesalius_assay(coordinates = query_coord,
    counts = query_counts)
query <- add_cells(query, cells = query_cells)


#use_cost <- c("feature", "niche", "cell_type", "composition", "territory")
use_cost <- c("feature", "niche", "composition", "territory")
#-----------------------------------------------------------------------------#
# Create embeddings 
#-----------------------------------------------------------------------------#
# REF

if (file.exists(paste0(output_data, stage_seed, "_arista_dev.rds"))) {
    seed <- readRDS(paste0(output_data, stage_seed, "_arista_dev.rds"))
} else {
    

    seed <- generate_embeddings(seed, normalization = "none", tensor_resolution = 1, filter_threshold = 1, filter_grid =1) %>%
    equalize_image(embedding = "PCA", dimensions = 1:20, sleft = 2, sright = 2) %>%
    smooth_image(dimensions = 1:20, method = c("iso", "box"), box = 10, sigma = 1, iter = 5) %>%
    segment_image(dimensions = 1:20, method = "kmeans", col_resolution = 20) %>%
    isolate_territories()


    file_out <- paste0(output_data, stage_seed, "_arista_dev.rds")
    saveRDS(seed, file = file_out)
}


# QUERY
if (file.exists(paste0(output_data, stage_query, "_arista_dev.rds"))) {
    query <- readRDS(paste0(output_data, stage_query, "_arista_dev.rds"))
} else {
    
    

    query <- generate_embeddings(query, normalization = "none",tensor_resolution = 1, filter_threshold = 1, filter_grid =1) %>%
    equalize_image(embedding = "PCA", dimensions = 1:20, sleft = 2, sright = 2) %>%
    smooth_image(dimensions = 1:20, method = c("iso", "box"), box = 10, sigma = 1, iter = 5) %>%
    segment_image(dimensions = 1:20, method = "kmeans", col_resolution = 20) %>%
    isolate_territories()
    file_out <- paste0(output_data, stage_query, "_arista_dev.rds")
    saveRDS(query, file = file_out)
}

#-----------------------------------------------------------------------------#
# Integrated
#-----------------------------------------------------------------------------#
matched <- map_assays(seed,
    query,
    signal = "variable_features",
    use_cost = use_cost,
    neighborhood = "graph",
    depth = 2,
    threshold = 0,
    batch_size = 2000,
    epochs = 20,
    use_norm = "raw",
    jitter = FALSE)

file_out <- paste0(output_data, stage_query,"_to_",stage_seed, "_arista_dev.rds")
saveRDS(matched, file = file_out)
