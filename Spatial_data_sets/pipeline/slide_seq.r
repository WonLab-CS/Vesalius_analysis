#-----------------------------------------------------------------------------#
# Loading Libraries and setting seed
#-----------------------------------------------------------------------------#
library(dplyr)
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
library(anndata, lib.loc = "/common/martinp4/R")
library(pwr, lib.loc = "/common/martinp4/R")
library(gsignal, lib.loc = "/common/martinp4/R")
library(kohonen, lib.loc = "/common/martinp4/R")
library(registry, lib.loc = "/common/martinp4/R")
library(rngtools, lib.loc = "/common/martinp4/R")
library(NMF, lib.loc = "/common/martinp4/R")
library(RcppHungarian, lib.loc = "/common/martinp4/R")
library(spatstat.utils, lib.loc = "/common/martinp4/R")
library(geometry, lib.loc = "/common/martinp4/R")
library(vesalius, lib.loc = "/common/martinp4/R")
library(RColorBrewer)
set.seed(1547)


#-----------------------------------------------------------------------------#
# Set future global for multicore processing
#-----------------------------------------------------------------------------#
plan(multicore, workers = 10)
max_size <- 10000 * 1024^2
options(future.globals.maxSize = max_size)

if (!dir.exists("/common/martinp4/benchmarking_out/Vesalius/bio_data/SSV2/")) {
    dir.create("/common/martinp4/benchmarking_out/Vesalius/bio_data/SSV2/")
}
output <- "/common/martinp4/benchmarking_out/Vesalius/bio_data/SSV2/"


#-----------------------------------------------------------------------------#
# Getting path to files
#-----------------------------------------------------------------------------#
input <- "/common/wonklab/SSv2"
coordinates <- list.files(path = input,
    pattern = "location|coord|Locations", full.names = TRUE)
counts <- list.files(path = input,
    pattern = "expression", full.names = TRUE)
tag <- list.files(path = input,
    pattern = "expression", full.names = FALSE)
tag <- gsub(".digital_expression.txt.gz|_expression_matrix.mtx.gz|.sparse_expression.txt",
    "", tag)
#-----------------------------------------------------------------------------#
# creating seed dataset
#-----------------------------------------------------------------------------#
f <- grep(pattern = "Puck_200115_08", tag)
seed_coord <- read.csv(coordinates[f], header = FALSE, skip = 1)
colnames(seed_coord) <- c("barcodes", "xcoord", "ycoord")
seed_counts <- read.table(counts[f], header = TRUE, row.names = 1)

seed_counts <- seed_counts[, apply(seed_counts, 2, sum) > 200]
seed_counts <- seed_counts[apply(seed_counts,1, sum) > 100, ]
seed_coord <- seed_coord[seed_coord$barcodes %in% colnames(seed_counts),]


seed <- build_vesalius_assay(seed_coord, seed_counts) %>%
    generate_embeddings(filter_threshold = 1,
        filter_grid = 1,
        tensor_resolution = 1,
        dim_reduction = "PCA",
        nfeatures = 2000) %>%
    equalize_image(dimensions = 1:30, sleft = 5, sright = 5) %>%
    smooth_image(dimensions = 1:30,
        method = c("iso", "box"),
        box = 10,
        sigma = 2,
        iter = 10) %>%
    segment_image(dimensions = 1:30,
        col_resolution = 25) %>%
    isolate_territories()

file_name <- paste0(output, "ssv2_Puck_200115_08_seed.rds")
saveRDS(seed, file = file_name)
#-----------------------------------------------------------------------------#
# creating query 
#-----------------------------------------------------------------------------#
f <- grep(pattern = "Puck_190921_21", tag)
query_coord <- read.csv(coordinates[f], header = FALSE, skip = 1)
colnames(query_coord) <- c("barcodes", "xcoord", "ycoord")
query_counts <- read.table(counts[f], header = TRUE, row.names = 1)
query_counts <- query_counts[, apply(query_counts, 2, sum) > 200]
query_counts <- query_counts[apply(query_counts,1, sum) > 100, ]
query_coord <- query_coord[query_coord$barcodes %in% colnames(query_counts),]


query <- build_vesalius_assay(query_coord, query_counts) %>%
    generate_embeddings(filter_threshold = 1,
        filter_grid = 1,
        tensor_resolution = 1,
        dim_reduction = "PCA",
        nfeatures = 2000) %>%
    equalize_image(dimensions = 1:30, sleft = 2.5, sright = 2.5) %>%
    smooth_image(dimensions = 1:30,
        method = c("iso", "box"),
        box = 10,
        sigma = 1,
        iter = 10) %>%
    segment_image(dimensions = 1:30,
        col_resolution = 20) %>%
    isolate_territories()

file_name <- paste0(output, "ssv2_Puck_190921_21_query.rds")
saveRDS(query, file = file_name)

#-----------------------------------------------------------------------------#
# Align 
#-----------------------------------------------------------------------------#
matched <- map_assays(seed,
    query,
    use_norm = "log_norm",
    k = 30,
    scale = FALSE,
    use_cost = c("feature", "niche", "territory"),
    batch_size = 10000,
    epochs = 10,
    jitter = FALSE)


matched <- generate_embeddings(matched, filter_threshold = 1,
        filter_grid = 1,
        tensor_resolution = 1,
        dim_reduction = "PCA",
        nfeatures = 2000) %>%
    equalize_image(dimensions = 1:30, sleft = 2.5, sright = 2.5) %>%
    smooth_image(dimensions = 1:30,
        method = c("iso", "box"),
        box = 15,
        sigma = 1.5,
        iter = 10) %>%
    segment_image(dimensions = 1:30,
        col_resolution = 20) %>%
    isolate_territories()

file_out <- paste0(output, "ssv2_matched.rds")
saveRDS(matched, file = file_out)
