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
#plan(multicore, workers = 2)
max_size <- 100000 * 1024^2
options(future.globals.maxSize = max_size)

#-----------------------------------------------------------------------------#
# directories
#-----------------------------------------------------------------------------#
args <- commandArgs(TRUE)
idx <- as.numeric(args[1])
input <- args[2]
output_plots <- output_data <- args[3]

#-----------------------------------------------------------------------------#
# load and prepare data
#-----------------------------------------------------------------------------#
h5 <- read_h5ad(input)
counts <- t(h5$X)
annot <- as.character(h5$obs$annotation)

if (!file.exists(paste0(output_data,"labels.csv"))){
    types <- unique(annot)
    n_colors <- length(types)
    base_colours <- c(
      "#E69F00",
      "#56B4E9",
      "#009E73",
      "#F0E442",
      "#0072B2",
      "#D55E00",
      "#CC79A7",
      "#999999")
    pal <- colorRampPalette(base_colours)
    colors <- pal(n_colors)[sample(seq(1, n_colors),size = n_colors)]
    labels <- data.frame("labels" = types, "colors" = colors)
    write.csv(labels, file = paste0(output_data,"labels.csv"),row.names = FALSE)
}


coord <- cbind(rownames(h5$obs), as.data.frame(h5$obsm$spatial), annot)
names(annot) <- rownames(h5$obs)
colnames(coord) <- c("barcodes", "x", "y", "trial")

stages <- sort(unique(h5$obs$timepoint))
stage_seed <- stages[idx + 1]
stage_query <- stages[idx]



locs_seed <- which(h5$obs$timepoint == stage_seed)
locs_seed <- sample(locs_seed, size = min(c(length(locs_seed), 50000)))
seed_coord <- coord[locs_seed, ]
seed_counts <- as(counts[, locs_seed], "CsparseMatrix")



locs_query <- which(h5$obs$timepoint == stage_query)
locs_query <- sample(locs_query, size = min(c(length(locs_query), 50000)))
query_coord <- coord[locs_query, ]
query_counts <- as(counts[, locs_query], "CsparseMatrix")




use_cost <- c("feature", "niche", "territory")
#-----------------------------------------------------------------------------#
# Create embeddings 
#-----------------------------------------------------------------------------#
# REF

if (file.exists(paste0(output_data, stage_seed, "_stereo.rds"))) {
    seed <- readRDS(paste0(output_data, stage_seed, "_stereo.rds"))
} else {
    seed <- build_vesalius_assay(coordinates = seed_coord,
    counts = seed_counts)

    seed_annot <- seed_coord$trial
    names(seed_annot) <- seed_coord$barcodes
    
    seed <- add_cells(seed, seed_annot, add_name = "Cells")

    seed <- generate_embeddings(seed, normalization = "none", tensor_resolution = 1, filter_threshold = 1, filter_grid =1) %>%
    equalize_image(embedding = "PCA", dimensions = 1:20, sleft = 2, sright = 2) %>%
    smooth_image(dimensions = 1:20, method = c("iso", "box"), box = 10, sigma = 1, iter = 5) %>%
    segment_image(dimensions = 1:20, method = "kmeans", col_resolution = 20) %>%
    isolate_territories()


    file_out <- paste0(output_data, stage_seed, "_stereo.rds")
    saveRDS(seed, file = file_out)
}


# QUERY
if (file.exists(paste0(output_data, stage_query, "_stereo.rds"))) {
    query <- readRDS(paste0(output_data, stage_query, "_stereo.rds"))
} else {
    query <- build_vesalius_assay(coordinates = query_coord,
        counts = query_counts)

    query_annot <- query_coord$trial
    names(query_annot) <- query_coord$barcodes
    query <- add_cells(query, query_annot, add_name = "Cells")

    query <- generate_embeddings(query, normalization = "none",tensor_resolution = 1, filter_threshold = 1, filter_grid =1) %>%
    equalize_image(embedding = "PCA", dimensions = 1:20, sleft = 2, sright = 2) %>%
    smooth_image(dimensions = 1:20, method = c("iso", "box"), box = 10, sigma = 1, iter = 5) %>%
    segment_image(dimensions = 1:20, method = "kmeans", col_resolution = 20) %>%
    isolate_territories()
    file_out <- paste0(output_data, stage_query, "_stereo.rds")
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
    depth = 3,
    threshold = 0,
    batch_size = 10000,
    epochs = 10,
    use_norm = "raw",
    jitter = 0)

file_out <- paste0(output_data, stage_query,"_to_",stage_seed, "_stereo.rds")
saveRDS(matched, file = file_out)
