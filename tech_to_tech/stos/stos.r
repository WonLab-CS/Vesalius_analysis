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

plan(multicore, workers = 2)
max_size <- 10000 * 1024^2
options(future.globals.maxSize = max_size)
#-----------------------------------------------------------------------------#
# Set future global for multicore processing
#-----------------------------------------------------------------------------#

if (!dir.exists("/common/martinp4/stos/report/stos/")) {
    dir.create("/common/martinp4/stos/report/stos/")
}
output <- "/common/martinp4/stos/report/stos/"
args <- commandArgs(TRUE)
idx <- as.numeric(args[1])

use_cost <- c("feature", "niche")
slices <- c( "embryo1","embryo2","embryo3")
#-----------------------------------------------------------------------------#
# Loading seqFISH
#-----------------------------------------------------------------------------#
coord <- list.files(path = "/common/wonklab/seqFISH/",
    pattern = "metadata", full.names = TRUE)
coord <- readRDS(coord)

coord <- coord[, c("uniqueID",
    "x_global",
    "y_global",
    "embryo",
    "celltype_mapped_refined")]
colnames(coord) <- gsub("x_global", "x", colnames(coord))
colnames(coord) <- gsub("y_global", "y", colnames(coord))
colnames(coord) <- gsub("uniqueID", "barcodes", colnames(coord))
colnames(coord) <- gsub("celltype_mapped_refined", "trial", colnames(coord))

coord <- coord[coord$trial != "Low quality", ]

counts <- list.files(path = "/common/wonklab/seqFISH/",
    pattern = "counts.Rds", full.names = TRUE)
counts <- readRDS(counts)

counts <- counts[, colnames(counts) %in% coord$barcodes]
nfeatures <- nrow(counts)

query_coord <- coord[grep(slices[idx], coord$barcodes), ]
query_counts <- counts[, grep(slices[idx], coord$barcodes)]
query_cells <- as.character(query_coord$trial)
names(query_cells) <- as.character(query_coord$barcodes)

query <- build_vesalius_assay(query_coord, query_counts) %>%
    generate_embeddings(filter_threshold = 1,
        filter_grid = 1,
        tensor_resolution = 1,
        nfeatures = nfeatures)
query <- add_cells(query, query_cells)
saveRDS(query, file = paste0(output,"seqFISH_",slices[idx],".rds"))
#-----------------------------------------------------------------------------#
# load and prepare data
#-----------------------------------------------------------------------------#
h5 <- "/common/wonklab/Stereo_seq/Mouse_embryo_all_stage.h5ad"
stage_query <- "E9.5"
h5 <- read_h5ad(h5)
counts <- t(h5$X)
annot <- as.character(h5$obs$annotation)
coord <- cbind(rownames(h5$obs), as.data.frame(h5$obsm$spatial), annot)
names(annot) <- rownames(h5$obs)
colnames(coord) <- c("barcodes", "x", "y", "trial")

locs_query <- which(h5$obs$timepoint == stage_query)
locs_query <- sample(locs_query, size = min(c(length(locs_query), 50000)))
seed_coord <- coord[locs_query, ]
seed_counts <- as(counts[, locs_query], "CsparseMatrix")
seed <- build_vesalius_assay(coordinates = seed_coord,
        counts = seed_counts)

seed_annot <- seed_coord$trial
names(seed_annot) <- seed_coord$barcodes
seed <- add_cells(seed, seed_annot, add_name = "Cells")

seed <- generate_embeddings(seed,
    tensor_resolution = 1,
    filter_threshold = 1,
    filter_grid = 1)
saveRDS(seed, file = paste0(output,"stereo_seq.rds"))
#-----------------------------------------------------------------------------#
# Integrated
#-----------------------------------------------------------------------------#
matched <- map_assays(seed,
    query,
    signal = "all_features",
    use_cost = use_cost,
    neighborhood = "graph",
    depth = 2,
    threshold = 0,
    batch_size = 10000,
    epochs = 10,
    use_norm = "log_norm",
    jitter = FALSE)

saveRDS(matched, file = paste0(output, "seqFISH_",slices[idx],"_to_stereo_matched.rds"))

matched <- map_assays(query,
    seed,
    signal = "all_features",
    use_cost = use_cost,
    neighborhood = "graph",
    depth = 2,
    threshold = 0,
    batch_size = 10000,
    epochs = 10,
    use_norm = "log_norm",
    jitter = FALSE)
saveRDS(matched, file = paste0(output, "stereo_to_seqFISH_",slices[idx],"_reverse.rds"))
