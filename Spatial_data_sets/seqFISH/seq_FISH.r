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
library(mcclust, lib.loc = "/common/martinp4/R")
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
args <- commandArgs(TRUE)
idx <- as.numeric(args[1])

if (!dir.exists("/common/martinp4/benchmarking_out/Vesalius/bio_data/seqFISH/")) {
    dir.create("/common/martinp4/benchmarking_out/Vesalius/bio_data/seqFISH/")
}
output <- "/common/martinp4/benchmarking_out/Vesalius/bio_data/seqFISH/"

file_name <- paste0(output, "combination_score_seqFISH.csv")
if (!file.exists(file_name)){
    score <- data.frame("Combination", "ARI", "VI")
    write.table(score,
        file = file_name,
        quote = FALSE,
        col.names = TRUE,
        row.names = FALSE,
        sep = ",")
}
seed_tag <- "embryo3"
query_tag <- "embryo1"
#-----------------------------------------------------------------------------#
# Loading seqFISH
#-----------------------------------------------------------------------------#
if (file.exists(paste0(output,"seqFISH_",seed_tag,"_seed.rds")) &
    file.exists(paste0(output,"seqFISH_",query_tag,"_query.rds"))) {
    seed <- readRDS(paste0(output,"seqFISH_",seed_tag,"_seed.rds"))
    query <- readRDS(paste0(output,"seqFISH_",query_tag,"_query.rds"))
} else {
    # original data
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
    # create seed
    
    seed_coord <- coord[grep(seed_tag, coord$barcodes), ]
    seed_counts <- counts[, grep(seed_tag, coord$barcodes)]
    seed_cells <- as.character(seed_coord$trial)
    names(seed_cells) <- as.character(seed_coord$barcodes)

    seed <- build_vesalius_assay(seed_coord, seed_counts) %>%
        generate_embeddings(filter_threshold = 1,
            filter_grid = 1,
            tensor_resolution = 1,
            nfeatures = nfeatures)
    seed <- add_cells(seed, seed_cells)
    saveRDS(seed, file  = paste0(output,"seqFISH_",seed_tag,"_seed.rds"))
    # create query
    
    query_coord <- coord[grep(query_tag, coord$barcodes), ]
    query_counts <- counts[, grep(query_tag, coord$barcodes)]
    query_cells <- as.character(query_coord$trial)
    names(query_cells) <- as.character(query_coord$barcodes)
    query <- build_vesalius_assay(query_coord, query_counts) %>%
        generate_embeddings(filter_threshold = 1,
            filter_grid = 1,
            tensor_resolution = 1,
            nfeatures = nfeatures) 
    query <- add_cells(query, query_cells)
    saveRDS(query, file  = paste0(output,"seqFISH_",query_tag,"_query.rds"))
}

#-----------------------------------------------------------------------------#
# Align 
#-----------------------------------------------------------------------------#
use_cost <- list(
    c("feature"),
    c("niche"),
    c("composition"),
    c("cell_type"),
    c("feature", "niche"),
    c("feature","niche","composition"),
    c("feature","niche","composition","cell_type"),
    c("feature", "composition"),
    c("niche","composition"))
matched <- map_assays(seed,
    query,
    neighborhood = "knn",
    k = 6,
    threshold = 0.3,
    use_norm = "log_norm",
    signal = "all_features",
    epochs = 10,
    batch_size = 5000,
    use_cost = use_cost[[idx]],
    jitter = FALSE)
file_name <- paste0(output,"seqFISH_",paste0(use_cost[[idx]],collapse = "_"),"_matched.rds")
saveRDS(matched, file  = file_name)

#-----------------------------------------------------------------------------#
# Quick score
#-----------------------------------------------------------------------------#
match_bar <- matched@map
query_cells <- matched@territories$Cells[match(match_bar$from, matched@territories$barcodes)]
ref_cells <- seed@territories$Cells[match(match_bar$to, seed@territories$barcodes)]
ari <- adjustedRandIndex(ref_cells, query_cells)
vi <- vi.dist(ref_cells, query_cells)

combi_score <- data.frame("Combination" = paste0(use_cost[[idx]],collapse = "_"),
    "ARI" = ari,
    "VI" = vi)
write.table(combi_score,
    file = paste0(output, "combination_score_seqFISH.csv"),
    quote = FALSE,
    col.names = FALSE,
    row.names = FALSE,
    sep = ",",
    append  = TRUE)

