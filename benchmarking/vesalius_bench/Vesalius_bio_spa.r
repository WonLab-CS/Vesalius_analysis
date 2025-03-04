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
cat("Libs: DONE\n")

#-----------------------------------------------------------------------------#
# Set future global for multicore processing
# Getting arguments 
# Setting use_cost for combinations runs if required 
#-----------------------------------------------------------------------------#
plan(multicore, workers = 12)
max_size <- 10000 * 1024^2
options(future.globals.maxSize = max_size)
args <- commandArgs(TRUE)
data_type <- args[1]
input <- args[2]
output <- args[3]
seed_tag <- args[4]
query_tag <- args[5]
cost <- as.numeric(args[6])

cat(paste(data_type,"\n"))

use_cost <- list(
    c("feature"),
    c("niche"),
    c("territory"),
    c("composition"),
    c("feature", "niche"),
    c("feature", "niche","cell_type"),
    c("feature","composition"),
    c("feature", "territory"),
    c("niche","territory"),
    c("niche","composition"),
    c("feature","niche","territory"),
    c("feature", "niche","composition"),
    c("feature", "niche","composition","territory"),
    c("feature", "niche","composition","territory","cell_type"))
use_cost <- use_cost[[cost]]
compressed_tag <- paste(sapply(gsub("cell_type", "y", use_cost),substr,1,1), collapse = "")



#-----------------------------------------------------------------------------#
# Getting path to files
#-----------------------------------------------------------------------------#
coordinates <- list.files(path = input,
    pattern = paste0(data_type,"_formatted_spatial_coordinates"), full.names = TRUE)
counts <- list.files(path = input,
    pattern = paste0(data_type,"_formatted_gene_counts"), full.names = TRUE)
cat("Set Up: DONE\n")
#-----------------------------------------------------------------------------#
# creating seed 
#-----------------------------------------------------------------------------#
seed_coord <- read.csv(grep(seed_tag, coordinates, value = TRUE))
seed_counts <- read.csv(grep(seed_tag, counts, value = TRUE), row.names = 1)
seed_counts$genes <- NULL 


seed <- build_vesalius_assay(seed_coord, seed_counts, assay = seed_tag)
cell_labels <- seed_coord$cell_labels
names(cell_labels) <- seed_coord$barcodes
seed <- add_cells(seed, cells = cell_labels, add_name = "cell_labels")
interactions <- seed_coord$interactions
names(interactions) <- seed_coord$barcodes
seed <- add_cells(seed, cells = interactions, add_name = "interactions")
seed <- generate_embeddings(seed,
    filter_threshold = 1,
    filter_grid = 1,
    tensor_resolution = 1,
    dim_reduction = "PCA",
    nfeatures = 2000)

if ("territory" %in% use_cost){
    seed <- seed %>%
        smooth_image(dimensions = 1:30,
            method = c("iso", "box"),
            box = 15,
            sigma = 1.5,
            iter = 10) %>%
        equalize_image(dimensions = 1:30, sleft = 2.5, sright = 2.5) %>%
        segment_image(dimensions = 1:30,
            col_resolution = 20) %>%
        isolate_territories()
}
cat("Ref: DONE\n")
#-----------------------------------------------------------------------------#
# creating query 
#-----------------------------------------------------------------------------#
query_coord <- read.csv(grep(query_tag, coordinates, value = TRUE))
query_counts <- read.csv(grep(query_tag, counts, value = TRUE), row.names = 1)
query_counts$genes <- NULL


query <- build_vesalius_assay(query_coord, query_counts, assay = query_tag)
cell_labels <- query_coord$cell_labels
names(cell_labels) <- query_coord$barcodes
query <- add_cells(query, cells = cell_labels, add_name = "cell_labels")
interactions <- query_coord$interactions
names(interactions) <- query_coord$barcodes
query <- add_cells(query, cells = interactions, add_name = "interactions")
query <- generate_embeddings(query,
    filter_threshold = 1,
    filter_grid = 1,
    tensor_resolution = 1,
    dim_reduction = "PCA",
    nfeatures = 2000) 

if ("territory" %in% use_cost){
    query <- query %>%
        smooth_image(dimensions = 1:30,
            method = c("iso", "box"),
            box = 15,
            sigma = 1.5,
            iter = 10) %>%
        equalize_image(dimensions = 1:30, sleft = 2.5, sright = 2.5) %>%
        segment_image(dimensions = 1:30,
            col_resolution = 20) %>%
        isolate_territories()
}
cat("Query : DONE\n")
#-----------------------------------------------------------------------------#
# Mapping - no need to rebuild embeddings here since that is downstream
# analysis
# Note that the cell type labels will be used for composition and cell type
# but cell_labels and interactions will be added to the final object
#-----------------------------------------------------------------------------#
matched <- map_assays(seed_assay = seed,
    query_assay = query,
    neighborhood = "knn",
    use_norm = "log_norm",
    method = "pearson",
    k = 20,
    epochs = 25,
    batch_size = 5000,
    allow_duplicates = TRUE,
    threshold = 0,
    filter_cells = TRUE,
    jitter = 0,
    use_cost = use_cost,
    seed_meta_labels = c("cell_labels","interactions"),
    query_meta_labels = c("cell_labels","interactions"),
    digits = 5)
cat("Mapping: DONE\n")
#-----------------------------------------------------------------------------#
# Adding and exporting scores for combinations
#-----------------------------------------------------------------------------#

contributions <- c("CV","IQR","POC")
for (i in contributions) {
    matched <- get_cost_contribution(matched,method = i)
}
contrib_loc <- grep("contribution_score", names(matched@cost))
contribution_list <- vector("list", length(contrib_loc))
for (i in seq_along(contrib_loc)) {
    contrib_local <- matched@cost[[contrib_loc[i]]]
    contrib_local$method <- contributions[i]
    colnames(contrib_local) <- gsub(contributions[i],"score",colnames(contrib_local))
    contribution_list[[i]] <- contrib_local
}
    
contribution_list <- do.call("rbind", contribution_list)
query_tag <- paste0(query_tag,"_",compressed_tag)
file_name <- paste0(output,"Vesalius_aligned_",data_type,"_",seed_tag,"_",query_tag,"_contribution_score.csv")
write.csv(contribution_list, file = file_name)
    


#-----------------------------------------------------------------------------#
# Export Mapping 
#-----------------------------------------------------------------------------#
export_match <- matched@territories

file_name <- paste0(output,"Vesalius_aligned_",data_type,"_",seed_tag,"_",query_tag,".csv")

write.csv(export_match, file = file_name) 
gc()
q("no")