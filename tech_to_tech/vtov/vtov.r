#-----------------------------------------------------------------------------#
# Loading Libraries and setting seed
#-----------------------------------------------------------------------------#
library(Seurat)
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
library(arrow, lib.loc = "/common/martinp4/R")
library(quadprog,lib.loc = "/common/martinp4/R")
library(spacexr,lib.loc = "/common/martinp4/R")
library(rjson)
library(lpSolve, lib.loc = "/common/martinp4/R")
library(TreeDist, lib.loc = "/common/martinp4/R")
library(mclust, lib.loc = "/common/martinp4/R")
library(anndata, lib.loc = "/common/martinp4/R")
library(pwr, lib.loc = "/common/martinp4/R")
library(gsignal, lib.loc = "/common/martinp4/R")
library(kohonen, lib.loc = "/common/martinp4/R")
library(registry, lib.loc = "/common/martinp4/R")
library(rngtools, lib.loc = "/common/martinp4/R")
library(NMF, lib.loc = "/common/martinp4/R")
library(vesalius, lib.loc = "/common/martinp4/R")
library(RColorBrewer)
set.seed(1547)

plan(multicore, workers = 1)
max_size <- 10000 * 1024^2
options(future.globals.maxSize = max_size)
#-----------------------------------------------------------------------------#
# Set future global for multicore processing
#-----------------------------------------------------------------------------#

cat("Output setup: DONE \n")

args <- commandArgs(TRUE)
output <- as.numeric(args[1])

use_cost <- c("niche", "territory")

#-----------------------------------------------------------------------------#
# Load HD data
#-----------------------------------------------------------------------------#
file_out <- paste0(output, "visiumHD_mouse_brain_query.rds")
query <- readRDS(file_out)
file_out <- paste0(output, "visiumHD_mouse_brain_cells.rds")
cells <- readRDS(file_out)


#-----------------------------------------------------------------------------#
# Load Visium data
#-----------------------------------------------------------------------------#
file_out <- paste0(output, "visium_mouse_brain_seed.rds")
seed <- readRDS(file_out)
file_out <- paste0(output, "visium_mouse_brain_prop.rds")
prop <- readRDS(file_out)
#-----------------------------------------------------------------------------#
# Process
#-----------------------------------------------------------------------------#
radius <- seed@meta$scale$scale / 2

convert_to_cont <- function(matched, prop, cells) {
    map <- matched@map
    prefix <- sapply(map$to, function(dat){
        prefix <- paste0(strsplit(x = dat, split = "-")[[1]][1:2], collapse = "-")
        return(prefix)
    })
    map$to <- prefix
    map <- split(map, map$to) 
    map <- map[prefix]
    map <- lapply(map, function(map, prop, cells) {
        local_prop <- prop[[map$to[1L]]]
        local_cells <- unlist(cells[names(cells) %in% map$from])
        if (length(local_prop) == 0 || length(local_cells) == 0) {
            return(0)
        }
        local_cells <- (table(local_cells) / length(local_cells)) * 100
        types <- intersect(names(local_prop), names(local_cells))
        if(length(types) == 0){
            freq <- 0
        } else{
            local_prop <- local_prop[types]
            local_cells <- local_cells[types]
            freq <- cbind(local_prop,local_cells)
            freq <- suppressWarnings(chisq.test(freq)$`p.value`)
        }
        
        return(freq)
    }, prop = prop, cells = cells)
    matched@map$prop <- unlist(map)
    return(matched)
}

#-----------------------------------------------------------------------------#
# Map - with jitter
#-----------------------------------------------------------------------------#
matched <- map_assays(seed,
    query,
    signal = "variable_features",
    use_cost = use_cost,
    neighborhood = "radius",
    radius = radius,
    threshold = 0,
    batch_size = sum(seed@tiles$origin ==1),
    epochs = 20,
    use_norm = "log_norm",
    jitter = radius)
cat("Mapping: DONE \n")
file_out <- paste0(output, "visiumHD_to_visium_mouse_brain_matched.rds")
saveRDS(matched, file = file_out)
cat("Mapping Saved: DONE \n")
matched <- convert_to_cont(matched, prop, cells)

matched <- matched %>%
    generate_embeddings(filter_threshold = 1,
        filter_grid = 1,
        tensor_resolution = 1,
        dim_reduction = "PCA",
        nfeatures = 2000) %>%
    equalize_image(dimensions = 1:30, sleft = 2.5, sright = 2.5) %>%
    smooth_image(dimensions = 1:30,
        method = c("iso","box"),
        box = 15,
        sigma = 3,
        iter = 20) %>%
    segment_image(dimensions = 1:30,
        col_resolution = 18) %>%
    isolate_territories(capture_radius = 0.01)
cells <- unlist(cells)
cells <- cells[match(matched@territories$barcodes, names(cells))]
names(cells)<-make.unique(names(cells), sep ="-")
matched <- add_cells(matched, cells)
file_out <- paste0(output, "visiumHD_to_visium_mouse_brain_matched.rds")
saveRDS(matched, file = file_out)
cat("Mapping Saved: DONE \n")

#-----------------------------------------------------------------------------#
# Map without jitter - mergung to look like original visium
#-----------------------------------------------------------------------------#

matched <- map_assays(seed,
    query,
    signal = "variable_features",
    use_cost = use_cost,
    neighborhood = "radius",
    radius = radius,
    threshold = 0,
    batch_size = sum(seed@tiles$origin ==1),
    epochs = 20,
    use_norm = "log_norm",
    jitter = 0)
cat("Mapping: DONE \n")
file_out <- paste0(output, "visiumHD_to_visium_mouse_brain_no_jitter.rds")
saveRDS(matched, file = file_out)
cat("Mapping Saved: DONE \n")
matched <- convert_to_cont(matched, prop, cells)

matched <- matched %>%
    generate_embeddings(filter_threshold = 1,
        filter_grid = 1,
        tensor_resolution = 0.95,
        dim_reduction = "PCA",
        nfeatures = 2000) %>%
    equalize_image(dimensions = 1:30, sleft = 5, sright = 5) %>%
    smooth_image(dimensions = 1:30,
        method = c("iso"),
        sigma = 1,
        iter = 5) %>%
    segment_image(dimensions = 1:30,
        col_resolution = 20) %>%
    isolate_territories(capture_radius = 0.01)

file_out <- paste0(output, "visiumHD_to_visium_mouse_brain_no_jitter.rds")
saveRDS(matched, file = file_out)
cat("Mapping Saved: DONE \n")