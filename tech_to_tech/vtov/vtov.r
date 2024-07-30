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

if (!dir.exists("/common/martinp4/stos/report/vtov/")) {
    dir.create("//common/martinp4/stos/report/vtov/")
}
output <- "/common/martinp4/stos/report/vtov/"

cat("Output setup: DONE \n")

args <- commandArgs(TRUE)
idx <- as.numeric(args[1])

use_cost <- c("niche", "territory")

#-----------------------------------------------------------------------------#
# Load HD data
#-----------------------------------------------------------------------------#
file_out <- paste0(output, "visiumHD_mouse_brain_query.rds")
if(file.exists(file_out)){
    query <- readRDS(file_out)
} else {
    hd_coord <- arrow::read_parquet("/common/wonklab/VisiumHD/Mouse_brain/square_008um/spatial/tissue_positions.parquet")
    hd_coord <- hd_coord[hd_coord$in_tissue != 0, c("barcode","pxl_col_in_fullres","pxl_row_in_fullres")]
    colnames(hd_coord) <- c("barcodes", "x","y")

    scale <- fromJSON(file = "/common/wonklab/VisiumHD/Mouse_brain/square_008um/spatial/scalefactors_json.json")

    image <- imager::load.image("/common/wonklab/VisiumHD/Mouse_brain/square_008um/spatial/tissue_hires_image.png")

    counts <- Seurat::Read10X_h5("/common/wonklab/VisiumHD/Mouse_brain/square_008um/filtered_feature_bc_matrix.h5")

    hd_coord <- hd_coord[sample(seq(1, nrow(hd_coord)), 100000),]
    counts <- counts[,colnames(counts) %in% hd_coord$barcodes]

    query <- build_vesalius_assay(hd_coord,
    counts,
    image = image,
    scale = "auto")
    cat("Query setup: DONE \n")

    query <- query %>%
    generate_embeddings(filter_threshold = 1,
        filter_grid = 1,
        tensor_resolution = 0.5,
        dim_reduction = "PCA",
        nfeatures = 2000) %>%
    equalize_image(dimensions = 1:30, sleft = 2.5, sright = 2.5) %>%
    smooth_image(dimensions = 1:30,
        method = c("iso","box"),
        box = 15,
        sigma = 3,
        iter = 20) %>%
    segment_image(dimensions = 1:30,
        col_resolution = 25) %>%
    isolate_territories()
    cat("Ouery Processing: DONE \n")
    file_out <- paste0(output, "visiumHD_mouse_brain_query.rds")
    saveRDS(query, file = file_out)
    cat("Ouery Saving: DONE \n")
}

#-----------------------------------------------------------------------------#
# Load Visium data
#-----------------------------------------------------------------------------#
file_out <- paste0(output, "visium_mouse_brain_seed.rds")
if (file.exists(file_out)) {
    seed <- readRDS(file_out)
} else {
    input <- "/common/wonklab/visium_brain/"
    coordinates <- paste0(input, "rep1/spatial/tissue_positions_list.csv")
    coord <- read.csv(coordinates, header = TRUE)
    coord <- coord[coord[, 2] == 1, ]
    coord <- coord[, c(1, 5, 6)]
    colnames(coord) <- c("barcodes", "x", "y")
    coord$x <- as.numeric(coord$x)
    coord$y <- as.numeric(coord$y)

    counts <- paste0(input, "rep1/CytAssist_FFPE_Mouse_Brain_Rep1_filtered_feature_bc_matrix.h5")
    counts <- Seurat::Read10X_h5(counts)

    img <- paste0(input,"rep1/spatial/tissue_hires_image.png")
    img <- imager::load.image(img)

    scale_2 <- paste0(input,"rep1/spatial/scalefactors_json.json")
    scale_2 <- fromJSON(file = scale_2)


    seed <- build_vesalius_assay(coord,
    counts,
    image = img,
    scale = "auto")

    cat("Seed setup: DONE \n")
    seed <- seed %>%
    generate_embeddings(filter_threshold = 1,
        filter_grid = 1,
        tensor_resolution = 0.5,
        dim_reduction = "PCA",
        nfeatures = 2000) %>%
    equalize_image(dimensions = 1:30, sleft = 5, sright = 5) %>%
    smooth_image(dimensions = 1:30,
        method = c("iso"),
        sigma = 1,
        iter = 5) %>%
    segment_image(dimensions = 1:30,
        col_resolution = 25) %>%
    isolate_territories()
    cat("Seed Processing: DONE \n")
    file_out <- paste0(output, "visium_mouse_brain_seed.rds")
    saveRDS(seed, file = file_out)
    cat("Seed Saving: DONE \n")
}


#-----------------------------------------------------------------------------#
# Process
#-----------------------------------------------------------------------------#

radius <- seed@meta$scale$scale / 2

#-----------------------------------------------------------------------------#
# Map 
#-----------------------------------------------------------------------------#
matched <- map_assays(seed,
    query,
    signal = "variable_features",
    use_cost = use_cost,
    neighborhood = "radius",
    radius = radius,
    threshold = 0,
    batch_size = sum(seed@tiles$origin ==1),
    epochs = 30,
    use_norm = "log_norm",
    jitter = radius)
cat("Mapping: DONE \n")


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
        col_resolution = 25) %>%
    isolate_territories()

file_out <- paste0(output, "visiumHD_to_visium_mouse_brain_matched.rds")
saveRDS(matched, file = file_out)
cat("Mapping Saved: DONE \n")