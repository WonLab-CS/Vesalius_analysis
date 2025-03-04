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
library(oneiric, lib.loc = "/common/martinp4/R")
library(RColorBrewer)
set.seed(1547)

plan(multicore, workers = 2)
max_size <- 10000 * 1024^2
options(future.globals.maxSize = max_size)
#-----------------------------------------------------------------------------#
# Get cli args
#-----------------------------------------------------------------------------#
args <- commandArgs(TRUE)
idx <- as.numeric(args[1])
input <- args[2]
output <- args[3]
file_tag <- args[4]
#-----------------------------------------------------------------------------#
# Create a common format to be used for all other tools
#-----------------------------------------------------------------------------#
coord <- list.files(path = input,
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
colnames(coord) <- gsub("embryo", "sample", colnames(coord))
colnames(coord) <- gsub("celltype_mapped_refined", "cell_labels", colnames(coord))
coord <- coord[coord$cell_labels != "Low quality", ]
spatial <- split(coord, coord$sample)
spatial <- lapply(spatial, function(spa){
    spa$x <- spa$x - min(spa$x) + 1
    spa$y <- spa$y - min(spa$y) + 1
    return(spa)
})
spatial <- oneiric::add_interactions(spatial, k = 10)


cells <- list.files(path = input,
    pattern = "counts.Rds", full.names = TRUE)
cells <- readRDS(cells)
cells <- lapply(spatial, function(spa, counts){
    return(counts[, colnames(counts) %in% spa$barcodes])
}, counts = cells)
names(cells) <- names(spatial)

#-----------------------------------------------------------------------------#
# export and save with common format
#-----------------------------------------------------------------------------#
for (i in seq_along(spatial)) {
    file_name <- paste0(output, file_tag, "_spatial_coordinates_sample_", names(spatial)[i],".csv")
    write.csv(spatial[[i]], file = file_name, quote = FALSE, row.names = FALSE)
    file_name <- paste0(output, file_tag, "_gene_counts_sample_", names(spatial)[i], ".csv")
    cell_tmp <- data.frame("genes" = rownames(cells[[i]]), cells[[i]])
    write.csv(cell_tmp, file = file_name, quote = FALSE, row.names = TRUE)
}
