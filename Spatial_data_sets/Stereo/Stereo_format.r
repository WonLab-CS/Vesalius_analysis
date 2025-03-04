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
library(oneiric, lib.loc = "/common/martinp4/R")
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
seed_input <- args[2]
query_input <- args[3]
output <- args[4]
file_tag <- args[5]
seed_tag <- args[6]
query_tag <- args[7]



#-----------------------------------------------------------------------------#
# load and prepare data
# Down sample for speed
#-----------------------------------------------------------------------------#
seed <- read_h5ad(seed_input)
seed_counts <- t(seed$X)
colnames(seed_counts) <- paste0("cell_", colnames(seed_counts))
seed_annot <- as.character(seed$obs$annotation)

seed_coord <- cbind(rownames(seed$obs), as.data.frame(seed$obsm$spatial), seed_annot)
colnames(seed_coord) <- c("barcodes", "x", "y", "cell_labels")
seed_coord$barcodes <- paste0("cell_", seed_coord$barcodes)
seed_coord$sample <- seed_tag
locs_seed <- sample(seed_coord$barcodes, size = 50000, replace = FALSE)
seed_coord <- seed_coord[seed_coord$barcodes %in% locs_seed, ]
seed_coord$x <- seed_coord$x - min(seed_coord$x) + 1
seed_coord$y <- seed_coord$y - min(seed_coord$y) + 1
seed_counts <- as(seed_counts[, colnames(seed_counts) %in% locs_seed], "CsparseMatrix")


query <- read_h5ad(query_input)
query_counts <- t(query$X)
colnames(query_counts) <- paste0("qcell_", colnames(query_counts))
query_annot <- as.character(query$obs$annotation)

query_coord <- cbind(rownames(query$obs), as.data.frame(query$obsm$spatial), query_annot)
colnames(query_coord) <- c("barcodes", "x", "y", "cell_labels")
query_coord$barcodes <- paste0("qcell_", query_coord$barcodes)
query_coord$sample <- query_tag
locs_query <- sample(query_coord$barcodes, size = 50000, replace = FALSE)
query_coord <- query_coord[query_coord$barcodes %in% locs_query, ]
query_coord$x <- query_coord$x - min(query_coord$x) + 1
query_coord$y <- query_coord$y - min(query_coord$y) + 1
query_counts <- as(query_counts[, colnames(query_counts) %in% locs_query], "CsparseMatrix")

#-----------------------------------------------------------------------------#
# add innteractions and export
#-----------------------------------------------------------------------------#
spatial <- list(seed_coord, query_coord)
names(spatial) <- c(seed_tag,query_tag)
spatial <- oneiric::add_interactions(spatial, k = 10)
spatial <- spatial[c(seed_tag,query_tag)]

counts <- list(seed_counts,query_counts)
names(counts) <- c(seed_tag,query_tag)

for (i in seq_along(spatial)) {
    file_name <- paste0(output, file_tag, "_spatial_coordinates_sample_", names(spatial)[i],".csv")
    write.csv(spatial[[i]], file = file_name, quote = FALSE, row.names = FALSE)
    file_name <- paste0(output, file_tag, "_gene_counts_sample_", names(spatial)[i], ".csv")
    cell_tmp <- data.frame("genes" = rownames(counts[[i]]), counts[[i]])
    write.csv(cell_tmp, file = file_name, quote = FALSE, row.names = TRUE)
}

