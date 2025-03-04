#-----------------------------------------------------------------------------#
# Loading Libraries and setting seed
#-----------------------------------------------------------------------------#
library(spatstat.utils, lib.loc = "/common/martinp4/R")
library(Seurat)
library(dplyr)
library(future)
library(future.apply)
library(ggplot2)
library(dplyr)
library(patchwork)
library(Matrix)
library(arrow, lib.loc = "/common/martinp4/R")
library(quadprog,lib.loc = "/common/martinp4/R")
library(spacexr,lib.loc = "/common/martinp4/R")
library(rjson)
library(RColorBrewer)
library(lpSolve, lib = "/common/martinp4/R")
library(TreeDist, lib.loc = "/common/martinp4/R")
library(mclust, lib.loc = "/common/martinp4/R")
library(pwr, lib.loc = "/common/martinp4/R")
library(kohonen, lib.loc = "/common/martinp4/R")
library(registry, lib.loc = "/common/martinp4/R")
library(rngtools, lib.loc = "/common/martinp4/R")
library(NMF, lib.loc = "/common/martinp4/R")
library(vesalius, lib.loc = "/common/martinp4/R")
library(oneiric, lib.loc = "/common/martinp4/R")
set.seed(1547)

plan(multicore, workers = 5)
max_size <- 10000 * 1024^2
options(future.globals.maxSize = max_size)
#-----------------------------------------------------------------------------#
# Get cli args
#-----------------------------------------------------------------------------#
args <- commandArgs(TRUE)
idx <- as.numeric(args[1])
input <- args[2]
input_ref <- args[3]
output <- args[4]
file_tag <- args[5]
seed <- args[6]
query <- args[7]

#-----------------------------------------------------------------------------#
# Getting path to files
#-----------------------------------------------------------------------------#
coordinates <- list.files(path = input,
    pattern = "location|coord|Locations", full.names = TRUE)
counts <- list.files(path = input,
    pattern = "expression", full.names = TRUE)
tag <- list.files(path = input,
    pattern = "expression", full.names = FALSE)
tag <- gsub(".digital_expression.txt.gz|_expression_matrix.mtx.gz|.sparse_expression.txt",
    "", tag)
#-----------------------------------------------------------------------------#
# Loading seed and query data sets 
#-----------------------------------------------------------------------------#
f <- grep(pattern = seed, tag)
seed_coord <- read.csv(coordinates[f], header = FALSE, skip = 1)
colnames(seed_coord) <- c("barcodes", "x", "y")
rownames(seed_coord) <- seed_coord$barcodes
seed_counts <- read.table(counts[f], header = TRUE, row.names = 1)
seed_counts <- seed_counts[, apply(seed_counts, 2, sum) > 200]
seed_counts <- seed_counts[apply(seed_counts,1, sum) > 100, ]
seed_coord <- seed_coord[seed_coord$barcodes %in% colnames(seed_counts),]

f <- grep(pattern = query, tag)
query_coord <- read.csv(coordinates[f], header = FALSE, skip = 1)
colnames(query_coord) <- c("barcodes", "x", "y")
rownames(query_coord) <- query_coord$barcodes
query_counts <- read.table(counts[f], header = TRUE, row.names = 1)
query_counts <- query_counts[, apply(query_counts, 2, sum) > 200]
query_counts <- query_counts[apply(query_counts,1, sum) > 100, ]
query_coord <- query_coord[query_coord$barcodes %in% colnames(query_counts),]

#-----------------------------------------------------------------------------#
# Loading ref scRNA
#-----------------------------------------------------------------------------#
rds <- readRDS(input_ref)
rds <- Seurat::UpdateSeuratObject(rds)

ref_counts <- as.matrix(Seurat::GetAssayData(rds, layer = "counts"))
ref_counts <- ref_counts[!duplicated(rownames(ref_counts)),]
cell_types <- rds@meta.data$liger_ident_coarse
names(cell_types) <- rownames(rds@meta.data)
cell_types <- as.factor(cell_types)
nUMI_ref <- rds@meta.data$nUMI
names(nUMI_ref) <- rownames(rds@meta.data)
reference <- Reference(ref_counts, cell_types, nUMI_ref)
cat("Ref build: DONE \n")
#-----------------------------------------------------------------------------#
# Annotating SEED
#-----------------------------------------------------------------------------#

nUMI <- colSums(seed_counts) 
puck <- SpatialRNA(seed_coord[,c("x","y")], seed_counts, nUMI)
barcodes <- colnames(puck@counts)

RCTD <- create.RCTD(puck, reference, max_cores = 1)
RCTD <- run.RCTD(RCTD, doublet_mode = 'doublet')


cells <- RCTD@results$results_df
print(head(cells))
cells <- cells[cells$spot_class != "reject" |
    cells$spot_class != "doublet_uncertain" |
    cells$spot_class != "doublet_certain", ]
seed_coord <- seed_coord[match(rownames(cells),seed_coord$barcodes),]
seed_coord$sample <- seed
seed_coord$cell_labels <- cells$first_type
seed_coord$x <- seed_coord$x - min(seed_coord$x) + 1
seed_coord$y <- seed_coord$y - min(seed_coord$y) + 1
seed_counts <- seed_counts[,colnames(seed_counts) %in% seed_coord$barcodes]
rm(RCTD,cells,nUMI,puck,barcodes)
#-----------------------------------------------------------------------------#
# Annotating QUERY
#-----------------------------------------------------------------------------#
nUMI <- colSums(query_counts)
puck <- SpatialRNA(query_coord[,c("x","y")], query_counts, nUMI)
barcodes <- colnames(puck@counts)

RCTD <- create.RCTD(puck, reference, max_cores = 1)
RCTD <- run.RCTD(RCTD, doublet_mode = 'doublet')

cells <- RCTD@results$results_df
cells <- cells[cells$spot_class != "reject" |
    cells$spot_class[idx] != "doublet_uncertain" |
    cells$spot_class[idx] != "doublet_certain", ]
query_coord <- query_coord[match(rownames(cells),query_coord$barcodes),]
query_coord$sample <- query
query_coord$cell_labels <- cells$first_type
query_coord$x <- query_coord$x - min(query_coord$x) + 1
query_coord$y <- query_coord$y - min(query_coord$y) + 1
query_counts <- query_counts[,colnames(query_counts) %in% query_coord$barcodes]
rm(RCTD,cells,nUMI,puck,barcodes)
#-----------------------------------------------------------------------------#
# add innteractions and export
#-----------------------------------------------------------------------------#
spatial <- list(seed_coord, query_coord)
names(spatial) <- c(seed,query)
spatial <- oneiric::add_interactions(spatial, k = 10)
spatial <- spatial[c(seed,query)]

counts <- list(seed_counts,query_counts)
names(counts) <- c(seed,query)

for (i in seq_along(spatial)) {
    file_name <- paste0(output, file_tag, "_spatial_coordinates_sample_", names(spatial)[i],".csv")
    write.csv(spatial[[i]], file = file_name, quote = FALSE, row.names = FALSE)
    file_name <- paste0(output, file_tag, "_gene_counts_sample_", names(spatial)[i], ".csv")
    cell_tmp <- data.frame("genes" = rownames(counts[[i]]), counts[[i]])
    write.csv(cell_tmp, file = file_name, quote = FALSE, row.names = TRUE)
}
