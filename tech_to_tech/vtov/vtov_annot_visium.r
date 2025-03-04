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
set.seed(1547)

plan(multicore, workers = 1)
max_size <- 10000 * 1024^2
options(future.globals.maxSize = max_size)
#-----------------------------------------------------------------------------#
# directory prep
#-----------------------------------------------------------------------------#


cat("Output setup: DONE \n")

args <- commandArgs(TRUE)

input <- args[1]
rds <- args[2]
output <- args[3]
#-----------------------------------------------------------------------------#
# build from annot
#-----------------------------------------------------------------------------#
rds <- paste0(rds,"SCRef_hippocampus.RDS")
rds <- readRDS(rds)
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


coordinates <- paste0(input, "rep1/spatial/tissue_positions_list.csv")
coord <- read.csv(coordinates, header = TRUE)
coord <- coord[coord[, 2] == 1, ]
coord <- coord_ves <-coord[, c(1, 5, 6)]
coord_ves <- coord
colnames(coord_ves) <- c("barcodes", "x","y")
rownames(coord) <- coord[,1]
coord <- coord[,2:3]
coord[,1] <- as.numeric(coord[,1])
coord[,2] <- as.numeric(coord[,2])

counts <- paste0(input, "rep1/CytAssist_FFPE_Mouse_Brain_Rep1_filtered_feature_bc_matrix.h5")
counts <- as.matrix(Seurat::Read10X_h5(counts))
cat("Visium Build: DONE \n")

nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI
puck <- SpatialRNA(coord, counts, nUMI)
barcodes <- colnames(puck@counts)

RCTD <- create.RCTD(puck, reference, max_cores = 1)
RCTD <- run.RCTD(RCTD, doublet_mode = 'multi')
file_out <- paste0(output, "visium_mouse_brain_RCTD.rds")
saveRDS(RCTD, file = file_out)

all_weights <- lapply(RCTD@results, "[[", "all_weights")
conf <- lapply(RCTD@results, "[[", "conf_list")
all_weights <- mapply(function(w,co){
        w <- w[names(w) %in% names(co)[as.numeric(co) ==1]]
        w <- w / sum(w)
        return(w * 100)
    }, all_weights, conf)

names(all_weights) <- barcodes

saveRDS(all_weights, file = paste0(output, "visium_mouse_brain_prop.rds"))
cat("RCTD: DONE \n")
#-----------------------------------------------------------------------------#
# Build ves
#-----------------------------------------------------------------------------#
seed <- build_vesalius_assay(coord_ves,
    counts,
    scale = "auto")

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
        col_resolution = 20) %>%
    isolate_territories(capture_radius = 0.2)
cat("Vesalius: DONE \n")  
file_out <- paste0(output, "visium_mouse_brain_seed.rds")
saveRDS(seed, file = file_out)
    