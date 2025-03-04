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
set.seed(1547)

plan(multicore, workers = 1)
max_size <- 10000 * 1024^2
options(future.globals.maxSize = max_size)

#-----------------------------------------------------------------------------#
# out dir
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

hd_coord <- arrow::read_parquet(paste0(input,"spatial/tissue_positions.parquet"))
hd_coord <- hd_coord[hd_coord$in_tissue != 0, c("barcode","pxl_col_in_fullres","pxl_row_in_fullres")]
barcodes <- as.character(hd_coord$barcode)
hd_coord <- as.data.frame(hd_coord[,2:3])
rownames(hd_coord) <- barcodes


counts <- as.matrix(Seurat::Read10X_h5(paste0(input,"filtered_feature_bc_matrix.h5")))
hd_coord <- hd_coord[sample(seq(1, nrow(hd_coord)), 200000),]
counts <- counts[,colnames(counts) %in% rownames(hd_coord)]

cat("Visium HD build: DONE \n")


nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI
puck <- SpatialRNA(hd_coord, counts, nUMI)
barcodes <- colnames(puck@counts)

RCTD <- create.RCTD(puck, reference, max_cores = 1)
RCTD <- run.RCTD(RCTD, doublet_mode = 'doublet')
file_out <- paste0(output, "visiumHD_mouse_brain_RCTD.rds")
saveRDS(RCTD, file = file_out)

cells <- RCTD@results$results_df
cells <- cells[cells$spot_class != "reject", ]
barcodes <- rownames(cells)
cells <- lapply(seq_len(nrow(cells)), function(idx, cells){
    if (cells$spot_class[idx] == "singlet" | cells$spot_class[idx] == "doublet_uncertain") {
        cell_type <- cells$first_type[idx]
    } else if (cells$spot_class[idx] == "doublet_certain") {
        cell_type <- c(cells$first_type[idx],cells$second_type[idx])
    }
    return(as.character(cell_type))
},cells = cells)
file_out <- paste0(output, "visiumHD_mouse_brain_cells.rds")
names(cells)<- barcodes
saveRDS(cells, file = file_out)

cat("RCTD: DONE \n")
#-----------------------------------------------------------------------------#
# Visium HD vesalius
#-----------------------------------------------------------------------------#
file_out <- paste0(output, "visiumHD_mouse_brain_cells.rds")
cells <- readRDS(file_out)
hd_coord <- hd_coord[rownames(hd_coord) %in% names(cells),]
hd_coord$barcodes <- rownames(hd_coord)
colnames(hd_coord) <- c("x", "y", "barcodes")
hd_coord <- hd_coord[, c("barcodes","x", "y" )]
counts <- counts[,colnames(counts) %in% hd_coord$barcodes]

query <- build_vesalius_assay(hd_coord,
    counts,
    scale = "auto")


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
        col_resolution = 20) %>%
    isolate_territories()

file_out <- paste0(output, "visiumHD_mouse_brain_query.rds")
saveRDS(query, file = file_out)
