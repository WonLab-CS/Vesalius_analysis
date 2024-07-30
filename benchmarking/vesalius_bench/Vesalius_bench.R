#-----------------------------------------------------------------------------#
# Loading Libraries and setting seed
#-----------------------------------------------------------------------------#
library(dplyr, lib.loc = "/common/martinp4/R")
library(igraph, lib.loc = "/common/martinp4/R")
library(future)
library(Matrix, lib.loc = "/common/martinp4/R")
library(ggplot2)
library(patchwork)
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
library(spatstat.utils, lib.loc = "/common/martinp4/R")
library(vesalius, lib.loc = "/common/martinp4/R")
set.seed(1547)
#-----------------------------------------------------------------------------#
# Output set up 
#-----------------------------------------------------------------------------#
input <- "/common/wonklab/synthetic_spatial"

if(!dir.exists("/common/martinp4/benchmarking_out/Vesalius/synthetic/")){
    dir.create("/common/martinp4/benchmarking_out/Vesalius/synthetic/")
}
output_data <- "/common/martinp4/benchmarking_out/Vesalius/synthetic/"
cat("Output setup: DONE \n")
#-----------------------------------------------------------------------------#
# Data set up 
#-----------------------------------------------------------------------------#
args <- commandArgs(TRUE)
idx <- as.numeric(args[1])
data_type <- args[2]

synth_coord <- paste0(data_type, "_spatial_territories_spatial_coordinates")
synth_coord <- list.files(input,
    pattern = synth_coord,
    full.names = TRUE)

synth_counts <- paste0(data_type, "_spatial_territories_gene_counts")
synth_counts <- list.files(input,
    pattern = synth_counts,
    full.names = TRUE)

tag <- paste0(data_type, "_spatial_territories_spatial_coordinates")
tag <- list.files(input,
    pattern = tag,
    full.names = FALSE)


i <- rep(seq_along(synth_coord), times = length(synth_coord))[idx]
j <- rep(seq_along(synth_coord), each = length(synth_coord))[idx]


use_cost <- c("feature","niche","territory")


cat("Data setup: DONE\n")

#-----------------------------------------------------------------------------#
# Start mapping
#-----------------------------------------------------------------------------#
ref_coord <- read.csv(synth_coord[i])
ref_counts <- read.csv(synth_counts[i], row.names = 1)
ref_counts$genes <- NULL
ref <- build_vesalius_assay(coordinates = ref_coord,
        counts = ref_counts,
        assay = paste0(data_type,"_",unique(ref_coord$sample)))
ref_cells <- ref_coord$cell_labels
names(ref_cells) <- ref_coord$barcodes
ref <- add_cells(ref, ref_cells)
ref <- ref %>%
    generate_embeddings() %>%
    smooth_image(dimensions = seq(1, 30), method =c("iso","box"),sigma=1, box = 10, iter = 20) %>%
    segment_image(dimensions = seq(1, 30), method = "kmeans", col_resolution  = 12)%>% 
    isolate_territories()
cat("Ref Procesing: DONE\n")

query_coord <- read.csv(synth_coord[j])
query_counts <- read.csv(synth_counts[j], row.names = 1)
query_counts$genes <- NULL

# to overcome issues with same names when comparing same seed and query
# the merging of the barcodes does not seem to be working very well
if (i == j) {
    query_coord$barcodes <- paste0("q_",query_coord$barcodes)
    colnames(query_counts) <- paste0("q_",colnames(query_counts))
}

query <- build_vesalius_assay(coordinates = query_coord,
        counts = query_counts,
        assay = paste0(data_type,"_",unique(query_coord$sample)))
query_cells <- query_coord$cell_labels
names(query_cells) <- query_coord$barcodes
query <- add_cells(query, query_cells)
query <- query %>%
    generate_embeddings() %>%
    smooth_image(dimensions = seq(1, 30), method =c("iso","box"),sigma=1, box = 10, iter = 20) %>%
    segment_image(dimensions = seq(1, 30), method = "kmeans", col_resolution  = 12) %>%
    isolate_territories()
cat("Query Procesing: DONE\n")




matched <- map_assays(seed_assay = ref,
    query_assay = query,
    neighborhood = "graph",
    use_norm = "raw",
    depth = 2,
    threshold = -1,
    epochs = 20,
    batch_size = 1000,
    jitter = FALSE,
    use_cost = use_cost)


cat("Mapping: DONE\n")

#-----------------------------------------------------------------------------#
# Export Mapping 
#-----------------------------------------------------------------------------#
export_match <- matched@territories
#export_match$sample <- query_coord$sample
colnames(export_match) <- c("barcodes","x","y","cell_labels")
export_match$barcodes <- gsub("q_","", export_match$barcodes)
ref_tag <- gsub(".csv","", tag[i])
ref_tag <- gsub("spatial_territories_spatial_coordinates_","",ref_tag)
query_tag <- gsub(".csv","", tag[j])
query_tag <- gsub("spatial_territories_spatial_coordinates_","",query_tag)
write.csv(export_match,
    file = paste0(output_data,"Vesalius_aligned_",ref_tag,"_",query_tag,".csv")) 

