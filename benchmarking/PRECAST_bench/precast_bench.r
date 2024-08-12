#-----------------------------------------------------------------------------#
# Loading Libraries and setting seed
#-----------------------------------------------------------------------------#
library(scater, lib.loc = "/common/martinp4/R")
library(Seurat, lib.loc = "/common/martinp4/R")
library(PRECAST, lib.loc = "/common/martinp4/R")
library(mclust,lib.loc = "/common/martinp4/R")

set.seed(1547)

#-----------------------------------------------------------------------------#
# Output set up 
#-----------------------------------------------------------------------------#
input <- "/common/wonklab/synthetic_spatial"

if(!dir.exists("/common/martinp4/benchmarking_out/precast/report/")){
    dir.create("/common/martinp4/benchmarking_out/precast/report/")
}
output_data <- "/common/martinp4/benchmarking_out/precast/report/"
cat("Output setup: DONE \n")
#-----------------------------------------------------------------------------#
# Utility functions
# scoring common latent space and clusters is a bit tricky since there 
# is no guarantee that cluster labels will recover cell type labels 
#-----------------------------------------------------------------------------#

assign_cluster <- function(reslist, query, seed) {
    query$cluster <- unlist(reslist[[1]]$cluster[2,])
    seed$cluster <- unlist(reslist[[1]]$cluster[1,])
    query <- split(query, query$cluster)
    seed <- split(seed, seed$cluster)
    clusters <- names(query)
    for (i in seq_along(clusters)){
        tmp <- seed[[clusters[i]]]
        query[[i]]$col <- sample(tmp$col, size = nrow(query[[i]]), replace = TRUE)
        query[[i]]$row <- sample(tmp$row, size = nrow(query[[i]]), replace = TRUE)
    }
    query <- do.call("rbind", query)
    query <- query[, c("barcodes", "col","row","cluster")]
    colnames(query) <- c("barcodes","x","y","cell_labels")
    return(query)
}


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


ref_coord <- read.csv(synth_coord[i])
colnames(ref_coord) <- gsub("x","col",colnames(ref_coord))
colnames(ref_coord) <- gsub("y","row",colnames(ref_coord))
rownames(ref_coord) <- ref_coord$barcodes
ref_counts <- read.csv(synth_counts[i], row.names = 1)
ref_counts$genes <- NULL

query_coord <- read.csv(synth_coord[j])
colnames(query_coord) <- gsub("x","col",colnames(query_coord))
colnames(query_coord) <- gsub("y","row",colnames(query_coord))
rownames(query_coord) <- query_coord$barcodes
query_counts <- read.csv(synth_counts[j], row.names = 1)
query_counts$genes <- NULL

#-----------------------------------------------------------------------------#
# Creating objects
#-----------------------------------------------------------------------------#
ref <- Seurat::CreateSeuratObject(counts = ref_counts,
    meta.data = ref_coord,
    project = "benchmarking")

query <- Seurat::CreateSeuratObject(counts = query_counts,
    meta.data = query_coord,
    project = "benchmarking")

object_list <- list(ref,query)

precast <- PRECAST::CreatePRECASTObject(seuList = object_list,
    selectGenesMethod = "HVGs",
    gene.number = 2000,
    rawData.preserve = TRUE)


#-----------------------------------------------------------------------------#
# Run PRECAST
#-----------------------------------------------------------------------------#
precast <-  PRECAST::AddAdjList(precast, platform = "Visium")
precast <- PRECAST::AddParSetting(precast,
    Sigma_equal = FALSE,
    maxIter = 20,
    verbose = TRUE,
    coreNum = 1)

precast <- PRECAST::PRECAST(precast, K = 12) 
reslist <- precast@resList
precast <- PRECAST::SelectModel(precast)

query_tmp <- assign_cluster(reslist,query_coord, ref_coord)



ref_tag <- gsub(".csv","", tag[i])
ref_tag <- gsub("spatial_territories_spatial_coordinates_","",ref_tag)
query_tag <- gsub(".csv","", tag[j])
query_tag <- gsub("spatial_territories_spatial_coordinates_","",query_tag)
write.csv(query_tmp,
    file = paste0(output_data,"PRECAST_aligned_",ref_tag,"_",query_tag,".csv")) 