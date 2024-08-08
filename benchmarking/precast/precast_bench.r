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

if(!dir.exists("/common/martinp4/benchmarking_out/Vesalius/synthetic/")){
    dir.create("/common/martinp4/benchmarking_out/Vesalius/synthetic/")
}
output_data <- "/common/martinp4/benchmarking_out/precast/report/"
cat("Output setup: DONE \n")
#-----------------------------------------------------------------------------#
# Utility functions
# scoring common latent space and clusters is a bit tricky since there 
# is no guarantee that cluster labels will recover cell type labels 
#-----------------------------------------------------------------------------#

assign_cluster <- function(clusters, query, seed, counts = NULL, method = "sort") {
    query <- switch(EXPR = method,
        "sort" = sort_cluster(clusters, query),
        "dist" = dist_cluster(clusters, query, counts))
    query <- assign_coord(query, seed)
    return(query)
}

sort_cluster <- function(clusters, query) {
    query_cluster <- sort(table(clusters))
    query_cell_labels <- sort(table(query$cell_labels))
    for (i in seq_along(query_cluster)) {
        locs_query <- query_cluster == names(query_cluster)[i]
        query$cell_labels[locs_query] <- names(query_cell_labels)[i]
    }
    return(query)

}

# Assigning random coordinates from the matched cell label
assign_coord <- function(query, seed) {
    locs <- rep(0, nrow(query))
    cell_labels <- unique(query$cell_labels)
    for (i in seq_along(cell_labels)) {
        q_locs <- which(query$cell_labels == cell_labels[i])
        s_locs <- which(seed$cell_labels == cell_labels[i])
        locs[q_locs] <- sample(s_locs, size = length(q_locs), replace = TRUE)
    }
    query$col <- seed$col[locs]
    query$row <- seed$row[locs]
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
    maxIter = 30,
    verbose = TRUE,
    coreNum = 1)

precast <- PRECAST::PRECAST(precast, K = length(unique(ref_coord$cell_labels)))
reslist <- precast@resList
precast <- PRECAST::SelectModel(precast)

query_tmp <- assign_cluster(unlist(reslist[[1]]$cluster[2,]),query_coord, ref_coord)
query_tmp <- query_tmp[,c("barcodes", "col","row","cell_labels")]
colnames(query_tmp) <- c("barcodes","x","y","cell_labels")


ref_tag <- gsub(".csv","", tag[i])
ref_tag <- gsub("spatial_territories_spatial_coordinates_","",ref_tag)
query_tag <- gsub(".csv","", tag[j])
query_tag <- gsub("spatial_territories_spatial_coordinates_","",query_tag)
write.csv(query_tmp,
    file = paste0(output_data,"PRECAST_aligned_",ref_tag,"_",query_tag,".csv")) 