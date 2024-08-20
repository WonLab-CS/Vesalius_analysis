#-----------------------------------------------------------------------------#
# Loading Libraries and setting seed
#-----------------------------------------------------------------------------#
library(ggplot2)
library(scater, lib.loc = "/common/martinp4/R")
library(Seurat, lib.loc = "/common/martinp4/R")
library(PRECAST, lib.loc = "/common/martinp4/R")
library(mclust,lib.loc = "/common/martinp4/R")

set.seed(1547)

#-----------------------------------------------------------------------------#
# Output set up 
#-----------------------------------------------------------------------------#
input <- "/common/wonklab/synthetic_spatial"

if(!dir.exists("/common/martinp4/benchmarking_out/precast/report/seqFISH/")){
    dir.create("/common/martinp4/benchmarking_out/precast/report/seqFISH/")
}
output <- "/common/martinp4/benchmarking_out/precast/report/seqFISH/"
cat("Output setup: DONE \n")
#-----------------------------------------------------------------------------#
# Utility functions
# scoring common latent space and clusters is a bit tricky since there 
# is no guarantee that cluster labels will recover cell type labels 
#-----------------------------------------------------------------------------#

assign_cluster <- function(reslist, query, seed) {
    query$cluster <- unlist(reslist$cluster[2,])
    seed$cluster <- unlist(reslist$cluster[1,])
    query <- split(query, query$cluster)
    seed <- split(seed, seed$cluster)
    clusters <- names(query)
    for (i in seq_along(clusters)){
        tmp <- seed[[clusters[i]]]
        locs <- sample(seq_len(nrow(tmp)),size = nrow(query[[i]]), replace = TRUE)
        query[[i]]$col <- tmp$col[locs]
        query[[i]]$row <- tmp$row[locs]
        query[[i]]$cell_labels <- tmp$cell_labels[locs]
    }
    query <- do.call("rbind", query)
    query <- query[, c("barcodes", "col","row","cell_labels","cluster")]
    colnames(query) <- c("barcodes","x","y","cell_labels","cluster")
    return(query)
}


#-----------------------------------------------------------------------------#
# Data set up 
#-----------------------------------------------------------------------------#
seed_tag <- "embryo3"
query_tag <- "embryo1"

if (file.exists(paste0(output,"seqFISH_",seed_tag,"_seed.rds")) &
    file.exists(paste0(output,"seqFISH_",query_tag,"_query.rds"))) {
    seed <- readRDS(paste0(output,"seqFISH_",seed_tag,"_seed.rds"))
    query <- readRDS(paste0(output,"seqFISH_",query_tag,"_query.rds"))
} else {
    # original data
    coord <- list.files(path = "/common/wonklab/seqFISH/",
        pattern = "metadata", full.names = TRUE)
    coord <- readRDS(coord)

    coord <- coord[, c("uniqueID",
        "x_global",
        "y_global",
        "embryo",
        "celltype_mapped_refined")]
    colnames(coord) <- gsub("x_global", "col", colnames(coord))
    colnames(coord) <- gsub("y_global", "row", colnames(coord))
    colnames(coord) <- gsub("uniqueID", "barcodes", colnames(coord))
    colnames(coord) <- gsub("celltype_mapped_refined", "trial", colnames(coord))
    coord <- coord[coord$trial != "Low quality", ]
    counts <- list.files(path = "/common/wonklab/seqFISH/",
        pattern = "counts.Rds", full.names = TRUE)
    counts <- readRDS(counts)
    counts <- counts[, colnames(counts) %in% coord$barcodes]
    nfeatures <- nrow(counts)
    # create seed

    seed_coord <- coord[grep(seed_tag, coord$barcodes), ]
    seed_counts <- counts[, grep(seed_tag, coord$barcodes)]
    seed <- Seurat::CreateSeuratObject(counts = seed_counts,
        meta.data = seed_coord,
        project = "benchmarking")

    saveRDS(seed, file  = paste0(output,"seqFISH_",seed_tag,"_seed.rds"))
    # create query
    
    query_coord <- coord[grep(query_tag, coord$barcodes), ]
    query_counts <- counts[, grep(query_tag, coord$barcodes)]
    query <- Seurat::CreateSeuratObject(counts = query_counts,
        meta.data = query_coord,
        project = "benchmarking")
    saveRDS(query, file  = paste0(output,"seqFISH_",query_tag,"_query.rds"))
}

#-----------------------------------------------------------------------------#
# Creating objects
#-----------------------------------------------------------------------------#


object_list <- list(seed,query)

precast <- PRECAST::CreatePRECASTObject(seuList = object_list,
    selectGenesMethod = "HVGs",
    gene.number = 400,
    rawData.preserve = TRUE)


#-----------------------------------------------------------------------------#
# Run PRECAST
#-----------------------------------------------------------------------------#
precast <-  PRECAST::AddAdjList(precast,
    platform = "Other_SRT",
    type = "fixed_number")
precast <- PRECAST::AddParSetting(precast,
    Sigma_equal = FALSE,
    maxIter = 20,
    verbose = TRUE,
    coreNum = 1)

buffer <- try(PRECAST::PRECAST(precast, K = 5:35), silent = TRUE)
if (is(buffer, "try-error")) {
    q("no")
}
precast <- PRECAST::SelectModel(buffer)
reslist <- precast@resList

query_coord <- query@meta.data
seed_coord <- seed@meta.data
save(query_coord, seed_coord, reslist,precast, file = paste0(output,"intermed_seqFISH.Rda"))
query_tmp <- assign_cluster(reslist,query_coord, seed_coord)

write.csv(query_tmp,
    file = paste0(output,"PRECAST_aligned_",seed_tag,"_",query_tag,".csv")) 

#-----------------------------------------------------------------------------#
# plotting
#-----------------------------------------------------------------------------#
n_colors <- length(unique(coord$trial))
base_colours <- c(
      "#E69F00",
      "#56B4E9",
      "#009E73",
      "#F0E442",
      "#0072B2",
      "#D55E00",
      "#CC79A7",
      "#999999")
pal <- colorRampPalette(base_colours)
colors <- pal(n_colors)

p <- ggplot(query_tmp, aes(x,y, col = as.factor(cell_lables))) +
    geom_point(size = 0.7) +
    scale_color_manual(values = colors) +
    theme_void() +
    theme(legend.title = element_text(size = 20),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.position = "left") +
    labs(color = "", title = "PRECAST - Query Cells")

file_name <- paste0(output,"PRECAST_aligned_",seed_tag,"_",query_tag,".pdf")
pdf(file_name, width = 12, height = 8)