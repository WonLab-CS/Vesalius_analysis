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
    }
    query <- do.call("rbind", query)
    query <- query[, c("barcodes", "col","row","cluster")]
    colnames(query) <- c("barcodes","x","y","cluster")
    return(query)
}

#-----------------------------------------------------------------------------#
# Getting path to files
#-----------------------------------------------------------------------------#
input <- "/common/wonklab/SSv2"
coordinates <- list.files(path = input,
    pattern = "location|coord|Locations", full.names = TRUE)
counts <- list.files(path = input,
    pattern = "expression", full.names = TRUE)
tag <- list.files(path = input,
    pattern = "expression", full.names = FALSE)
tag <- gsub(".digital_expression.txt.gz|_expression_matrix.mtx.gz|.sparse_expression.txt",
    "", tag)
#-----------------------------------------------------------------------------#
# creating seed dataset
#-----------------------------------------------------------------------------#
f <- grep(pattern = "Puck_200115_08", tag)
seed_coord <- read.csv(coordinates[f], header = FALSE, skip = 1)
colnames(seed_coord) <- c("barcodes", "col", "row")
seed_counts <- read.table(counts[f], header = TRUE, row.names = 1)

seed_counts <- seed_counts[, apply(seed_counts, 2, sum) > 200]
seed_counts <- seed_counts[apply(seed_counts,1, sum) > 100, ]
seed_coord <- seed_coord[seed_coord$barcodes %in% colnames(seed_counts),]
seed_tag <-"Puck_200115_08"
seed <- Seurat::CreateSeuratObject(counts = seed_counts,
        meta.data = seed_coord,
        project = "benchmarking")
#-----------------------------------------------------------------------------#
# creating query 
#-----------------------------------------------------------------------------#
f <- grep(pattern = "Puck_190921_21", tag)
query_coord <- read.csv(coordinates[f], header = FALSE, skip = 1)
colnames(query_coord) <- c("barcodes", "col", "row")
query_counts <- read.table(counts[f], header = TRUE, row.names = 1)
query_counts <- query_counts[, apply(query_counts, 2, sum) > 200]
query_counts <- query_counts[apply(query_counts,1, sum) > 100, ]
query_coord <- query_coord[query_coord$barcodes %in% colnames(query_counts),]
query_tag <- "Puck_190921_21"
query <- Seurat::CreateSeuratObject(counts = query_counts,
        meta.data = query_coord,
        project = "benchmarking")
#-----------------------------------------------------------------------------#
# Creating objects
#-----------------------------------------------------------------------------#


object_list <- list(seed,query)

precast <- PRECAST::CreatePRECASTObject(seuList = object_list,
    selectGenesMethod = "HVGs",
    gene.number = 2000,
    rawData.preserve = TRUE)


#-----------------------------------------------------------------------------#
# Run PRECAST
#-----------------------------------------------------------------------------#
precast <-  PRECAST::AddAdjList(precast,
    platform = "Other_SRT",
    type = "fixed_number",
    number = 30)
precast <- PRECAST::AddParSetting(precast,
    Sigma_equal = FALSE,
    maxIter = 20,
    verbose = TRUE,
    coreNum = 1)

buffer <- try(PRECAST::PRECAST(precast, K = 5:35), silent = TRUE)
if (is(buffer, "try-error")){
    q("no")
}
precast <- PRECAST::SelectModel(buffer)
reslist <- precast@resList
save(query_coord, seed_coord, reslist,precast, file = paste0(output,"intermed_slide.Rda"))
query_coord$cluster <-  unlist(reslist$cluster[2,])
seed_coord$cluster <-  unlist(reslist$cluster[1,])



#-----------------------------------------------------------------------------#
# plotting
#-----------------------------------------------------------------------------#
n_colors <- length(unique(query_coord$cluster))
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


p <- ggplot(seed_coord, aes(col,row, col = as.factor(cluster))) +
    geom_point(size = 0.7) +
    scale_color_manual(values = colors) +
    theme_void() +
    theme(legend.title = element_text(size = 20),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.position = "left") +
    labs(color = "", title = "PRECAST - Reference")

p1 <- ggplot(query_coord, aes(col,row, col = as.factor(cluster))) +
    geom_point(size = 0.7) +
    scale_color_manual(values = colors) +
    theme_void() +
    theme(legend.title = element_text(size = 20),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.position = "left") +
    labs(color = "", title = "PRECAST - Query")



file_name <- paste0(output,"PRECAST_aligned_",seed_tag,"_",query_tag,".pdf")
pdf(file_name, width = 22, height = 8)