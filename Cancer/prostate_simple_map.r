#-----------------------------------------------------------------------------#
# Loading Libraries and setting seed
#-----------------------------------------------------------------------------#
library(dplyr, lib.loc = "/common/martinp4/R")
library(igraph, lib.loc = "/common/martinp4/R")
library(future)
library(Matrix, lib.loc = "/common/martinp4/R")
library(ggplot2)
library(rjson)
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

library(ggplot2)
library(patchwork)
library(tidyr)
library(RColorBrewer)
library(ggnewscale, lib.loc = "/common/martinp4/R")
library(ggpubr)
library(ggtext,  lib.loc = "/common/martinp4/R")
set.seed(1547)


#-----------------------------------------------------------------------------#
# Get args
#-----------------------------------------------------------------------------#
args <- commandArgs(TRUE)
counts <- args[1]
coordinates <- args[2]
annotations <- args[3]
output <- args[4]
ref_tag <- args[5]
query_tag <- args[6]
cores <- as.numeric(args[7])


plan(multicore, workers = cores)
#-----------------------------------------------------------------------------#
# Utilities
#-----------------------------------------------------------------------------#

reformat <- function(counts, coordinates, annotations, data_set){
    counts <- readRDS(counts)
    coordinates <- readRDS(coordinates)
    annotations <- readRDS(annotations)
    counts <- counts[[data_set]]
    coordinates <- coordinates[[data_set]]
    coordinates$barcodes <- rownames(coordinates)
    coordinates <- coordinates[,c("barcodes","xcoord","ycoord")]
    colnames(coordinates) <- c("barcodes" , "x", "y")
    cells <- annotations$cell1[match(coordinates$barcodes,rownames(annotations))]
    names(cells) <- coordinates$barcodes
    ves <- build_vesalius_assay(coordinates, counts)
    ves <- add_cells(ves, cells = cells, add_name = "cell_labels")
    return(ves)
}




#-----------------------------------------------------------------------------#
# Load and process Refs
#-----------------------------------------------------------------------------#
ref <- reformat(counts, coordinates, annotations, ref_tag)

ref <- ref %>%
    generate_embeddings(tensor_resolution = 1, verbose = FALSE) %>%
    smooth_image(dimensions = seq(1, 30), sigma = 2, iter = 10,verbose = FALSE) %>%
    equalize_image(dimensions = seq(1, 30), sleft = 2.5, sright = 2.5, verbose = FALSE) %>%
    segment_image(dimensions = seq(1, 30), method = "kmeans", col_resolution  = 15, verbose = FALSE) %>%
    isolate_territories(verbose = FALSE)

cat("Ref Intra Procesing: DONE\n")


#-----------------------------------------------------------------------------#
# Load and process Query
#-----------------------------------------------------------------------------#
query <- reformat(counts, coordinates, annotations, query_tag)
query <- query %>%
    generate_embeddings(tensor_resolution = 1, verbose = FALSE) %>%
    smooth_image(dimensions = seq(1, 30), sigma = 2, iter = 10,verbose = FALSE) %>%
    equalize_image(dimensions = seq(1, 30), sleft = 2.5, sright = 2.5, verbose = FALSE) %>%
    segment_image(dimensions = seq(1, 30), method = "kmeans", col_resolution  = 15, verbose = FALSE) %>%
    isolate_territories(verbose = FALSE)
cat("Query Procesing: DONE\n")

#-----------------------------------------------------------------------------#
# map integrate and process
#-----------------------------------------------------------------------------#

matched <- map_assays(
    seed_assay = ref,
    query_assay = query,
    use_norm = "log_norm",
    neighborhood = "knn",
    k = 15,
    threshold = 0,
    batch_size = 5000,
    epochs = 20,
    jitter = 0,
    use_cost = c("feature","niche","territory","cell_type"),
    seed_meta_labels = c("cell_labels"),
    query_meta_labels = c("cell_labels"),
    digits = 5)




#-----------------------------------------------------------------------------#
# Prepare data for plotting
#-----------------------------------------------------------------------------#
ref_ter <- ref@territories[,c("barcodes","x","y","cell_labels")]
ref_ter$sample <- "Reference"
query_ter <- query@territories[,c("barcodes","x","y","cell_labels")]
query_ter$sample <- "Query"
matched_ter <- matched@territories[,c("barcodes","x","y","cell_labels")]
matched_ter$sample <- "Mapped"
merged <- rbind(ref_ter, query_ter,matched_ter)
#-----------------------------------------------------------------------------#
# plot and export
#-----------------------------------------------------------------------------#
base_colours <- c(
        "#E69F00",
        "#56B4E9",
        "#009E73",
        "#F0E442",
        "#0072B2",
        "#D55E00",
        "#CC79A7",
        "#999999")
    
ter_pal <- colorRampPalette(base_colours)

merged$cell_labels <- as.factor(merged$cell_labels)

cols  <- ter_pal(length(levels(merged$cell_labels)))

g1 <- ggplot(merged, aes(x,y,col = cell_labels)) +
    geom_point(size = 0.8, alpha = 1) +
    scale_color_manual(values = cols) +
    theme_bw() + 
    labs(title = "", color = "") +
    theme(strip.background =element_rect(fill="#082233ff"),
        strip.text = element_text(colour = 'white', size = 15),
        legend.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.position = "bottom") +
    guides(colour = guide_legend(
        override.aes = list(size = 5))) +
    facet_wrap(~sample)

file_name <- paste0(output,"htan_map_cells.pdf")
pdf(file_name, height = 8 , width = 16)
print(g1)
dev.off()
  
