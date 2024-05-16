#-----------------------------------------------------------------------------#
# Loading Libraries and setting seed
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
library(dplyr,lib.loc = "/common/martinp4/R")
library(ggplot2)
library(patchwork)
library(ggpubr)
library(lpSolve, lib.loc = "/common/martinp4/R")
library(mclust, lib.loc = "/common/martinp4/R")
library(mcclust, lib.loc = "/common/martinp4/R")


set.seed(1547)
#-----------------------------------------------------------------------------#
# set directories
#-----------------------------------------------------------------------------#
if (!dir.exists("/common/martinp4/benchmarking_out/CytoSpace/")) {
    dir.create("/common/martinp4/benchmarking_out/CytoSpace/")
}
output_data <- "/common/martinp4/benchmarking_out/CytoSpace/"

ref <- list.files("/common/wonklab/synthetic_spatial",
    pattern = "spatial_coordinates",
    full.names = TRUE)

query <- list.dirs(output_data,
    full.names = TRUE,
    recursive = FALSE)


#-----------------------------------------------------------------------------#
# plots
#-----------------------------------------------------------------------------#
file_name <- paste0(output_data,"CytoSpace_output_plots_no_lab.pdf")
pdf(file_name, width = 8, height = 4)

for (i in seq_along(query)) {
    files_used <- paste0(query[i],"/files_used.txt")
    files <- readLines(files_used)
    ref_tag <- gsub("Reference: ","", files[1])
    ref_tag <- gsub(" ","",ref_tag)
    ref_tag <- gsub("gene_counts","spatial_coordinates", ref_tag)
    local_ref_data <- read.csv(paste0("/common/wonklab/synthetic_spatial/",ref_tag))
    ref_tag <- gsub("_spatial_territories_spatial_coordinates_"," ", ref_tag)
    ref_tag <- gsub("_"," ", ref_tag)
    ref_tag <- gsub(".csv", "", ref_tag)

    query_tag <- gsub("Query: ","", files[2])
    query_tag <- gsub(" ","",query_tag)
    query_tag <- gsub("_spatial_territories_gene_counts_"," ", query_tag)
    query_tag <- gsub("_"," ", query_tag)
    query_tag <- gsub(".csv", "", query_tag)

    local_query <- list.files(query[i], pattern = "assigned_locations.csv", full.names = TRUE)
    if (length(local_query) == 0) {
        next    
    }
    local_query <- read.csv(local_query)
    og_query <- gsub("Query: ","", files[2])
    og_query <- gsub(" ","",og_query)
    og_query <- gsub("gene_counts","spatial_coordinates", og_query)
    local_query_data <- read.csv(paste0("/common/wonklab/synthetic_spatial/",og_query))

    local_ref_data <- local_ref_data[match(local_query$SpotID,local_ref_data$barcodes),]
    local_query_data <- local_query_data[match(local_query$OriginalCID,local_query_data$barcodes),]
    local_query$cell_labels <- local_query_data$cell_labels

    cols <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))
   cols <- cols(max(c(length(unique(local_ref_data$cell_labels)),length(unique(local_query$cell_labels)))))
    g_ref <- ggplot(local_ref_data, aes(x = x, y = y, col = as.factor(cell_labels))) +
    geom_point(size = 0.5) +
    theme_bw() +
    labs(title = paste("Reference", ref_tag),color = "Cell Labels") +
    scale_color_manual(values = cols) +
    guides(colour = guide_legend(
        override.aes = list(size =  5)))


    g_query <- ggplot(local_query, aes(x = col, y = row, col = as.factor(cell_labels))) +
    geom_point(size = 0.5) +
    theme_bw() +
    labs(title = paste("Query", query_tag), color = "Cell Labels")+
    scale_color_manual(values = cols) +
    guides(colour = guide_legend(
        override.aes = list(size =  5)))

    g_bind <- g_ref + g_query
    print(g_bind)


}
dev.off()



#-----------------------------------------------------------------------------#
# Accuracy table
#-----------------------------------------------------------------------------#

file_name <- paste0(output_data,"CytoSpace_mapping_accuracy.csv")
map_score <- data.frame("from" = character(),
    "to" = character(),
    "ari" = numeric(),
    "vi" = numeric(),
    "accuracy" = numeric())
for (i in seq_along(query)) {
    files_used <- paste0(query[i],"/files_used.txt")
    files <- readLines(files_used)
    ref_tag <- gsub("Reference: ","", files[1])
    ref_tag <- gsub(" ","",ref_tag)
    ref_tag <- gsub("gene_counts","spatial_coordinates", ref_tag)
    local_ref_data <- read.csv(paste0("/common/wonklab/synthetic_spatial/",ref_tag))
    ref_tag <- gsub("_spatial_territories_spatial_coordinates_"," ", ref_tag)
    ref_tag <- gsub("_", " ", ref_tag)
    ref_tag <- gsub(".csv", "", ref_tag)

    query_tag <- gsub("Query: ","", files[2])
    query_tag <- gsub(" ","",query_tag)
    query_tag <- gsub("_spatial_territories_gene_counts_"," ", query_tag)
    query_tag <- gsub("_"," ", query_tag)
    query_tag <- gsub(".csv", "", query_tag)

    local_query <- list.files(query[i], pattern = "assigned_locations.csv", full.names = TRUE)
    if (length(local_query) == 0) {
        next    
    }
    local_query <- read.csv(local_query)

    local_ref_data <- local_ref_data[match(local_query$SpotID,local_ref_data$barcodes),]
    og_query <- gsub("Query: ","", files[2])
    og_query <- gsub(" ","",og_query)
    og_query <- gsub("gene_counts","spatial_coordinates", og_query)
    local_query_data <- read.csv(paste0("/common/wonklab/synthetic_spatial/",og_query))

    local_ref_data <- local_ref_data[match(local_query$SpotID,local_ref_data$barcodes),]
    local_query_data <- local_query_data[match(local_query$OriginalCID,local_query_data$barcodes),]
    #local_query$cell_labels <- local_query_data$cell_labels

    ari <- adjustedRandIndex(local_ref_data$cell_labels, local_query_data$cell_labels)
    vi <- vi.dist(local_ref_data$cell_labels, local_query_data$cell_labels)
    accuracy <- sum(local_ref_data$cell_labels == local_query_data$cell_labels) / nrow(local_query_data)
    local_score <- data.frame("from" = paste("Reference", ref_tag),
        "to" = paste("Query", query_tag),
        "ari" = ari,
        "vi" = vi,
        "accuracy" = accuracy)
    map_score <- rbind(map_score, local_score)

}

map_score$method <- "CytoSpace"

write.csv(map_score, file = file_name)