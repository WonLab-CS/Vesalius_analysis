#-----------------------------------------------------------------------------#
# Loading Libraries and setting seed
#-----------------------------------------------------------------------------#
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(lpSolve, lib = "/common/martinp4/R")
library(mclust, lib.loc = "/common/martinp4/R")
library(mcclust, lib.loc = "/common/martinp4/R")


set.seed(1547)

#-----------------------------------------------------------------------------#
# set directories
#-----------------------------------------------------------------------------#
if (!dir.exists("/common/martinp4/benchmarking_out/Tangram/")) {
    dir.create("/common/martinp4/benchmarking_out/Tangram/")
}
output_data <- "/common/martinp4/benchmarking_out/Tangram/"

ref <- list.files("/common/wonklab/synthetic_spatial",
    pattern = "spatial_coordinates",
    full.names = TRUE)

query <- list.files(output_data,
    pattern = ".csv",
    full.names = TRUE)
query <- query[!grepl("mapping_accuracy",query)]

#-----------------------------------------------------------------------------#
# plots
#-----------------------------------------------------------------------------#
file_name <- paste0(output_data,"Tangram_output_plots.pdf")
pdf(file_name, width = 9, height = 4)

for (i in seq_along(query)) {
    local_ref <- sapply(strsplit(query[i],"sample_"),"[[",2)
    local_ref <- sapply(strsplit(local_ref,"_"),"[[",1)
    if (grepl("tinker", query[i])) {
        ref_tag <- "tinker_spatial"
    } else if (grepl("circle", query[i])) {
        ref_tag <- "circle_spatial"
    } else if (grepl("layered", query[i])) {
        ref_tag <- "layered_spatial"
    }
    query_tag <- sapply(strsplit(query[i],"sample_"),"[[",3)
    query_tag <- gsub(".csv","",query_tag)


    local_ref_data <- ref[grepl(ref_tag, ref) & grepl(paste0("sample_",local_ref,".csv"), ref)]
    local_ref_data <- read.csv(local_ref_data, header = TRUE)
    local_ref_data$Territory <- gsub("\\.","_",local_ref_data$Territory)
    local_query <- read.csv(query[i], header = TRUE)

    cols <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))
    cols <- cols(max(c(length(unique(local_ref_data$Territory)),length(unique(local_query$Territory)))))
    g_ref <- ggplot(local_ref_data, aes(x = x, y = y, col = as.factor(Territory))) +
    geom_point(size = 0.5) +
    theme_bw() +
    labs(title = paste("Reference", gsub("_spatial","",ref_tag), "sample",local_ref),color = "Cell Labels") +
    scale_color_manual(values = cols) +
    guides(colour = guide_legend(
        override.aes = list(size =  5)))

    
    g_query <- ggplot(local_query, aes(x = x, y = y, col = as.factor(Territory))) +
    geom_point(size = 0.5) +
    theme_bw() +
    labs(title = paste("Query", gsub("_spatial","",ref_tag), "sample",query_tag),color = "Cell Labels")+
    scale_color_manual(values = cols ) +
    guides(colour = guide_legend(
        override.aes = list(size =  5)))

    g_bind <- g_ref + g_query
    print(g_bind)


}
dev.off()



#-----------------------------------------------------------------------------#
# Accuracy table
#-----------------------------------------------------------------------------#

file_name <- paste0(output_data,"Tangram_mapping_accuracy.csv")
map_score <- data.frame("from" = character(),
    "to" = character(),
    "ari" = numeric(),
    "vi" = numeric(),
    "accuracy" = numeric())
for (i in seq_along(query)) {
    local_ref <- sapply(strsplit(query[i],"sample_"),"[[",2)
    local_ref <- sapply(strsplit(local_ref,"_"),"[[",1)
    if (grepl("tinker", query[i])) {
        ref_tag <- "tinker_spatial"
    } else if (grepl("circle", query[i])) {
        ref_tag <- "circle_spatial"
    } else if (grepl("layered", query[i])) {
        ref_tag <- "layered_spatial"
    }
    query_tag <- sapply(strsplit(query[i],"sample_"),"[[",3)
    query_tag <- gsub(".csv","",query_tag)


    local_ref_data <- ref[grepl(ref_tag, ref) & grepl(paste0("sample_",local_ref,".csv"), ref)]
    local_ref_data <- read.csv(local_ref_data, header = TRUE)
    local_query <- read.csv(query[i], header = TRUE)
    matched <- RANN::nn2(local_query[,c("x","y")], query = local_ref_data[,c("x","y")], k = 1)$nn.idx
    local_query <- local_query[matched, ]
    local_ref_data$Territory <- gsub("\\.","_",local_ref_data$Territory)
    ari <- adjustedRandIndex(local_ref_data$Territory, local_query$Territory)
    vi <- vi.dist(local_ref_data$Territory, local_query$Territory)
    accuracy <- sum(local_ref_data$Territory == local_query$Territory) / nrow(local_query)
    local_score <- data.frame("from" = paste("Reference", gsub("_spatial","",ref_tag), "sample",local_ref),
        "to" = paste("Query", gsub("_spatial","",ref_tag), "sample",query_tag),
        "ari" = ari,
        "vi" = vi,
        "accuracy" = accuracy)
    map_score <- rbind(map_score, local_score)

}


map_score$method <- "Tangram"

write.csv(map_score, file = file_name)