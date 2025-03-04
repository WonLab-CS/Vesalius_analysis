#-----------------------------------------------------------------------------#
# Loading Libraries and setting seed
#-----------------------------------------------------------------------------#
library(dplyr, lib.loc = "/common/martinp4/R")
library(ggplot2)
library(patchwork)
library(ggpubr)
library(RColorBrewer)


set.seed(1547)
max_size <- 10000 * 1024^2
options(future.globals.maxSize = max_size)
#-----------------------------------------------------------------------------#
# Output set up 
#-----------------------------------------------------------------------------#
args <- commandArgs(TRUE)
input_matched <- args[1]
input_ref <- args[2]
data_type <- args[3]
output <- args[4]
ref <- args[5]
query <- args[6]


#-----------------------------------------------------------------------------#
# utils
#-----------------------------------------------------------------------------#
get_files <- function(
    input,
    data_type) {
    dirs <- list.dirs(input, recursive = FALSE)
    dirs <- grep("report", dirs, value = TRUE, invert = TRUE)
    files <- list.files(dirs,
        pattern = ".csv",
        recursive = TRUE,
        full.names = TRUE)
    file_tags <- list.files(dirs,
        pattern = ".csv",
        recursive = TRUE,
        full.names = FALSE)

    files <- grep(data_type, files, value = TRUE)
    file_tags <- grep(data_type, file_tags, value = TRUE)
    files <- grep("contribution_score", files, value = TRUE, invert = TRUE)
    file_tags <- grep("contribution_score", file_tags, value = TRUE, invert = TRUE)
    files <- grep("computational_performance", files, value = TRUE, invert = TRUE)
    file_tags <- grep("computational_performance", file_tags, value = TRUE, invert = TRUE)

    file_tags <- gsub(".csv","",file_tags)
    file_tags <- gsub("_aligned","",file_tags)
    file_tags <- gsub("report/","", file_tags)
    parts <- strsplit(file_tags,"_")
    file_tags <- sapply(parts, function(p){
        method <- p[1]
        if (method == "Vesalius" && !"contribution_score" %in% p) {
            combi <- p[length(p)]
            method <- paste0(method,"_", combi)
            return(method)
        } else if (method == "Vesalius" && "contribution_score" %in% p) {
            combi <- p[length(p) - 1]
            method <- paste0(method,"_", combi)
            return(method)
        } else if (method == "CytoSpace") {
            if (grepl("noLab" , p[2])) {
                method <- paste0(method,"_noLab")
            }
            return(method)
        } else {
            return(method) 
        }
        
    })
    names(files) <- file_tags
    return(files)
}

scale_coordinates <- function(ref,query){
    ref$x <- ref$x - min(ref$x)
    ref$y <- ref$y - min(ref$y)
    x <- max(query$x) / max(ref$x)
    y <- max(query$y) / max(ref$y)
    ref$x <- ref$x * x
    ref$y <- ref$y * y
    return(ref)
}
load_files <- function(mapped) {
    tags <- names(mapped)
    mapped <- lapply(mapped,read.csv)
    mapped <- mapply(function(m,t){
        m$sample <- t
        m <- m[c("barcodes","x","y","sample","cell_labels","interactions")]
        return(m)    
    }, m = mapped, t = tags, SIMPLIFY = FALSE)
    return(mapped)
}


#-----------------------------------------------------------------------------#
# Get ref files
#-----------------------------------------------------------------------------#
ref <- list.files(
    input_ref,
    ref,
    full.names = TRUE)
ref <- grep(x = ref, pattern = data_type, value = TRUE)
ref <- grep(x = ref, pattern = "spatial_coordinates", value = TRUE)
ref <- read.csv(ref)
ref$sample <- "Reference"

query <- list.files(
    input_ref,
    query,
    full.names = TRUE)
query <- grep(x = query, pattern = data_type, value = TRUE)
query <- grep(x = query, pattern = "spatial_coordinates", value = TRUE)
query <- read.csv(query)
query$sample <- "Query"


ref <- scale_coordinates(ref,query)

#-----------------------------------------------------------------------------#
# Get mapped files
#-----------------------------------------------------------------------------#
mapped <- get_files(
    input_matched,
    data_type)
mapped <- load_files(mapped)
mapped <- lapply(mapped, scale_coordinates, query)
mapped <- do.call("rbind",mapped)
#-----------------------------------------------------------------------------#
# Mergde and prepare for plotting
#-----------------------------------------------------------------------------#
scores <- list.files(
    input_matched,
    recursive = TRUE,
    pattern = "benchmarking_scores.csv",
    full.names = TRUE)
scores <- grep(data_type,scores, value = TRUE)
scores <- as.data.frame(read.csv(scores))
scores <- scores[grep(x = scores$Method,pattern ="Vesalius_fncty|Vesalius_fnct|Vesalius_fnt"),]
scores <- scores$Method[order(scores$ARI_cell, decreasing = TRUE)[1:3]]

all_methods <- rbind(ref, query, mapped)
levels <- c("Reference", "Query",scores,"CytoSpace","CytoSpace_noLab","Tangram","Scanorama","SLAT","GPSA","PASTE")
all_methods <- all_methods[all_methods$sample %in% levels,]
all_methods$sample <- as.factor(all_methods$sample)
all_methods$sample <- factor(all_methods$sample,levels = levels)

#-----------------------------------------------------------------------------#
# Plotting
#-----------------------------------------------------------------------------#
max_cols <- length(unique(all_methods$cell_labels))
base_colours <- c(
        "#E69F00",
        "#56B4E9",
        "#009E73",
        "#F0E442",
        "#0072B2",
        "#D55E00",
        "#CC79A7",
        "#999999")
if (max_cols < length(base_colours)) {
    ter_pal <- colorRampPalette(base_colours[seq(1, max_cols)])
} else {
    ter_pal <- colorRampPalette(base_colours)
}
cols <- ter_pal(max_cols)

p <- ggplot(all_methods, aes(x,y, col = cell_labels)) +
    geom_point(size = 0.2) +
    scale_color_manual(values = cols) +
    theme_bw() + 
    theme(strip.background =element_rect(fill="#082233ff"),
        strip.text = element_text(colour = 'white', size = 15),
        legend.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.position = "bottom") +
    labs(color = "") + 
    guides(colour = guide_legend(
            override.aes = list(size = 5)))+
    facet_wrap(~sample)

file_name <- paste0(output,data_type, "_mapped_events.pdf")
pdf(file_name, width = 14.5, height = 14)
print(p)
dev.off()
#q("no")

