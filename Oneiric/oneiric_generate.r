library(ggplot2)
library(RColorBrewer)
library(oneiric, lib.loc = "/common/martinp4/R")
library(scater,lib.loc = "/common/martinp4/R")
args <- commandArgs(TRUE)
output <- args[1]
seed <- as.numeric(args[2])
generate_sim_data(output = output,
    seed = seed,
    plot = TRUE,
    run_mem = TRUE,
    simple = TRUE)


### Plot

files <- list.files(path = output,
    pattern = "spatial_coordinates",
    full.names = TRUE)
tags <- list.files(path = output,
    pattern = "spatial_coordinates",
    full.names = FALSE)
tags <- sapply(strsplit(tags, "_spatial_coordinates_sample_"),"[[",1)

files <- mapply(function(f,t){
    df <- read.csv(f)
    df$Data_Type <- t
    return(df)
}, files, tags, SIMPLIFY = FALSE)

files <- do.call("rbind", files)
files <- split(files, files$Data_Type)

for (i in seq_along(files)) {
    local <- files[[i]]
    local$sample <- as.factor(local$sample)
    local$sample <- factor(local$sample, levels = as.character(1:12))
    max_cols <- length(unique(local$cell_labels))
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
     

    p <- ggplot(local, aes(x,y, col = as.factor(cell_labels))) +
        geom_point(size = 1, alpha=0.75) +
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
    data <- unique(local$Data_Type)
    file_name <- paste0(output,data, "_synthetic_territorries.pdf")
    pdf(file_name, width = 10.5, height = 8)
    print(p)
    dev.off()

}


   