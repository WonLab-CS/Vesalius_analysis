#-----------------------------------------------------------------------------#
# Loading Libraries and setting seed
#-----------------------------------------------------------------------------#
library(dplyr, lib.loc = "/common/martinp4/R")
library(tidyr)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(RColorBrewer)
set.seed(1547)

max_size <- 10000 * 1024^2
options(future.globals.maxSize = max_size)
#-----------------------------------------------------------------------------#
# Get args
#-----------------------------------------------------------------------------#
args <- commandArgs(TRUE)
input_matched <- args[1]
input_ref <- args[2]
data_type <- args[3]
output <- args[4]
ref <- args[5]
query <- args[6]

#-----------------------------------------------------------------------------#
# Get contribution files
#-----------------------------------------------------------------------------#
files <- list.files(
    path = input_matched,
    recursive = TRUE,
    pattern = "contribution_score.csv",
    full.names = TRUE)
tags <- list.files(
    path = input_matched,
    recursive = TRUE,
    pattern = "contribution_score.csv",
    full.names = FALSE)

files <- grep(
    pattern = data_type,
    x = files,
    value = TRUE)

tags <- grep(
    pattern = data_type,
    x = tags,
    value = TRUE)
tags <- gsub(
    paste0("Vesalius/report/Vesalius_aligned_",data_type,"_",ref,"_",query,"_"),
    "",
    tags)

tags <- gsub(
    "_contribution_score.csv",
    "",
    tags)

#-----------------------------------------------------------------------------#
# Load and plot contribution scores
#-----------------------------------------------------------------------------#
contributions <- vector("list", length(tags))
names(contributions) <- tags
for (i in seq_along(files)){
    local <- read.csv(files[i])
    local <- local[
        grep("total_cost",local$metric, invert = TRUE),
        c("metric","score","method")]
    local$combination <- tags[i]
    contributions[[i]] <- local
}

contributions <- do.call("rbind",contributions)
contributions <- split(contributions,contributions$method)

#-----------------------------------------------------------------------------#
# Plot 
#-----------------------------------------------------------------------------#
for (i in seq_along(contributions)){
    local <- contributions[[i]]
    levels <- c("f","n","c","t","fn","fc","ft","nc","nt","fnt","fnc","fny","fnct","fncty")
    local$combination <- as.factor(local$combination)
    local$combination <- factor(local$combination, levels = levels)
    if (names(contributions)[i] == "POC"){
        local$score <- local$score *100
        y <- "% of contribution"
    } else {
        y <- "Coef. of Variation"
    }
    
    max_cols <- length(unique(local$metric))
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
    g <- ggplot(local, aes(fill = metric,y = score, x=combination)) +
        geom_bar(position = "stack",stat = "identity") +
        scale_fill_manual(values = cols)+
        theme_bw()+
        labs(fill = "", y = y, x = "Combinations") +
        theme(strip.background =element_rect(fill="#082233ff"),
        strip.text = element_text(colour = 'white', size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size=15),
        legend.text = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1,vjust =1)) +
    facet_wrap(~method)
    file_name <- paste0(output, data_type,"_",names(contributions)[i],"_scores.pdf")
    pdf(file_name, width = 6, height = 5)
    print(g)
    dev.off()

}