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
# Get args
#-----------------------------------------------------------------------------#
args <- commandArgs(TRUE)
input <- output <- args[1]
data_type <- args[2]

#-----------------------------------------------------------------------------#
# Load files
#-----------------------------------------------------------------------------#
scores <- list.files(
    input,
    pattern = "benchmarking_scores.csv",
    full.names = TRUE)
scores <- grep(data_type,scores, value = TRUE)
scores <- as.data.frame(read.csv(scores))
scores$Data_Type <- data_type


ter_col <- length(unique(scores$Method))
#base_colors <- c("#08324e","#1a6ea7ff","#146269ff","#51b090","#de9c71","#b75b1cff","#8e0703ff")
base_colours <- c(
      "#E69F00",
      "#56B4E9",
      "#009E73",
      "#F0E442",
      "#0072B2",
      "#D55E00",
      "#CC79A7",
      "#999999")
if (ter_col < length(base_colours)) {
      ter_pal <- colorRampPalette(base_colours[seq(1, ter_col)])
} else {
      ter_pal <- colorRampPalette(base_colours)
}
cols <- ter_pal(ter_col)
names(cols) <- scores$Method
#-----------------------------------------------------------------------------#
# Showing all scores
#-----------------------------------------------------------------------------#
ord <- order(scores$ARI_cell, decreasing = TRUE)
new_order <- scores$Method[ord]
scores$Method <- as.factor(scores$Method)
scores$Method <- factor(scores$Method, levels = new_order)
cols <- cols[new_order]

ari_cell <- ggplot(scores, aes(x = Method, y = ARI_cell,fill = Method)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values = cols) +
    theme_bw() +
    theme(strip.background =element_rect(fill="#082233ff"),
        strip.text = element_text(colour = 'white', size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size=15),
        legend.text = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1,vjust =1)) +
    facet_wrap(~Data_Type)
file_out <- paste0(output,data_type, "_benchmarking_ARI_cell.pdf")
pdf(file_out, width = 9, height = 6.5)
print(ari_cell)
dev.off()


ord <- order(scores$ARI_interactions, decreasing = TRUE)
new_order <- scores$Method[ord]
scores$Method <- as.factor(scores$Method)
scores$Method <- factor(scores$Method, levels = new_order)
cols <- cols[new_order]

ari_inter <- ggplot(scores, aes(x = Method, y = ARI_interactions, fill = Method)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values = cols) +
    theme_bw() +
    theme(strip.background =element_rect(fill="#082233ff"),
        strip.text = element_text(colour = 'white', size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size=15),
        legend.text = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1,vjust =1)) +
    facet_wrap(~Data_Type) 
    
file_out <- paste0(output,data_type, "_benchmarking_ARI_inter.pdf")
pdf(file_out, width = 9, height = 6.5)
print(ari_inter)
dev.off()

ord <- order(scores$VI_cell, decreasing = FALSE)
new_order <- scores$Method[ord]
scores$Method <- as.factor(scores$Method)
scores$Method <- factor(scores$Method, levels = new_order)
cols <- cols[new_order]
vi_cell <- ggplot(scores, aes(x = Method, y = VI_cell, fill = Method)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values = cols) +
    theme_bw() +
    theme(strip.background =element_rect(fill="#082233ff"),
        strip.text = element_text(colour = 'white', size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size=15),
        legend.text = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1,vjust =1))+
    facet_wrap(~Data_Type) 
file_out <- paste0(output,data_type, "_benchmarking_VI_cell.pdf")
pdf(file_out, width = 9, height = 6.5)
print(vi_cell)
dev.off()


ord <- order(scores$VI_interactions, decreasing = FALSE)
new_order <- scores$Method[ord]
scores$Method <- as.factor(scores$Method)
scores$Method <- factor(scores$Method, levels = new_order)
cols <- cols[new_order]
vi_inter <- ggplot(scores, aes(x = Method, y = VI_interactions, fill = Method)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values = cols) +
    theme_bw() +
    theme(strip.background =element_rect(fill="#082233ff"),
        strip.text = element_text(colour = 'white', size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size=15),
        legend.text = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1,vjust =1)) +
    facet_wrap(~Data_Type) 
file_out <- paste0(output,data_type, "_benchmarking_VI_interactions.pdf")
pdf(file_out, width = 9, height = 6.5)
print(vi_cell)
dev.off()


ord <- order(scores$JI_interactions, decreasing = TRUE)
new_order <- scores$Method[ord]
scores$Method <- as.factor(scores$Method)
scores$Method <- factor(scores$Method, levels = new_order)
cols <- cols[new_order]
ji_inter <- ggplot(scores, aes(x = Method, y = JI_interactions, fill = Method)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values = cols) +
    theme_bw() +
    theme(strip.background =element_rect(fill="#082233ff"),
        strip.text = element_text(colour = 'white', size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size=15),
        legend.text = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1,vjust =1)) +
    facet_wrap(~Data_Type)
file_out <- paste0(output,data_type, "_benchmarking_JI_interactions.pdf")
pdf(file_out, width = 9, height = 6.5)
print(ji_inter)
dev.off()



