library(dplyr)
library(future)
library(future.apply)
library(ggplot2)
library(dplyr)
library(patchwork)
library(Matrix)
library(deldir)
library(imager)
library(imagerExtra)
library(lpSolve, lib = "/common/martinp4/R")
library(TreeDist, lib.loc = "/common/martinp4/R")
library(mclust, lib.loc = "/common/martinp4/R")
library(pwr, lib.loc = "/common/martinp4/R")
library(kohonen, lib.loc = "/common/martinp4/R")
library(registry, lib.loc = "/common/martinp4/R")
library(rngtools, lib.loc = "/common/martinp4/R")
library(NMF, lib.loc = "/common/martinp4/R")
library(spatstat.utils, lib.loc = "/common/martinp4/R")
library(vesalius, lib.loc = "/common/martinp4/R")
library(RColorBrewer)
library(ggalluvial, lib.loc = "/common/martinp4/R")
set.seed(1547)
#-----------------------------------------------------------------------------#
# Set future global for multicore processing
#-----------------------------------------------------------------------------#

args <- commandArgs(TRUE)
idx <- as.numeric(args[1])
input <- args[2]
output <- args[3]


slices <- c( "embryo1","embryo2","embryo3")
#-----------------------------------------------------------------------------#
# Utility functions
#-----------------------------------------------------------------------------#

get_acc <- function(rib_df, seqFISH, stereo) {
    seq <- rib_df$Stereo_seq_Cells[
            rib_df$seqFISH_Cells %in% seqFISH] %in% stereo
    stereo <- rib_df$seqFISH_Cells[
            rib_df$Stereo_seq_Cells %in% stereo] %in% seqFISH

    seq <- sum(seq) / (sum(seq) + sum(!seq))
    stereo <- sum(stereo) / (sum(stereo) + sum(!stereo))
    return(data.frame("seq_acc" = seq, "stereo_acc" = stereo))
}
#-----------------------------------------------------------------------------#
# Load data
#-----------------------------------------------------------------------------#
files <- list.files(output, pattern = ".rds", full.names = TRUE)
stereo <- readRDS(grep("stereo_seq", files, value = TRUE))
fish <- files[grepl("seqFISH", files) & !grepl("_to_", files)]
fish <- readRDS(grep(slices[idx],fish, value = TRUE))
matched <- grep("matched", files, value = TRUE)
matched <- readRDS(grep(slices[idx], matched, value = TRUE))
reverse <- grep("reverse", files, value = TRUE)
reverse <- readRDS(grep(slices[idx], reverse, value = TRUE))


#-----------------------------------------------------------------------------#
# plotting
#-----------------------------------------------------------------------------#

g1 <- view_mapping_metrics(matched, cex_pt = 1.5, trial = "feature") +
    theme_void() +
    labs(title = "Feature")
g2 <- view_mapping_metrics(matched, cex_pt = 1.5, trial = "niche") +
    theme_void()+
    labs(title = "Niche")
g3 <- view_mapping_metrics(matched, cex_pt = 1.5, trial = "cost") + 
    theme_void()+
    labs(title = "Cost")

file_out <- paste0(output, "mapping_scores_seqFISH_",slices[idx],"_to_stereo.pdf")
pdf(file_out, width = 24, height = 8)
all <- g1 + g2 + g3
print(all)
dev.off()


ter_1 <- territory_plot(stereo, trial = "Cells",randomise = FALSE, cex_pt = 3, alpha = 1) +
    #scale_color_manual(values = colors) +
    theme_void() +
    theme(legend.title = element_text(size = 20),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.position = "left") +
    labs(color = "", title = "Stereo-Seq")

ter_2 <- territory_plot(fish, trial = "Cells",randomise = FALSE, cex_pt = 1, alpha = 1) +
    #scale_color_manual(values = colors) +
    theme_void() +
    theme(legend.title = element_text(size = 20),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.position = "left") +
    labs(color = "", title = paste0("seqFISH ", slices[idx]))

ter_3 <- territory_plot(matched, trial = "Cells",randomise = FALSE, cex_pt = 1.3, alpha = 1) +
    #scale_color_manual(values = colors) +
    theme_void() +
    theme(legend.title = element_text(size = 20),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.position = "left") +
    labs(color = "", title = paste0("seqFISH ", slices[idx], " to Stereo-Seq"))

file_out <- paste0(output, "mapped_cells_seqFISH_",slices[idx],"_to_stereo.pdf")
pdf(file_out, width = 24, height = 8)
print(ter_1 + ter_2 + ter_3)
dev.off()


#-----------------------------------------------------------------------------#
# indiv plots
#-----------------------------------------------------------------------------#


file_out <- paste0(output, "stereo_cells.pdf")
pdf(file_out, width = 12, height = 8)
print(ter_1)
dev.off()
file_out <- paste0(output, "seqFISH_",slices[idx],"_cells.pdf")
pdf(file_out, width = 12, height = 8)
print(ter_2)
dev.off()
file_out <- paste0(output, "seqFISH_",slices[idx],"_to_stereo_cells.pdf")
pdf(file_out, width = 12, height = 8)
print(ter_3)
dev.off()


file_out <- paste0(output, "seqFISH_",slices[idx],"_to_stereo_feature_score.pdf")
pdf(file_out, width = 8, height = 7)
print(g1)
dev.off()
file_out <- paste0(output, "seqFISH_",slices[idx],"_to_stereo_niche_score.pdf")
pdf(file_out, width = 8, height = 7)
print(g2)
dev.off()
file_out <- paste0(output, "seqFISH_",slices[idx],"_to_stereo_cost_score.pdf")
pdf(file_out, width = 8, height = 7)
print(g3)
dev.off()


#-----------------------------------------------------------------------------#
# plotting - reverse
#-----------------------------------------------------------------------------#

g1 <- view_mapping_metrics(matched, cex_pt = 1.5, trial = "feature") +
    theme_void() +
    labs(title = "Feature")
g2 <- view_mapping_metrics(matched, cex_pt = 1.5, trial = "niche") +
    theme_void()+
    labs(title = "Niche")
g3 <- view_mapping_metrics(matched, cex_pt = 1.5, trial = "cost") + 
    theme_void()+
    labs(title = "Cost")

file_out <- paste0(output, "mapping_scores_stereo_to_seqFISH_",slices[idx],".pdf")
pdf(file_out, width = 24, height = 8)
all <- g1 + g2 + g3
print(all)
dev.off()


ter_1 <- territory_plot(stereo, trial = "Cells",randomise = FALSE, cex_pt = 1.5, alpha = 1) +
    #scale_color_manual(values = colors) +
    theme_void() +
    theme(legend.title = element_text(size = 20),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.position = "left") +
    labs(color = "", title = "Stereo-Seq")

ter_2 <- territory_plot(fish, trial = "Cells",randomise = FALSE, cex_pt = 1, alpha = 1) +
    #scale_color_manual(values = colors) +
    theme_void() +
    theme(legend.title = element_text(size = 20),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.position = "left") +
    labs(color = "", title = paste0("seqFISH ", slices[idx]))

ter_3 <- territory_plot(reverse, trial = "Cells",randomise = FALSE, cex_pt = 1.3, alpha = 1) +
    #scale_color_manual(values = colors) +
    theme_void() +
    theme(legend.title = element_text(size = 20),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.position = "left") +
    labs(color = "", title = paste0("Stereo-seq to seqFISH ", slices[idx]))

file_out <- paste0(output, "mapped_cells_stereo_to_seqFISH_",slices[idx],".pdf")
pdf(file_out, width = 24, height = 8)
print(ter_2 + ter_1 + ter_3)
dev.off()


#-----------------------------------------------------------------------------#
# indiv plots
#-----------------------------------------------------------------------------#


file_out <- paste0(output, "stereo_cells.pdf")
pdf(file_out, width = 12, height = 8)
print(ter_1)
dev.off()
file_out <- paste0(output, "seqFISH_",slices[idx],"_cells.pdf")
pdf(file_out, width = 12, height = 8)
print(ter_2)
dev.off()
file_out <- paste0(output, "stereo_to_seqFISH_",slices[idx],"_cells.pdf")
pdf(file_out, width = 12, height = 8)
print(ter_3)
dev.off()


file_out <- paste0(output, "stereo_to_seqFISH_",slices[idx],"_feature_score.pdf")
pdf(file_out, width = 8, height = 7)
print(g1)
dev.off()
file_out <- paste0(output, "stereo_to_seqFISH_",slices[idx],"_niche_score.pdf")
pdf(file_out, width = 8, height = 7)
print(g2)
dev.off()
file_out <- paste0(output, "stereo_to_seqFISH_",slices[idx],"_cost_score.pdf")
pdf(file_out, width = 8, height = 7)
print(g3)
dev.off()
#-----------------------------------------------------------------------------#
# Set brain cell labels for accuracy 
#-----------------------------------------------------------------------------#
seq_brain <- c("Forebrain/Midbrain/Hindbrain","Spinal cord")
stereo_brain <- c("Brain","Notochord")

#-----------------------------------------------------------------------------#
# Alluvial plot 
#-----------------------------------------------------------------------------#
from <- matched@map$from 
sf_cells <- fish@territories$Cells[match(from, fish@territories$barcodes)]
to <- sapply(sapply(strsplit(matched@map$to, "-"),"[",1:2,simplify = F),paste0, collapse ="-")
st_cells <- stereo@territories$Cells[match(to, stereo@territories$barcodes)]
rib_df <- data.frame("seqFISH_Cells" = sf_cells,
    "Stereo_seq_Cells" = st_cells,  
    Freq = rep(1, length(st_cells)))

sf_tabs <- table(sf_cells) > 200
rib_df <- rib_df[rib_df$seqFISH_Cells %in% names(sf_tabs)[sf_tabs],]

acc_brain <- get_acc(rib_df, seq_brain, stereo_brain)
acc_brain$from <- slices[idx]
acc_brain$to <- "Stereo"
print(acc_brain)

labels <- sort(unique(rib_df$Stereo_seq_Cells))
ter_col <- length(labels)
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

colors <- ter_pal(ter_col)

rib <- ggplot(rib_df, aes(axis1 = seqFISH_Cells, axis2 = Stereo_seq_Cells, y = Freq)) +
    geom_alluvium(aes(fill = Stereo_seq_Cells),
        curve_type = "sigmoid") +
    geom_stratum() +
    geom_text(stat = "stratum", aes(label = paste(after_stat(stratum))), size = 4)+
    theme_void() + 
    theme(legend.text = element_text(size=20),
        legend.title = element_text(size = 20),
        legend.key.size = unit(1.5, "cm")) +
    scale_fill_manual(values = colors) +
    labs(fill = "Stereo-seq Cell Types") +
    guides(fill = guide_legend(
            override.aes = list(size = 2)))

file_out <- paste0(output, "seqFISH_",slices[idx],"_to_stereo_alluvial_cells.pdf")
pdf(file_out, width = 14, height = 10)
print(rib)
dev.off()

file_out <- paste0(output, "seqFISH_",slices[idx],"_to_stereo_acc.csv")
write.csv(acc_brain, file = file_out)



#-----------------------------------------------------------------------------#
# Alluvial plot rev
#-----------------------------------------------------------------------------#
from <- reverse@map$to 
sf_cells <- fish@territories$Cells[match(from, fish@territories$barcodes)]
to <- sapply(sapply(strsplit(reverse@map$from, "-"),"[",1:2,simplify = F),paste0, collapse ="-")
st_cells <- stereo@territories$Cells[match(to, stereo@territories$barcodes)]
rib_df <- data.frame("seqFISH_Cells" = sf_cells,
    "Stereo_seq_Cells" = st_cells,  
    Freq = rep(1, length(st_cells)))

sf_tabs <- table(sf_cells) > 200
rib_df <- rib_df[rib_df$seqFISH_Cells %in% names(sf_tabs)[sf_tabs],]

acc_brain <- get_acc(rib_df, seq_brain, stereo_brain)
acc_brain$from <- "Stereo"
acc_brain$to <- slices[idx]
print(acc_brain)

labels <- sort(unique(rib_df$seqFISH_Cells))
ter_col <- length(labels)
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

colors <- ter_pal(ter_col)

rib <- ggplot(rib_df, aes(axis1 = Stereo_seq_Cells, axis2 = seqFISH_Cells, y = Freq)) +
    geom_alluvium(aes(fill = seqFISH_Cells),
        curve_type = "sigmoid") +
    geom_stratum() +
    geom_text(stat = "stratum", aes(label = paste(after_stat(stratum))), size = 4)+
    theme_void() + 
    theme(legend.text = element_text(size=20),
        legend.title = element_text(size = 20),
        legend.key.size = unit(1.1, "cm")) +
    scale_fill_manual(values = colors) +
    labs(fill = "SeqFISH Cell Types") +
    guides(fill = guide_legend(
            override.aes = list(size = 1.5)))

file_out <- paste0(output, "stereo_to_seqFISH_",slices[idx],"_alluvial_cells.pdf")
pdf(file_out, width = 14, height = 10)
print(rib)
dev.off()
file_out <- paste0(output, "stereo_to_seqFISH_",slices[idx],"_acc.csv")
write.csv(acc_brain, file = file_out)

#-----------------------------------------------------------------------------#
# Get accuracy of broad cell types
#-----------------------------------------------------------------------------#