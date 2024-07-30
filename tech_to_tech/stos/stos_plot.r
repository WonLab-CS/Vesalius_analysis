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
set.seed(1547)
#-----------------------------------------------------------------------------#
# Set future global for multicore processing
#-----------------------------------------------------------------------------#
if (!dir.exists("/common/martinp4/stos/report/stos/")) {
    dir.create("/common/martinp4/stos/report/stos/")
}
output <- "/common/martinp4/stos/report/stos/"
args <- commandArgs(TRUE)
idx <- as.numeric(args[1])

use_cost <- c("feature", "niche")
slices <- c( "embryo1","embryo2","embryo3")
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

