#-----------------------------------------------------------------------------#
# Loading Libraries and setting seed
#-----------------------------------------------------------------------------#
library(dplyr, lib.loc = "/common/martinp4/R")
library(tidyr)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(RColorBrewer)
library(vesalius, lib.loc = "/common/martinp4/R")
set.seed(1547)

#-----------------------------------------------------------------------------#
# set and load
#-----------------------------------------------------------------------------#
max_size <- 10000 * 1024^2
output <- "/common/martinp4/benchmarking_out/Vesalius/report/"
score <- read.csv("/common/martinp4/benchmarking_out/Vesalius/bio_data/seqFISH/combination_score_seqFISH.csv",header= T, skip=1)
score <- pivot_longer(score, cols = c("ARI","VI"),names_to = "Method", values_to = "Score")

use_cost <- c(
    c("feature"),
    c("niche"),
    c("composition"),
    c("cell_type"),
    c("feature_niche"),
    c("feature_niche_composition"),
    c("feature_niche_composition_cell_type"),
    c("feature_composition"),
    c("niche_composition"))

score$Combination <- as.factor(score$Combination)
score$Combination <- factor(score$Combination ,
    levels = use_cost)


#-----------------------------------------------------------------------------#
# plot 
#-----------------------------------------------------------------------------#
cols <- c("#082233ff","#1a6ea7ff","#146269ff","#fdcf72ff","#b75b1cff","#8e0703ff","#bfdcc3ff","#25b4c1ff","#bf4c40ff","#dd9a94ff")
ar <- score[score$Method == "ARI" ,]
ari <- ggplot(ar, aes(x = Combination, y = Score, fill = Combination)) +
    geom_bar() +
    scale_fill_manual(values = rep(cols[1],9)) +
    theme_bw() +
    theme(strip.background =element_rect(fill="#082233ff"),
        strip.text = element_text(colour = 'white', size = 20),
        axis.title = element_text(size = 15),
        legend.title = element_text(size=15),
        legend.text = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 20, angle = 45, hjust = 1,vjust =1),
        legend.position = "none") +
    facet_wrap(~Method)
file_out <- paste0(output, "combination_ARI_score_seqFISH.pdf")
pdf(file_out, width = 6, height = 6)
print(ari)
dev.off()

vi <- score[score$Method == "VI",]
vis <- ggplot(vi, aes(x = Combination, y = Score, fill = Combination)) +
    geom_violin() +
    scale_fill_manual(values = rep(cols[1],9)) +
    theme_bw() +
    theme(strip.background =element_rect(fill="#082233ff"),
        strip.text = element_text(colour = 'white', size = 20),
        axis.title = element_text(size = 15),
        legend.title = element_text(size=15),
        legend.text = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 20, angle = 45, hjust = 1,vjust =1),
        legend.position = "none") +
    facet_wrap(~Method)
file_out <- paste0(output, "combination_VI_score_seqFISH.pdf")
pdf(file_out, width = 6, height = 6)
print(vis)
dev.off()
