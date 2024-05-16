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
score <- read.csv("/common/martinp4/benchmarking_out/Vesalius/report/combination_score.csv",header= T, skip=1)
score <- pivot_longer(score, cols = c("ARI","VI"),names_to = "Method", values_to = "Score")
locs <- score$ref_sample != score$query_sample
score <- score[locs,]


score$Combination <- as.factor(score$Combination)
score$Combination <- factor(score$Combination ,
    levels = c("feature",
        "niche",
        "territory",
        "feature_niche",
        "feature_territory",
        "niche_territory",
        "feature_niche_territory"))
#-----------------------------------------------------------------------------#
# plot 
#-----------------------------------------------------------------------------#
cols <- c("#082233ff","#1a6ea7ff","#146269ff","#fdcf72ff","#b75b1cff","#8e0703ff","#bfdcc3ff","#25b4c1ff","#bf4c40ff","#dd9a94ff")
ar <- score[score$Method == "ARI" & score$Regime == "circle",]
ari <- ggplot(ar, aes(x = Combination, y = Score, fill = Combination)) +
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
file_out <- paste0(output, "combination_ARI_score_circle.pdf")
pdf(file_out, width = 6, height = 6)
print(ari)
dev.off()

vi <- score[score$Method == "VI"& score$Regime == "circle",]
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
file_out <- paste0(output, "combination_VI_score_circle.pdf")
pdf(file_out, width = 6, height = 6)
print(vis)
dev.off()

cols <- c("#082233ff","#1a6ea7ff","#146269ff","#fdcf72ff","#b75b1cff","#8e0703ff","#bfdcc3ff","#25b4c1ff","#bf4c40ff","#dd9a94ff")
ar <- score[score$Method == "ARI" & score$Regime == "layered",]
ari <- ggplot(ar, aes(x = Combination, y = Score, fill = Combination)) +
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
file_out <- paste0(output, "combination_ARI_score_layered.pdf")
pdf(file_out, width = 6, height = 6)
print(ari)
dev.off()

vi <- score[score$Method == "VI"& score$Regime == "layered",]
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
file_out <- paste0(output, "combination_VI_score_layered.pdf")
pdf(file_out, width = 6, height = 6)
print(vis)
dev.off()