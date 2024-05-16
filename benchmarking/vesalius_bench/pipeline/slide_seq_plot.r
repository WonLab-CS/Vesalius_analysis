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
max_size <- 10000 * 1024^2
options(future.globals.maxSize = max_size)

if (!dir.exists("/common/martinp4/benchmarking_out/Vesalius/bio_data/SSV2/")) {
    dir.create("/common/martinp4/benchmarking_out/Vesalius/bio_data/SSV2/")
}
output <- "/common/martinp4/benchmarking_out/Vesalius/bio_data/SSV2/"
#-----------------------------------------------------------------------------#
# Load data
#-----------------------------------------------------------------------------#
files <- list.files(output, pattern = ".rds", full.names = TRUE)
seed <- readRDS(grep("seed", files, value = TRUE))
query <- readRDS(grep("query", files, value = TRUE))
matched <- readRDS(grep("matched", files, value = TRUE))

tags <- list.files(output, pattern = ".rds", full.names = FALSE)
seed_tag <- grep("_seed.rds",tags, value = TRUE)
seed_tag <- gsub("_seed.rds","", seed_tag)

query_tag <- grep("_query.rds",tags, value = TRUE)
query_tag <- gsub("_query.rds","", query_tag)

#-----------------------------------------------------------------------------#
# plotting
#-----------------------------------------------------------------------------#
g1 <- image_plot(seed, embedding = "PCA") +
    theme_void() +
    theme(legend.title = element_text(size = 20),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.position = "none") +
    labs(color = "", title = paste0("Reference - ", seed_tag))
g2 <- image_plot(query, embedding = "PCA") +
    theme_void() +
    theme(legend.title = element_text(size = 20),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.position = "none") +
    labs(color = "", title =  paste0("Query - ", query_tag))
g3 <- image_plot(matched, embedding = "PCA") +
    theme_void() +
    theme(legend.title = element_text(size = 20),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.position = "none") +
    labs(color = "", title =  paste0(query_tag, " mapped to ", seed_tag))

file_out <- paste0(output, "mapped_embeddings_",query_tag,"_to_",seed_tag,".pdf")
pdf(file_out, width = 24, height = 8)
print(g1 + g2 + g3)
dev.off()




g4 <- view_mapping_metrics(matched, cex_pt = 1.5, trial = "feature") +
    theme_void() +
    labs(title = "Feature")
g5 <- view_mapping_metrics(matched, cex_pt = 1.5, trial = "niche") +
    theme_void()+
    labs(title = "Niche")
g6 <- view_mapping_metrics(matched, cex_pt = 1.5, trial = "territory") + 
    theme_void()+
    labs(title = "Territory")
g7 <- view_mapping_metrics(matched, cex_pt = 1.5, trial = "cost") + 
    theme_void()+
    labs(title = "Cost")

file_out <- paste0(output, "mapping_scores_",query_tag,"_to_",seed_tag,".pdf")
pdf(file_out, width = 24, height = 24)
all <- (g7 + g4) / (g5 + g6)
print(all)
dev.off()


ter_1 <- territory_plot(seed, trial = "Territory", cex_pt = 1, alpha = 1) +
    #scale_color_manual(values = colors) +
    theme_void() +
    theme(legend.title = element_text(size = 20),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.position = "none") +
    labs(color = "", title = paste0("Reference - ", seed_tag))

ter_2 <- territory_plot(query,trial = "Territory",  cex_pt = 1, alpha=1) +
    #scale_color_manual(values = colors) +
    theme_void() +
    theme(legend.title = element_text(size = 20),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.position = "none") +
    labs(color = "", title = paste0("Query - ", query_tag))

ter_3 <- territory_plot(matched, trial = "Territory",  cex_pt = 2.1, alpha=1) +
    #scale_color_manual(values = colors) +
    theme_void() +
    theme(legend.title = element_text(size = 20),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.position = "none") +
    labs(color = "", title = paste0(query_tag, " mapped to ", seed_tag))

file_out <- paste0(output, "mapped_territories_",query_tag,"_to_",seed_tag,".pdf")
pdf(file_out, width = 24, height = 8)
print(ter_1 + ter_2 + ter_3)
dev.off()


#-----------------------------------------------------------------------------#
# indiv plots
#-----------------------------------------------------------------------------#

file_out <- paste0(output, seed_tag,"_embbeding.pdf")
pdf(file_out, width = 12, height = 12)
print(g1)
dev.off()
file_out <- paste0(output, query_tag,"_embbeding.pdf")
pdf(file_out, width = 12, height = 12)
print(g2)
dev.off()
file_out <- paste0(output, query_tag,"_to_",seed_tag,"_embbeding.pdf")
pdf(file_out, width = 12, height = 12)
print(g3)
dev.off()

file_out <- paste0(output, seed_tag,"_territories.pdf")
pdf(file_out, width = 8, height = 8)
print(ter_1)
dev.off()
file_out <- paste0(output, query_tag,"_territories.pdf")
pdf(file_out, width = 8, height = 8)
print(ter_2)
dev.off()
file_out <- paste0(output, query_tag,"_to_",seed_tag,"_territories.pdf")
pdf(file_out, width = 8, height = 8)
print(ter_3)
dev.off()


file_out <- paste0(output, query_tag,"_to_",seed_tag,"_feature_score.pdf")
pdf(file_out, width = 8, height = 6)
print(g4)
dev.off()
file_out <- paste0(output, query_tag,"_to_",seed_tag,"_niche_score.pdf")
pdf(file_out, width = 8, height = 6)
print(g5)
dev.off()
file_out <- paste0(output, query_tag,"_to_",seed_tag,"_territory_score.pdf")
pdf(file_out, width = 8, height = 6)
print(g6)
dev.off()
file_out <- paste0(output, query_tag,"_to_",seed_tag,"_cost_score.pdf")
pdf(file_out, width = 8, height = 6)
print(g7)
dev.off()
