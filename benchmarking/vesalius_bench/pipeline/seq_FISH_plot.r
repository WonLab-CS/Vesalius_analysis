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
library(ggpubr)
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
args <- commandArgs(TRUE)
idx <- as.numeric(args[1])
if (!dir.exists("/common/martinp4/benchmarking_out/Vesalius/bio_data/seqFISH/")) {
    dir.create("/common/martinp4/benchmarking_out/Vesalius/bio_data/seqFISH/")
}
output <- "/common/martinp4/benchmarking_out/Vesalius/bio_data/seqFISH/"
#-----------------------------------------------------------------------------#
# Load data
#-----------------------------------------------------------------------------#
files <- list.files(output, pattern = ".rds", full.names = TRUE)
seed <- readRDS(grep("seed", files, value = TRUE))
query <- readRDS(grep("query", files, value = TRUE))
matched <- readRDS(grep("matched", files, value = TRUE)[idx])

tags <- list.files(output, pattern = ".rds", full.names = FALSE)
seed_tag <- grep("_seed.rds",tags, value = TRUE)
seed_tag <- gsub("_seed.rds","", seed_tag)

query_tag <- grep("_query.rds",tags, value = TRUE)
query_tag <- gsub("_query.rds","", query_tag)

matched_tag <- grep("_matched.rds", tags, value = TRUE)[idx]
matched_tag <- gsub("_matched.rds","", matched_tag)
matched_tag <- gsub("seqFISH_", "", matched_tag)


#-----------------------------------------------------------------------------#
# plotting
#-----------------------------------------------------------------------------#
combi <- colnames(matched@map[,seq(6, ncol(matched@map))])
mapping_scores <- vector("list", length(combi))

for (i in seq_along(mapping_scores)) {
    mapping_scores[[i]] <- view_mapping_metrics(matched, cex_pt = 1.5, trial = combi[i]) +
        theme_void() +
        labs(title = combi[i])
}

mapping_scores <- ggarrange(plotlist = mapping_scores,
    ncol = ceiling(sqrt(length(combi))),
    nrow = ceiling(sqrt(length(combi))))

file_out <- paste0(output, matched_tag,"_",query_tag,"_to_",seed_tag,".pdf")
pdf(file_out, width = 6 * ceiling(sqrt(length(combi))), height = 6 * ceiling(sqrt(length(combi))))
print(mapping_scores)
dev.off()


ter_1 <- territory_plot(seed, trial = "Cells",randomise = FALSE, cex_pt = 0.7, alpha = 1) +
    #scale_color_manual(values = colors) +
    theme_void() +
    theme(legend.title = element_text(size = 20),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.position = "left") +
    labs(color = "", title = "Reference")

ter_2 <- territory_plot(query, trial = "Cells",randomise = FALSE, cex_pt = 0.7, alpha = 1) +
    #scale_color_manual(values = colors) +
    theme_void() +
    theme(legend.title = element_text(size = 20),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.position = "left") +
    labs(color = "", title = "Query")

ter_3 <- territory_plot(matched, trial = "Cells",randomise = FALSE, cex_pt = 0.7, alpha = 1) +
    #scale_color_manual(values = colors) +
    theme_void() +
    theme(legend.title = element_text(size = 20),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.position = "left") +
    labs(color = "", title = "Mapped")

file_out <- paste0(output, matched_tag,"_",query_tag,"_to_",seed_tag,".pdf")
pdf(file_out, width = 24, height = 8)
print(ter_1 + ter_2 + ter_3)
dev.off()


#-----------------------------------------------------------------------------#
# indiv plots
#-----------------------------------------------------------------------------#


# file_out <- paste0(output, seed_tag,"_cells.pdf")
# pdf(file_out, width = 12, height = 8)
# print(ter_1)
# dev.off()
# file_out <- paste0(output, query_tag,"_cells.pdf")
# pdf(file_out, width = 12, height = 8)
# print(ter_2)
# dev.off()
# file_out <- paste0(output, query_tag,"_to_",seed_tag,"_cells.pdf")
# pdf(file_out, width = 12, height = 8)
# print(ter_3)
# dev.off()


# file_out <- paste0(output, query_tag,"_to_",seed_tag,"_feature_score.pdf")
# pdf(file_out, width = 8, height = 8)
# print(g1)
# dev.off()
# file_out <- paste0(output, query_tag,"_to_",seed_tag,"_niche_score.pdf")
# pdf(file_out, width = 8, height = 8)
# print(g2)
# dev.off()
# file_out <- paste0(output, query_tag,"_to_",seed_tag,"_composition_score.pdf")
# pdf(file_out, width = 8, height = 8)
# print(g3)
# dev.off()
# file_out <- paste0(output, query_tag,"_to_",seed_tag,"_cell_type_score.pdf")
# pdf(file_out, width = 8, height = 8)
# print(g4)
# dev.off()
# file_out <- paste0(output, query_tag,"_to_",seed_tag,"_cost_score.pdf")
# pdf(file_out, width = 8, height = 8)
# print(g5)
# dev.off()
