#-----------------------------------------------------------------------------#
# Loading Libraries and setting seed
#-----------------------------------------------------------------------------#
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
library(Morpho)
library(ggpubr)
library(lpSolve, lib = "/common/martinp4/R")
library(TreeDist, lib.loc = "/common/martinp4/R")
library(mclust, lib.loc = "/common/martinp4/R")
library(anndata, lib.loc = "/common/martinp4/R")
library(pwr, lib.loc = "/common/martinp4/R")
library(gsignal, lib.loc = "/common/martinp4/R")
library(kohonen, lib.loc = "/common/martinp4/R")
library(registry, lib.loc = "/common/martinp4/R")
library(rngtools, lib.loc = "/common/martinp4/R")
library(NMF, lib.loc = "/common/martinp4/R")
library(RcppHungarian, lib.loc = "/common/martinp4/R")
library(spatstat.utils, lib.loc = "/common/martinp4/R")
library(geometry, lib.loc = "/common/martinp4/R")
library(vesalius, lib.loc = "/common/martinp4/R")
library(RColorBrewer)

set.seed(1547)
#plan(multicore, workers = 5)
max_size <- 100000 * 1024^2
options(future.globals.maxSize = max_size)

#-----------------------------------------------------------------------------#
# directories
#-----------------------------------------------------------------------------#
args <- commandArgs(TRUE)
idx <- as.numeric(args[1])

if (!dir.exists("/common/wonklab/Stereo_seq/report/")) {
    dir.create("/common/wonklab/Stereo_seq/report/")
}
output_plots <- "/common/wonklab/Stereo_seq/report/"
output_data <- "/common/wonklab/Stereo_seq/report/"

#-----------------------------------------------------------------------------#
# Get Single Data sets
#-----------------------------------------------------------------------------#

files <- list.files(output_data, pattern = "_to_", full.names = TRUE)
files <- grep(".rds", files, value = TRUE)[idx]
tags <- list.files(output_data, pattern = "_to_", full.names = FALSE)
tags <- grep(".rds", tags, value = TRUE)[idx]
tags <- gsub("_stereo.rds","", tags)
matched <- readRDS(files)
matched_cells <- matched@territories
stage_seed <- sapply(strsplit(tags, "_to_"),"[[", 2)
stage_query <- sapply(strsplit(tags, "_to_"),"[[", 1)


ref_files <-  list.files(output_data, pattern = ".rds", full.names = TRUE)
ref_files <- ref_files[grepl(pattern = stage_seed, x = ref_files) & !grepl(pattern = "_to_", x = ref_files)]
seed <- readRDS(ref_files)
seed <- seed@territories

query_files <-  list.files(output_data, pattern = ".rds", full.names = TRUE)
query_files <- query_files[grepl(pattern = stage_query, x = query_files) & !grepl(pattern = "_to_", x = query_files)]
query <- readRDS(query_files)
query <- query@territories


#-----------------------------------------------------------------------------#
# Get full cell type labels
#-----------------------------------------------------------------------------#
# label_file <- paste0(output_data,"labels.csv")
# labels <- read.csv(label_file, header = TRUE)
# colors <- labels$colors
# labels <- labels$labels
labels <- sort(union(unique(seed$Cells), unique(query$Cells)))
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
#colors <- sample(ter_pal(ter_col), ter_col)
colors <- ter_pal(ter_col)

seed$Cells <- as.factor(seed$Cells)
locs <- match(levels(seed$Cells),labels)
seed_colors <- colors[locs]
#seed$Cells <- factor(seed$Cells, levels = labels)

query$Cells <- as.factor(query$Cells)
locs <- match( levels(query$Cells),labels)
query_colors <- colors[locs]
#query$Cells <- factor(query$Cells, levels = labels)

matched_cells$Cells <- as.factor(matched_cells$Cells)
locs <- match( levels(matched_cells$Cells),labels)
matched_colors <- colors[locs]
#matched_cells$Cells <- factor(matched_cells$Cells, levels = labels)

#-----------------------------------------------------------------------------#
# Plot cells
#-----------------------------------------------------------------------------#
pt <- ifelse(stage_seed == "E.95", 2, 0.5)
seed_plot <- ggplot(seed, aes(x,y, col = Cells)) +
        geom_point(size = pt, alpha = 1) + 
        scale_color_manual(values = seed_colors) +
        theme_void() +
        theme(legend.title = element_text(size = 20),
            plot.title = element_text(size = 20),
            legend.text = element_text(size = 15),
            legend.position = "left") +
        labs(color = "", title = stage_seed) +
        guides(colour = guide_legend(
            override.aes = list(size = 5)))

pt <- ifelse(stage_query == "E.95", 2, 0.5)
query_plot <- ggplot(query, aes(x,y, col = Cells)) +
        geom_point(size = pt, alpha = 1) + 
        scale_color_manual(values = query_colors) +
        theme_void() +
        theme(legend.title = element_text(size = 20),
            plot.title = element_text(size = 20),
            legend.text = element_text(size = 15),
            legend.position = "left") +
        labs(color = "", title = stage_query) +
        guides(colour = guide_legend(
            override.aes = list(size = 5)))

pt <- ifelse(stage_query == "E.95", 2, 0.5)
matched_plot <- ggplot(matched_cells, aes(x,y, col = Cells)) +
        geom_point(size = pt, alpha = 1) + 
        scale_color_manual(values = matched_colors) +
        theme_void() +
        theme(legend.title = element_text(size = 20),
            plot.title = element_text(size = 20),
            legend.text = element_text(size = 15),
            legend.position = "left") +
        labs(color = "", title = paste(stage_query,"to",stage_seed)) +
        guides(colour = guide_legend(
            override.aes = list(size = 5)))

file_name <- paste0(output_plots, stage_query,"_to_",stage_seed,"_reference_MOSTA.pdf")
pdf(file_name, width = 10, height = 8)
print(seed_plot)
dev.off()

file_name <- paste0(output_plots, stage_query,"_to_",stage_seed,"_query_MOSTA.pdf")
pdf(file_name, width = 10, height = 8)
print(query_plot)
dev.off()

file_name <- paste0(output_plots, stage_query,"_to_",stage_seed,"_mapped_MOSTA.pdf")
pdf(file_name, width = 10, height = 8)
print(matched_plot)
dev.off()



#-----------------------------------------------------------------------------#
# plot
#-----------------------------------------------------------------------------#

scores <- c("cost", "feature","niche","territory")
count <- 1
for (j in seq_along(scores)) {
    map_scores <- view_mapping_metrics(matched, scores[j], cex_pt = 0.5)+ 
        theme_void() +
        theme(plot.title = element_text(size = 20),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 20)) +
        labs(title = scores[j])
    file_name <- paste0(output_plots, tags,"_",scores[j],"_scores_MOSTA.pdf")
    pdf(file_name, width = 10,height = 8)
    print(map_scores)
    dev.off()
    file_name <- paste0(output_plots, tags,"_",scores[j],"_scores_MOSTA.png")
    png(file_name, width = 1000,height = 800, type = "cairo")
    print(map_scores)
    dev.off()
}

#-----------------------------------------------------------------------------#
# embeds 
#-----------------------------------------------------------------------------#



