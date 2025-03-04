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
input <- args[2]
output_plots <- output_data <- args[3]

#-----------------------------------------------------------------------------#
# Get all Data sets
#-----------------------------------------------------------------------------#
files <- list.files(input, pattern = "_to_", full.names = TRUE)
files <- grep("arista_regen.rds", files, value = TRUE)

tags <- list.files(outpinputut_data, pattern = "_to_", full.names = FALSE)
tags <- grep("arista_regen.rds", tags, value = TRUE)
tags <- gsub("arista_regen.rds","", tags)
cat("DATA FIND: DONE \n")



files <- files[c(4,6,1,2,3,5)]
tags <- tags[c(4,6,1,2,3,5)]
matched <- lapply(files, readRDS)
names(matched) <- tags
cat("MATCHED DATA LOAD: DONE \n")


ref_files <-  list.files(input, pattern = ".rds", full.names = TRUE)
ref_tags <-  list.files(input, pattern = ".rds", full.names = FALSE)
ref_files <- ref_files[!grepl(pattern = "_to_", x = ref_files) & grepl(pattern = "arista_regen", x = ref_files)]
ref_tags <- ref_tags[!grepl(pattern = "_to_", x = ref_tags) & grepl(pattern = "arista_regen", x = ref_tags)]
ref_tags <- gsub("arista_regen.rds","", ref_tags)

ref_files <- ref_files[c(4,6,1,2,3,5,7)]
ref_tags <- ref_tags[c(4,6,1,2,3,5,7)]
ref_files <- lapply(ref_files, readRDS)
names(ref_files) <- ref_tags
cat("REF DATA LOAD: DONE \n")
#-----------------------------------------------------------------------------#
# Cell type 
#-----------------------------------------------------------------------------#
n_colors <- length(unique(unlist(lapply(ref_files, function(x){
        return(unique(x@territories$Cells))
    }))))
base_colours <- c(
      "#E69F00",
      "#56B4E9",
      "#009E73",
      "#F0E442",
      "#0072B2",
      "#D55E00",
      "#CC79A7",
      "#999999")
pal <- colorRampPalette(base_colours)
colors <- pal(n_colors)#[sample(seq(1,n_colors),size=n_colors)]
labels <- unique(unlist(lapply(ref_files, function(x){
        return(unique(x@territories$Cells))
    })))

cell_plots <- vector("list", length(ref_files))
for (i in seq_along(cell_plots)) {
    pt <- 1.3
    tmp <- ref_files[[i]]@territories
    tmp$Cells <- as.factor(tmp$Cells)
    locs <- match(levels(tmp$Cells), labels)
    sub_cols <- colors[locs]
    tmp$Cells <- factor(tmp$Cells, levels = labels)
    cell_plots[[i]] <- ggplot(tmp, aes(x,y, col = Cells)) +
        geom_point(size = pt, alpha = 1) + 
        scale_color_manual(values = sub_cols) +
        theme_void() +
        theme(legend.title = element_text(size = 20),
            plot.title = element_text(size = 45),
            legend.text = element_text(size = 15),
            legend.position = "none") +
        labs(color = "", title = ref_tags[i]) +
        guides(colour = guide_legend(
            override.aes = list(size = 5)))
}

cell_plots_out <- ggarrange(plotlist = cell_plots,
    nrow = 3,
    ncol = ceiling(length(matched)/2),
    widths = rep(1, ceiling(length(matched)/2)))
file_name <- paste0(output_plots, "reference_ARISTA_regen.pdf")
pdf(file_name, width = 7.5 * length(ref_files)/2, height = 16)
print(cell_plots_out)
dev.off()


map_plots <- vector("list", length(matched))
for (i in seq_along(map_plots)) {
    pt <- 1.3
    tmp <- matched[[i]]@territories
    tmp$Cells <- as.factor(tmp$Cells)
    locs <- match(levels(tmp$Cells), labels)
    sub_cols <- colors[locs]
    tmp$Cells <- factor(tmp$Cells, levels = labels)
    map_plots[[i]] <- ggplot(tmp, aes(x,y, col = Cells)) +
        geom_point(size = pt, alpha = 1) + 
        scale_color_manual(values = sub_cols) +
        theme_void() +
        theme(legend.title = element_text(size = 20),
            plot.title = element_text(size = 25),
            legend.text = element_text(size = 15),
            legend.position = "none") +
        labs(color = "", title = tags[i]) +
        guides(colour = guide_legend(
            override.aes = list(size = 5)))
}

map_plots_out <- ggarrange(plotlist = map_plots,
    nrow = 2,
    ncol = ceiling(length(matched)/2),
    widths = rep(1, ceiling(length(matched)/2)))
file_name <- paste0(output_plots, "mapped_ARISTA_regen.pdf")
pdf(file_name, width = 7.5 * length(matched) /2, height = 16)
print(map_plots_out)
dev.off()

sq <- ceiling(sqrt(n_colors))
x <- rep(seq(1, 30, l = sq), sq)[seq(1,n_colors)]
y <- rep(seq(sq,1), each = sq)[seq(1,n_colors)]
legend <- data.frame(x =x, y = y,labels = labels)
legend$labels <- as.factor(legend$labels)
legend$labels <- factor(legend$labels , levels = labels)
leg_plot <- ggplot(legend, aes(x,y, col = labels)) +
    geom_point(size = 12, alpha = 1) +
    xlim(0,32) +
    ylim(0,32) +
    geom_text(x = legend$x +0.75,hjust=0, y = legend$y, label = levels(legend$labels), size = 5, color = "black") + 
    scale_color_manual(values = colors)+
    theme_void()+
    theme(legend.position = "none",
        plot.margin = margin(2, 2, 2, 2, "cm"))
file_name <- paste0(output_plots, "legend_ARISTA_regen.pdf")
pdf(file_name, width = 15,height = 15)
print(leg_plot)
dev.off()


#-----------------------------------------------------------------------------#
# plot
#-----------------------------------------------------------------------------#
scores <- c("cost", "feature","niche","territory")
map_scores <- vector("list", length(matched) * length(scores))

count <- 1
for (i in seq_along(matched)) {
    tmp <- matched[[i]]
    
    for (j in seq_along(scores)) {
        map_scores[[count]] <- view_mapping_metrics(tmp, scores[j]) + 
            theme_void() +
            theme(plot.title = element_text(size = 40),
                legend.text = element_text(size = 15),
                legend.title = element_text(size = 20)) +
            labs(title = scores[j])

        count <- count + 1
    }
}

map_score_out <- ggarrange(plotlist = map_scores, ncol = length(scores), nrow = length(matched))
file_name <- paste0(output_plots, "scores_ARISTA_regen.pdf")
pdf(file_name, width = 35,height = 35)
print(map_score_out)
dev.off()

file_name <- paste0(output_plots, "scores_ARISTA_regen.png")
png(file_name, width = 3500,height = 3500, type = "cairo")
print(map_score_out)
dev.off()