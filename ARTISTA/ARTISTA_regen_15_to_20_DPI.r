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
# Loading
#-----------------------------------------------------------------------------#


if (!dir.exists("/common/wonklab/Stereo_seq_arista/report/")) {
    dir.create("/common/wonklab/Stereo_seq_arista/report/")
}
output_plots <- "/common/wonklab/Stereo_seq_arista/report/"
output_data <- "/common/wonklab/Stereo_seq_arista/report/"

idx <-""
#-----------------------------------------------------------------------------#
# Get Single Data sets
#-----------------------------------------------------------------------------#
files <- list.files(output_data, pattern = "_to_", full.names = TRUE)
files <-files[grepl("_arista_regen.rds", files) &
    !grepl("integrated", files) &
    !grepl("metric_clusters", files)]
files <- grep(pattern = "15DPI_1_to_20DPI_1",files, value = TRUE)
tags <- list.files(output_data, pattern = "_to_", full.names = FALSE)
tags <- grep(".rds", tags, value = TRUE)
tags <- grep("_arista_regen.rds", tags, value = TRUE)
tags <-tags[grepl(".rds", tags) &
    !grepl("integrated", tags) &
    !grepl("metric_clusters", tags)]
tags <- grep(pattern = "15DPI_1_to_20DPI_1",tags, value = TRUE)
tags <- gsub("_arista_regen.rds","", tags)



matched <- readRDS(files)
matched_cells <- matched@territories
stage_seed <- sapply(strsplit(tags, "_to_"),"[[", 2)
stage_query <- sapply(strsplit(tags, "_to_"),"[[", 1)


ref_files <-  list.files(output_data, pattern = ".rds", full.names = TRUE)
ref_files <- ref_files[grepl(pattern = stage_seed, x = ref_files) & !grepl(pattern = "_to_", x = ref_files)]
seed <- readRDS(ref_files)
seed_cells <- seed@territories

query_files <-  list.files(output_data, pattern = ".rds", full.names = TRUE)
query_files <- query_files[grepl(pattern = stage_query, x = query_files) & !grepl(pattern = "_to_", x = query_files)]
query <- readRDS(query_files)
query_cells <- query@territories


file_name <- paste0(output_data,tags,"_DEGs_inter_sample.csv")
# if (!file.exists(file_name)) {
#     degs <- data.frame("genes","p_value","p_value_adj","seed_pct","query_pct","fold_change","effect_size","seed","query","Cells")
#     write.table(degs,
#         file = file_name,
#         quote = FALSE,
#         col.names = TRUE,
#         row.names = FALSE,
#         sep = ",")
# }
#-----------------------------------------------------------------------------#
# Integrate
#-----------------------------------------------------------------------------#
inter <- integrate_assays(matched, seed, infer = TRUE,
    labels_mapped = c("Cells","Territory"),
    labels_reference = c("Cells","Territory"))
inter <- generate_tiles(inter, filter_threshold = 1, filter_grid =1) %>%
    equalize_image(embedding = "integrated", dimensions = 1:20, sleft = 2.5, sright = 2.5) %>%
    smooth_image(embedding = "integrated",dimensions = 1:20, method = c("iso", "box"), box = 10, sigma = 1, iter = 15) %>%
    segment_image(dimensions = 1:20, method = "kmeans", col_resolution = 10) %>%
    isolate_territories(capture_radius = 0.01, min_spatial_index = 20)
inter_cells <- inter@territories
file_name <- paste0(output_data,tags,"_integrated.rds")
saveRDS(inter, file = file_name)
cat("Integration: done\n")
#-----------------------------------------------------------------------------#
# Cells 
#-----------------------------------------------------------------------------#
labels <- sort(union(unique(seed_cells$Cells), unique(query_cells$Cells)))
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

seed_cells$Cells <- as.factor(seed_cells$Cells)
locs <- match(levels(seed_cells$Cells),labels)
seed_colors <- colors[locs]
#seed$Cells <- factor(seed$Cells, levels = labels)

query_cells$Cells <- as.factor(query_cells$Cells)
locs <- match( levels(query_cells$Cells),labels)
query_colors <- colors[locs]
#query$Cells <- factor(query$Cells, levels = labels)

matched_cells$Cells <- as.factor(matched_cells$Cells)
locs <- match( levels(matched_cells$Cells),labels)
matched_colors <- colors[locs]


#-----------------------------------------------------------------------------#
# Plot cells
#-----------------------------------------------------------------------------#
pt <- 1
seed_plot <- ggplot(seed_cells, aes(x,y, col = Cells)) +
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

pt <- 1
query_plot <- ggplot(query_cells, aes(x,y, col = Cells)) +
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

pt <- 1
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

pt <- 1
inter_plot <- ggplot(inter_cells, aes(x,y, col = Cells_Cells)) +
        geom_point(size = pt, alpha = 1) + 
        scale_color_manual(values = colors) +
        facet_wrap(~sample)+
        theme_void() +
        theme(legend.title = element_text(size = 20),
            plot.title = element_text(size = 20),
            legend.text = element_text(size = 15),
            legend.position = "left") +
        labs(color = "", title = paste(stage_query,"to",stage_seed)) +
        guides(colour = guide_legend(
            override.aes = list(size = 5)))

file_name <- paste0(output_plots, stage_seed,"_",idx,"_reference_ARISTA_regen.pdf")
pdf(file_name, width = 10, height = 8)
print(seed_plot)
dev.off()

file_name <- paste0(output_plots, stage_query,"_",idx,"_reference_ARISTA_regen.pdf")
pdf(file_name, width = 10, height = 8)
print(query_plot)
dev.off()

file_name <- paste0(output_plots, stage_query,"_to_",stage_seed,"_reference_ARISTA_regen.pdf")
pdf(file_name, width = 10, height = 8)
print(matched_plot)
dev.off()

file_name <- paste0(output_plots, stage_query,"_to_",stage_seed,"_reference_ARISTA_regen_inter.pdf")
pdf(file_name, width = 18, height = 8)
print(inter_plot)
dev.off()

#-----------------------------------------------------------------------------#
# Territories
#-----------------------------------------------------------------------------#
labels <- sort(union(unique(seed_cells$Territory), unique(query_cells$Territory)))
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
colors_match <- ter_pal(length(unique(matched_cells$Territory))) 
colors_inter <- ter_pal(length(unique(inter_cells$Territory)))

seed_cells$Territory <- as.factor(seed_cells$Territory)
query_cells$Territory <- as.factor(query_cells$Territory)
matched_cells$Territory <- as.factor(matched_cells$Territory)

#-----------------------------------------------------------------------------#
# Plot territories
#-----------------------------------------------------------------------------#
pt <- 1
seed_plot <- ggplot() +
        geom_point(data = subset(seed_cells,Territory != "isolated"), aes(x,y, col = Territory),size = pt, alpha = 1) + 
        scale_color_manual(values = colors) +
        theme_void() +
        theme(legend.title = element_text(size = 20),
            plot.title = element_text(size = 20),
            legend.text = element_text(size = 15),
            legend.position = "left") +
        labs(color = "", title = stage_seed) +
        guides(colour = guide_legend(
            override.aes = list(size = 5)))

pt <- 1
query_plot <- ggplot() +
        geom_point(data = subset(query_cells,Territory != "isolated"), aes(x,y, col = Territory),size = pt, alpha = 1) + 
        scale_color_manual(values = colors) +
        theme_void() +
        theme(legend.title = element_text(size = 20),
            plot.title = element_text(size = 20),
            legend.text = element_text(size = 15),
            legend.position = "left") +
        labs(color = "", title = stage_query) +
        guides(colour = guide_legend(
            override.aes = list(size = 5)))

pt <- 1
matched_plot <- ggplot() +
        geom_point(data = subset(matched_cells,Territory != "isolated"), aes(x,y, col = Territory),size = pt, alpha = 1) + 
        scale_color_manual(values = colors_match) +
        theme_void() +
        theme(legend.title = element_text(size = 20),
            plot.title = element_text(size = 20),
            legend.text = element_text(size = 15),
            legend.position = "left") +
        labs(color = "", title = paste(stage_query,"to",stage_seed)) +
        guides(colour = guide_legend(
            override.aes = list(size = 5)))

pt <- 1
inter_plot <- ggplot() +
        geom_point(data = subset(inter_cells,Territory != "isolated"), aes(x,y, col = Territory),size = pt, alpha = 1) + 
        scale_color_manual(values = sample(colors_inter, length(colors_inter))) +
        facet_wrap(~sample) +
        theme_void() +
        theme(legend.title = element_text(size = 20),
            plot.title = element_text(size = 20),
            legend.text = element_text(size = 15),
            legend.position = "left") +
        labs(color = "", title = paste(stage_query,"to",stage_seed)) +
        guides(colour = guide_legend(
            override.aes = list(size = 5)))

file_name <- paste0(output_plots, stage_seed,"_",idx,"_reference_ARISTA__terregen.pdf")
pdf(file_name, width = 10, height = 8)
print(seed_plot)
dev.off()

file_name <- paste0(output_plots, stage_query,"_",idx,"_reference_ARISTA_terregen.pdf")
pdf(file_name, width = 10, height = 8)
print(query_plot)
dev.off()

file_name <- paste0(output_plots, stage_query,"_to_",stage_seed,"_reference_ARISTA_terregen.pdf")
pdf(file_name, width = 10, height = 8)
print(matched_plot)
dev.off()

file_name <- paste0(output_plots, stage_query,"_to_",stage_seed,"_reference_ARISTA_terregen_inter.pdf")
pdf(file_name, width = 18, height = 8)
print(inter_plot)
dev.off()

#-----------------------------------------------------------------------------#
# split territory
#-----------------------------------------------------------------------------#
pt <- 1
seed_plot <- ggplot() +
        geom_point(data = subset(seed_cells,Territory != "isolated"), aes(x,y, col = Territory),size = pt, alpha = 1) + 
        facet_wrap(~Territory)+
        scale_color_manual(values = colors) +
        theme_void() +
        theme(legend.title = element_text(size = 20),
            plot.title = element_text(size = 20),
            legend.text = element_text(size = 15),
            legend.position = "left") +
        labs(color = "", title = stage_seed) +
        guides(colour = guide_legend(
            override.aes = list(size = 5)))

pt <- 1
query_plot <- ggplot() +
        geom_point(data = subset(query_cells,Territory != "isolated"), aes(x,y, col = Territory),size = pt, alpha = 1) + 
         facet_wrap(~Territory)+        
        scale_color_manual(values = colors) +
        theme_void() +
        theme(legend.title = element_text(size = 20),
            plot.title = element_text(size = 20),
            legend.text = element_text(size = 15),
            legend.position = "left") +
        labs(color = "", title = stage_query) +
        guides(colour = guide_legend(
            override.aes = list(size = 5)))

pt <- 1
matched_plot <- ggplot() +
        geom_point(data = subset(matched_cells,Territory != "isolated"), aes(x,y, col = Territory),size = pt, alpha = 1) + 
         facet_wrap(~Territory)+
        scale_color_manual(values = colors_match) +
        theme_void() +
        theme(legend.title = element_text(size = 20),
            plot.title = element_text(size = 20),
            legend.text = element_text(size = 15),
            legend.position = "left") +
        labs(color = "", title = paste(stage_query,"to",stage_seed)) +
        guides(colour = guide_legend(
            override.aes = list(size = 5)))
pt <- 1
inter_plot <- ggplot() +
        geom_point(data = subset(inter_cells,Territory != "isolated"), aes(x,y, col = Territory),size = pt, alpha = 1) + 
         facet_wrap(~Territory)+
        scale_color_manual(values = colors_inter) +
        theme_void() +
        theme(legend.title = element_text(size = 20),
            plot.title = element_text(size = 20),
            legend.text = element_text(size = 15),
            legend.position = "left") +
        labs(color = "", title = paste(stage_query,"to",stage_seed)) +
        guides(colour = guide_legend(
            override.aes = list(size = 5)))

file_name <- paste0(output_plots, stage_seed,"_",idx,"_reference_ARISTA__terregen_split.pdf")
pdf(file_name, width = 14, height = 14)
print(seed_plot)
dev.off()

file_name <- paste0(output_plots, stage_query,"_",idx,"_reference_ARISTA_terregen_split.pdf")
pdf(file_name, width = 14, height = 14)
print(query_plot)
dev.off()

file_name <- paste0(output_plots, stage_query,"_to_",stage_seed,"_reference_ARISTA_terregen_split.pdf")
pdf(file_name, width = 14, height = 14)
print(matched_plot)
dev.off()

file_name <- paste0(output_plots, stage_query,"_to_",stage_seed,"_reference_ARISTA_terregen_split_inter.pdf")
pdf(file_name, width = 14, height = 14)
print(inter_plot)
dev.off()




#-----------------------------------------------------------------------------#
# compare mapped cells with cell types
#-----------------------------------------------------------------------------#
inter <- identify_markers(inter,
    trial = "Territory",
    norm_method = "scaled",
    seed = c(8),
    query = c(8),
    sample = TRUE)

degs <- get_markers(inter) %>% arrange(desc(fold_change))

plot_degs <- degs
file_out <- paste0(output_plots,tags,"_integrated_diff_ter8_.pdf")
pdf(file_out, width = 24, height = 8)
for (i in seq_len(nrow(plot_degs))) {
    g <- view_gene_expression(matched, norm_method = "log_norm", genes = plot_degs$genes[i], cex_pt = 0.35) +
        theme_void() +
        labs(title = paste0(plot_degs$genes[i]))
    g1 <- view_gene_expression(seed, norm_method = "raw", genes = plot_degs$genes[i], cex_pt = 0.35) +
        theme_void() +
        labs(title = paste0("Reference ", plot_degs$genes[i]))
    g2 <- view_gene_expression(query, norm_method = "last", genes = plot_degs$genes[i], cex_pt = 0.35) +
        theme_void() +
        labs(title = paste0("Query ", plot_degs$genes[i]))
    g_all <- g1 + g + g2
    
    
    print(g_all)
    
}
dev.off()

file_name <- paste0(output_data,tags,"_DEGs_inter_sample.csv")
write.csv(degs, file = file_name,
        quote = FALSE,
        col.names = TRUE,
        row.names = FALSE)