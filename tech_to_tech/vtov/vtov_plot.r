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
args <- commandArgs(TRUE)
output <- args[1]
#-----------------------------------------------------------------------------#
# Load data
#-----------------------------------------------------------------------------#
files <- list.files(output, pattern = ".rds", full.names = TRUE)
seed <- readRDS(grep("seed", files, value = TRUE))
query <- readRDS(grep("query", files, value = TRUE))
matched <- readRDS(grep("matched", files, value = TRUE))


tags <- list.files(output, pattern = ".rds", full.names = FALSE)
seed_tag <- gsub("_seed.rds","", grep("seed",tags, value = TRUE))
query_tag <- gsub("_query.rds","", grep("query",tags, value = TRUE))


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
g2 <- image_plot(query, embedding = "last") +
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





g4 <- view_mapping_metrics(matched, cex_pt = 1.5, trial = "niche") +
    theme_void()+
    labs(title = "Niche")
g5 <- view_mapping_metrics(matched, cex_pt = 1.5, trial = "territory") + 
    theme_void()+
    labs(title = "Territory")
g6 <- view_mapping_metrics(matched, cex_pt = 1.5, trial = "cost") + 
    theme_void()+
    labs(title = "Cost")

file_out <- paste0(output, "mapping_scores_",query_tag,"_to_",seed_tag,".pdf")
pdf(file_out, width = 24, height = 8)
all <- g4 + g5 + g6
print(all)
dev.off()


ter_1 <- territory_plot(seed, trial = "Territory", cex_pt = 3.5, alpha = 1,use_image = FALSE) +
    #scale_color_manual(values = colors) +
    theme_void() +
    theme(legend.title = element_text(size = 20),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.position = "none") +
    labs(color = "", title = paste0("Reference - ", seed_tag))

ter_2 <- territory_plot(query,trial = "Territory",  cex_pt = 0.75, alpha=1, use_image = FALSE) +
    #scale_color_manual(values = colors) +
    theme_void() +
    theme(legend.title = element_text(size = 20),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.position = "none") +
    labs(color = "", title = paste0("Query - ", query_tag))

ter_3 <- territory_plot(matched, trial = "Territory",  cex_pt = 0.75, alpha=1, use_image = FALSE) +
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

cell_dense <- matched@map$to
cell_dense <- sapply(cell_dense, function(dat){
    prefix <- paste0(strsplit(x = dat, split = "-")[[1]][1:2], collapse = "-")
    return(prefix)
})
cell_dense <- as.data.frame(table(cell_dense))
colnames(cell_dense) <- c("barcodes","Freq")
seed_coord <- seed@territories[,c("barcodes","x","y")]
cell_dense <- right_join(cell_dense, seed_coord, by = "barcodes")
cell_dense_plot <- ggplot(cell_dense, aes(x,y,col = Freq)) +
    geom_point(size = 3.5)+
    scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral"))) +
    theme_void()+
    labs(title = "Mapping Density",color = "Density")


file_out <- paste0(output, seed_tag,"_density.pdf")
pdf(file_out, width = 8, height = 8)
print(cell_dense_plot)
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
pdf(file_out, width = 12, height = 12)
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



file_out <- paste0(output, query_tag,"_to_",seed_tag,"_niche_score.pdf")
pdf(file_out, width = 8, height = 6)
print(g4)
dev.off()
file_out <- paste0(output, query_tag,"_to_",seed_tag,"_territory_score.pdf")
pdf(file_out, width = 8, height = 6)
print(g5)
dev.off()
file_out <- paste0(output, query_tag,"_to_",seed_tag,"_cost_score.pdf")
pdf(file_out, width = 8, height = 6)
print(g6)
dev.off()

#-----------------------------------------------------------------------------#
# Jaccard function
# a bit hardcoded for now not worth the effort or making it work for all cases
#-----------------------------------------------------------------------------#
jaccard_cells <- function(viz, hd, matched) {
    to <- paste0(sapply(strsplit(matched@map$to,"-"),"[[",1),"-1")
    from <- matched@map$from
    viz <- viz[unique(to)]
    names(hd) <- to
    hd <- split(hd, names(hd))
    hd <- lapply(hd, unlist)
    hd <- hd[unique(to)]
    jacques <- mapply(function(h,v){
        jack <- length(intersect(h, names(v))) / length(union(h, names(v)))
        return(jack)
    },hd,viz)
    jacques <- jacques[to]
    ter <- matched@territories
    ter$Jaccard <- jacques
    ter$prop <- matched@map$prop
    ter$signif <- matched@map$prop > 0.05
    return(ter)
}

#-----------------------------------------------------------------------------#
# load and plot
#-----------------------------------------------------------------------------#

prop_viz <- readRDS(paste0(output, "visium_mouse_brain_prop.rds"))
prop_hd <- readRDS(paste0(output, "visiumHD_mouse_brain_cells.rds" ))

jack <- jaccard_cells(prop_viz, prop_hd, matched)
jack$Sample <- "Visium to VisiumHD - Jaccard I."


h <- ggplot(jack,aes(x = Jaccard)) +
    geom_histogram(
        breaks=seq(0,max(jack$Jaccard)+0.1,0.1),
        fill="#082233ff",
        colour="#082233ff") +
    geom_vline(xintercept = quantile(jack$Jaccard, 0.25),color = "#E69F00",linetype ="solid", linewidth = 1.2)+
    geom_vline(xintercept = quantile(jack$Jaccard, 0.5),color = "#56B4E9",linetype ="solid",linewidth = 1.2)+
    geom_vline(xintercept = quantile(jack$Jaccard, 0.75),color = "#009E73",linetype ="solid", linewidth = 1.2)+
    theme_bw() +
    theme(strip.background =element_rect(fill="#082233ff"),
        strip.text = element_text(colour = 'white', size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size=15),
        legend.text = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1,vjust =1)) +
    facet_wrap(~Sample)
file_out <- paste0(output,"Visium_to_VisiumHD_Jaccard.pdf")
pdf(file_out, width = 5, height = 4)
print(h)
dev.off()

p <- ggplot(jack, aes(x = x,  y = y, col = Jaccard)) +
    geom_point() +
    scale_color_viridis_c(option = "mako") +
    theme_void() +
    labs(title = "Jaccard Index",color = "Score")

file_out <- paste0(output,"Visium_to_VisiumHD_Jaccard_mapped.pdf")
pdf(file_out, width = 8, height = 6)
print(p)
dev.off()

p <- ggplot(jack, aes(x = x,  y = y, col = prop)) +
    geom_point() +
    scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral"))) +
    theme_void() +
    labs(title = "Cell Type Proportions",color = "p-val")

file_out <- paste0(output,"Visium_to_VisiumHD_Prop_mapped.pdf")
pdf(file_out, width = 8, height = 6)
print(p)
dev.off()



p <- ggplot(jack, aes(x = x,  y = y, col = signif)) +
    geom_point(size = 0.85) +
    scale_color_manual(values = c("#E69F00", "#56B4E9")) +
    theme_void() +
    labs(title = "Cell Type Proportions Discretized",color = "pval < 0.05") +
    guides(colour = guide_legend(
            override.aes = list(size = 5)))

file_out <- paste0(output,"Visium_to_VisiumHD_PropDisc_mapped.pdf")
pdf(file_out, width = 8, height = 6)
print(p)
dev.off()
