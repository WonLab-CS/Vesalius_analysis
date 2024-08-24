#-----------------------------------------------------------------------------#
# Loading Libraries and setting seed
#-----------------------------------------------------------------------------#
library(dplyr, lib.loc = "/common/martinp4/R")
library(ggplot2)
library(patchwork)
library(ggpubr)
library(RColorBrewer)


library(vesalius, lib.loc = "/common/martinp4/R")
set.seed(1547)

max_size <- 10000 * 1024^2
options(future.globals.maxSize = max_size)
#-----------------------------------------------------------------------------#
# Output set up 
#-----------------------------------------------------------------------------#
input <- "/common/martinp4/benchmarking_out/"

if(!dir.exists("/common/martinp4/benchmarking_out/report/")){
    dir.create("/common/martinp4/benchmarking_out/report/")
}
output <- "/common/martinp4/benchmarking_out/report/"
#output <- "/Users/martinp4/Documents/Cedars/Vesalius/Scenes/benchmarking/report/"

#-----------------------------------------------------------------------------#
# load scores
#-----------------------------------------------------------------------------#
file_name <- paste0(output,"benchmarking_scores.csv")
scores <- read.csv(file_name)
scores$Method <- gsub("synthetic/","", scores$Method)
scores$X <- NULL
scores$Regime <- factor(scores$Regime, levels = c("circle","layered","dropped"))
#-----------------------------------------------------------------------------#
# Plot - scores
#-----------------------------------------------------------------------------#
ari <- ggplot(scores, aes(x = Method, y = ARI, fill = Method)) +
    geom_boxplot() +
    #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1) +
    scale_fill_manual(values =c("#08324e","#1a6ea7ff","#146269ff","#51b090","#de9c71","#b75b1cff","#8e0703ff")) +
    theme_bw() +
    theme(strip.background =element_rect(fill="#082233ff"),
        strip.text = element_text(colour = 'white', size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size=15),
        legend.text = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1,vjust =1)) +
    facet_wrap(~Regime)
file_out <- paste0(output, "benchmarking_ARI.pdf")
pdf(file_out, width = 8, height = 4)
print(ari)
dev.off()


vi <- ggplot(scores, aes(x = Method, y = VI, fill = Method)) +
    geom_boxplot() +
    #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1) +
    scale_fill_manual(values =c("#08324e","#1a6ea7ff","#146269ff","#51b090","#de9c71","#b75b1cff","#8e0703ff")) +
    theme_bw() +
    theme(strip.background =element_rect(fill="#082233ff"),
        strip.text = element_text(colour = 'white', size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size=15),
        legend.text = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1,vjust =1)) +
    facet_wrap(~Regime)
file_out <- paste0(output, "benchmarking_VI.pdf")
pdf(file_out, width = 8, height = 4)
print(vi)
dev.off()



#-----------------------------------------------------------------------------#
# Get best matching
#-----------------------------------------------------------------------------#

best_circle <- scores[scores$Method == "Vesalius" & scores$Regime == "circle",]
best_circle <- best_circle[best_circle$ARI == max(best_circle$ARI), ]


best_layer <- scores[scores$Method == "Vesalius" & scores$Regime == "layered",]
best_layer <- best_layer[best_layer$ARI == sort(best_layer$ARI,decreasing = TRUE)[3], ]

best_dropped <- scores[scores$Method == "Vesalius" & scores$Regime == "dropped",]
best_dropped <- best_drop[best_drop$ARI == sort(best_drop$ARI,decreasing = TRUE)[1], ]


#-----------------------------------------------------------------------------#
# Get best matching and ref files
#-----------------------------------------------------------------------------#
matched_path <- "/common/martinp4/benchmarking_out/"
ref_path <- "/common/wonklab/synthetic_spatial/"

files <- list.files(list.dirs(matched_path, recursive = FALSE),
    pattern = ".csv",
    recursive = TRUE,
    full.names = TRUE)
files <- sort(grep("aligned",files, value = TRUE))
ref_path <- list.files(ref_path, pattern = "spatial_coordinates", full.names = TRUE)

### Circle
ref_circle <- ref_path[grepl("circle",ref_path) &
    grepl(paste0("sample_",best_circle$ref_sample), ref_path)]
ref_circle <- read.csv(ref_circle)
ref_circle$Method <- "Reference"
ref_circle <- ref_circle[, c("barcodes","x","y","cell_labels","Method")]
query_circle <- ref_path[grepl("circle",ref_path) &
    grepl(paste0("sample_",best_circle$query_sample), ref_path)]
query_circle <- read.csv(query_circle)
query_circle$Method <- "Query"
query_circle <- query_circle[, c("barcodes","x","y","cell_labels","Method")]

matched <- grep(paste0("circle_sample_",best_circle$ref_sample,"_circle_sample_",best_circle$query_sample),files, value = TRUE)
methods <- gsub("/common/martinp4/benchmarking_out//","", matched)
methods <- sapply(strsplit(split = "/report/",x = methods),"[[",1)
matched <- mapply(function(f,tag){
    f <- read.csv(f)
    f <- f[,-1]
    f$Method <- tag
    f$barcodes <- gsub("q_","", f$barcodes)
    f$cell_labels <- gsub("celltype_","",f$cell_labels)
    f$cluster <- NULL
    rownames(f) <- NULL
    f <- f[,1:5]
    colnames(f) <- c("barcodes","x","y","cell_labels","Method")
    return(f)
}, matched, methods, SIMPLIFY = FALSE)
matched <- do.call("rbind", matched)
rownames(matched) <- NULL

all_data_circle <- rbind(ref_circle, query_circle, matched)
all_data_circle$Method[all_data_circle$Method == "precast"] <- "PRECAST"

### Layered
ref_layered <- ref_path[grepl("layered",ref_path) &
    grepl(paste0("sample_",best_layer$ref_sample), ref_path)]
ref_layered <- read.csv(ref_layered)
ref_layered$Method <- "Reference"
ref_layered <- ref_layered[, c("barcodes","x","y","cell_labels","Method")]
query_layered <- ref_path[grepl("layered",ref_path) &
    grepl(paste0("sample_",best_layer$query_sample), ref_path)]
query_layered <- read.csv(query_layered)
query_layered$Method <- "Query"
query_layered <- query_layered[, c("barcodes","x","y","cell_labels","Method")]

matched <- grep(paste0("layered_sample_",best_layer$ref_sample,"_layered_sample_",best_layer$query_sample),files, value = TRUE)
methods <- gsub("/common/martinp4/benchmarking_out//","", matched)
methods <- sapply(strsplit(split = "/report/",x = methods),"[[",1)
matched <- mapply(function(f,tag){
    f <- read.csv(f)
    f <- f[,-1]
    f$Method <- tag
    f$barcodes <- gsub("q_","", f$barcodes)
    f$cell_labels <- gsub("celltype_","",f$cell_labels)
    f$cluster <- NULL
    rownames(f) <- NULL
    f <- f[,1:5]
    colnames(f) <- c("barcodes","x","y","cell_labels","Method")
    return(f)
}, matched, methods, SIMPLIFY = FALSE)
matched <- do.call("rbind", matched)
rownames(matched) <- NULL

all_data_layered <- rbind(ref_layered, query_layered, matched)
all_data_layered$Method[all_data_layered$Method == "precast"] <- "PRECAST"



### dropped
ref_dropped <- ref_path[grepl("dropped",ref_path) &
    grepl(paste0("sample_",best_dropped$ref_sample), ref_path)]
ref_dropped <- read.csv(ref_dropped)
ref_dropped$Method <- "Reference"
ref_dropped <- ref_dropped[, c("barcodes","x","y","cell_labels","Method")]
query_dropped <- ref_path[grepl("dropped",ref_path) &
    grepl(paste0("sample_",best_dropped$query_sample), ref_path)]
query_dropped <- read.csv(query_dropped)
query_dropped$Method <- "Query"
query_dropped <- query_dropped[, c("barcodes","x","y","cell_labels","Method")]

matched <- grep(paste0("dropped_sample_",best_layer$ref_sample,"_dropped_sample_",best_layer$query_sample),files, value = TRUE)
methods <- gsub("/common/martinp4/benchmarking_out//","", matched)
methods <- sapply(strsplit(split = "/report/",x = methods),"[[",1)
matched <- mapply(function(f,tag){
    f <- read.csv(f)
    f <- f[,-1]
    f$Method <- tag
    f$barcodes <- gsub("q_","", f$barcodes)
    f$cell_labels <- gsub("celltype_","",f$cell_labels)
    f$cluster <- NULL
    rownames(f) <- NULL
    f <- f[,1:5]
    colnames(f) <- c("barcodes","x","y","cell_labels","Method")
    return(f)
}, matched, methods, SIMPLIFY = FALSE)
matched <- do.call("rbind", matched)
rownames(matched) <- NULL

all_data_dropped<- rbind(ref_dropped, query_dropped, matched)
all_data_dropped$Method[all_data_dropped$Method == "precast"] <- "PRECAST"

#-----------------------------------------------------------------------------#
# Plot - best matched
#-----------------------------------------------------------------------------#
all_data_circle$cell_labels <- as.factor(all_data_circle$cell_labels)
all_data_circle$Method <- as.factor(all_data_circle$Method)
all_data_circle$Method <- factor(all_data_circle$Method,
    levels = c("Reference","Query","Vesalius","CytoSpace","Tangram","SLAT","PASTE","GPSA","PRECAST"))
cols <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))
cols <- cols(length(levels(all_data_circle$cell_labels)))
g1 <- ggplot(all_data_circle, aes(x = x, y = y, col = cell_labels)) +
    geom_point(size = 0.8) +
    theme_bw() +
    theme(strip.background =element_rect(fill="#082233ff"),
        strip.text = element_text(colour = 'white', size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size=15),
        legend.text = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
    scale_color_manual(values = cols) +
    facet_wrap(~Method) +
    guides(colour = guide_legend(
        override.aes = list(size =  5)))+
    labs(col = "Cell Labels")
file_name <- paste0(output, "benchmarking_best_match_circle.pdf")
pdf(file_name, width = 10, height = 8)
print(g1)
dev.off()

all_data_layered$cell_labels <- as.factor(all_data_layered$cell_labels)
all_data_layered$Method <- as.factor(all_data_layered$Method)
all_data_layered$Method <- factor(all_data_layered$Method,
    levels = c("Reference","Query","Vesalius","CytoSpace","Tangram","SLAT","PASTE","GPSA","PRECAST"))
cols <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))
cols <- cols(length(levels(all_data_layered$cell_labels)))
g1 <- ggplot(all_data_layered, aes(x = x, y = y, col = cell_labels)) +
    geom_point(size = 0.8) +
    theme_bw() +
    theme(strip.background =element_rect(fill="#082233ff"),
        strip.text = element_text(colour = 'white', size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size=15),
        legend.text = element_text(size = 12),
        axis.text.y = element_text(size = 12))+
    scale_color_manual(values = cols) +
    facet_wrap(~Method) +
    guides(colour = guide_legend(
        override.aes = list(size =  5)))+
    labs(col = "Cell Labels")
file_name <- paste0(output, "benchmarking_best_match_layered.pdf")
pdf(file_name, width = 10, height = 8)
print(g1)
dev.off()


all_data_dropped$cell_labels <- as.factor(all_data_dropped$cell_labels)
all_data_dropped$Method <- as.factor(all_data_dropped$Method)
all_data_dropped$Method <- factor(all_data_dropped$Method,
    levels = c("Reference","Query","Vesalius","CytoSpace","Tangram","SLAT","PASTE","GPSA","PRECAST"))
cols <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))
cols <- cols(length(levels(all_data_dropped$cell_labels)))
g1 <- ggplot(all_data_dropped, aes(x = x, y = y, col = cell_labels)) +
    geom_point(size = 0.8) +
    theme_bw() +
    theme(strip.background =element_rect(fill="#082233ff"),
        strip.text = element_text(colour = 'white', size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size=15),
        legend.text = element_text(size = 12),
        axis.text.y = element_text(size = 12))+
    scale_color_manual(values = cols) +
    facet_wrap(~Method) +
    guides(colour = guide_legend(
        override.aes = list(size =  5)))+
    labs(col = "Cell Labels")
file_name <- paste0(output, "benchmarking_best_match_dropped.pdf")
pdf(file_name, width = 10, height = 8)
print(g1)
dev.off()


ref_circle$regime <- "Circle"
ref_circle$cell_labels <- as.factor(ref_circle$cell_labels)
cols <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))
cols <- cols(length(levels(ref_circle$cell_labels)))
g1 <- ggplot(ref_circle, aes(x = x, y = y, col = cell_labels)) +
    geom_point(size = 1.2) +
    theme_bw() +
    theme(strip.background =element_rect(fill="#082233ff"),
        strip.text = element_text(colour = 'white', size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size=15),
        legend.text = element_text(size = 12),
        axis.text.y = element_text(size = 12))+
    scale_color_manual(values = cols) +
    facet_wrap(~regime) +
    guides(colour = guide_legend(
        override.aes = list(size =  5)))+
    labs(col = "Cell Labels")
file_name <- paste0(output, "benchmarking_reference_circle.pdf")
pdf(file_name, width = 6, height = 5)
print(g1)
dev.off()


ref_layered$regime <- "Layer"
ref_layered$cell_labels <- as.factor(ref_layered$cell_labels)
cols <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))
cols <- cols(length(levels(ref_layered$cell_labels)))
g1 <- ggplot(ref_layered, aes(x = x, y = y, col = cell_labels)) +
    geom_point(size = 1.2) +
    theme_bw() +
    theme(strip.background =element_rect(fill="#082233ff"),
        strip.text = element_text(colour = 'white', size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size=15),
        legend.text = element_text(size = 12),
        axis.text.y = element_text(size = 12))+
    scale_color_manual(values = cols) +
    facet_wrap(~regime) +
    guides(colour = guide_legend(
        override.aes = list(size =  5)))+
    labs(col = "Cell Labels")
file_name <- paste0(output, "benchmarking_reference_layer.pdf")
pdf(file_name, width = 6, height = 5)
print(g1)
dev.off()


ref_dropped$regime <- "Dropped"
ref_dropped$cell_labels <- as.factor(ref_dropped$cell_labels)
cols <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))
cols <- cols(length(levels(ref_dropped$cell_labels)))
g1 <- ggplot(ref_dropped, aes(x = x, y = y, col = cell_labels)) +
    geom_point(size = 1.2) +
    theme_bw() +
    theme(strip.background =element_rect(fill="#082233ff"),
        strip.text = element_text(colour = 'white', size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size=15),
        legend.text = element_text(size = 12),
        axis.text.y = element_text(size = 12))+
    scale_color_manual(values = cols) +
    facet_wrap(~regime) +
    guides(colour = guide_legend(
        override.aes = list(size =  5)))+
    labs(col = "Cell Labels")
file_name <- paste0(output, "benchmarking_reference_layer.pdf")
pdf(file_name, width = 6, height = 5)
print(g1)
dev.off()