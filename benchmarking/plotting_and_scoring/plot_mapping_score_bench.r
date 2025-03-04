#-----------------------------------------------------------------------------#
# Loading Libraries and setting seed
#-----------------------------------------------------------------------------#
library(dplyr, lib.loc = "/common/martinp4/R")
library(tidyr)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(RColorBrewer)
set.seed(1547)

max_size <- 10000 * 1024^2
options(future.globals.maxSize = max_size)
#-----------------------------------------------------------------------------#
# Get args
#-----------------------------------------------------------------------------#
args <- commandArgs(TRUE)
input <- output <- args[1]


#-----------------------------------------------------------------------------#
# Load files
#-----------------------------------------------------------------------------#
scores <- list.files(
    input,
    pattern = "synthetic_benchmarking_scores.csv",
    full.names = TRUE)

scores <- as.data.frame(read.csv(scores))


ter_col <- length(unique(scores$Method))
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
cols <- ter_pal(ter_col)
names(cols) <- unique(scores$Method)

scores <- scores %>%
    pivot_longer(
        cols = c("ARI_cell","ARI_interactions","VI_cell","VI_interactions","JI_interactions"),
        names_to = "Score_method",
        values_to = "Score")
scores <- split(scores, list(scores$Score_method, scores$Data_Type))
#-----------------------------------------------------------------------------#
# Showing all scores
#-----------------------------------------------------------------------------#

for (i in seq_along(scores)) {
    local_scores <- scores[[i]]
    cols_local <- cols
    sub_scores <- split(local_scores, local_scores$Method)
    pvals <- lapply(sub_scores, function(s,a){
        local <- lapply(a, function(a,s){
                pval <- wilcox.test(a$Score,s$Score)$`p.value`
                df <- data.frame(
                    "target" = unique(a$Method),
                    "source" = unique(s$Method),
                    "Data_Type" = unique(a$Data_Type),
                    "pval" = pval)
                return(df)
            }, s = s)
        local <- do.call("rbind", local)
        return(local)
    }, a = sub_scores)
    pvals <- do.call("rbind", pvals)
    
    ns <- sapply(sub_scores, nrow) 
    sub_scores <- sapply(sub_scores,function(x){return(median(x$Score))})
    ord <- order(sub_scores, decreasing = TRUE)
    new_order <- names(sub_scores)[ord]
    cols_local <- cols_local[new_order]
    ns <- ns[new_order]
    ns_all <- ns[match(as.character(local_scores$Method),names(ns))]
    local_scores$Method <- paste0(local_scores$Method," - (n = ", ns_all,")")
    local_scores$Method <- as.factor(local_scores$Method)
    local_scores$Method  <- factor(local_scores$Method , levels = paste0(new_order," - (n = ", ns,")"))
    names(cols_local) <- paste0(new_order," - (n = ", ns,")")

    pvals$target <- as.factor(pvals$target)
    pvals$target <- factor(pvals$target, levels = new_order)
    pvals$source <- as.factor(pvals$source)
    pvals$source <- factor(pvals$source, levels = new_order)
    

    box <- ggplot(local_scores, aes(x = Method, y = Score, fill = Method)) +
        geom_boxplot() +
        scale_fill_manual(values = cols_local) +
        theme_bw() +
        theme(strip.background =element_rect(fill="#082233ff"),
            strip.text = element_text(colour = 'white', size = 15),
            axis.title = element_text(size = 15),
            legend.title = element_text(size=15),
            legend.text = element_text(size = 12),
            axis.text.y = element_text(size = 12),
            axis.text.x = element_text(size = 12, angle = 75, hjust = 1,vjust =1))+
        facet_wrap(~Data_Type)

    heat_cont <- ggplot() +
        geom_tile(data = pvals, 
                  aes(x = target, y = source, fill = pval)) +
        scale_fill_viridis_c(option = "mako") + 
        theme_bw() +
        theme(strip.background =element_rect(fill="#082233ff"),
            strip.text = element_text(colour = 'white', size = 15),
            axis.text.x = element_text(angle = 75, vjust = 1, hjust=1, size = 12),
            legend.title = element_text(size = 15),
            axis.text.y = element_text(size = 12),
            axis.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            plot.title  = element_text(size = 18)) +
        guides(colour = guide_legend(
                override.aes = list(linewidth = 5))) +
        labs(color = "pval < 0.05") +
        facet_wrap(~Data_Type)

    pvals$signif <- pvals$pval < 0.05
    heat_disc <- ggplot() +
        geom_tile(data = pvals, 
                  aes(x = target, y = source, fill = signif)) +
        scale_fill_manual(values = c("#E69F00", "#56B4E9")) + 
        theme_bw() +
        theme(strip.background =element_rect(fill="#082233ff"),
            strip.text = element_text(colour = 'white', size = 15),
            axis.text.x = element_text(angle = 75, vjust = 1, hjust=1, size = 12),
            legend.title = element_text(size = 15),
            axis.text.y = element_text(size = 12),
            axis.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            plot.title  = element_text(size = 18)) +
        guides(colour = guide_legend(
                override.aes = list(linewidth = 5)))+
        labs(fill = "pval < 0.05")+
        facet_wrap(~Data_Type)
    
    file_out <- paste0(
        output,
        "benchmarking_box_",
        unique(local_scores$Score_method),
        "_",
        unique(local_scores$Data_Type),
        ".pdf")
    pdf(file_out, width = 14, height = 6.5)
    print(box)
    dev.off()
    file_out <- paste0(
        output,
        "benchmarking_pval_",
        unique(local_scores$Score_method),
        "_",
        unique(local_scores$Data_Type),
        ".pdf")
    pdf(file_out, width = 8, height = 6.5)
    print(heat_cont)
    dev.off()
    file_out <- paste0(
        output,
        "benchmarking_pvalDisc_",
        unique(local_scores$Score_method),
        "_",
        unique(local_scores$Data_Type),
        ".pdf")
    pdf(file_out, width = 8, height = 6.5)
    print(heat_disc)
    dev.off()
}







