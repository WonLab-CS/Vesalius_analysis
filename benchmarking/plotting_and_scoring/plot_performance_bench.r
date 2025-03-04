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
input_matched <- args[1]
input_ref <- args[2]
data_type <- args[3]
output <- args[4]
ref <- args[5]
query <- args[6]


#-----------------------------------------------------------------------------#
# Load files
#-----------------------------------------------------------------------------#
scores <- list.files(
    input_matched,
    pattern = "comp_performance_metrics",
    recursive = TRUE,
    full.names = TRUE)
tags <- list.files(
    input_matched,
    pattern = "comp_performance_metrics",
    recursive = TRUE,
    full.names = FALSE)
tags <- sapply(strsplit(tags,"/"),"[[",1)
scores <- lapply(scores,read.csv)
scores <- lapply(scores, function(x){
    local <- split(x, x$NumericValue)
    values <- names(local)
    local <- lapply(local,function(l){
        buffer <- l[,c("MaxMemoryGB","RuntimeSeconds")]

        if ("Exceed_limit" %in% buffer$MaxMemoryGB | "Exceed_limit" %in% buffer$RuntimeSeconds) {
            
            return(c("Exceed_limit","Exceed_limit"))
        } else {
            buffer$MaxMemoryGB <- as.numeric(buffer$MaxMemoryGB)
            buffer$RuntimeSeconds <- as.numeric(buffer$RuntimeSeconds)
            tmp <- apply(buffer,2,median)
            return(tmp)
        }
       
    })
    df <- data.frame(
        "Method" = unique(x$Method),
        "N_cells" = values,
        do.call("rbind",local))
    return(df)
})
scores <- do.call("rbind",scores)
scores <- scores[scores$MaxMemoryGB != "Exceed_limit",]
scores$run_time <- "Mean Run Time"
scores$mem <- "Max Memory Usage"
#-----------------------------------------------------------------------------#
# Plot
#-----------------------------------------------------------------------------#

cell_levels <- as.character(c(100, 500, 1000, 2500, 5000, 10000, 20000, 50000, 100000, 200000))
scores$N_cells <- as.factor(scores$N_cells)
scores$N_cells <- factor(scores$N_cells, levels = cell_levels)


run <- ggplot(scores, aes(x = N_cells, y = Method , fill = log(as.numeric(RuntimeSeconds)))) +
    geom_tile()+
    theme_bw() +
    scale_fill_viridis_c(option = "mako")+
    theme(strip.background =element_rect(fill="#082233ff"),
        strip.text = element_text(colour = 'white', size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size=15),
        legend.text = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 35, hjust = 1,vjust =1))+
    labs(y = "Method",x = "Number of Cells",fill = "Log(sec)") +
    guides(colour = guide_legend(
                override.aes = list(linewidth = 5))) +
    facet_wrap(~run_time)


mem <- ggplot(scores, aes(x = N_cells, y = Method, fill = as.numeric(MaxMemoryGB))) +
    geom_tile()+
    theme_bw() +
    scale_fill_viridis_c(option = "mako")+
    theme(strip.background =element_rect(fill="#082233ff"),
        strip.text = element_text(colour = 'white', size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size=15),
        legend.text = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 35, hjust = 1,vjust =1))+
    labs(y = "Method",x = "Number of Cells",fill = "GB") +
    guides(colour = guide_legend(
                override.aes = list(linewidth = 5))) +
    facet_wrap(~mem)

    
file_out <- paste0(
    output,"run_time.pdf")
pdf(file_out, width = 6, height = 5)
print(run)
dev.off()

file_out <- paste0(
    output,"mem_usage.pdf")
pdf(file_out, width = 6, height = 5)
print(mem)
dev.off()
    
    


#-----------------------------------------------------------------------------#
# Load files
#-----------------------------------------------------------------------------#
scores <- list.files(
    input_matched,
    pattern = "Vesalius_comp_performance_combinations.csv",
    recursive = TRUE,
    full.names = TRUE)
tags <- list.files(
    input_matched,
    pattern = "Vesalius_comp_performance_combinations.csv",
    recursive = TRUE,
    full.names = FALSE)
tags <- sapply(strsplit(tags,"/"),"[[",1)
scores <- read.csv(scores)
# removing failed runs - forgot to remove the master file where verything was being dumped
scores <- split(scores, scores$TaskID)[3:4]
scores <- do.call("rbind",scores)
scores <- split(scores, list(scores$Combination,scores$NumericValue))
scores <- lapply(scores, function(s){
    mem <- mean(s$MaxMemoryGB)
    run_time <- mean(s$RuntimeSeconds)
    df <- data.frame(
        "Combination" = unique(s$Combination),
        "N_cells" = unique(s$NumericValue),
        "mem" = mem,
        "run_time" = run_time)
    return(df)


})
scores <- do.call("rbind", scores)




levels <- c("f","n","c","t","fn","fc","ft","nc","nt","fnt","fnc","fny","fnct","fncty")
scores$Combination <- as.factor(scores$Combination)
scores$Combination <- factor(scores$Combination, levels = levels)

run <- ggplot(scores, aes(x = Combination, y = run_time , fill = as.numeric(run_time))) +
    geom_bar(position = "stack",stat = "identity")+
    theme_bw() +
    scale_fill_viridis_c(option = "mako")+
    theme(strip.background =element_rect(fill="#082233ff"),
        strip.text = element_text(colour = 'white', size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size=15),
        legend.text = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 65, hjust = 1,vjust =1))+
    labs(y = "Run Time in Seconds",x = "Combinations",fill = "Seconds") +
    guides(colour = guide_legend(
                override.aes = list(linewidth = 5))) +
    facet_wrap(~N_cells)




mem <- ggplot(scores, aes(x = Combination, y = mem, fill = as.numeric(mem))) +
    geom_bar(position = "stack",stat = "identity")+
    theme_bw() +
    scale_fill_viridis_c(option = "mako")+
    theme(strip.background =element_rect(fill="#082233ff"),
        strip.text = element_text(colour = 'white', size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size=15),
        legend.text = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 65, hjust = 1,vjust =1))+
    labs(y = "Max Memory Usage in GB",x = "Combinations",fill = "GB") +
    guides(colour = guide_legend(
                override.aes = list(linewidth = 5))) +
    facet_wrap(~N_cells)

    
file_out <- paste0(
    output,"run_time_combi.pdf")
pdf(file_out, width = 8, height = 5)
print(run)
dev.off()

file_out <- paste0(
    output,"mem_usage_combi.pdf")
pdf(file_out, width = 8, height = 5)
print(mem)
dev.off()
    

