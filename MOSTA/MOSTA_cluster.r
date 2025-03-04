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
output_data <- output_plots <- args[3]

if (!dir.exists("/common/wonklab/Stereo_seq/report/")) {
    dir.create("/common/wonklab/Stereo_seq/report/")
}
output_plots <- "/common/wonklab/Stereo_seq/report/"
output_data <- "/common/wonklab/Stereo_seq/report/"

cat("Directory setup: done \n")
#-----------------------------------------------------------------------------#
# Get Single Data sets and preparing output files
#-----------------------------------------------------------------------------#
files <- list.files(output_data, pattern = "_to_", full.names = TRUE)
files <-files[grepl(".rds", files) &
    !grepl("integrated", files) &
    !grepl("metric_clusters", files)]
files <- files[idx]
tags <- list.files(output_data, pattern = "_to_", full.names = FALSE)
tags <- grep(".rds", tags, value = TRUE)
tags <-tags[grepl(".rds", tags) &
    !grepl("integrated", tags) &
    !grepl("metric_clusters", tags)]
tags <- tags[idx]
tags <- gsub("_stereo.rds","", tags)
# dirty skip
if (tags != "E11.5_to_E12.5"){
    q("no")
}
matched <- readRDS(files)
matched_cells <- matched@territories
stage_seed <- sapply(strsplit(tags, "_to_"),"[[", 2)
stage_query <- sapply(strsplit(tags, "_to_"),"[[", 1)


ref_files <-  list.files(output_data, pattern = ".rds", full.names = TRUE)
ref_files <- ref_files[grepl(pattern = stage_seed, x = ref_files) & !grepl(pattern = "_to_", x = ref_files)]
seed <- readRDS(ref_files)

query_files <-  list.files(output_data, pattern = ".rds", full.names = TRUE)
query_files <- query_files[grepl(pattern = stage_query, x = query_files) & !grepl(pattern = "_to_", x = query_files)]
query <- readRDS(query_files)


cat("Input load: done\n")


file_name <- paste0(output_data,tags,"_DEGs_intra_cluster.csv")
if (!file.exists(file_name)) {
    degs <- data.frame("genes","p_value","p_value_adj","seed_pct","query_pct","fold_change","effect_size","seed","query","Cells","cost")
    write.table(degs,
        file = file_name,
        quote = FALSE,
        col.names = TRUE,
        row.names = FALSE,
        sep = ",")
}

file_name <- paste0(output_data,tags,"_DEGs_inter_sample_inferred.csv")
if (!file.exists(file_name)) {
    degs <- data.frame("genes","p_value","p_value_adj","seed_pct","query_pct","fold_change","effect_size","seed","query","Cells")
    write.table(degs,
        file = file_name,
        quote = FALSE,
        col.names = TRUE,
        row.names = FALSE,
        sep = ",")
}
file_name <- paste0(output_data,tags,"_DEGs_inter_sample_scaled.csv")
if (!file.exists(file_name)) {
    degs <- data.frame("genes","p_value","p_value_adj","seed_pct","query_pct","fold_change","effect_size","seed","query","Cells")
    write.table(degs,
        file = file_name,
        quote = FALSE,
        col.names = TRUE,
        row.names = FALSE,
        sep = ",")
}
cat("Output setup: done\n")
#-----------------------------------------------------------------------------#
# Run clustering for each cell
#-----------------------------------------------------------------------------#
cells <- unique(matched_cells$Cells)
use_cost <- list("feature", "niche","territory")#, c("feature", "niche","territory"))
for (i in seq_along(cells)){
    for (j in seq_along(use_cost)) {
        cat(paste0("Cost: ",paste0(use_cost[[j]],collapse = " "))," Cells: ",cells[i],"\n")
        matched <- get_metric_clusters(matched,
            use_cost = use_cost[[j]],
            trial = "Cells",
            group_identity = cells[i],
            k = 5)
        clusters <- unique(matched@territories$Map_cluster)
        clusters <- clusters[clusters != "Not Selected"]
        degs <- vector("list", length(clusters))
        for (k in seq_along(clusters)){
            matched <- identify_markers(matched,
                trial = "Map_cluster",
                seed = clusters[k],
                query = clusters[clusters != clusters[k]],
                min_spatial_index = 10)
            df <- get_markers(matched)
            if (is.null(df)){
                degs[[k]] <- NULL
            } else if (nrow(df) == 0){
                degs[[k]] <- df
            } else {
                df$Cells <- cells[i]
                df$cost <- paste0(use_cost[[j]], collapse = "_")
                degs[[k]] <- df
            }
        }
        degs <- do.call("rbind", degs)
        file_name <- paste0(output_data,tags,"_DEGs_intra_cluster.csv")
        write.table(degs,
            file = file_name,
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            sep = ",",
            append = TRUE)

        cost_tag <- ifelse(length(use_cost[[j]])>1,"all",use_cost[[j]])
        colnames(matched@territories) <- gsub("Map_cluster",
            paste0(cost_tag,"_",cells[i]),
            colnames(matched@territories))
        file_out <- paste0(output_plots,tags,"_",paste0(cost_tag,"_",cells[i]), ".pdf")
        pdf(file_out, width = 10, height = 8)
        print(territory_plot(matched, trial = paste0(cost_tag,"_",cells[i])),cex_pt = 0.2) +
            labs(title = paste0(cost_tag,"_",cells[i]))
        dev.off()
    }
}
file_name <- paste0(output_data,tags,"_metric_clusters.rds")
saveRDS(matched, file = file_name)
cat("Clustering and DEG by cell: done\n")
#-----------------------------------------------------------------------------#
# integrate and compare mapped cells with cell types
#-----------------------------------------------------------------------------#

inter <- integrate_assays(matched, seed, infer = TRUE)
file_name <- paste0(output_data,tags,"_integrated.rds")
saveRDS(inter, file = file_name)
cat("Integration: done\n")

cells <- intersect(unique(matched@territories$Cells),unique(seed@territories$Cells))
degs <- vector("list", length(cells))
for (i in seq_along(cells)){
    tmp <- identify_markers(inter,
        norm_method = "inferred",
        trial = "Cells",
        seed = cells[i],
        query = cells[i],
        sample = TRUE)
    df <- get_markers(tmp)
    if (is.null(df)){
        degs[[i]] <- NULL
    } else if (nrow(df) == 0){
        degs[[i]] <- df
    } else {
        df$Cells <- cells[i]
        degs[[i]] <- df
    }
}
    
degs <- do.call("rbind", degs)
file_name <- paste0(output_data,tags,"_DEGs_inter_sample_inferred.csv")
write.table(degs,
    file = file_name,
    quote = FALSE,
    col.names = TRUE,
    row.names = FALSE,
    sep = ",",
    append = TRUE)


cells <- intersect(unique(matched@territories$Cells),unique(seed@territories$Cells))
degs <- vector("list", length(cells))
for (i in seq_along(cells)){
    tmp <- identify_markers(inter,
        norm_method = "scaled",
        trial = "Cells",
        seed = cells[i],
        query = cells[i],
        sample = TRUE)
    df <- get_markers(tmp)
    if (is.null(df)){
        degs[[i]] <- NULL
    } else if (nrow(df) == 0){
        degs[[i]] <- df
    } else {
        df$Cells <- cells[i]
        degs[[i]] <- df
    }
}
    
degs <- do.call("rbind", degs)
file_name <- paste0(output_data,tags,"_DEGs_inter_sample_scaled.csv")
write.table(degs,
    file = file_name,
    quote = FALSE,
    col.names = TRUE,
    row.names = FALSE,
    sep = ",",
    append = TRUE)

cat("DEG by sample: done\n")


#-----------------------------------------------------------------------------#
# Plotting genes and clusters
#-----------------------------------------------------------------------------#

file_name <- paste0(output_data,tags,"_metric_clusters.rds")
matched <- readRDS(file_name)


cluster_col <- matched@territories
for (i in seq(5, ncol(cluster_col))) {
    g <- territory_plot(matched, trial = colnames(cluster_col)[i], randomise = FALSE, cex_pt = 0.5,alpha =1, highlight = 1:5) +
    theme_void()
    file_out <- paste0(output_plots,tags,"_",colnames(cluster_col)[i], ".pdf")
    pdf(file_out, width = 9, height = 8)
    print(g)
    dev.off()
}

file_name <- paste0(output_data,tags,"_DEGs_intra_cluster.csv")
degs <- read.csv(file_name, header = TRUE, skip=2)

brain <- degs[degs$Cells == "Brain", ]
brain <- brain %>%
    group_by(seed,cost) %>%
    slice_max(fold_change,n=10)
for (i in seq_len(nrow(brain))){
    g <- view_gene_expression(matched, norm_method = "raw", genes = brain$genes[i], cex_pt = 0.35) +
        theme_void() +
        labs(title = paste0(brain$genes[i]," ", brain$cost[i]))
    g1 <- view_gene_expression(seed, norm_method = "raw", genes = brain$genes[i], cex_pt = 0.35) +
        theme_void() +
        labs(title = paste0("Reference ", brain$genes[i]," ", brain$cost[i]))
    g2 <- view_gene_expression(query, norm_method = "raw", genes = brain$genes[i], cex_pt = 0.35) +
        theme_void() +
        labs(title = paste0("Query ", brain$genes[i]," ", brain$cost[i]))
    g_all <- g1 + g + g2
    file_out <- paste0(output_plots,tags,"_brain_",paste0(brain$genes[i],"_", brain$cost[i],"_",brain$seed[i]), ".pdf")
    pdf(file_out, width = 21, height = 8)
    print(g_all)
    dev.off()

}



file_name <- paste0(output_data,tags,"_DEGs_inter_sample.csv")
degs_inter <- read.csv(file_name, header = TRUE, skip =2)

brain <- degs_inter[degs_inter$seed == "Brain_matched", ]
brain <- brain %>%
    slice_max(fold_change,n=13)
for (i in seq_len(nrow(brain))){
    g <- view_gene_expression(matched, norm_method = "raw", genes = brain$genes[i], cex_pt = 0.35) +
        theme_void() +
        labs(title = brain$genes[i])
    g1 <- view_gene_expression(seed, norm_method = "raw", genes = brain$genes[i], cex_pt = 0.35) +
        theme_void() +
        labs(title = paste0("Reference ", brain$genes[i]))
    g2 <- view_gene_expression(query, norm_method = "raw", genes = brain$genes[i], cex_pt = 0.35) +
        theme_void() +
        labs(title = paste0("Query ", brain$genes[i]))
    g_all <- g1 + g + g2
    file_out <- paste0(output_plots,tags,"_brain_inter_",paste0(brain$genes[i],"_", brain$cost[i],"_",brain$seed[i]), ".pdf")
    pdf(file_out, width = 21, height = 8)
    print(g_all)
    dev.off()

}


q("no")