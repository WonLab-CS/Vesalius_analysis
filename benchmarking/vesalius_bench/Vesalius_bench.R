#-----------------------------------------------------------------------------#
# Loading Libraries and setting seed
#-----------------------------------------------------------------------------#
library(dplyr, lib.loc = "/common/martinp4/R")
library(igraph, lib.loc = "/common/martinp4/R")
library(future)
library(Matrix, lib.loc = "/common/martinp4/R")
library(ggplot2)
library(patchwork)
library(deldir)
library(imager)
library(imagerExtra)
library(Morpho)
library(lpSolve, lib = "/common/martinp4/R")
library(TreeDist, lib.loc = "/common/martinp4/R")
library(mclust, lib.loc = "/common/martinp4/R")
library(mcclust, lib.loc = "/common/martinp4/R")
library(pwr, lib.loc = "/common/martinp4/R")
library(gsignal, lib.loc = "/common/martinp4/R")
library(kohonen, lib.loc = "/common/martinp4/R")
library(registry, lib.loc = "/common/martinp4/R")
library(rngtools, lib.loc = "/common/martinp4/R")
library(NMF, lib.loc = "/common/martinp4/R")
library(spatstat.utils, lib.loc = "/common/martinp4/R")
library(vesalius, lib.loc = "/common/martinp4/R")
set.seed(1547)
#-----------------------------------------------------------------------------#
# Output set up 
#-----------------------------------------------------------------------------#
args <- commandArgs(TRUE)
idx <- as.numeric(args[1])
data_type <- args[2]
input <- args[3]
output_data <- args[4]
cost <- as.numeric(args[5])
cat("Output setup: DONE \n")


max_size <- 10000 * 1024^2
options(future.globals.maxSize = max_size)
#-----------------------------------------------------------------------------#
# Data set up
#-----------------------------------------------------------------------------#
synth_coord <- list.files(input,
    pattern = "_spatial_coordinates_",
    full.names = TRUE)
synth_coord <- grep(data_type, synth_coord, value = TRUE, fixed = TRUE)


synth_counts <- list.files(input,
    pattern = "_gene_counts_",
    full.names = TRUE)
synth_counts <- grep(data_type, synth_counts, value = TRUE, fixed = TRUE)


tag <- list.files(input,
    pattern = "_spatial_coordinates_",
    full.names = FALSE)
tag <- grep(data_type, tag, value = TRUE, fixed = TRUE)

i <- rep(seq_along(synth_coord), times = length(synth_coord))[idx]
j <- rep(seq_along(synth_coord), each = length(synth_coord))[idx]



if (i == j ) {
    q("no")
}


ref_tag <- gsub(".csv","", tag[i])
ref_tag <- gsub("spatial_coordinates_","",ref_tag)
query_tag <- gsub(".csv","", tag[j])
query_tag <- gsub("spatial_coordinates_","",query_tag)

#-----------------------------------------------------------------------------#
# Setting up parameters based on bash input
#-----------------------------------------------------------------------------#


if (grepl(x = data_type, pattern = "computational_performance")){
    batch_size <- 5000
} else {
    batch_size <- 1000
}

if (!is.na(cost)) {
    use_cost <- list(
        c("feature"),
        c("niche"),
        c("territory"),
        c("composition"),
        c("feature", "niche"),
        c("feature", "niche","cell_type"),
        c("feature","composition"),
        c("feature", "territory"),
        c("niche","territory"),
        c("niche","composition"),
        c("feature","niche","territory"),
        c("feature", "niche","composition"),
        c("feature", "niche","composition","territory"),
        c("feature", "niche","composition","territory","cell_type"))
    use_cost <- use_cost[[cost]]
    compressed_tag <- paste(sapply(gsub("cell_type", "y", use_cost),substr,1,1), collapse = "")
    query_tag <- paste0(query_tag,"_",compressed_tag)
} else {
    use_cost <- c("feature", "cell_type")
}
cat(paste(data_type, "\n"))
cat(paste(query_tag, "\n"))


cat("Data setup: DONE\n")


#-----------------------------------------------------------------------------#
# Prepare seed
# Only run territory isolation if required
#-----------------------------------------------------------------------------#
ref_coord <- read.csv(synth_coord[i])
ref_counts <- read.csv(synth_counts[i], row.names = 1)
ref_counts$genes <- NULL
ref <- build_vesalius_assay(coordinates = ref_coord,
        counts = ref_counts,
        assay = paste0(data_type,"_",unique(ref_coord$sample)))
ref_cells <- ref_coord$cell_labels
names(ref_cells) <- ref_coord$barcodes
ref <- add_cells(ref, ref_cells, add_name = "cell_labels")
ref_interactions <- ref_coord$interactions
names(ref_interactions) <- ref_coord$barcodes
ref <- add_cells(ref,ref_interactions, add_name = "interactions")

ref <- ref %>%
     generate_embeddings()

if ("territory" %in% use_cost){
    ref <- ref %>%
        smooth_image(sigma = 2, iter = 15) %>%
        equalize_image(sleft = 2.5, sright = 2.5) %>%
        segment_image(col_resolution = 12) %>%
        isolate_territories()
}
cat("Ref Procesing: DONE\n")

#-----------------------------------------------------------------------------#
# Prepare query
# Only run territory isolation if required
#-----------------------------------------------------------------------------#
query_coord <- read.csv(synth_coord[j])
query_counts <- read.csv(synth_counts[j], row.names = 1)
query_counts$genes <- NULL
query <- build_vesalius_assay(coordinates = query_coord,
        counts = query_counts,
        assay = paste0(data_type,"_",unique(query_coord$sample)))
query_cells <- query_coord$cell_labels
names(query_cells) <- query_coord$barcodes
query <- add_cells(query, query_cells, add_name = "cell_labels")
query_interactions <- query_coord$interactions
names(query_interactions) <- query_coord$barcodes
query <- add_cells(query, query_interactions, add_name = "interactions")


query <- query %>%
    generate_embeddings()

if ("territory" %in% use_cost){
    query <- query %>%
        smooth_image(sigma = 2, iter = 10) %>%
        equalize_image(sleft = 2.5, sright = 2.5) %>%
        segment_image(col_resolution = 12) %>%
        isolate_territories()
}
cat("Query Procesing: DONE\n")


#-----------------------------------------------------------------------------#
# Mapping
#-----------------------------------------------------------------------------#
matched <- map_assays(seed_assay = ref,
    query_assay = query,
    neighborhood = "knn",
    use_norm = "log_norm",
    method = "pearson",
    k = 6,
    epochs = 25,
    batch_size = batch_size,
    allow_duplicates = TRUE,
    threshold = 0.9,
    #threshold = 0,
    filter_cells = TRUE,
    jitter = 0,
    use_cost = use_cost,
    seed_meta_labels = c("cell_labels","interactions"),
    query_meta_labels = c("cell_labels","interactions"),
    digits = 5)

cat("Mapping: DONE\n")
#-----------------------------------------------------------------------------#
# Adding and exporting scores for combinations
#-----------------------------------------------------------------------------#
if (!grepl(x = data_type, pattern = "computational_performance")){
    contributions <- c("CV","IQR","POC")
    for (i in contributions) {
        matched <- get_cost_contribution(matched,method = i)
    }
    contrib_loc <- grep("contribution_score", names(matched@cost))
    contribution_list <- vector("list", length(contrib_loc))
    for (i in seq_along(contrib_loc)) {
        contrib_local <- matched@cost[[contrib_loc[i]]]
        contrib_local$method <- contributions[i]
        colnames(contrib_local) <- gsub(contributions[i],"score",colnames(contrib_local))
        contribution_list[[i]] <- contrib_local
    }
    
    contribution_list <- do.call("rbind", contribution_list)
    file_name <- paste0(output_data,"Vesalius_aligned_",data_type,"_",ref_tag,"_",query_tag,"_contribution_score.csv")
    write.csv(contribution_list, file = file_name)
    
} else {
    q("no")
}

#-----------------------------------------------------------------------------#
# Export Mapping 
#-----------------------------------------------------------------------------#
export_match <- matched@territories
file_name <- paste0(output_data,"Vesalius_aligned_",data_type,"_",ref_tag,"_",query_tag,".csv")

write.csv(export_match, file = file_name) 

gc()
q("no")