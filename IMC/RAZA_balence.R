#-----------------------------------------------------------------------------#
# Loading Libraries and setting seed
#-----------------------------------------------------------------------------#
library(dplyr)
library(igraph)
library(future)
library(Matrix)
library(ggplot2)
library(patchwork)
library(deldir)
library(imager)
library(imagerExtra)
library(Morpho)
library(lpSolve)
library(TreeDist)
library(mclust)
library(mcclust)
library(pwr)
library(gsignal)
library(kohonen)
library(registry)
library(rngtools)
library(NMF)
library(spatstat.utils)
library(vesalius)
set.seed(1547)

plan(multicore, workers = 2)
max_size <- 10000 * 1024^2
options(future.globals.maxSize = max_size)
#-----------------------------------------------------------------------------#
# Set future global for multicore processing
#-----------------------------------------------------------------------------#

if (!dir.exists("/common/wonklab/RAZA/output_plots/")) {
    dir.create("/common/wonklab/RAZA/output_plots/")
}
output_plots <- "/common/wonklab/RAZA/output_plots/"

if (!dir.exists("/common/wonklab/RAZA/matched/")) {
    dir.create("/common/wonklab/RAZA/matched/")
}
output_data <- "/common/wonklab/RAZA/matched/"
cat("Output setup: DONE \n")

#-----------------------------------------------------------------------------#
# Data set up 
#-----------------------------------------------------------------------------#
args <- commandArgs(TRUE)
idx <- as.numeric(args[1])
train <- list.files("/common/wonklab/RAZA/split_data",
    pattern = "_balence.Rda",
    full.names = TRUE)

test <- list.files("/common/wonklab/RAZA/split_data",
    pattern = "_balence.Rda",
    full.names = TRUE)

i <- rep(seq_along(train), times = length(train))[idx]
j <- rep(seq_along(test), each = length(test))[idx]

if (i == j) {
    q("no")
}

use_cost <- c("feature","niche","composition","territory","cell_type")

cat("Data setup: DONE\n")

#-----------------------------------------------------------------------------#
# Initialize output
#-----------------------------------------------------------------------------#
file_name <- paste0(output_data, "feature_score_matrix_balence.csv")
if (!file.exists(file_name)){
    feature_score <- data.frame("from", "to", "score")
    write.table(feature_score,
        file = file_name,
        quote = FALSE,
        col.names = TRUE,
        row.names = FALSE,
        sep = ",")
}

file_name <- paste0(output_data, "niche_score_matrix_balence.csv")
if (!file.exists(file_name)){
    niche_score <- data.frame("from", "to", "score")
    write.table(niche_score,
        file = file_name,
        quote = FALSE,
        col.names = TRUE,
        row.names = FALSE,
        sep = ",")
}


file_name <- paste0(output_data, "cell_label_score_matrix_balence.csv")
if (!file.exists(file_name)){
    label_score <- data.frame("from", "to", "score")
    write.table(label_score,
        file = file_name,
        quote = FALSE,
        col.names = TRUE,
        row.names = FALSE,
        sep = ",")
}

file_name <- paste0(output_data, "territory_score_matrix_balence.csv")
if (!file.exists(file_name)){
    territory_score <- data.frame("from", "to", "score")
    write.table(territory_score,
        file = file_name,
        quote = FALSE,
        col.names = TRUE,
        row.names = FALSE,
        sep = ",")
}

file_name <- paste0(output_data, "composition_score_matrix_balence.csv")
if (!file.exists(file_name)){
    composition_score <- data.frame("from", "to", "score")
    write.table(composition_score,
        file = file_name,
        quote = FALSE,
        col.names = TRUE,
        row.names = FALSE,
        sep = ",")
}


file_name <- paste0(output_data, "cost_score_matrix_balence.csv")
if (!file.exists(file_name)){
    cost_score <- data.frame("from", "to", "score")
    write.table(cost_score,
        file = file_name,
        quote = FALSE,
        col.names = TRUE,
        row.names = FALSE,
        sep = ",")
}
cat("Output init: DONE\n")

#-----------------------------------------------------------------------------#
# Start mapping
#-----------------------------------------------------------------------------#
ref <- get(load(train[i]))
ref <- ref %>%
    generate_embeddings() %>%
    smooth_image(dimensions = seq(1, 30), sigma = 1, iter = 10) %>%
    segment_image(dimensions = seq(1, 30), method = "kmeans", col_resolution  = 10) %>%
    isolate_territories()
cat("Ref Procesing: DONE\n")


query <- get(load(test[j]))
query <- query %>%
    generate_embeddings() %>%
    smooth_image(dimensions = seq(1, 30), sigma = 1, iter = 10) %>%
    segment_image(dimensions = seq(1, 30), method = "kmeans", col_resolution  = 10) %>%
    isolate_territories()
cat("Query Procesing: DONE\n")

no_batch <- min(c(nrow(ref@territories),nrow(query@territories), 1000)) - 1

matched <- map_assays(seed_assay = ref,
    query_assay = query,
    neighborhood = "graph",
    depth = 1,
    threshold = -1,
    batch_size = no_batch,
    epochs = 10,
    use_cost = use_cost)

file_name <- paste0(output_data, "matched_RAZA_ref_",get_assay_names(ref),"_query_",get_assay_names(query),"_balence.rda")
#save(matched, file = file_name)
cat("Mapping: DONE\n")

#-----------------------------------------------------------------------------#
# Output scores
#-----------------------------------------------------------------------------#
feature_sum <- mean(matched@map$feature)
feature_sum <- data.frame("from" = get_assay_names(query),
            "to" = get_assay_names(ref),
            "score" = feature_sum)
write.table(feature_sum,
    file = paste0(output_data, "feature_score_matrix_balence.csv"),
    quote = FALSE,
    col.names = FALSE,
    row.names = FALSE,
    sep = ",",
    append  = TRUE)
        
niche_sum <- mean(matched@map$niche)
niche_sum <- data.frame("from" = get_assay_names(query),
    "to" = get_assay_names(ref),
    "score" = niche_sum)
write.table(niche_sum,
    file = paste0(output_data,"niche_score_matrix_balence.csv"),
    quote = FALSE,
    col.names = FALSE,
    row.names = FALSE,
    sep = ",",
    append  = TRUE)

label_sum <- mean(matched@map$cell_type)
label_sum <- data.frame("from" = get_assay_names(query),
    "to" = get_assay_names(ref),
    "score" = label_sum)
write.table(label_sum,
    file = paste0(output_data,"cell_label_score_matrix_balence.csv"),
    quote = FALSE,
    col.names = FALSE,
    row.names = FALSE,
    sep = ",",
    append  = TRUE)

territory_sum <- mean(matched@map$territory)
territory_sum <- data.frame("from" = get_assay_names(query),
    "to" = get_assay_names(ref),
     "score" = territory_sum)
write.table(territory_sum,
    file = paste0(output_data,"territory_score_matrix_balence.csv"),
    quote = FALSE,
    col.names = FALSE,
    row.names = FALSE,
    sep = ",",
    append  = TRUE)


composition_sum <- mean(matched@map$composition)
composition_sum <- data.frame("from" = get_assay_names(query),
    "to" = get_assay_names(ref),
     "score" = composition_sum)
write.table(composition_sum,
    file = paste0(output_data,"composition_score_matrix_balence.csv"),
    quote = FALSE,
    col.names = FALSE,
    row.names = FALSE,
    sep = ",",
    append  = TRUE)

cost_sum <- mean(matched@map$cost)
cost_sum <- data.frame("from" = get_assay_names(query),
    "to" = get_assay_names(ref),
     "score" = cost_sum)
write.table(cost_sum,
    file = paste0(output_data,"cost_score_matrix_balence.csv"),
    quote = FALSE,
    col.names = FALSE,
    row.names = FALSE,
    sep = ",",
    append  = TRUE)
cat("Export scores: DONE\n")