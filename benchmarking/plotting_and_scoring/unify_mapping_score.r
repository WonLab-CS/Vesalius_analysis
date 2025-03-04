#-----------------------------------------------------------------------------#
# Loading Libraries and setting seed
#-----------------------------------------------------------------------------#
library(dplyr, lib.loc = "/common/martinp4/R")
library(igraph, lib.loc = "/common/martinp4/R")
library(future)
library(future.apply)
library(Matrix, lib.loc = "/common/martinp4/R")
library(ggplot2)
library(patchwork)
library(RANN)
library(lpSolve,lib.loc = "/common/martinp4/R")
library(mclust, lib.loc = "/common/martinp4/R")
library(mcclust, lib.loc = "/common/martinp4/R")
set.seed(1547)
plan(multicore, workers = 1)
max_size <- 10000 * 1024^2
options(future.globals.maxSize = max_size)
#-----------------------------------------------------------------------------#
# Output set up 
#-----------------------------------------------------------------------------#
args <- commandArgs(TRUE)
input_matched <- args[1]
input_ref <- args[2]
data_type <- args[3]
output <- args[4]
ref <- args[5]
query <- args[6]
#-----------------------------------------------------------------------------#
# Funks
#-----------------------------------------------------------------------------#


jaccard <- function(ref, query) {
    ref_split <- strsplit(ref, "@")
    query_split <- strsplit(query, "@")
    all_elements <- unique(unlist(c(ref_split, query_split)))
    ref_mat <- sapply(ref_split, function(x) as.integer(all_elements %in% x))
    query_mat <- sapply(query_split, function(x) as.integer(all_elements %in% x))
    intersection <- colSums(ref_mat & query_mat)
    union <- colSums(ref_mat | query_mat)
    jacc <- intersection / union
    return(mean(jacc))
}


get_files <- function(
    input,
    data_type,
    ref,
    query) {
    dirs <- list.dirs(input, recursive = FALSE)
    dirs <- grep("report", dirs, value = TRUE, invert = TRUE)
    files <- list.files(dirs,
        pattern = ".csv",
        recursive = TRUE,
        full.names = TRUE)
    file_tags <- list.files(dirs,
        pattern = ".csv",
        recursive = TRUE,
        full.names = FALSE)
    if (data_type == "synthetic"){
        tag <- "one_cell|two_cell|contact_one|contact_two|circle|layered|dropped"
    } else {
        tag <- data_type
    }
    files <- grep(tag, files, value = TRUE)
    file_tags <- grep(tag, file_tags, value = TRUE)
    files <- grep("contribution_score", files, value = TRUE, invert = TRUE)
    file_tags <- grep("contribution_score", file_tags, value = TRUE, invert = TRUE)
    files <- grep("computational_performance", files, value = TRUE, invert = TRUE)
    file_tags <- grep("computational_performance", file_tags, value = TRUE, invert = TRUE)
    
    files_list <- mapply(
        decompose_files,
        files,
        file_tags,
        MoreArgs = list(data_type, ref,query),
        SIMPLIFY = FALSE)
    names(files_list) <- NULL
    return(files_list)
}

decompose_files <- function(
    file,
    file_tag,
    data_type,
    ref,
    query){

    file_tag <- gsub(".csv","",file_tag)
    file_tag <- gsub("_aligned","",file_tag)
    file_tag <- gsub("report/","", file_tag)
    parts <- strsplit(file_tag,"_")[[1]]
    method <- parts[1]
    if (data_type == "synthetic"){
        types <- c("one_cell","two_cell","contact_one",
            "contact_two","circle","layered","dropped")
        loc <- sapply(types, grepl, file_tag)
        data_type <- types[loc]
        ref <- parts[min(grep(ref,parts)):min(grep(ref,parts)) + 1]
        query <- parts[max(grep(query,parts)):max(grep(query,parts)) + 1]
    }
    if (method == "Vesalius" && !"contribution_score" %in% parts) {
        combi <- parts[length(parts)]
        method <- paste0(method,"_", combi)
    } else if (method == "Vesalius" && "contribution_score" %in% parts) {
        combi <- parts[length(parts) - 1]
        method <- paste0(method,"_", combi)
    }
    if (method == "CytoSpace" && grepl(pattern = "noLab", x = file_tag)) {
        method <- paste0(method,"_","noLab")
    } 
    
    file_list <- list(
        "Method" = method,
        "Data_Type" = data_type,
        "Ref" = ref,
        "Query" = query,
        "file" = file)
    return(file_list)
}

load_files <- function(file_list, input_ref) {
    ref_file <- list.files(
        input_ref,
        pattern = ".csv",
        full.names = TRUE)
    ref_file <- grep(file_list$Ref, ref_file, value = TRUE)
    ref_file <- grep(file_list$Data_Type, ref_file, value = TRUE)
    ref_file <- grep("spatial_coordinates", ref_file, value = TRUE)
    ref_file <- read.csv(ref_file)
    query_file <- read.csv(file_list$file)
    loaded_files <- c(file_list, list("ref" = ref_file, "query" = query_file))
    return(loaded_files)
}

compute_score <- function(data) {
    print(data$Method)
    
    ref <- data$ref
    query <- data$query
    matched_locs <- RANN::nn2(
        data = ref[,c("x", "y")],
        query = query[,c("x","y")],
        k = 1)$nn.idx[,1]
    ref <- ref[matched_locs,]
    ari_cell <- adjustedRandIndex(ref$cell_labels, query$cell_labels)
    vi_cell <- vi.dist(ref$cell_labels, query$cell_labels)
    ari_inter <- adjustedRandIndex(ref$interactions, query$interactions)
    vi_inter <- vi.dist(ref$interactions, query$interactions)
    jaccard_index <- jaccard(ref$interactions,query$interactions)
    data <- c(
        data[!names(data) %in% c("file","ref","query")],
        setNames(list(ari_cell),"ARI_cell"),
        setNames(list(vi_cell),"VI_cell"),
        setNames(list(ari_inter),"ARI_interactions"),
        setNames(list(vi_inter),"VI_interactions"),
        setNames(list(jaccard_index),"JI_interactions"))
    return(data)
}

sanitize_labels <- function(labels) {
    labels <- strsplit(labels, "-")
    labels <- sapply(labels, "[[",1)
    return(unlist(labels))
}

score_cell_match <- function(file_list, input_ref) {
    files <- load_files(file_list, input_ref)
    score <- compute_score(files)
    return(score)
}

#-----------------------------------------------------------------------------#
# Get Files 
#-----------------------------------------------------------------------------#
file_list <- get_files(
    input_matched,
    data_type,
    ref,
    query)
#-----------------------------------------------------------------------------#
# build files
#-----------------------------------------------------------------------------#
scores <- lapply(file_list, score_cell_match, input_ref)
scores <- do.call("rbind", scores)
file_name <- paste0(output,data_type,"_benchmarking_scores.csv")
write.csv(scores,file = file_name)

