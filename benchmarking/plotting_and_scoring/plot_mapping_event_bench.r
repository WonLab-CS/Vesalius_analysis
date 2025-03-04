#-----------------------------------------------------------------------------#
# Loading Libraries and setting seed
#-----------------------------------------------------------------------------#
library(dplyr, lib.loc = "/common/martinp4/R")
library(ggplot2)
library(patchwork)
library(ggpubr)
library(RColorBrewer)


set.seed(1547)
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
# utils
#-----------------------------------------------------------------------------#

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


sanitize_labels <- function(labels) {
    labels <- strsplit(labels, "-")
    labels <- sapply(labels, "[[",1)
    return(unlist(labels))
}


prune_file_list <- function(scores, file_list, data_type) {
    ref <- scores$Ref[1L]
    query <- scores$Query[1L]
    data_loc <- sapply(file_list, function(f,d){
        return(f$Data_Type == d)
    }, d = data_type)
    file_list <- file_list[data_loc]
    ves_locs <- sapply(file_list, function(f,s){
            
            if (grepl("Vesalius",f$Method)){
                
                if (f$Method %in% s$Method){
                    return(TRUE)
                } else {
                    return(FALSE)
                }
            } else {
                return(TRUE)
            }
        }, s = scores)
    file_list <- file_list[ves_locs]
    ref_locs <- sapply(file_list, function(f, ref){
        return(as.character(f$Ref) == as.character(ref))
            
    },ref)
    query_locs <- sapply(file_list, function(f, query){
        return(as.character(f$Query) == as.character(query))
    },query)
    file_list <- file_list[ref_locs & query_locs]
    
    
    return(file_list)
}


scale_coordinates <- function(ref,query){
    ref$x <- ref$x - min(ref$x)
    ref$y <- ref$y - min(ref$y)
    x <- max(query$x) / max(ref$x)
    y <- max(query$y) / max(ref$y)
    ref$x <- ref$x * x
    ref$y <- ref$y * y
    return(ref)
}

load_files <- function(file_list, input_ref) {
    input_files <- list.files(
        input_ref,
        pattern = ".csv",
        full.names = TRUE)
    
    ref_file <- grep(paste0("sample_",file_list[[1]]$Ref,"\\.csv$"), input_files, value = TRUE)
    ref_file <- grep(file_list[[1]]$Data_Type, ref_file, value = TRUE)
    ref_file <- grep("spatial_coordinates", ref_file, value = TRUE)
    ref_file <- read.csv(ref_file)
    ref_file$Method <- "Reference"
    ref_file$Data_Type <- file_list[[1]]$Data_Type
    ref_file <- ref_file[,c("x","y","cell_labels","Method","Data_Type")]
    query_file <- grep(paste0("sample_",file_list[[1]]$Query,"\\.csv$"), input_files, value = TRUE)
    query_file <- grep(file_list[[1]]$Data_Type, query_file, value = TRUE)
    query_file <- grep("spatial_coordinates", query_file, value = TRUE)
    query_file <- read.csv(query_file)
    query_file$Method <- "Query"
    query_file$Data_Type <- file_list[[1]]$Data_Type
    query_file <- query_file[,c("x","y","cell_labels","Method","Data_Type")]
    mapped_files <- vector("list", length(file_list))
    for (i in seq_along(file_list)) {
        local <- read.csv(file_list[[i]]$file)
        local$Method <- file_list[[i]]$Method
        local$Data_Type <- file_list[[i]]$Data_Type
        local <- scale_coordinates(local,query_file)
        mapped_files[[i]] <- local[,c("x","y","cell_labels","Method","Data_Type")]
    }
    
    loaded_files <- c(list(ref_file), list(query_file),mapped_files)
    loaded_files <- do.call("rbind",loaded_files)
    return(loaded_files)
}




#-----------------------------------------------------------------------------#
# Mergde and prepare for plotting
#-----------------------------------------------------------------------------#
scores <- list.files(
    input_matched,
    recursive = TRUE,
    pattern = "benchmarking_scores.csv",
    full.names = TRUE)
scores <- grep(data_type,scores, value = TRUE)
scores <- as.data.frame(read.csv(scores))
scores$X <- NULL
scores <- scores[grep(x = scores$Method,pattern ="Vesalius"),]
scores <- split(scores, scores$Data_Type)

for (i in seq_along(scores)) {
    local <- scores[[i]]
    local <- split(local, local$Method)
    local <- local[match(c("Vesalius_fncty","Vesalius_fnct","Vesalius_fnt"), names(local))]
    local <- lapply(local, function(l){
        l <- l[which(l$ARI_cell == max(l$ARI_cell))[1L], c("Method","Ref","Query")]
        return(l)
    })
    scores[[i]] <- do.call("rbind", local)
}


#-----------------------------------------------------------------------------#
# Get mapped files and plot
#-----------------------------------------------------------------------------#
file_list <- get_files(
    input_matched,
    data_type,
    ref,
    query)


for (i in seq_along(scores)) {
    local_files <- prune_file_list(scores[[i]],file_list, names(scores)[i])
    local_files <- load_files(local_files, input_ref)
    levels <- c("Reference", "Query",scores[[i]]$Method,"CytoSpace","CytoSpace_noLab","Tangram","Scanorama","SLAT","GPSA","PASTE")
    local_files <- local_files[local_files$Method %in% levels,]
    local_files$Method <- as.factor(local_files$Method)
    local_files$Method <- factor(local_files$Method,levels = levels)
    max_cols <- length(unique(local_files$cell_labels))
    base_colours <- c(
        "#E69F00",
        "#56B4E9",
        "#009E73",
        "#F0E442",
        "#0072B2",
        "#D55E00",
        "#CC79A7",
        "#999999")
    if (max_cols < length(base_colours)) {
        ter_pal <- colorRampPalette(base_colours[seq(1, max_cols)])
    } else {
        ter_pal <- colorRampPalette(base_colours)
    }
    cols <- ter_pal(max_cols)

    p <- ggplot(local_files, aes(x,y, col = as.factor(cell_labels))) +
        geom_point(size = 1, alpha=0.75) +
        scale_color_manual(values = cols) +
        theme_bw() + 
        theme(strip.background =element_rect(fill="#082233ff"),
            strip.text = element_text(colour = 'white', size = 15),
            legend.title = element_text(size = 20),
            plot.title = element_text(size = 20),
            legend.text = element_text(size = 15),
            legend.position = "bottom") +
        labs(color = "") + 
        guides(colour = guide_legend(
            override.aes = list(size = 5)))+
        facet_wrap(~Method)
    data <- unique(local_files$Data_Type)
    file_name <- paste0(output,data, "_mapped_events.pdf")
    pdf(file_name, width = 10.5, height = 8)
    print(p)
    dev.off()
}
