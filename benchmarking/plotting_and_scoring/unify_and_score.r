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
library(deldir)
library(imager)
library(imagerExtra)
library(Morpho)
library(RANN)
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
plan(multicore, workers = 10)
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
data_path <- "/common/wonklab/synthetic_spatial/"
#-----------------------------------------------------------------------------#
# Funks
#-----------------------------------------------------------------------------#
rmse <- function(ref, matched) {
    ref <- as.vector(as.matrix(ref))
    matched <- as.vector(as.matrix(matched))
    mse <- mean(ref - matched) ^ 2
    rmse <- sqrt(mse)
    return(rmse)
}

score_cell_match <- function(file, file_tag, data_path) {
    cat(file_tag)
    cat("\n")
    dat <- unlist(strsplit(file_tag ,"_"))
    method <- dat[1L]
    regime <- dat[3L]
    ref <- dat[5L]
    query <- gsub(".csv","",dat[8L])
    if (ref == query){
        # no self comparison
        return(NULL)
    }
    path <- list.files(data_path, full.names = TRUE)
    ref_path <- path[grepl(paste0("sample_",ref,".csv"),path) &
        grepl(regime,path)]
    query_path <- path[grepl(paste0("sample_",query,".csv"),path) &
        grepl(regime,path)]
    ref_coord <- read.csv(ref_path[2L])
    ref_counts <- read.csv(ref_path[1L], row.names = 1)
    ref_counts$genes <- NULL
    query_coord <- read.csv(query_path[2L])
    query_counts <- read.csv(query_path[1L],row.names =1)
    query_counts$genes <- NULL
    matched <- read.csv(file)
    matched <- matched[,-1]
    if (any(colnames(matched) %in% c("row","col"))){
        colnames(matched) <- c("barcodes","x","y","cluster","cell_labels")
    }
    matched$barcodes <- gsub("q_","", matched$barcodes)
    matched$cell_labels <- gsub("celltype_","",matched$cell_labels)
    matched_locs <- RANN::nn2(data = ref_coord[,c("x", "y")],
        query = matched[,c("x","y")],
        k = 1)$nn.idx[,1]
    ref_coord <- ref_coord[matched_locs,]
    ari <- adjustedRandIndex(ref_coord$cell_labels, matched$cell_labels)
    vi <- vi.dist(ref_coord$cell_labels, matched$cell_labels)
    root_mse <- rmse(ref_counts[,matched_locs],query_counts[,matched$barcodes])
    score <- data.frame("Method" = method,
        "Regime" = regime,
        "ref_sample" = ref,
        "query_sample" = query,
        "ARI" = ari,
        "VI" = vi,
        "RMSE" = root_mse)
    return(score)
}

#-----------------------------------------------------------------------------#
# Get Files 
#-----------------------------------------------------------------------------#
files <- list.files(list.dirs(input, recursive = FALSE),
    pattern = ".csv",
    recursive = TRUE,
    full.names = TRUE)
files <- sort(grep("aligned",files, value = TRUE))

file_tags <- list.files(list.dirs(input, recursive = FALSE),
    pattern = ".csv",
    recursive = TRUE,
    full.names = FALSE)
file_tags <- grep("aligned",file_tags, value = TRUE)
file_tags <- sort(gsub("report/","", file_tags))

#-----------------------------------------------------------------------------#
# build files
#-----------------------------------------------------------------------------#
scores <- mapply(score_cell_match, files, file_tags, MoreArgs = list(data_path))
scores <- do.call("rbind", scores)
file_name <- paste0(output,"benchmarking_scores.csv")
write.csv(scores,file = file_name)

