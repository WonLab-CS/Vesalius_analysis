#-----------------------------------------------------------------------------#
# Loading Libraries and setting seed
#-----------------------------------------------------------------------------#
library(dplyr)
library(future)
library(ggplot2)
library(dplyr)
library(patchwork)
library(Matrix)
library(deldir)
library(imager)
library(imagerExtra)
library(vesalius)

set.seed(1547)


#-----------------------------------------------------------------------------#
# Set future global for multicore processing 
# Adding data paths 
#-----------------------------------------------------------------------------#
plan(multicore, workers = 4)
max_size <- 10000 * 1024^2
options(future.globals.maxSize = max_size)

input <-"/common/wonklab/RAZA/MBTMEIMCPublic/"

if (!dir.exists("/common/wonklab/RAZA/split_data/")) {
    dir.create("/common/wonklab/RAZA/split_data/")
}
output_data <- "/common/wonklab/RAZA/split_data/"


if (!dir.exists("/common/wonklab/RAZA/output_plots/")) {
    dir.create("/common/wonklab/RAZA/output_plots/")
}
output_plots <- "/common/wonklab/RAZA/output_plots/"
#-----------------------------------------------------------------------------#
# Load data
#-----------------------------------------------------------------------------#

dat <- read.csv(paste0(input,"SingleCells.csv"))

dat <- split(dat, dat$metabric_id)


dat <- dat[sapply(dat, nrow) > 1000]

cli <- read.csv(paste0(input, "IMCClinical.csv"))
cli <- cli[cli$metabric_id %in% names(dat), ]
status <- apply(cli[,c("ERStatus","ERBB2_pos","Grade","DeathBreast","yearsToStatus")],1,
    function(x){
        return(sum(is.na(x))==0)
    })
retain <- cli$metabric_id[status]


dat <- dat[names(dat) %in% retain]

cli <- cli[cli$metabric_id %in% retain,]


pos <- sample(which(cli$ERStatus == "pos"), size = 50, replace = FALSE)
neg <- sample(which(cli$ERStatus == "neg"), size = 50, replace = FALSE)
cli <- cli[c(pos,neg), ]

dat <- dat[names(dat) %in% cli$metabric_id]

#-----------------------------------------------------------------------------#
# out coord data 
#-----------------------------------------------------------------------------#
coord <- vector("list", length(dat))
counts <- vector("list", length(dat))
cells <- vector("list", length(dat))

for (i in seq_along(dat)) {
    tmp <- dat[[i]][, c("Location_Center_X", "Location_Center_Y")]
    tmp <- data.frame(
        "barcodes" = paste0(dat[[i]]$metabric_id,"_", seq_len(nrow(dat[[i]]))),
        tmp)
    colnames(tmp) <- c("barcodes", "x", "y")
    coord[[i]] <- tmp
    tmp <- dat[[i]][, seq(grep("Histone.H3", colnames(dat[[i]])),
        grep("c.Caspase3",colnames(dat[[i]])))]
    rownames(tmp) <- paste0(dat[[i]]$metabric_id, "_", seq_len(nrow(dat[[i]])))
    counts[[i]] <- t(tmp)
    tmp<- dat[[i]]$cellPhenotype
    names(tmp) <- paste0(dat[[i]]$metabric_id, "_", seq_len(nrow(dat[[i]])))
    cells[[i]] <- tmp
    ves <- build_vesalius_assay(coord[[i]],
        counts[[i]],
        assay = unique(dat[[i]]$metabric_id))
    ves <- add_cells(ves, cells[[i]])
    file_name <- paste0(output_data, names(dat)[i], "_balence.Rda")
    save(ves, file = file_name)
}



