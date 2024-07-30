#-----------------------------------------------------------------------------#
# Loading Libraries and setting seed
#-----------------------------------------------------------------------------#
library(dplyr)
library(future)
library(ggplot2)
library(patchwork)
library(tidyr)
library(RColorBrewer)
library(ggnewscale)
library(ggpubr)
library(viridis)
library(ggtext)
set.seed(1547)

#-----------------------------------------------------------------------------#
# Output set up 
#-----------------------------------------------------------------------------#
input <- "/Users/martinp4/Documents/Cedars/Vesalius/Scenes/RAZA/"
output_plots <- "/Users/martinp4/Documents/Cedars/Vesalius/Scenes/RAZA/"


score_matrix <- list.files(input, pattern = "score_matrix_balence.csv", full.names = TRUE)
tag_list <- list.files(input, pattern = "score_matrix_balence.csv", full.names = FALSE)
status <- read.csv(paste0(input, "MBTMEIMCPublic/IMCClinical.csv"))

#-----------------------------------------------------------------------------#
# Unified clustering using cost only for clusters
# Not sure how I want to order thing here
#-----------------------------------------------------------------------------#
cost <- read.csv(grep("cost", score_matrix,value = TRUE),
    header = TRUE,
    skip = 1)

cost$from <- as.factor(cost$from)
cost$to <- as.factor(cost$to)


g <- ggplot(cost, aes(x = from, y = to)) +
        geom_tile(data = cost, aes(fill = score)) +
        scale_fill_gradientn(colors = rev(brewer.pal(11, "Spectral"))) +
        labs(title = "Cost", x = "Query", y = "Reference", fill = "Score") +
        theme_minimal() + 
        theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),
            legend.title = element_text(size = 15),
            axis.text = element_text(size = 8),
            axis.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            plot.title  = element_text(size = 18)) +
        guides(colour = guide_legend(
        override.aes = list(linewidth = 5)))

#-----------------------------------------------------------------------------#
# Clinical matching - cost
#-----------------------------------------------------------------------------#
cost <- read.csv(grep("cost", score_matrix,value = TRUE),
    header = TRUE,
    skip = 1)

min_cost <- split(cost, cost$from)
min_cost <- lapply(min_cost, function(x){
      return(x[x$score == min(x$score),])
    })
min_cost <- do.call("rbind", min_cost)
min_cost$matched <- paste0(min_cost$from,"_",min_cost$to)


min_cost$ERStatus <- status$ERStatus[match(min_cost$from, status$metabric_id)] ==
    status$ERStatus[match(min_cost$to, status$metabric_id)]
min_cost$Grade <- status$Grade[match(min_cost$from, status$metabric_id)] ==
    status$Grade[match(min_cost$to, status$metabric_id)]
min_cost$ERBB2_pos <- status$ERBB2_pos[match(min_cost$from, status$metabric_id)] ==
    status$ERBB2_pos[match(min_cost$to, status$metabric_id)]
min_cost$DeathBreast <- status$DeathBreast[match(min_cost$from, status$metabric_id)] ==
    status$DeathBreast[match(min_cost$to, status$metabric_id)]

min_cost$yearsToStatus <- abs(status$yearsToStatus[match(min_cost$from, status$metabric_id)] -
    status$yearsToStatus[match(min_cost$to, status$metabric_id)])

min_cost$match_cost <- apply(min_cost[,5:8],1,sum, na.rm =T) / sum(!apply(min_cost[,5:8],1,is.na))

#-----------------------------------------------------------------------------#
# Clinical matching - Feature
#-----------------------------------------------------------------------------#
feature <- read.csv(grep("feature", score_matrix,value = TRUE),
    header = TRUE,
    skip = 1)


feature$matched <- paste0(feature$from,"_",feature$to)
min_feature <- feature[match(min_cost$matched,feature$matched), ]
min_feature$ERStatus <- status$ERStatus[match(min_feature$from, status$metabric_id)] ==
    status$ERStatus[match(min_feature$to, status$metabric_id)]
min_feature$Grade <- status$Grade[match(min_feature$from, status$metabric_id)] ==
    status$Grade[match(min_cost$to, status$metabric_id)]
min_feature$ERBB2_pos <- status$ERBB2_pos[match(min_feature$from, status$metabric_id)] ==
    status$ERBB2_pos[match(min_feature$to, status$metabric_id)]
min_feature$DeathBreast <- status$DeathBreast[match(min_feature$from, status$metabric_id)] ==
    status$DeathBreast[match(min_feature$to, status$metabric_id)]

min_feature$yearsToStatus <- abs(status$yearsToStatus[match(min_feature$from, status$metabric_id)] -
    status$yearsToStatus[match(min_feature$to, status$metabric_id)])
min_feature$match_feature <- apply(min_feature[,5:8],1,sum, na.rm=T) / 
    apply(min_feature[,5:8],1,function(r){return(sum(!is.na(r)))})
#-----------------------------------------------------------------------------#
# Clinical matching - niche
#-----------------------------------------------------------------------------#
niche <- read.csv(grep("niche", score_matrix,value = TRUE),
    header = TRUE,
    skip = 1)

niche$matched <- paste0(niche$from,"_",niche$to)
min_niche <- niche[match(min_cost$matched,niche$matched), ]
min_niche$ERStatus <- status$ERStatus[match(min_niche$from, status$metabric_id)] ==
    status$ERStatus[match(min_niche$to, status$metabric_id)]
min_niche$Grade <- status$Grade[match(min_niche$from, status$metabric_id)] ==
    status$Grade[match(min_niche$to, status$metabric_id)]
min_niche$ERBB2_pos <- status$ERBB2_pos[match(min_niche$from, status$metabric_id)] ==
    status$ERBB2_pos[match(min_niche$to, status$metabric_id)]
min_niche$DeathBreast <- status$DeathBreast[match(min_niche$from, status$metabric_id)] ==
    status$DeathBreast[match(min_niche$to, status$metabric_id)]

min_niche$yearsToStatus <- abs(status$yearsToStatus[match(min_niche$from, status$metabric_id)] -
    status$yearsToStatus[match(min_niche$to, status$metabric_id)])
min_niche$match_niche <- apply(min_niche[,5:8],1,sum, na.rm=T) / 
    apply(min_niche[,5:8],1,function(r){return(sum(!is.na(r)))})
#-----------------------------------------------------------------------------#
# Clinical matching - territory
#-----------------------------------------------------------------------------#
territory <- read.csv(grep("territory", score_matrix,value = TRUE),
    header = TRUE,
    skip = 1)


territory$matched <- paste0(territory$from,"_",territory$to)
min_territory <- territory[match(min_cost$matched,territory$matched), ]
min_territory$ERStatus <- status$ERStatus[match(min_territory$from, status$metabric_id)] ==
    status$ERStatus[match(min_territory$to, status$metabric_id)]
min_territory$Grade <- status$Grade[match(min_territory$from, status$metabric_id)] ==
    status$Grade[match(min_territory$to, status$metabric_id)]
min_territory$ERBB2_pos <- status$ERBB2_pos[match(min_territory$from, status$metabric_id)] ==
    status$ERBB2_pos[match(min_territory$to, status$metabric_id)]
min_territory$DeathBreast <- status$DeathBreast[match(min_territory$from, status$metabric_id)] ==
    status$DeathBreast[match(min_territory$to, status$metabric_id)]

min_territory$yearsToStatus <- abs(status$yearsToStatus[match(min_territory$from, status$metabric_id)] -
    status$yearsToStatus[match(min_territory$to, status$metabric_id)])
min_territory$match_territory <- apply(min_territory[,5:8],1,sum, na.rm=T) / 
    apply(min_territory[,5:8],1,function(r){return(sum(!is.na(r)))})
#-----------------------------------------------------------------------------#
# Clinical matching - composition
#-----------------------------------------------------------------------------#
composition <- read.csv(grep("composition", score_matrix,value = TRUE),
    header = TRUE,
    skip = 1)

composition$matched <- paste0(composition$from,"_",composition$to)
min_composition <- composition[match(min_cost$matched,composition$matched), ]
min_composition$ERStatus <- status$ERStatus[match(min_composition$from, status$metabric_id)] ==
    status$ERStatus[match(min_composition$to, status$metabric_id)]
min_composition$Grade <- status$Grade[match(min_composition$from, status$metabric_id)] ==
    status$Grade[match(min_composition$to, status$metabric_id)]
min_composition$ERBB2_pos <- status$ERBB2_pos[match(min_composition$from, status$metabric_id)] ==
    status$ERBB2_pos[match(min_composition$to, status$metabric_id)]
min_composition$DeathBreast <- status$DeathBreast[match(min_composition$from, status$metabric_id)] ==
    status$DeathBreast[match(min_composition$to, status$metabric_id)]

min_composition$yearsToStatus <- abs(status$yearsToStatus[match(min_composition$from, status$metabric_id)] -
    status$yearsToStatus[match(min_composition$to, status$metabric_id)])

min_composition$match_composition <- apply(min_composition[,5:8],1,sum, na.rm=T) / 
    apply(min_composition[,5:8],1,function(r){return(sum(!is.na(r)))})


#-----------------------------------------------------------------------------#
# Clinical matching - cell label
#-----------------------------------------------------------------------------#
cell_label <- read.csv(grep("cell_label", score_matrix,value = TRUE),
    header = TRUE,
    skip = 1)


cell_label$matched <- paste0(cell_label$from,"_",cell_label$to)
min_cell_label <- cell_label[match(min_cost$matched,cell_label$matched), ]
min_cell_label$ERStatus <- status$ERStatus[match(min_cell_label$from, status$metabric_id)] ==
    status$ERStatus[match(min_cell_label$to, status$metabric_id)]
min_cell_label$Grade <- status$Grade[match(min_cell_label$from, status$metabric_id)] ==
    status$Grade[match(min_cell_label$to, status$metabric_id)]
min_cell_label$ERBB2_pos <- status$ERBB2_pos[match(min_cell_label$from, status$metabric_id)] ==
    status$ERBB2_pos[match(min_cell_label$to, status$metabric_id)]
min_cell_label$DeathBreast <- status$DeathBreast[match(min_cell_label$from, status$metabric_id)] ==
    status$DeathBreast[match(min_cell_label$to, status$metabric_id)]

min_cell_label$yearsToStatus <- abs(status$yearsToStatus[match(min_cell_label$from, status$metabric_id)] -
    status$yearsToStatus[match(min_cell_label$to, status$metabric_id)])

min_cell_label$match_cell_label <- apply(min_cell_label[,5:8],1,sum, na.rm=T) / 
    apply(min_cell_label[,5:8],1,function(r){return(sum(!is.na(r)))})

#-----------------------------------------------------------------------------#
# Combined status
#-----------------------------------------------------------------------------#
combined <- data.frame("matched" = min_cost$matched,
    "cost" = min_cost$score,
    "feature" = min_feature$score,
    "niche" = min_niche$score,
    "territory" = min_territory$score,
    "composition" = min_composition$score,
    "cell_label" = min_cell_label$score,
    "cost_match" = min_cost$match_cost,
    "feature_match" = min_feature$match_feature,
    "niche_match" = min_niche$match_niche,
    "territory_match" = min_territory$match_territory,
    "composition_match" = min_composition$match_composition,
    "cell_label_match" = min_cell_label$match_cell_label)
combined$cost_match <- combined$cost_match * 100
score_long <- combined %>%
    select(c("matched","cost","feature","niche","territory","composition","cell_label","cost_match"))
ord <- order(score_long$cost_match)
score_long <- score_long %>% pivot_longer(cols = c("cost","feature","territory" ,"niche","composition","cell_label","cost_match"))
score_long$name <- as.factor(score_long$name)
score_long$name <- factor(score_long$name, levels = c("cost","feature","niche","territory","composition","cell_label","cost_match"))
score_long$matched <- as.factor(score_long$matched)
score_long$matched <- factor(score_long$matched, levels = levels(score_long$matched)[ord])

g <- ggplot(score_long, aes(x = matched, y = name)) +
    geom_tile(data = filter(score_long, name == "cost"), aes(fill = value)) +
    scale_fill_gradientn(colors = rev(brewer.pal(9, "Blues"))) +
    labs(fill = "Cost") +
    new_scale("fill") +
    geom_tile(data = filter(score_long, name %in% c("feature","niche","territory","composition","cell_label")), aes(fill = value)) +
    scale_fill_gradientn(colors = brewer.pal(9, "Oranges")) +
    labs(title = "", y = "", x = "Matched Data sets", fill = "Score") +
    new_scale("fill") +
    geom_tile(data = filter(score_long, name == "cost_match"), aes(fill = value)) +
    labs(fill = "Clinical Out.") +
    scale_fill_gradientn(colors = brewer.pal(9, "Greens")) +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 9),
            legend.title = element_text(size = 15),
            axis.text.y = element_text(size = 18),
            axis.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            plot.title  = element_text(size = 18)) +
    guides(colour = guide_legend(
        override.aes = list(linewidth = 5)))
file_out <- paste0(output_plots,"best_matching_scores_full.pdf")
pdf(file_out, height = 6, width=18)
print(g)
dev.off()



#-----------------------------------------------------------------------------#
# Clinical outcomes for min cost 
#-----------------------------------------------------------------------------#
clinical_cost <- min_cost %>%
    select(c("matched","ERStatus","Grade","ERBB2_pos","DeathBreast","yearsToStatus"))
n_pred_values <- apply(clinical_cost[, 2:5], 2, sum, na.rm = TRUE)
n_preds <- c(order(apply(clinical_cost[, 2:5], 2, sum, na.rm = TRUE), decreasing = TRUE),5)
n_preds <-c("ERStatus","Grade","ERBB2_pos","DeathBreast","yearsToStatus")[n_preds]
n_pred_values <- n_pred_values[n_preds]
ord <- order(clinical_cost[, n_preds[1]])
ord <- order(apply(clinical_cost[,2:5],1,sum), decreasing = FALSE)
clinical_cost <- clinical_cost %>%
    pivot_longer(cols = c("ERStatus","Grade","ERBB2_pos","DeathBreast","yearsToStatus"))
clinical_cost$name <- as.factor(clinical_cost$name)
clinical_cost$name <- factor(clinical_cost$name, levels = n_preds)
clinical_cost$matched <- as.factor(clinical_cost$matched)
clinical_cost$matched <- factor(clinical_cost$matched, levels = levels(clinical_cost$matched)[ord])

text_format <- paste0(names(n_pred_values)," = " ,n_pred_values)[-5]
text_format <- paste0(text_format, collapse = "\n")
#text_format <- data.frame("x" = 6, "y" = 40, labels = text_format)

g <- ggplot(clinical_cost, aes(x = matched, y = name)) +
    geom_tile(data = filter(clinical_cost, name!="yearsToStatus"), aes(fill = as.factor(value)), color = "white") +
    scale_fill_manual(values =  brewer.pal(3, "Greens")[c(1,3)]) +
    labs(title = "", y = "", x = "Matched Data sets", fill = "Matched") +
    new_scale("fill") +
    geom_tile(data = filter(clinical_cost, name == "yearsToStatus"), aes(fill = value)) +
    labs(fill = "Delta Years To Status") +
    labs(tag = text_format) +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 9),
            legend.title = element_text(size = 15),
            axis.text.y = element_text(size = 18),
            axis.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            plot.title  = element_text(size = 18),
            plot.tag.position = c(0.92, 0.2)) +
    guides(colour = guide_legend(
        override.aes = list(linewidth = 5))) 

file_out <- paste0(output_plots,"clinical_outcomes.pdf")
pdf(file_out, height = 6, width=18)
print(g)
dev.off()

#-----------------------------------------------------------------------------#
# ER Pos
#-----------------------------------------------------------------------------#
df <- cost
df$erpos_to <- status$ERStatus[match(as.character(df$to), status$metabric_id)]
df$erpos_from <- status$ERStatus[match(as.character(df$from), status$metabric_id)]

df <- df  %>% filter(erpos_from == "pos" & erpos_to == "pos")
# get positive clustering
heatmap_data <- df %>%
        select(from, to, score) %>%
        spread(to, score)
fr <- as.character(heatmap_data$from)
to <- colnames(heatmap_data)[-1]
heatmap_matrix <- as.matrix(heatmap_data[, -1])
colnames(heatmap_matrix) <- to
rownames(heatmap_matrix) <- fr


#heatmap_matrix <- vesalius:::overlap_distance_matrix(heatmap_matrix, 20, FALSE)
#
d <- as.dist(heatmap_matrix)
d <- (d - min(d)) / (max(d) - min(d))

row_order_cluster <- hclust(as.dist(d))
row_order_cut <- cutree(row_order_cluster, h = 0.5)
row_order <- names(row_order_cut)[row_order_cluster$order]
col_order_cluster <- hclust(as.dist(t(d)))
col_order_cut <- cutree(col_order_cluster, h = 0.5)
col_order <- names(col_order_cut)[col_order_cluster$order]
df <- reshape2::melt(as.matrix(heatmap_matrix), c("row", "col"), value.name = "scores")
df$row <- factor(df$row, levels = c(row_order))
df$col <- factor(df$col, levels = c(col_order))
df$scores <- df$scores
df$trial[match(df$row, names(row_order_cut))] <- row_order_cut[match(df$row, names(row_order_cut))]

sub_df <- df
sub_df$clust_2 <- sub_df$trial[match(sub_df$col, sub_df$row)]
sub_df <- sub_df[sub_df$trial == sub_df$clust_2,]
sub_df$trial <- as.factor(sub_df$trial)
cols <- vesalius:::create_palette(sub_df, FALSE)


g <- ggplot(df, aes(x = row, y = col)) +
        geom_tile(data = df, aes(fill = scores)) +
        #scale_fill_gradientn(colors = brewer.pal(9, "RdBu")) +
        scale_fill_viridis_c(option = "magma")+
        labs(title = "Cost Score", x = "", y = "", fill = "Score") +
        new_scale("color") +
        geom_tile(data = sub_df, mapping = aes(col = trial), alpha = 0,linewidth = 0.75) +
        labs(color = "Cluster") +
        scale_color_manual(values = cols) +
        theme_minimal() + 
        theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),
            legend.title = element_text(size = 15),
            axis.text = element_text(size = 8),
            axis.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            plot.title  = element_text(size = 18)) +
        guides(colour = guide_legend(
        override.aes = list(linewidth = 5)))
   

file_out <- paste0(output_plots,"ERStatus_pos.pdf")
pdf(file_out, width = 12, height = 10)
print(g)
dev.off()

feat <- feature
feat$erpos_to <- status$ERStatus[match(as.character(feat$to), status$metabric_id)]
feat$erpos_from <- status$ERStatus[match(as.character(feat$from), status$metabric_id)]

feat <- feat  %>% filter(erpos_from == "pos" & erpos_to == "pos")
feat$row <- factor(feat$from, levels = c(row_order))
feat$col <- factor(feat$to, levels = c(col_order))
feat$trial[match(feat$from, names(row_order_cut))] <- row_order_cut[match(feat$from, names(row_order_cut))]

g <- ggplot(feat, aes(x = row, y = col)) +
        geom_tile(data = feat, aes(fill = score)) +
        scale_fill_viridis_c(option = "mako")+
        labs(title = "Feature Score", x = "", y = "", fill = "Score") +
        new_scale("color") +
        geom_tile(data = sub_df, mapping = aes(col = trial), alpha = 0,linewidth = 0.75) +
        labs(color = "Cluster") +
        scale_color_manual(values = cols) +
        theme_minimal() + 
        theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),
            legend.title = element_text(size = 15),
            axis.text = element_text(size = 8),
            axis.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            plot.title  = element_text(size = 18)) +
        guides(colour = guide_legend(
        override.aes = list(linewidth = 5)))
   

file_out <- paste0(output_plots,"ERStatus_pos_feat.pdf")
pdf(file_out, width = 12, height = 10)
print(g)
dev.off()


ni <- niche
ni$erpos_to <- status$ERStatus[match(as.character(ni$to), status$metabric_id)]
ni$erpos_from <- status$ERStatus[match(as.character(ni$from), status$metabric_id)]

ni <- ni  %>% filter(erpos_from == "pos" & erpos_to == "pos")
ni$row <- factor(ni$from, levels = c(row_order))
ni$col <- factor(ni$to, levels = c(col_order))
ni$trial[match(ni$from, names(row_order_cut))] <- row_order_cut[match(ni$from, names(row_order_cut))]

g <- ggplot(ni, aes(x = row, y = col)) +
        geom_tile(data = ni, aes(fill = score)) +
        scale_fill_viridis_c(option = "mako")+
        labs(title = "Niche Score", x = "", y = "", fill = "Score") +
        new_scale("color") +
        geom_tile(data = sub_df, mapping = aes(col = trial), alpha = 0,linewidth = 0.75) +
        labs(color = "Cluster") +
        scale_color_manual(values = cols) +
        theme_minimal() + 
        theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),
            legend.title = element_text(size = 15),
            axis.text = element_text(size = 8),
            axis.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            plot.title  = element_text(size = 18)) +
        guides(colour = guide_legend(
        override.aes = list(linewidth = 5)))
   

file_out <- paste0(output_plots,"ERStatus_pos_niche.pdf")
pdf(file_out, width = 12, height = 10)
print(g)
dev.off()


comp <- composition
comp$erpos_to <- status$ERStatus[match(as.character(comp$to), status$metabric_id)]
comp$erpos_from <- status$ERStatus[match(as.character(comp$from), status$metabric_id)]

comp <- comp  %>% filter(erpos_from == "pos" & erpos_to == "pos")
comp$row <- factor(comp$from, levels = c(row_order))
comp$col <- factor(comp$to, levels = c(col_order))
comp$trial[match(comp$from, names(row_order_cut))] <- row_order_cut[match(comp$from, names(row_order_cut))]

g <- ggplot(comp, aes(x = row, y = col)) +
        geom_tile(data = comp, aes(fill = score)) +
        scale_fill_viridis_c(option = "mako")+
        labs(title = "Composition Score", x = "", y = "", fill = "Score") +
        new_scale("color") +
        geom_tile(data = sub_df, mapping = aes(col = trial), alpha = 0,linewidth = 0.75) +
        labs(color = "Cluster") +
        scale_color_manual(values = cols) +
        theme_minimal() + 
        theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),
            legend.title = element_text(size = 15),
            axis.text = element_text(size = 8),
            axis.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            plot.title  = element_text(size = 18)) +
        guides(colour = guide_legend(
        override.aes = list(linewidth = 5)))
   

file_out <- paste0(output_plots,"ERStatus_pos_comp.pdf")
pdf(file_out, width = 12, height = 10)
print(g)
dev.off()


ter <- territory
ter$erpos_to <- status$ERStatus[match(as.character(ter$to), status$metabric_id)]
ter$erpos_from <- status$ERStatus[match(as.character(ter$from), status$metabric_id)]

ter <- ter  %>% filter(erpos_from == "pos" & erpos_to == "pos")
ter$row <- factor(ter$from, levels = c(row_order))
ter$col <- factor(ter$to, levels = c(col_order))
ter$trial[match(ter$from, names(row_order_cut))] <- row_order_cut[match(ter$from, names(row_order_cut))]

g <- ggplot(ter, aes(x = row, y = col)) +
        geom_tile(data = ter, aes(fill = score)) +
        scale_fill_viridis_c(option = "mako")+
        labs(title = "Territory Score", x = "", y = "", fill = "Score") +
        new_scale("color") +
        geom_tile(data = sub_df, mapping = aes(col = trial), alpha = 0,linewidth = 0.75) +
        labs(color = "Cluster") +
        scale_color_manual(values = cols) +
        theme_minimal() + 
        theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),
            legend.title = element_text(size = 15),
            axis.text = element_text(size = 8),
            axis.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            plot.title  = element_text(size = 18)) +
        guides(colour = guide_legend(
        override.aes = list(linewidth = 5)))
   

file_out <- paste0(output_plots,"ERStatus_pos_ter.pdf")
pdf(file_out, width = 12, height = 10)
print(g)
dev.off()


cell <- cell_label
cell$erpos_to <- status$ERStatus[match(as.character(cell$to), status$metabric_id)]
cell$erpos_from <- status$ERStatus[match(as.character(cell$from), status$metabric_id)]

cell <- cell  %>% filter(erpos_from == "pos" & erpos_to == "pos")
cell$row <- factor(cell$from, levels = c(row_order))
cell$col <- factor(cell$to, levels = c(col_order))
cell$trial[match(cell$from, names(row_order_cut))] <- row_order_cut[match(cell$from, names(row_order_cut))]

g <- ggplot(cell, aes(x = row, y = col)) +
        geom_tile(data = cell, aes(fill = score)) +
        scale_fill_viridis_c(option = "mako")+
        labs(title = "Cell Label Score", x = "", y = "", fill = "Score") +
        new_scale("color") +
        geom_tile(data = sub_df, mapping = aes(col = trial), alpha = 0,linewidth = 0.75) +
        labs(color = "Cluster") +
        scale_color_manual(values = cols) +
        theme_minimal() + 
        theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),
            legend.title = element_text(size = 15),
            axis.text = element_text(size = 8),
            axis.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            plot.title  = element_text(size = 18)) +
        guides(colour = guide_legend(
        override.aes = list(linewidth = 5)))
   

file_out <- paste0(output_plots,"ERStatus_pos_cell_lab.pdf")
pdf(file_out, width = 12, height = 10)
print(g)
dev.off()


co <- cost
co$erpos_to <- status$ERStatus[match(as.character(co$to), status$metabric_id)]
co$erpos_from <- status$ERStatus[match(as.character(co$from), status$metabric_id)]

co <- co  %>% filter(erpos_from == "pos" & erpos_to == "pos")
co$row <- factor(co$from, levels = c(row_order))
co$col <- factor(co$to, levels = c(col_order))
co$trial[match(co$from, names(row_order_cut))] <- row_order_cut[match(co$from, names(row_order_cut))]

g <- ggplot(co, aes(x = row, y = col)) +
        geom_tile(data = co, aes(fill = score)) +
        scale_fill_viridis_c(option = "mako")+
        labs(title = "Cost Score", x = "", y = "", fill = "Score") +
        new_scale("color") +
        geom_tile(data = sub_df, mapping = aes(col = trial), alpha = 0,linewidth = 0.75) +
        labs(color = "Cluster") +
        scale_color_manual(values = cols) +
        theme_minimal() + 
        theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),
            legend.title = element_text(size = 15),
            axis.text = element_text(size = 8),
            axis.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            plot.title  = element_text(size = 18)) +
        guides(colour = guide_legend(
        override.aes = list(linewidth = 5)))
   

file_out <- paste0(output_plots,"ERStatus_pos_cost.pdf")
pdf(file_out, width = 12, height = 10)
print(g)
dev.off()


    
#-----------------------------------------------------------------------------#
# ER NEG
#-----------------------------------------------------------------------------#

df <- cost
df$erpos_to <- status$ERStatus[match(as.character(df$to), status$metabric_id)]
df$erpos_from <- status$ERStatus[match(as.character(df$from), status$metabric_id)]

df <- df  %>% filter(erpos_from == "neg" & erpos_to == "neg")
# get positive clustering
heatmap_data <- df %>%
        select(from, to, score) %>%
        spread(to, score)
fr <- as.character(heatmap_data$from)
to <- colnames(heatmap_data)[-1]
heatmap_matrix <- as.matrix(heatmap_data[, -1])
colnames(heatmap_matrix) <- to
rownames(heatmap_matrix) <- fr


#heatmap_matrix <- vesalius:::overlap_distance_matrix(heatmap_matrix, 20, FALSE)
#
d <- as.dist(heatmap_matrix)
d <- (d - min(d)) / (max(d) - min(d))

row_order_cluster <- hclust(as.dist(d))
row_order_cut <- cutree(row_order_cluster, h = 0.35)
row_order <- names(row_order_cut)[row_order_cluster$order]
col_order_cluster <- hclust(as.dist(t(d)))
col_order_cut <- cutree(col_order_cluster, h = 0.35)
col_order <- names(col_order_cut)[col_order_cluster$order]
df <- reshape2::melt(as.matrix(heatmap_matrix), c("row", "col"), value.name = "scores")
df$row <- factor(df$row, levels = c(row_order))
df$col <- factor(df$col, levels = c(col_order))
df$scores <- df$scores
df$trial[match(df$row, names(row_order_cut))] <- row_order_cut[match(df$row, names(row_order_cut))]

sub_df <- df
sub_df$clust_2 <- sub_df$trial[match(sub_df$col, sub_df$row)]
sub_df <- sub_df[sub_df$trial == sub_df$clust_2,]
sub_df$trial <- as.factor(sub_df$trial)
cols <- vesalius:::create_palette(sub_df, FALSE)


g <- ggplot(df, aes(x = row, y = col)) +
        geom_tile(data = df, aes(fill = scores)) +
        #scale_fill_gradientn(colors = brewer.pal(9, "RdBu")) +
        scale_fill_viridis_c(option = "magma")+
        labs(title = "Cost Score", x = "", y = "", fill = "Score") +
        new_scale("color") +
        geom_tile(data = sub_df, mapping = aes(col = trial), alpha = 0,linewidth = 0.75) +
        labs(color = "Cluster") +
        scale_color_manual(values = cols) +
        theme_minimal() + 
        theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),
            legend.title = element_text(size = 15),
            axis.text = element_text(size = 8),
            axis.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            plot.title  = element_text(size = 18)) +
        guides(colour = guide_legend(
        override.aes = list(linewidth = 5)))
   

file_out <- paste0(output_plots,"ERStatus_neg.pdf")
pdf(file_out, width = 12, height = 10)
print(g)
dev.off()

feat <- feature
feat$erpos_to <- status$ERStatus[match(as.character(feat$to), status$metabric_id)]
feat$erpos_from <- status$ERStatus[match(as.character(feat$from), status$metabric_id)]

feat <- feat  %>% filter(erpos_from == "neg" & erpos_to == "neg")
feat$row <- factor(feat$from, levels = c(row_order))
feat$col <- factor(feat$to, levels = c(col_order))
feat$trial[match(feat$from, names(row_order_cut))] <- row_order_cut[match(feat$from, names(row_order_cut))]

g <- ggplot(feat, aes(x = row, y = col)) +
        geom_tile(data = feat, aes(fill = score)) +
        scale_fill_viridis_c(option = "mako")+
        labs(title = "Feature Score", x = "", y = "", fill = "Score") +
        new_scale("color") +
        geom_tile(data = sub_df, mapping = aes(col = trial), alpha = 0,linewidth = 0.75) +
        labs(color = "Cluster") +
        scale_color_manual(values = cols) +
        theme_minimal() + 
        theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),
            legend.title = element_text(size = 15),
            axis.text = element_text(size = 8),
            axis.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            plot.title  = element_text(size = 18)) +
        guides(colour = guide_legend(
        override.aes = list(linewidth = 5)))
   

file_out <- paste0(output_plots,"ERStatus_neg_feat.pdf")
pdf(file_out, width = 12, height = 10)
print(g)
dev.off()


ni <- niche
ni$erpos_to <- status$ERStatus[match(as.character(ni$to), status$metabric_id)]
ni$erpos_from <- status$ERStatus[match(as.character(ni$from), status$metabric_id)]

ni <- ni  %>% filter(erpos_from == "neg" & erpos_to == "neg")
ni$row <- factor(ni$from, levels = c(row_order))
ni$col <- factor(ni$to, levels = c(col_order))
ni$trial[match(ni$from, names(row_order_cut))] <- row_order_cut[match(ni$from, names(row_order_cut))]

g <- ggplot(ni, aes(x = row, y = col)) +
        geom_tile(data = ni, aes(fill = score)) +
        scale_fill_viridis_c(option = "mako")+
        labs(title = "Niche Score", x = "", y = "", fill = "Score") +
        new_scale("color") +
        geom_tile(data = sub_df, mapping = aes(col = trial), alpha = 0,linewidth = 0.75) +
        labs(color = "Cluster") +
        scale_color_manual(values = cols) +
        theme_minimal() + 
        theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),
            legend.title = element_text(size = 15),
            axis.text = element_text(size = 8),
            axis.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            plot.title  = element_text(size = 18)) +
        guides(colour = guide_legend(
        override.aes = list(linewidth = 5)))
   

file_out <- paste0(output_plots,"ERStatus_neg_niche.pdf")
pdf(file_out, width = 12, height = 10)
print(g)
dev.off()


comp <- composition
comp$erpos_to <- status$ERStatus[match(as.character(comp$to), status$metabric_id)]
comp$erpos_from <- status$ERStatus[match(as.character(comp$from), status$metabric_id)]

comp <- comp  %>% filter(erpos_from == "neg" & erpos_to == "neg")
comp$row <- factor(comp$from, levels = c(row_order))
comp$col <- factor(comp$to, levels = c(col_order))
comp$trial[match(comp$from, names(row_order_cut))] <- row_order_cut[match(comp$from, names(row_order_cut))]

g <- ggplot(comp, aes(x = row, y = col)) +
        geom_tile(data = comp, aes(fill = score)) +
        scale_fill_viridis_c(option = "mako")+
        labs(title = "Composition Score", x = "", y = "", fill = "Score") +
        new_scale("color") +
        geom_tile(data = sub_df, mapping = aes(col = trial), alpha = 0,linewidth = 0.75) +
        labs(color = "Cluster") +
        scale_color_manual(values = cols) +
        theme_minimal() + 
        theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),
            legend.title = element_text(size = 15),
            axis.text = element_text(size = 8),
            axis.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            plot.title  = element_text(size = 18)) +
        guides(colour = guide_legend(
        override.aes = list(linewidth = 5)))
   

file_out <- paste0(output_plots,"ERStatus_neg_comp.pdf")
pdf(file_out, width = 12, height = 10)
print(g)
dev.off()


ter <- territory
ter$erpos_to <- status$ERStatus[match(as.character(ter$to), status$metabric_id)]
ter$erpos_from <- status$ERStatus[match(as.character(ter$from), status$metabric_id)]

ter <- ter  %>% filter(erpos_from == "neg" & erpos_to == "neg")
ter$row <- factor(ter$from, levels = c(row_order))
ter$col <- factor(ter$to, levels = c(col_order))
ter$trial[match(ter$from, names(row_order_cut))] <- row_order_cut[match(ter$from, names(row_order_cut))]

g <- ggplot(ter, aes(x = row, y = col)) +
        geom_tile(data = ter, aes(fill = score)) +
        scale_fill_viridis_c(option = "mako")+
        labs(title = "Territory Score", x = "", y = "", fill = "Score") +
        new_scale("color") +
        geom_tile(data = sub_df, mapping = aes(col = trial), alpha = 0,linewidth = 0.75) +
        labs(color = "Cluster") +
        scale_color_manual(values = cols) +
        theme_minimal() + 
        theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),
            legend.title = element_text(size = 15),
            axis.text = element_text(size = 8),
            axis.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            plot.title  = element_text(size = 18)) +
        guides(colour = guide_legend(
        override.aes = list(linewidth = 5)))
   

file_out <- paste0(output_plots,"ERStatus_neg_ter.pdf")
pdf(file_out, width = 12, height = 10)
print(g)
dev.off()


cell <- cell_label
cell$erpos_to <- status$ERStatus[match(as.character(cell$to), status$metabric_id)]
cell$erpos_from <- status$ERStatus[match(as.character(cell$from), status$metabric_id)]

cell <- cell  %>% filter(erpos_from == "neg" & erpos_to == "neg")
cell$row <- factor(cell$from, levels = c(row_order))
cell$col <- factor(cell$to, levels = c(col_order))
cell$trial[match(cell$from, names(row_order_cut))] <- row_order_cut[match(cell$from, names(row_order_cut))]

g <- ggplot(cell, aes(x = row, y = col)) +
        geom_tile(data = cell, aes(fill = score)) +
        scale_fill_viridis_c(option = "mako")+
        labs(title = "Cell Label Score", x = "", y = "", fill = "Score") +
        new_scale("color") +
        geom_tile(data = sub_df, mapping = aes(col = trial), alpha = 0,linewidth = 0.75) +
        labs(color = "Cluster") +
        scale_color_manual(values = cols) +
        theme_minimal() + 
        theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),
            legend.title = element_text(size = 15),
            axis.text = element_text(size = 8),
            axis.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            plot.title  = element_text(size = 18)) +
        guides(colour = guide_legend(
        override.aes = list(linewidth = 5)))
   

file_out <- paste0(output_plots,"ERStatus_neg_cell_lab.pdf")
pdf(file_out, width = 12, height = 10)
print(g)
dev.off()


co <- cost
co$erpos_to <- status$ERStatus[match(as.character(co$to), status$metabric_id)]
co$erpos_from <- status$ERStatus[match(as.character(co$from), status$metabric_id)]

co <- co  %>% filter(erpos_from == "neg" & erpos_to == "neg")
co$row <- factor(co$from, levels = c(row_order))
co$col <- factor(co$to, levels = c(col_order))
co$trial[match(co$from, names(row_order_cut))] <- row_order_cut[match(co$from, names(row_order_cut))]

g <- ggplot(co, aes(x = row, y = col)) +
        geom_tile(data = co, aes(fill = score)) +
        scale_fill_viridis_c(option = "mako")+
        labs(title = "Cost Score", x = "", y = "", fill = "Score") +
        new_scale("color") +
        geom_tile(data = sub_df, mapping = aes(col = trial), alpha = 0,linewidth = 0.75) +
        labs(color = "Cluster") +
        scale_color_manual(values = cols) +
        theme_minimal() + 
        theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),
            legend.title = element_text(size = 15),
            axis.text = element_text(size = 8),
            axis.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            plot.title  = element_text(size = 18)) +
        guides(colour = guide_legend(
        override.aes = list(linewidth = 5)))
   

file_out <- paste0(output_plots,"ERStatus_neg_cost.pdf")
pdf(file_out, width = 12, height = 10)
print(g)
dev.off()



# #-----------------------------------------------------------------------------#
# # Showing status
# #-----------------------------------------------------------------------------#
# full_df <- do.call("rbind", full_df)
# full_df$metabric_id <- full_df$from
# status <- left_join(full_df, status, by = "metabric_id")

# er <- split(status, status$type)
# bar <- vector("list", length(er))
# for (i in seq_along(er)){
#     tmp <- er[[i]]
#     clust <- unique(tmp$clust)
#     tmp$status_tab <- 0
#     tmp$q <- 0
#     tmp$ord <- 0
#     for (j in seq_along(clust)){
#         tot <- sum(tmp$clust == clust[j])
#         pos <- which(tmp$clust == clust[j] & tmp$ERStatus == "pos")
#         tmp$status_tab[pos] <- (length(pos) / tot) * 100
#         tmp$q[pos[1]] <- 1 
#         neg <- which(tmp$clust == clust[j] & tmp$ERStatus == "neg")
#         tmp$status_tab[neg] <- (length(neg) / tot) * 100
#         tmp$q[neg[1]] <- 1
#         tmp$ord[tmp$clust == clust[j]]<- (length(pos) / tot) - (length(neg) / tot)
        
#     }
#     tmp <- tmp[tmp$q == 1,]
#     tmp$clust <- factor(tmp$clust, levels = unique(tmp$clust[order(tmp$ord, decreasing = TRUE)]))
#     bar[[i]] <- ggplot(tmp, aes(fill=ERStatus, y=status_tab, x=clust)) + 
#     geom_bar(position = "stack", stat = "identity") +
#     scale_fill_manual(values = c("#D55E00", "#0072B2")) + 
#     theme_bw() +
#     labs(x = "Cluster", y = "ER Status Distribution", title = unique(tmp$type)) +
#      theme(legend.title = element_text(size = 15),
#             axis.text = element_text(size = 12),
#             axis.title = element_text(size = 15),
#             legend.text = element_text(size = 15),
#             plot.title  = element_text(size = 18))
# }

# file_name <- paste0(output_plots, "cluster_distribution.pdf")
# pdf(file_name, width = 16, height = 12)
# ggarrange(plotlist = bar, ncol = 3, nrow=2)
# dev.off()

