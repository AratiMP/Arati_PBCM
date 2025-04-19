################################################################################
######################## Breast Cancer Subtypes - DIABLO #######################
########################## PLOT INFO EXTRACTION FILE ###########################
################################################################################

# Setwd
setwd('C:\\Users\\arati\\OneDrive - FCT NOVA\\PBCM\\plotsData') #FILL

############################### IMPORTS ########################################
install.packages("devtools")
install.packages("remotes")
remotes::install_github("msxakk89/dat")

library(mixOmics)
library(RColorBrewer)
library(ggtext)
library(ggplot2)
library(dplyr)
############################## CONSTANTS #######################################

# Model Name
MODEL_PATH = 'C:\\users\\arati\\OneDrive - FCT NOVA\\PBCM\\Model-Ari' 

################################## CODE ########################################

#### Get the model ####
model <- readRDS(MODEL_PATH)

#### Data for PCA-like plot ####

# See the distribution of patients and if the model managed to segregate the subtypes
plot_obj <- plotIndiv(model,
                      comp = c(1,2),
                      ind.names = FALSE, 
                      legend = TRUE,
                      pch = 19,
                      cex = 1,
                      style = 'lattice')

# Save the data
write.csv(plot_obj$df, 
          'Plot_Indiv_Data.csv', 
          row.names = T)

#### Data for ROC curves ####

# See the performance inside each dataset in distinguishing between subtypes
dfs_aurocs <- list()

for (dataset in names(model$keepX)) {
  
  plot_auroc_obj <- auroc(model, roc.block = dataset,
                          roc.comp = 2, print = FALSE)
  graph_name <- paste('graph.', dataset, sep = '')
  
  df <- plot_auroc_obj[[graph_name]]$comp2$data
  
  df$Dataset <- dataset
  
  dfs_aurocs[[dataset]] <- df
}

aurocs <- do.call(rbind, dfs_aurocs)

write.csv(aurocs, 
          'Plot_AUROC_Data.csv', 
          row.names = T)

#### Get the Features Selected ####

for (n in 1: 2) {
  feature_selection_info <- selectVar(model, comp = n)
  features_mRNA <- feature_selection_info$mRNA$name
  features_miRNA <- feature_selection_info$miRNA$name
  max.len = max(length(features_mRNA), length(features_miRNA))
  x = c(features_mRNA, rep(NA, max.len - length(features_mRNA)))
  y = c(features_miRNA, rep(NA, max.len - length(features_miRNA)))
  
  features_selected <- list(mRNA = x, miRNA = y)
  
  write.csv(features_selected, 
            file = paste('Features_Selected_Comp', n, '.csv', sep = ''), 
            row.names = FALSE)
}

#Export the loadings data

loadings_graph_comp1_mRNA <- plotLoadings(model, comp = 1, contrib = 'max', method = 'median')YT

write.csv(loadings_graph_comp1$mRNA, 'Loadings_Plot_mRNA_Comp1.csv', row.names = T)

loadings_graph_comp1_miRNA <- plotLoadings(model, comp = 1, contrib = 'max', method = 'median')
write.csv(loadings_graph_comp1$miRNA, 'Loadings_Plot_miRNA_Comp1.csv', row.names = T)

loadings_graph_comp2_mRNA <- plotLoadings(model, comp = 2, contrib = 'max', method = 'median')
write.csv(loadings_graph_comp2$mRNA, 'Loadings_Plot_mRNA_Comp2.csv', row.names = T)

loadings_graph_comp2_miRNA <- plotLoadings(model, comp = 2, contrib = 'max', method = 'median')
write.csv(loadings_graph_comp2$miRNA, 'Loadings_Plot_miRNA_Comp2.csv', row.names = T)

#Circos Plot image export

png('cicrosPlot, cutoff 07.png', width = 10, height =20, units = 'in', res= 300)
circosPlot(Model, cutoff = 0.7,
           color.blocks= c('yellow2', 'magenta3'),
           color.cor = c("red2","dodgerblue1"), size.labels = 1)
dev.off()

# Heatmap image export

pdf('Heatmap.pdf')
cimDiablo(model, size.legend = 0.8)
dev.off()