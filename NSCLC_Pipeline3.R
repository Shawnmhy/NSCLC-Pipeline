library(spatstat); library(Rfast2); library(data.table); library(Rfast)
library(dplyr); library(tidyr); library(reshape2); library(tibble); library(tidyverse)
library(ggplot2); library(umap); library(ggvoronoi); library(Rtsne); library(ComplexHeatmap)
library(ggforce); library(UpSetR); library(rstatix); library(dendextend)
library(tripack); library(igraph); library(RANN)
library(readxl); library(stringr); library(rjson); library(pracma)
library(survival); library(survminer)
library(ggh4x)
# working directory
current_path <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
setwd('..')
source('./Codes/Function.R')


#---------------------- DATA INPUT ----------------#
# read single-cell data
NSCLCdata <- readRDS('./NSCLCdataset.rds')
# read patient data
patient_data <- read.csv('./NSCLC Clinical Data.csv')
# core demographic data
Core_demo <- read.csv('./Clinicopathology parameters.csv')
# polygon data

#--------------------------------------------------#
#


# PREREQUISITE DATA: INHERITED FROM PIPELINE 1
#-------------- Cell type fraction per image with hierarchical clustering ------------#
# Create table with celltype fractions

cur_df <- NSCLCdata %>%
  # remove control cores
  dplyr::filter(core %in% Core_demo$Core) %>%
  group_by(core, Phenotype) %>%
  dplyr::summarise(n=n()) %>% 
  group_by(core) %>%
  mutate(fraction = n / sum(n)) %>%
  reshape2::dcast(core ~ Phenotype, value.var = "fraction", fill=0)


matrixrownames <- cur_df$core 

# now we create a matrix from the data and cluster the data based on the cell fractions
hm_dat = as.matrix(cur_df[,-1])
rownames(hm_dat) <- as.character(matrixrownames)

# calculate distance and then cluster images based on cluster fraction
dd <- dist((hm_dat), method = "euclidean")
hc <- hclust(dd, method = "ward.D2")




# see cluster tree
cl_tree <- cutree(hc, k = 3) %>%
  data.frame() %>%
  rownames_to_column() %>%
  `colnames<-` (c('Core', 'Cluster')) %>%
  merge(Core_demo, by = 'Core') %>%
  dplyr::filter(Region == 'Central')

table_tree <- table(cl_tree[, c(2, 5)]) %>%
  data.frame()
clus_dist <- reshape(table_tree, idvar = "PatientID", timevar = "Cluster", direction = "wide")


# method: find the maximum
group <- apply(clus_dist[, -1], 1, which.max)
clus_dist$group <- as.numeric(group)
patient_data_clus <- merge(clus_dist, patient_data, by = 'PatientID')
patient_data_clus <- patient_data_clus %>%
  dplyr::filter(group != 2)


mean(patient_data_clus[patient_data_clus$group == 1,]$Time)
mean(patient_data_clus[patient_data_clus$group == 3,]$Time)

# CONSTRUCT SPATIAL TOPOLOGICAL GRAPH
#-------------- Cell type fraction per image with hierarchical clustering ------------#
# Create table for data storage



# index cell type
index_types <- unique(NSCLCdata$Phenotype)


# assortativity data frame
assort_df <- data.frame(matrix(nrow = 0, ncol = length(index_types) + 2))
colnames(assort_df) <- c(index_types, 'Group', 'Core')

assort_df2 <- data.frame(matrix(nrow = 0, ncol = 8))
colnames(assort_df2) <- c(unique(NSCLCdata$Phenotype), 'Group', 'Core')


# centrality coefficients data frame
degrees_all <- data.frame(matrix(nrow = 0, ncol = 0))
betweenness_all <- data.frame(matrix(nrow = 0, ncol = 0))
closeness_all <- data.frame(matrix(nrow = 0, ncol = 0))
transitivities_all <- data.frame(matrix(nrow = 0, ncol = 0))
colnames(NSCLCdata)[1] <- 'Core'



for(file in list.files('./polygonList')){
  
  #file <- list.files('./polygonList')[4]
  spName <- strsplit(file, '\\]|,') # regex to extract cores
  
  core <- paste0(LETTERS[as.numeric(spName[[1]][2])], spName[[1]][3]) # pieces to cores
  
  # cell position
  CellPos <- NSCLCdata %>%
    dplyr::filter(Core == core)
  
  # cell boundary
  CellBound = fread(paste0('./polygonList/', file)) %>%
    data.frame() %>%
    merge(CellPos, by = 'CellID') 
  
  
  
  tryCatch(
    expr = {
      
      # patient data
      location <- Core_demo[Core_demo$Core == core, 'Region']
      
      if(location == 'Central'){
        
        
        PatientID <- Core_demo[Core_demo$Core == core, 'PatientID']
        
        # group number
        Group <- patient_data_clus[patient_data_clus$PatientID == PatientID, 'group']
        #--------- Cell Map -----------#
        
        
        
        
        # save image
        
        #--------- Cell Map -----------#
        
        p <- ggplot(CellBound) +
          theme_void() +
          geom_polygon(aes(PosX, PosY, group = CellID, fill = Phenotype)) +
          scale_fill_manual(values = c('CD8' = '#f8f8f8', 'FoxP3' = '#33a02c', 'CD163' = '#f8f8f8', 'Other' = '#f8f8f8',
                                       'Tumor' = '#e31b1c', 'FoxP3CD8' = '#f8f8f8')) +
          theme(legend.position="none") +
          #ylim(700, 1428) + # 670 - 1398 (Patient 16)
          #xlim(1200, 1928) + # 830 - 1558 (Patient 16)
          coord_fixed(ratio = 1)
        p
        
        #ggsave(p, file=paste0("./Figures/Spatial Clustering/Patient_", PatientID, '_', core, '.jpeg'), width = 8, height = 8, units = "in", dpi = 300)
        
        
        
        #--------- Spatial Graph -------#
        
        #-------------- USE FUNCTION HERE ---------------#
        
        lists <- Network(CellPos[, c('CellXPos', 'CellYPos', 'Phenotype', 'ExprPhenotype', 'CellID')], 60)
        Dual_EdgeList <- lists[[1]]
        Dual_NodeList <- lists[[2]]
        
        Dual_EdgeList$from <- Dual_NodeList[Dual_EdgeList$from, 'nodes']
        Dual_EdgeList$to <- Dual_NodeList[Dual_EdgeList$to, 'nodes']
        
        # split Dual_Edgelist to two part, part 1 contains connected nodes of same kinds
        
        # ------ part 1 -----#
        
        Dual_EdgeList_types <- data.frame(from = CellPos[match(Dual_EdgeList$from, CellPos$CellID), 'Phenotype'],
                                          to = CellPos[match(Dual_EdgeList$to, CellPos$CellID), 'Phenotype'])
        
        
        Dual_EdgeList_pt1 <- Dual_EdgeList[which(Dual_EdgeList_types$from == Dual_EdgeList_types$to), ]
        Dual_EdgeList_pt1_types <- Dual_EdgeList_types[which(Dual_EdgeList_types$from == Dual_EdgeList_types$to), ]
        
        # ------ part 2 -----#
        
        Dual_EdgeList_pt2 <- Dual_EdgeList[which(Dual_EdgeList_types$from != Dual_EdgeList_types$to), ]
        Dual_EdgeList_pt2_types <- Dual_EdgeList_types[which(Dual_EdgeList_types$from == Dual_EdgeList_types$to), ]
        
        p <- ggplot() +
          #theme_bw() +
          theme_void() +
          #geom_polygon(data = mask, aes(Y, X, group = cellLabelInImage, fill = immuneGroup)) +
          geom_segment(aes(x = Dual_NodeList[match(Dual_EdgeList_pt1$from, Dual_NodeList$nodes), 2],
                           xend = Dual_NodeList[match(Dual_EdgeList_pt1$to, Dual_NodeList$nodes), 2],
                           y = Dual_NodeList[match(Dual_EdgeList_pt1$from, Dual_NodeList$nodes), 3],
                           yend = Dual_NodeList[match(Dual_EdgeList_pt1$to, Dual_NodeList$nodes), 3],
                           color = Dual_EdgeList_pt1_types$from), size = 0.4) +
          geom_segment(aes(x = Dual_NodeList[match(Dual_EdgeList_pt2$from, Dual_NodeList$nodes), 2],
                           xend = Dual_NodeList[match(Dual_EdgeList_pt2$to, Dual_NodeList$nodes), 2],
                           y = Dual_NodeList[match(Dual_EdgeList_pt2$from, Dual_NodeList$nodes), 3],
                           yend = Dual_NodeList[match(Dual_EdgeList_pt2$to, Dual_NodeList$nodes), 3],
                           color = '#b4b3b3'), size = 0.4) +
          geom_point(data = CellPos, aes(CellXPos, CellYPos, color = Phenotype), size = 1) +
          scale_color_manual(values = c('CD8' = '#a6cee3', 'FoxP3' = '#33a02c', 'CD163' = '#cab2d6', 'Other' = '#1f78b4',
                                        'Tumor' = '#e31b1c', 'FoxP3CD8' = '#fdbf6f')) + #f8f8f8
          theme(legend.position = 'none') +
          #ylim(670, 1398) +
          #xlim(830, 1558) +
          #ylim(700, 1428) + # 670 - 1398 (Patient 16)
          #xlim(1200, 1928) + # 830 - 1558 (Patient 16)
          coord_fixed(ratio = 1)
        #p
        #ggsave(p, file=paste0("./Figures/Spatial_Clustering/Patient_", PatientID, '_', core, '.jpeg'), width = 8, height = 8, units = "in", dpi = 300)
        
        # ------------------ #
        #      iGraph        #
        # ------------------ #
        
        
        #---------------- Compute centrality coefficients -------------#
        ig <- graph_from_data_frame(vertices = Dual_NodeList, d = Dual_EdgeList, directed = FALSE)
        
        ig <- set_vertex_attr(ig, 'cell_type', value = CellPos$Phenotype)  
        
        
        for(ctype in c('CD8', 'CD163', 'FoxP3', 'Tumor')){
          #ctype <- 'CD163'
          degrees <- igraph::degree(ig, v = V(ig), normalized = TRUE)%>%
            as.data.frame() %>%
            cbind(CellPos$Phenotype) %>%
            `colnames<-` (c('degrees', 'Phenotype')) %>%
            dplyr::filter(Phenotype == ctype)
          
          
          closenesses <- closeness(ig, vids = V(ig), weights = NULL, normalized = TRUE)%>%
            as.data.frame() %>%
            cbind(CellPos$Phenotype) %>%
            `colnames<-` (c('closenesses', 'Phenotype')) %>%
            dplyr::filter(Phenotype == ctype)
          
          
          betweenesses <- betweenness(ig, v = V(ig), normalized = TRUE)%>%
            as.data.frame() %>%
            cbind(CellPos$Phenotype) %>%
            `colnames<-` (c('betweennesses', 'Phenotype')) %>%
            dplyr::filter(Phenotype == ctype)
          
          
          #transitivities <- transitivity(ig, type = 'localundirected', vids = V(ig)[V(ig)$`cell_type` == ctype])
          transitivities <- transitivity(ig, type = 'localundirected', vids = V(ig)) %>%
            as.data.frame() %>%
            cbind(CellPos$Phenotype) %>%
            `colnames<-` (c('transitivities', 'Phenotype')) %>%
            dplyr::filter(Phenotype == ctype)
          
          degrees_all <- rbind.data.frame(degrees_all, cbind.data.frame(degrees, core, Group, PatientID))
          closeness_all <- rbind.data.frame(closeness_all, cbind.data.frame(closenesses, core, Group, PatientID))
          betweenness_all <- rbind.data.frame(betweenness_all, cbind.data.frame(betweenesses, core, Group, PatientID))
          transitivities_all <- rbind.data.frame(transitivities_all, cbind.data.frame(transitivities, core, Group, PatientID))
          
        }
        
        
        # construct iGraph object
        ig <- graph_from_data_frame(vertices = Dual_NodeList, d = Dual_EdgeList, directed = FALSE)
        
        # assign vertex attributes
        # For this part, we set the objective cell type as type 1 and all other cell types as type 2
        assort_vec <- data.frame(matrix(nrow = 0, ncol = length(index_types)))
        colnames(assort_vec) <- index_types
        
        for(type1 in unique(CellPos$Phenotype)){
          
          type1_id <- which(CellPos$Phenotype == type1)
          vertex_attr <- rep(2, nrow(CellPos))
          
          vertex_attr[type1_id] <- 1
          ig <- set_vertex_attr(ig, 'cell_type', value = vertex_attr)  
          
          print(assortativity_nominal(ig, types = V(ig)$cell_type))
          
          assort_vec[1, type1] <- assortativity_nominal(ig, types = V(ig)$cell_type)
        }
        
        
        assort_vec$Group <- Group
        assort_vec$Core <- core
        assort_df <- rbind.data.frame(assort_df, assort_vec)
        
        #lt_assort_df <- assort_df[assort_df$Group == 'long-term',]
        #st_assort_df <- assort_df[assort_df$Group == 'short-term',]
        
        
        # the index of the second minimum distance
        #index2 = apply(dist_mutual, 1, function(x) names(sort(x)[2]))
        
        # the cell type pair of index cell and the top closest cell
        #cell_pair <- data.frame(index_type = CellPos$Phenotype, CellPos[index2, 'Phenotype'])
        
        #cell_pair_Mat <- as.data.frame.matrix(table(cell_pair))
        
        
      }
    },
    error = function(e){ 
      # (Optional)
      # Do this if an error is caught...
    }
  )  
  
}

#--------------------------------------------------------------------------#
#-------------- Visualizations for network measurements -------------------#
#--------------------------------------------------------------------------#
library(hrbrthemes)

degrees_all_vis <- closeness_all %>%
  dplyr::filter(ctype != 'FoxP3CD8') %>%
  dplyr::filter(ctype != 'Other') 




pd = position_dodge(width = 0.75)

signif_pairs <- degrees_all_vis %>%
  select(-c('core', 'PatientID')) %>%
  melt(id.vars = c('Phenotype', 'Group')) #%>%
#select(-'core')

stat.test <- signif_pairs %>%
  group_by(Phenotype) %>%
  wilcox_test(value ~ Group) %>%
  add_significance() %>%
  add_xy_position(x = "Phenotype")

stat.test$y.position <- log(stat.test$y.position) - 2.5 # - closeness
 
p <- ggplot(degrees_all_vis, aes(x= Phenotype, y= log(closenesses + 0.000001), colour=as.factor(Group))) + 
  stat_boxplot(geom = "errorbar", width = 0.25, position=pd) + 
  scale_color_manual(values = c('1' = '#0a377b', '3' = "#d63434")) +
  geom_boxplot(outlier.color = NA) +
  #ylim(-9, -4) + # degree centrality
  #ylim(2, 20) + # betweenness centrality
  #ylim(-2, 1) + # transitivity centrality
  ylim(-5, -2) + # transitivity centrality
  theme_classic() +
  stat_pvalue_manual(stat.test, size = 5) + rremove("x.title") + rremove('y.title') +
  theme(strip.text = element_text(size = 15)) +
  #ylim(800, 1300) +
  font("xy.text", size = 20) +
  font("xy.text", size = 20) +
  rremove('legend') 
p
ggsave(p, file=paste0("./Figures/closenesses_centrality.jpeg"), width = 4, height = 5, units = "in", dpi = 300)



#------------------- Heatmap similar to Won et al ----------------------------#

nn_DF_LT <- data.frame(matrix(nrow = 4, ncol = 4))
nn_DF_LT[is.na(nn_DF_LT)] <- 0
colnames(nn_DF_LT) <- c('CD163', 'CD8', 'FoxP3', 'Tumor')
rownames(nn_DF_LT) <- c('CD163', 'CD8', 'FoxP3', 'Tumor')

nn_DF_ST <- nn_DF_LT

for(file in list.files('./polygonList')){
  
  #file <- list.files('./polygonList')[4]
  spName <- strsplit(file, '\\]|,') # regex to extract cores
  
  core <- paste0(LETTERS[as.numeric(spName[[1]][2])], spName[[1]][3]) # pieces to cores
  
  # cell position
  CellPos <- NSCLCdata %>%
    dplyr::filter(Core == core) %>%
    dplyr::filter(Phenotype != 'Other') %>%
    dplyr::filter(Phenotype != 'FoxP3CD8') 
  
  
  tryCatch(
    expr = {
      
      # patient data
      location <- Core_demo[Core_demo$Core == core, 'Region']
      
      if(location == 'Central'){
        
        print(core)
        
        PatientID <- Core_demo[Core_demo$Core == core, 'PatientID']
        
        # group number
        Group <- patient_data_clus[patient_data_clus$PatientID == PatientID, 'group']
        
        
        
        
        
        res <- nn2(CellPos[, c('CellXPos', 'CellYPos')], query = CellPos[, c('CellXPos', 'CellYPos')], 
                   k = 2, treetype = 'kd', searchtype = 'priority')[['nn.idx']] %>%
          data.frame() %>%
          `colnames<-` (c('index', 'nn1'))
        
        # find the nearest neighbors of type i
        
        res$index <- CellPos[res$index, 'Phenotype']
        res$nn1 <- CellPos[res$nn1, 'Phenotype']
        #res$nn2 <- CellPos[res$nn2, 'Phenotype']
        #res$nn3 <- CellPos[res$nn3, 'Phenotype']
        
        res_mlt <- melt(res, id.vars=c("index"))
        
        nn_DF <- table(res_mlt[, c('index', 'value')]) %>%
          as.data.frame.matrix() 
        
       for(rid in 1:4){
         
         rn <- nn_DF[rid,] %>%
           rownames()
         
         # total number
         tn <- CellPos %>%
           dplyr::filter(Phenotype == rn) %>%
           nrow()
        
         # normalization
         nn_DF[rid,] <- nn_DF[rid,] / tn
         
       }
       
        if(Group == 1){
          nn_DF_LT <- nn_DF_LT + nn_DF  
        }
        if(Group == 3){
          nn_DF_ST <- nn_DF_ST + nn_DF  
          
        }
        
        
        

      }
    },
    error = function(e){ 
      # (Optional)
      # Do this if an error is caught...
    }
  )
}


heatmap(as.matrix(nn_DF_ST - nn_DF_LT))
heatmap(as.matrix(nn_DF_LT))



#------------------- Mimimum spanning tree ----------------------------#

miniDist_DF_LT <- data.frame(matrix(nrow = 4, ncol = 4))
miniDist_DF_LT[is.na(miniDist_DF_LT)] <- 0
colnames(miniDist_DF_LT) <- c('CD163', 'CD8', 'FoxP3', 'Tumor')
rownames(miniDist_DF_LT) <- c('CD163', 'CD8', 'FoxP3', 'Tumor')
miniDist_DF_all <- data.frame(matrix(nrow = 0, ncol = 0))
res_all <- data.frame(matrix(nrow = 0, ncol = 0))
miniDist_DF_ST <- miniDist_DF_LT

for(file in list.files('./polygonList')){
  
  #file <- list.files('./polygonList')[1]
  spName <- strsplit(file, '\\]|,') # regex to extract cores
  
  core <- paste0(LETTERS[as.numeric(spName[[1]][2])], spName[[1]][3]) # pieces to cores
  
  # cell position
  CellPos <- NSCLCdata %>%
    dplyr::filter(Core == core) %>%
    dplyr::filter(Phenotype != 'Other') %>%
    dplyr::filter(Phenotype != 'FoxP3CD8')%>%
    dplyr::filter(ExprPhenotype != '68') %>%
    dplyr::filter(!(Phenotype == 'Tumor' & ExprPhenotype == '64')) %>%
    dplyr::filter(!(Phenotype == 'CD163' & ExprPhenotype == '64'))
  
  # 
  CellPos$TruePheno <- as.factor(paste0(CellPos$Phenotype, '_', CellPos$ExprPhenotype)) %>%
    as.numeric
  
  
  unique(CellPos$TruePheno)
  tryCatch(
    expr = {
      
      # patient data
      location <- Core_demo[Core_demo$Core == core, 'Region']
      
      if(location == 'Central'){
        
        print(core)
        
        PatientID <- Core_demo[Core_demo$Core == core, 'PatientID']
        
        # group number
        Group <- patient_data_clus[patient_data_clus$PatientID == PatientID, 'group']
        
        
        
        
        
        res <- nn2(CellPos[, c('CellXPos', 'CellYPos')], query = CellPos[, c('CellXPos', 'CellYPos')], 
                   k = 2, treetype = 'kd', searchtype = 'priority')[['nn.idx']] %>%
          data.frame() %>%
          `colnames<-` (c('index', 'nn1'))
        
        res$dist <- nn2(CellPos[, c('CellXPos', 'CellYPos')], query = CellPos[, c('CellXPos', 'CellYPos')], 
                        k = 2, treetype = 'kd', searchtype = 'priority')[['nn.dists']][,2]
          
        
        # find the nearest neighbors of type i

        res$index <- CellPos[res$index, 'TruePheno']
        res$nn1 <- CellPos[res$nn1, 'TruePheno']
        
        res$combination <- paste0(res$index, '_', res$nn1)
        res_aggregated <- aggregate(res, list(res$combination), FUN = 'mean')
        
        res_aggregated$index <- strsplit(res_aggregated[, 'Group.1'], '_') %>%
          sapply("[[",1)
        
        res_aggregated$nn1 <- strsplit(res_aggregated[, 'Group.1'], '_') %>%
          sapply("[[",2)
        
        colnames(res_aggregated)[1] <- 'combination'
        res_all <- rbind.data.frame(res_all, cbind.data.frame(res_aggregated[, 1:4], Group, core))
        
        
        
      }
    },
    error = function(e){ 
      # (Optional)
      # Do this if an error is caught...
    }
  )
}



res_aggregated <- aggregate(res_all, list(res_all$combination, res_all$Group), FUN = 'mean')

res_aggregated[res_aggregated$Group == 1, 'index']<- strsplit(res_aggregated[res_aggregated$Group == 1, 'Group.1'], '_') %>%
  sapply("[[",1)

res_aggregated[res_aggregated$Group == 1, 'nn1']<- strsplit(res_aggregated[res_aggregated$Group == 1, 'Group.1'], '_') %>%
  sapply("[[",2)

res_aggregated[res_aggregated$Group == 3, 'index']<- strsplit(res_aggregated[res_aggregated$Group == 3, 'Group.1'], '_') %>%
  sapply("[[",1)

res_aggregated[res_aggregated$Group == 3, 'nn1']<- strsplit(res_aggregated[res_aggregated$Group == 3, 'Group.1'], '_') %>%
  sapply("[[",2)

res_aggregated_LT <- res_aggregated[res_aggregated$Group == 1,]
res_aggregated_LT <-  xtabs(dist ~index + nn1, res_aggregated_LT)
res_aggregated_ST <- res_aggregated[res_aggregated$Group == 3,]
res_aggregated_ST <-  xtabs(dist ~index + nn1, res_aggregated_ST)




#Mutual.dist <- data.frame(matrix(nrow = 0, ncol = 0))



cols <- c('#cab2d6', '#e31b1c', 
          '#cab2d6', '#a6cee3', 
          '#a6cee3', '#a6cee3', '#33a02c',
          '#33a02c', '#33a02c', '#e31b1c')

pdf("Figures/network_mindist.pdf",width=7,height=5)
set.seed(4)
g<-graph.adjacency(-res_aggregated_ST,weighted=T)
g_mst<-mst(g)
network5<-as_adjacency_matrix(g_mst)
network5<-graph_from_adjacency_matrix(network5, mode="undirected", diag=F)

plot.igraph(network5,
            vertex.color=cols,
            vertex.size=12,
            edge.arrow.size = 0.5,
            edge.arrow.width = 1,
            vertex.label.dist=0,
            vertex.label.family="sans", #sans = Arial, serif= Times
            vertex.label.font=2,
            vertex.label.cex=0.75,
            vertex.label.color="black",
            vertex.frame.color="transparent")
title("short-term survivor",cex.main= 1,col.main="black")



g<-graph.adjacency(-res_aggregated_LT,weighted=T)
g_mst<-mst(g)
network5<-as_adjacency_matrix(g_mst)
network5<-graph_from_adjacency_matrix(network5, mode="undirected", diag=F)
plot.igraph(network5,
            vertex.color=cols,
            vertex.size=12,
            edge.arrow.size = 0.5,
            edge.arrow.width = 1,
            vertex.label.dist=0,
            vertex.label.family="sans", #sans = Arial, serif= Times
            vertex.label.font=2,
            vertex.label.cex=0.75,
            vertex.label.color="black",
            vertex.frame.color="transparent")
title("long-term survivor",cex.main=1,col.main="black")


dev.off()
#cols <- c('CD8' = '#a6cee3', 'FoxP3' = '#33a02c', 'CD163' = '#cab2d6','Tumor' = '#e31b1c')






#---------------- pairwise Gcross function ------------------#

group1_all <- data.frame(matrix(nrow = 0, ncol = 0))
group2_all <- data.frame(matrix(nrow = 0, ncol = 0))
group3_all <- data.frame(matrix(nrow = 0, ncol = 0))
Mutual.dist <- data.frame(matrix(nrow = 0, ncol = 0))

for(file in list.files('./boundary')) {
  
  #file <- list.files('./polygonList')[80] # This is M9
  spName <- strsplit(file, '\\]|,') # regex to extract cores
  
  Sample <- paste0(LETTERS[as.numeric(spName[[1]][2])], spName[[1]][3]) # pieces to cores
  
  #ID <- 2089
  Location <- Core_demo[Core_demo$Core == Sample, 'Region']
  
  
  PatientID <- Core_demo[Core_demo$Core == Sample, 'PatientID']
  Group <- patient_data_clus[patient_data_clus$PatientID == PatientID, 'group']
  
  coreDat <- NSCLCdata %>%
    filter(Core == Sample)
  
  # ROI data
  #ROI <- 'ROI03'
  # CD163 cells coords
  posCD163 <- coreDat %>%
    filter(Phenotype == 'CD163')
  
  # FoxP3 cells coords
  posFoxP3 <- coreDat %>%
    filter(Phenotype == 'FoxP3')
  
  # CD8 cells coords
  posCD8 <- coreDat %>%
    filter(Phenotype == 'CD8')
  
  # Tumor cell coords
  posTumor <- coreDat %>%
    filter(Phenotype == 'Tumor')
  
  
  # define region to compute Kcross function
  # Kcross functions
  # read boundaries
  tisse_files <- list.files(paste0('./QuPath/Tissue_Boundary/', Sample))
  
  
  tryCatch(
    expr = {
      
      if(Location == 'Central'){
        Tissue_all <- data.frame(matrix(nrow = 0, ncol = 0))
        for(piece in tisse_files){
          Tissue <- fromJSON(file = paste0('./QuPath/Tissue_Boundary/', Sample, '/', piece))%>%
            with(geometry) %>%
            with(coordinates) %>%
            `[[`(1) %>%
            as.data.frame() %>%
            t() %>%
            as.data.frame()
          
          
          flag <- clockwise(data.frame(x = Tissue$V1, y = Tissue$V2))
          
          # if clockwise, counter-clockwise it
          if(flag){
            Tissue$V1 <- rev(Tissue$V1)
            Tissue$V2 <- rev(Tissue$V2)
          }
          
          
          
          Tissue_all <- rbind(Tissue_all, cbind(Tissue, piece))
        }
        
        colnames(Tissue_all) <- c('x', 'y', 'group')
        
        #ggplot() +
        #  geom_polygon(aes(Tissue_all[,1], Tissue_all[,2], group = Tissue_all[,3]), fill = NA, color = 'black') +
        #  geom_point(data = coreDat, aes(CellXPos, CellYPos))
        
        
        # CD8 + FoxP3
        list1 <- bivarAnalysis.Kcross('posCD8', 'posFoxP3', Region = list(Tissue_all), neigbor_thresh = 60, step_thresh = 1)
        
        list1_subgroups1 <- split(list1[[1]], list1[[1]]$interval)
        list1_subgroups1_g1 <- list1_subgroups1$`1` %>%
          mutate(category = 'CD8-FoxP3') %>%
          mutate(Group = Group)
        list1_subgroups1_g2 <- list1_subgroups1$`2` %>%
          mutate(category = 'CD8-FoxP3') %>%
          mutate(Group = Group)
        list1_subgroups1_g3 <- list1_subgroups1$`3` %>%
          mutate(category = 'CD8-FoxP3') %>%
          mutate(Group = Group)
        
        list1_subgroups2 <- split(list1[[2]], list1[[2]]$interval)
        list1_subgroups2_g1 <- list1_subgroups2$`1` %>%
          mutate(category = 'FoxP3-CD8') %>%
          mutate(Group = Group)
        list1_subgroups2_g2 <- list1_subgroups2$`2` %>%
          mutate(category = 'FoxP3-CD8') %>%
          mutate(Group = Group)
        list1_subgroups2_g3 <- list1_subgroups2$`3` %>%
          mutate(category = 'FoxP3-CD8') %>%
          mutate(Group = Group)
        
        
        list1_AUCs <- bivarAnalysis.Kcross_AUC('posCD8', 'posFoxP3', Region = list(Tissue_all), neigbor_thresh = 60)
        
        
        
        # CD8 + CD163
        list2 <- bivarAnalysis.Kcross('posCD8', 'posCD163', Region = list(Tissue_all), neigbor_thresh = 60, step_thresh = 1)
        list2_subgroups1 <- split(list2[[1]], list2[[1]]$interval)
        list2_subgroups1_g1 <- list2_subgroups1$`1` %>%
          mutate(category = 'CD8-CD163') %>%
          mutate(Group = Group)
        list2_subgroups1_g2 <- list2_subgroups1$`2` %>%
          mutate(category = 'CD8-CD163') %>%
          mutate(Group = Group)
        list2_subgroups1_g3 <- list2_subgroups1$`3` %>%
          mutate(category = 'CD8-CD163') %>%
          mutate(Group = Group)
        
        list2_subgroups2 <- split(list2[[2]], list2[[2]]$interval)
        list2_subgroups2_g1 <- list2_subgroups2$`1` %>%
          mutate(category = 'CD163-CD8') %>%
          mutate(Group = Group)
        list2_subgroups2_g2 <- list2_subgroups2$`2` %>%
          mutate(category = 'CD163-CD8') %>%
          mutate(Group = Group)
        list2_subgroups2_g3 <- list2_subgroups2$`3` %>%
          mutate(category = 'CD163-CD8') %>%
          mutate(Group = Group)
        
        
        list2_AUCs <- bivarAnalysis.Kcross_AUC('posCD8', 'posCD163', Region = list(Tissue_all), neigbor_thresh = 60)
        
        
        
        # CD8 + Tumor
        list3 <- bivarAnalysis.Kcross('posCD8', 'posTumor', Region = list(Tissue_all), neigbor_thresh = 60, step_thresh = 1)
        list3_subgroups1 <- split(list3[[1]], list3[[1]]$interval)
        list3_subgroups1_g1 <- list3_subgroups1$`1` %>%
          mutate(category = 'CD8-Tumor') %>%
          mutate(Group = Group)
        list3_subgroups1_g2 <- list3_subgroups1$`2` %>%
          mutate(category = 'CD8-Tumor') %>%
          mutate(Group = Group)
        list3_subgroups1_g3 <- list3_subgroups1$`3` %>%
          mutate(category = 'CD8-Tumor') %>%
          mutate(Group = Group)
        
        list3_subgroups2 <- split(list3[[2]], list3[[2]]$interval)
        list3_subgroups2_g1 <- list3_subgroups2$`1` %>%
          mutate(category = 'Tumor-CD8') %>%
          mutate(Group = Group)
        list3_subgroups2_g2 <- list3_subgroups2$`2` %>%
          mutate(category = 'Tumor-CD8') %>%
          mutate(Group = Group)
        list3_subgroups2_g3 <- list3_subgroups2$`3` %>%
          mutate(category = 'Tumor-CD8') %>%
          mutate(Group = Group)
        
        list3_AUCs <- bivarAnalysis.Kcross_AUC('posCD8', 'posTumor', Region = list(Tissue_all), neigbor_thresh = 60)
        
        
        # FoxP3 + CD163
        list4 <- bivarAnalysis.Kcross('posFoxP3', 'posCD163', Region = list(Tissue_all), neigbor_thresh = 60, step_thresh = 1)
        list4_subgroups1 <- split(list4[[1]], list4[[1]]$interval)
        list4_subgroups1_g1 <- list4_subgroups1$`1` %>%
          mutate(category = 'FoxP3-CD163') %>%
          mutate(Group = Group)
        list4_subgroups1_g2 <- list4_subgroups1$`2` %>%
          mutate(category = 'FoxP3-CD163') %>%
          mutate(Group = Group)
        list4_subgroups1_g3 <- list4_subgroups1$`3` %>%
          mutate(category = 'FoxP3-CD163') %>%
          mutate(Group = Group)
        
        list4_subgroups2 <- split(list4[[2]], list4[[2]]$interval)
        list4_subgroups2_g1 <- list4_subgroups2$`1` %>%
          mutate(category = 'CD163-FoxP3') %>%
          mutate(Group = Group)
        list4_subgroups2_g2 <- list4_subgroups2$`2` %>%
          mutate(category = 'CD163-FoxP3') %>%
          mutate(Group = Group)
        list4_subgroups2_g3 <- list4_subgroups2$`3` %>%
          mutate(category = 'CD163-FoxP3') %>%
          mutate(Group = Group)
        
        list4_AUCs <- bivarAnalysis.Kcross_AUC('posFoxP3', 'posCD163', Region = list(Tissue_all), neigbor_thresh = 60)
        
        
        # FoxP3 + Tumor
        list5 <- bivarAnalysis.Kcross('posFoxP3', 'posTumor', Region = list(Tissue_all), neigbor_thresh = 60, step_thresh = 1)
        list5_subgroups1 <- split(list5[[1]], list5[[1]]$interval)
        list5_subgroups1_g1 <- list5_subgroups1$`1` %>%
          mutate(category = 'FoxP3-Tumor') %>%
          mutate(Group = Group)
        list5_subgroups1_g2 <- list5_subgroups1$`2` %>%
          mutate(category = 'FoxP3-Tumor') %>%
          mutate(Group = Group)
        list5_subgroups1_g3 <- list5_subgroups1$`3` %>%
          mutate(category = 'FoxP3-Tumor') %>%
          mutate(Group = Group)
        
        list5_subgroups2 <- split(list5[[2]], list5[[2]]$interval)
        list5_subgroups2_g1 <- list5_subgroups2$`1` %>%
          mutate(category = 'Tumor-FoxP3') %>%
          mutate(Group = Group)
        list5_subgroups2_g2 <- list5_subgroups2$`2` %>%
          mutate(category = 'Tumor-FoxP3') %>%
          mutate(Group = Group)
        list5_subgroups2_g3 <- list5_subgroups2$`3` %>%
          mutate(category = 'Tumor-FoxP3') %>%
          mutate(Group = Group)
        
        list5_AUCs <- bivarAnalysis.Kcross_AUC('posFoxP3', 'posTumor', Region = list(Tissue_all), neigbor_thresh = 60)
        
        
        # Tumor + CD163
        list6 <- bivarAnalysis.Kcross('posTumor', 'posCD163', Region = list(Tissue_all), neigbor_thresh = 60, step_thresh = 1)
        list6_subgroups1 <- split(list6[[1]], list6[[1]]$interval)
        list6_subgroups1_g1 <- list1_subgroups1$`1` %>%
          mutate(category = 'Tumor-CD163') %>%
          mutate(Group = Group)
        list6_subgroups1_g2 <- list1_subgroups1$`2` %>%
          mutate(category = 'Tumor-CD163') %>%
          mutate(Group = Group)
        list6_subgroups1_g3 <- list1_subgroups1$`3` %>%
          mutate(category = 'Tumor-CD163') %>%
          mutate(Group = Group)
        
        list6_subgroups2 <- split(list6[[2]], list6[[2]]$interval)
        list6_subgroups2_g1 <- list1_subgroups2$`1` %>%
          mutate(category = 'CD163-Tumor') %>%
          mutate(Group = Group)
        list6_subgroups2_g2 <- list1_subgroups2$`2` %>%
          mutate(category = 'CD163-Tumor') %>%
          mutate(Group = Group)
        list6_subgroups2_g3 <- list1_subgroups2$`3` %>%
          mutate(category = 'CD163-Tumor') %>%
          mutate(Group = Group)
        
        list6_AUCs <- bivarAnalysis.Kcross_AUC('posTumor', 'posCD163', Region = list(Tissue_all), neigbor_thresh = 60)
        
        
        #------------compute AUC -----------#
        
        group1_all <- rbind.data.frame(group1_all, list1_subgroups1_g1, list1_subgroups2_g1, list2_subgroups1_g1, list2_subgroups2_g1, list3_subgroups1_g1,
                                       list3_subgroups2_g1, list4_subgroups1_g1, list4_subgroups2_g1, list5_subgroups1_g1, list5_subgroups2_g1,
                                       list6_subgroups1_g1, list6_subgroups2_g1)
        
        
        group2_all <- rbind.data.frame(group2_all, list1_subgroups1_g2, list1_subgroups2_g2, list2_subgroups1_g2, list2_subgroups2_g2, list3_subgroups1_g2,
                                       list3_subgroups2_g2, list4_subgroups1_g2, list4_subgroups2_g2, list5_subgroups1_g2, list5_subgroups2_g2,
                                       list6_subgroups1_g2, list6_subgroups2_g2)
        
        
        group3_all <- rbind.data.frame(group3_all, list1_subgroups1_g3, list1_subgroups2_g3, list2_subgroups1_g3, list2_subgroups2_g3, list3_subgroups1_g3,
                                       list3_subgroups2_g3, list4_subgroups1_g3, list4_subgroups2_g3, list5_subgroups1_g3, list5_subgroups2_g3,
                                       list6_subgroups1_g3, list6_subgroups2_g3)
        
        Mutual.dist <- rbind.data.frame(Mutual.dist, cbind.data.frame(Sample, list1_AUCs[[1]], list1_AUCs[[2]],  list2_AUCs[[1]], list2_AUCs[[2]], list3_AUCs[[1]], list3_AUCs[[2]],
                                                                      list4_AUCs[[1]], list4_AUCs[[2]], list5_AUCs[[1]], list5_AUCs[[2]], list6_AUCs[[1]], list6_AUCs[[2]], Sample, PatientID, Group))
        
        print(Sample)
        
      }
    },
    error = function(e){ 
      # (Optional)
      # Do this if an error is caught...
    }
  )  
  
  
  
  
}

#------------------ Voronoi Tesselation -------------------#
colnames(Mutual.dist)[1] <- 'core'

newdata <- Mutual.dist[order(-Mutual.dist[,8]),]

test <- merge(newdata, cur_df, by = 'core')


cellPos <- NSCLCdata %>%
  filter(Core == 'E6') %>%
  #filter(Phenotype != 'Other') %>%
  select('CellXPos', 'CellYPos', 'Phenotype', 'ExprPhenotype', 'CellID') #%>%


tisse_files <- list.files(paste0('./QuPath/Tissue_Boundary/E6'))

Tissue_all <- data.frame(matrix(nrow = 0, ncol = 0))
for(piece in tisse_files){
  Tissue <- fromJSON(file = paste0('./QuPath/Tissue_Boundary/E6', '/', piece))%>%
    with(geometry) %>%
    with(coordinates) %>%
    `[[`(1) %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame()
  
  
  flag <- clockwise(data.frame(x = Tissue$V1, y = Tissue$V2))
  
  # if clockwise, counter-clockwise it
  if(flag){
    Tissue$V1 <- rev(Tissue$V1)
    Tissue$V2 <- rev(Tissue$V2)
  }
  
  
  
  Tissue_all <- rbind(Tissue_all, cbind(Tissue, piece))
}

p <- ggplot(cellPos, aes(CellXPos, CellYPos, fill = Phenotype)) +
  geom_voronoi(outline = Tissue_all[,1:2],
               color = 1, size = 0.1) +
  scale_fill_manual(values = c('CD8' = '#f8f8f8', 'FoxP3' = '#33a02c', 'CD163' = '#cab2d6', 'Other' = '#f8f8f8',
                               'Tumor' = '#f8f8f8', 'FoxP3CD8' = '#f8f8f8')) + #f8f8f8
  #scale_fill_gradient(low = "#B9DDF1",
  #                    high = "#2A5783",
  #                    guide = FALSE) +
  theme_void() +
  theme(legend.position = 'none') +
  coord_fixed()
p
ggsave(p, file=paste0("./Figures/E6_CD163_FoxP3.jpeg"), width = 6, height = 4, units = "in", dpi = 300)



# statistical test to find significant pairs
p_all <- c()
for(id in 2:13){
  
  id <- 10
  g1 <- Mutual.dist[Mutual.dist$Group == 1, id]
  g2 <- Mutual.dist[Mutual.dist$Group == 3, id]
  
  mean(g1)
  mean(g2)
  test <- wilcox.test(g1, g2)  
  p_all <- c(p_all, test$p.value)
}

p.adjust(p_all)
#---------- FoxP3 -> CD163 ---------#
signif_pairs <- Mutual.dist[, c(1,8:10, 13, 15:16)]
colnames(signif_pairs)[2:4] <- c('FoxP3 -> CD163', 'CD163 -> FoxP3', 'FoxP3 -> Tumor')

signif_pairs <- signif_pairs %>%
  select(-'PatientID') %>%
  melt(id.vars = c('core', 'Group')) %>%
  select(-'core')

stat.test <- signif_pairs %>%
  group_by(variable) %>%
  wilcox_test(value ~ Group) %>%
  add_significance()

stat.test 


p <- ggplot(data = signif_pairs, aes(x = as.factor(Group), y = value)) +
  stat_boxplot(geom ='errorbar', width = 0.3) +
  geom_boxplot(width = 0.6, fill = "white") +
  ylim(-40, 25) +
  #stat_pvalue_manual(stat.test, size = 6) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        axis.text.y = element_text(angle = 90),
        strip.text = element_text(size = 20),
        strip.background = element_rect(fill = "white")) +
  ylab('') +
  ggh4x::facet_wrap2(.~variable, scales = 'free', nrow = 1) 
p


ggsave(p, file=paste0("./Figures/Gcross_comparison.jpeg"), width = 12, height = 5, units = "in", dpi = 300)





#-- Archives ----#

ls <- group3_all[group3_all$Group == 1 & group3_all$category == 'FoxP3-CD163', 'rs']

ss <- group3_all[group3_all$Group == 3 & group3_all$category == 'FoxP3-CD163', 'rs']

boxplot(ls, ss)



# combine selected data into one data frame

Group123_combined <- rbind.data.frame(group1_all[group1_all$category == 'FoxP3-CD163', ],
                                      group2_all[group2_all$category == 'FoxP3-CD163', ],
                                      group3_all[group3_all$category == 'FoxP3-CD163', ],
                                      group1_all[group1_all$category == 'FoxP3-Tumor', ],
                                      group2_all[group2_all$category == 'FoxP3-Tumor', ],
                                      group3_all[group3_all$category == 'FoxP3-Tumor', ])



stat.test <- Group123_combined %>% 
  group_by(category, interval) %>%
  wilcox_test(rs ~ Group) %>% 
  add_significance() %>%
  add_xy_position(x = "Group")



# New facet label names for supp variable



p <- ggplot(data = Group123_combined, aes(x = as.factor(Group), y = rs, color = as.factor(Group))) +
  stat_boxplot(geom ='errorbar', width = 0.4) +
  geom_boxplot(width = 0.6, fill = "white") +
  theme_bw() +
  stat_pvalue_manual(stat.test, size = 6) +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        #axis.text.y = element_text(angle = 90),
        strip.text = element_text(size = 15),
        legend.position = 'none',
        strip.background = element_rect(fill ='white')) +
  ylab('') +
  ylim(0, 1.5) +
  scale_y_continuous(breaks = seq(0, 1.5, by = 0.5)) +
  facet_grid(category ~ interval) +
  scale_color_manual(values = c('1' = '#0a377b', '3' = "#d63434")) 
p

ggsave(p, file=paste0("./Figures/Gfunction_value.jpeg"), width = 10, height = 6, units = "in", dpi = 300)



group1_all[group1_all$category == 'FoxP3-CD163', 'rs']

#-----------------------------------------------#
#--- GGBoxplot for assortativity coefficient ---#
#-----------------------------------------------#

st_assort_df <- assort_df[assort_df$Group == 3,]
lt_assort_df <- assort_df[assort_df$Group == 1,]

#------------- facet boxplot -------------#
a <- cbind.data.frame(stack(st_assort_df[, c(1,3,4,5)]), '3') %>%
  'colnames<-' (c('values', 'ind', 'group'))
b <- cbind.data.frame(stack(lt_assort_df[, c(1,3,4,5)]), '1') %>%
  'colnames<-' (c('values', 'ind', 'group'))

assort_stacked_comb <- rbind.data.frame(a, b)


stat.test <- assort_stacked_comb %>%
  group_by(ind) %>%
  wilcox_test(values ~ group) %>%
  adjust_pvalue(method = "bonferroni") %>%
  dplyr::filter(ind == 'CD8' | ind == 'FoxP3') %>%
  #dplyr::filter(ind == 'Tumor' | ind == 'CD163') %>%
  
  add_significance()
stat.test 


# CD8 + FoxP3 part (because their scales are different, so split to two parts)
pt1 <- assort_stacked_comb %>%
  dplyr::filter(ind == 'CD8' | ind == 'FoxP3') #%>%
  #dplyr::filter(ind == 'Tumor' | ind == 'CD163') #%>%

pt1$ind <- as.character(pt1$ind)

pt1$group <- factor(pt1$group, levels = c('1', '3'))
bxp <- ggboxplot(
  pt1, x = "group", y = "values", 
  facet.by = "ind", #alpha = 0.6, 
  add = "jitter", add.params = list(size = 1),
  panel.labs.font = list(size = 14),
  color = "group", palette = c("#72b6bb", "#bf584f"),
  fill = 'white',
  nrow = 2,
  ncol = 2,
  bxp.errorbar = TRUE
)

bxp

# Make facet and add p-values
stat.test <- stat.test %>% add_xy_position(x = "group")
p <- bxp + stat_pvalue_manual(stat.test, size = 6) + rremove("x.title") + rremove('y.title') +
  scale_color_manual(values = c('1' = '#0a377b', '3' = "#d63434")) +
  theme(strip.text = element_text(size = 15),
        strip.background = element_rect(fill ='white')) +
  #ylim(-0.1, 1.1) +
  font("xy.text", size = 15) +
  font("xy.text", size = 15) +
  rremove('legend')

p$layers[[2]]$aes_params$textsize <- 6

p

ggsave(p, file=paste0("./Figures/Assortativity_CD8FoxP3.jpeg"), width = 4, height = 5 , units = "in", dpi = 300)

#------------------------------------------------#

st_assort_stacked <- stack(st_assort_df[, 1:15])

p <- st_assort_stacked %>% filter(is.na(st_assort_stacked) == FALSE) %>% 
  ggplot(aes(x=values, y=ind, fill = ind)) +
  geom_boxplot(outlier.shape = NA, lwd = 1) +
  geom_jitter(alpha = 0.5) +
  theme_classic() +
  scale_fill_manual(values = c('B' = '#82b5e5', 'CD4 T' = '#639d77', 'CD8 T' = '#4975b1', 'CD3 T' = '#f0e8c8', 
                               'NK' = '#ff9955', 'Neutrophil' = '#d38d5f', 'Macrophage' = '#b9631b', 
                               'DC' = '#ebe15a', 'DC/Mono' = '#eac6af', 'Mono/Neu' = '#a9eeff', 'Endothelial' = '#e78ac3',
                               'Treg' = '#8da1cc', 'Tumor' = '#a6d854', 'Keratin-positive tumor' = '#b4b3b3',
                               'Mesenchymal-like' = '#66c2a4')) +
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size = 35),
        #axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        #axis.ticks.y = element_blank(),
        legend.position = 'NA') +
  xlab('Assortativity Coefficient') +
  ylab('Cell type') 

ggsave(p, file=paste0("./StanfordTMA/Figures/Assortativity/short_term_assortativity.jpeg"), width = 10, height = 15, units = "in", dpi = 300)



lt_assort_stacked <- stack(lt_assort_df[, 1:15])

p <- lt_assort_stacked %>% filter(is.na(lt_assort_stacked) == FALSE) %>% 
  ggplot(aes(x=values, y=ind, fill = ind)) +
  geom_boxplot(outlier.shape = NA, lwd = 1) +
  geom_jitter(alpha = 0.5) +
  theme_classic() +
  scale_fill_manual(values = c('B' = '#82b5e5', 'CD4 T' = '#639d77', 'CD8 T' = '#4975b1', 'CD3 T' = '#f0e8c8', 
                               'NK' = '#ff9955', 'Neutrophil' = '#d38d5f', 'Macrophage' = '#b9631b', 
                               'DC' = '#ebe15a', 'DC/Mono' = '#eac6af', 'Mono/Neu' = '#a9eeff', 'Endothelial' = '#e78ac3',
                               'Treg' = '#8da1cc', 'Tumor' = '#a6d854', 'Keratin-positive tumor' = '#b4b3b3',
                               'Mesenchymal-like' = '#66c2a4')) +
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size = 35),
        #axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        #axis.ticks.y = element_blank(),
        legend.position = 'NA') +
  xlab('Assortativity Coefficient') +
  ylab('Cell type') 
p
ggsave(p, file=paste0("./StanfordTMA/Figures/Assortativity/long_term_assortativity.jpeg"), width = 10, height = 15, units = "in", dpi = 300)



#-----------------------------------------------#
#--- GGBoxplot for assortativity coefficient ---#
#-----------------------------------------------#



for(ct in unique(transitivities_all$ctype)){
  
  #ct <- 'Tumor'
  #ct <- 'Endothelial'
  #--------- aggregate to core level for comparison ----------#
  lt_val <- transitivities_all %>%
    dplyr::filter(core %in% cl_profile$Core) %>%
    dplyr::filter(ctype == ct) %>%
    dplyr::filter(Group == '1') %>%
    group_by(core) %>%
    summarise_at(vars(transitivities), list(var = ~mean(., na.rm = T))) %>%
    as.data.frame()
  
  
  st_val <- transitivities_all %>%
    dplyr::filter(ctype == ct) %>%
    dplyr::filter(core %in% cl_profile$Core) %>%
    dplyr::filter(Group == '3') %>%
    group_by(core) %>%
    summarise_at(vars(transitivities), list(var = ~mean(., na.rm = T))) %>%
    as.data.frame()
  
  #lt_val <- degrees_all[degrees_all$ctype == ct & degrees_all$category == 'long-term', 1] %>%
  #  as.numeric()
  #st_val <- degrees_all[degrees_all$ctype == ct & degrees_all$category == 'short-term', 1] %>%
  #  as.numeric()
  
  
  test_res <- wilcox.test(lt_val$var, st_val$var)
  print(paste0('Cell type: ', ct, '; value: ', test_res$p.value))
}

