library(spatstat); library(Rfast2); library(data.table); library(Rfast)
library(dplyr); library(tidyr); library(reshape2); library(tibble); library(tidyverse)
library(ggplot2); library(umap); library(ggvoronoi); library(Rtsne); library(ComplexHeatmap)
library(ggforce); library(UpSetR); library(rstatix); library(dendextend)
library(tripack); library(igraph)
library(readxl); library(stringr); library(rjson); library(pracma); library(Rmisc)
library(survival); library(survminer)
# working directory
current_path <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
setwd('..')
source('./Codes/Function.R')


#---------------------- DATA INPUT ----------------#
# read single-cell data
NSCLCdata <- readRDS('./NSCLCdataset.rds')
colnames(NSCLCdata)[1] <- 'core'

NSCLCdata %>%
  group_by(Phenotype) %>%
  tally()
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

colnames(NSCLCdata)[1] <- 'Core'





#---------- CONSTRUCT RISK SCORES ---------#
#-------- FIND OUT WHICH MARKER CONTRIBUTE TO THE PROXIMITY BETWEEN FoxP3, Tumor, and CD163


es.CD163_All <- data.frame(matrix(nrow = 0, ncol = 0))
es.FoxP3_All <- data.frame(matrix(nrow = 0, ncol = 0))
es.FoxP3_All_st <- data.frame(matrix(nrow = 0, ncol = 0))
es.Tumor_All <- data.frame(matrix(nrow = 0, ncol = 0))

for(Core in unique(NSCLCdata$Core)){
  
 # Core <- 'E6'
  PatientID <- Core_demo[Core_demo$Core == Core, 'PatientID']
  
  #PatientID <- '13'
  # group number
  Group <- patient_data_clus[patient_data_clus$PatientID == PatientID, 'group']
  # Region
  Region <- Core_demo[Core_demo$Core == Core, 'Region']
  
  if(is.empty(Group)){
    Group <- '0'
  }
  
  if(is.empty(PatientID)){
    PatientID <- '0'
  }
  
  if(is.empty(Region)){
    Region <- 'Zero'
  }
  
  # retrieve the core data
  core_data <- NSCLCdata[NSCLCdata$Core == Core, ]
  core_data$ExprPhenotype <- factor(core_data$ExprPhenotype, levels = c("0", 
                                                            "4",
                                                            "64",
                                                            "68"))  # FoxP3 data
  posFoxP3 <- core_data %>%
    dplyr::filter(Phenotype == 'FoxP3') %>%
    select(CellID, CellXPos, CellYPos, ExprPhenotype)
  
  # CD163 data
  posCD163 <- core_data %>%
    dplyr::filter(Phenotype == 'CD163') %>%
    select(CellID, CellXPos, CellYPos, ExprPhenotype)
  
  
  # FoxP3 data
  posTumor <- core_data %>%
    dplyr::filter(Phenotype == 'Tumor') %>%
    select(CellID, CellXPos, CellYPos, ExprPhenotype)
  
  
  #------- For Long-Term survivors, why CD163 and FoxP3 are close?
  
  
  if(nrow(posCD163) >= 10 & nrow(posFoxP3) >= 10 & Region == 'Central' & Group == '1'){
    
    print(Region)
    
    #Core <- 'M11'
    
    #--------- Get the tissue boundary to do the simulation-------#
    Tissue_all <- data.frame(matrix(nrow = 0, ncol = 0))
    tisse_files <- list.files(paste0('./QuPath/Tissue_Boundary/', Core))
    
    for(piece in tisse_files){
      Tissue <- fromJSON(file = paste0('./QuPath/Tissue_Boundary/', Core, '/', piece))%>%
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
      
      colnames(Tissue) <- c('x', 'y')
      
      # single polgon defined by owin
      
      Tissue_all <- rbind(Tissue_all, cbind(Tissue, piece))
    }
    
    colnames(Tissue_all) <- c('x', 'y', 'group')
    
    # convert polygon data frames to owin objects
    poly_list <- split(Tissue_all, Tissue_all$group)
    # only want lon-lats in the list, not the names
    poly_list <- lapply(poly_list, function(x) { x["group"] <- NULL; x })    
    ps <- lapply(poly_list, Polygon)
    p1 <- lapply(seq_along(ps), function(i) Polygons(list(ps[[i]]), 
                                                     ID = names(poly_list)[i]  ))
    # create SpatialPolygons object
    spatial_polys <- SpatialPolygons(p1) 
      
    owins <- maptools::as.owin.SpatialPolygons(spatial_polys)

    # 
    #ggplot() +
    #  geom_polygon(aes(Tissue_all[,1], Tissue_all[,2], group = Tissue_all[,3]), fill = NA, color = 'black') #+
      #geom_point(data = coreDat, aes(CellXPos, CellYPos))
    

    
    ### Case 1: CD4 T cell to Macrophage
    distMatrix <- flexclust::dist2(posCD163[, c('CellXPos', 'CellYPos')], posFoxP3[, c('CellXPos', 'CellYPos')]) %>%
      as.matrix()
    
    # These CD4 T cells has at least 1 adjacent Macrophages
    adj.CD163.id <- which(apply(distMatrix, 1, FUN = min) < 60)
    
    # Get the single-cell data for these CD4 T cells
    adj.CD163 <- posCD163[adj.CD163.id, ]
    
    if(!(is.empty(adj.CD163.id))){
      
      # Real ratio
     
      realRatio_CD163 <- table(adj.CD163$ExprPhenotype) %>%
        data.frame() %>%
        `colnames<-` (c('ExprPhenotype','count'))
      
      #realRatio_CD163[realRatio_CD163$ExprPhenotype == '4', 'count'] <- realRatio_CD163[realRatio_CD163$ExprPhenotype == '4', 'count'] #+
        #realRatio_CD163[realRatio_CD163$ExprPhenotype == '68', 'count']
      
      #realRatio_CD163[realRatio_CD163$ExprPhenotype == '64', 'count'] <- realRatio_CD163[realRatio_CD163$ExprPhenotype == '64', 'count'] #+
        #realRatio_CD163[realRatio_CD163$ExprPhenotype == '68', 'count']
      
      realRatio_CD163 <- realRatio_CD163 %>%
        mutate(ratio = count / nrow(adj.CD163))
    }
    
    
    ### Case 2: Macrophage to CD4 T
    
    distMatrix <- flexclust::dist2(posFoxP3[, c('CellXPos', 'CellYPos')], posCD163[, c('CellXPos', 'CellYPos')]) %>%
      as.matrix()
    # These macrophages have at least 1 adjacent CD4T
    adj.FoxP3.id <- which(apply(distMatrix, 1, FUN = min) < 60)
    
    # Get the single-cell data for these Macrophages
    adj.FoxP3 <- posFoxP3[adj.FoxP3.id, ]
    
    
    
    
    if(!(is.empty(adj.FoxP3.id))){
      # Real ratio
      realRatio_FoxP3 <- table(adj.FoxP3$ExprPhenotype) %>%
        data.frame() %>%
        `colnames<-` (c('ExprPhenotype', 'count'))
      
      #realRatio_FoxP3[realRatio_FoxP3$ExprPhenotype == '4', 'count'] <- realRatio_FoxP3[realRatio_FoxP3$ExprPhenotype == '4', 'count'] #+
        #realRatio_FoxP3[realRatio_FoxP3$ExprPhenotype == '68', 'count']
      
      #realRatio_FoxP3[realRatio_FoxP3$ExprPhenotype == '64', 'count'] <- realRatio_FoxP3[realRatio_FoxP3$ExprPhenotype == '64', 'count'] #+
        #realRatio_FoxP3[realRatio_FoxP3$ExprPhenotype == '68', 'count']
      
      realRatio_FoxP3 <- realRatio_FoxP3 %>%
        mutate(ratio = count / nrow(adj.FoxP3))
      
    }
    
    
    
    
    # simulated ratios
    simRatios <- ratio.sim(posCD163, posFoxP3, spatial_polys, 500)
    
    # enrichment score for CD4 T cells clusters
    es.CD163_all <- data.frame(matrix(nrow = 1, ncol = 0))
    for(id in c(2,3)){
      
      #id <- 2
      curExprPhenotype <- realRatio_CD163[id, 'ExprPhenotype'] %>%
        as.character() %>%
        as.numeric()
      # enrichment score
      #es.CD163 <- (length(which(realRatio_CD163[id,'ratio'] > simRatios[[1]][simRatios[[1]]$ExprPhenotype == curExprPhenotype, 'ratio'])) 
      #             + 500 - nrow(simRatios[[1]][simRatios[[1]]$ExprPhenotype == curExprPhenotype,]))/ 500
      
      #es.CD163_all <- cbind.data.frame(es.CD163_all, es.CD163)
      
      # alternative: z-score and p value
      zscore <- (realRatio_CD163[id,'ratio'] - mean(simRatios[[1]][simRatios[[1]]$ExprPhenotype == curExprPhenotype, 'ratio'])) / sd(simRatios[[1]][simRatios[[1]]$ExprPhenotype == curExprPhenotype, 'ratio'])
      
      # p value
      pval <-  2*pnorm(-abs(zscore))
      es.CD163_all <- cbind.data.frame(es.CD163_all, pval)
      
    }
    
    es.CD163_all <- es.CD163_all %>%
      cbind.data.frame(Core, PatientID) %>%
      `colnames<-` (c('es.CD163_PDL1', 'es.CD163_PD1', 'Core', 'PatientID'))
    
    es.CD163_All <- rbind.data.frame(es.CD163_All, es.CD163_all)
    
    # enrichment score for CD4 T cells clusters
    es.FoxP3_all <- data.frame(matrix(nrow = 1, ncol = 0))
    for(id in c(2,3)){
      
      #id <- 2
      curExprPhenotype <- realRatio_FoxP3[id, 'ExprPhenotype'] %>%
        as.character() %>%
        as.numeric()
      # enrichment score
      #es.FoxP3 <- (length(which(realRatio_FoxP3[id,'ratio'] > simRatios[[2]][simRatios[[2]]$ExprPhenotype == curExprPhenotype,'ratio'])) 
      #             + 500 - nrow(simRatios[[2]][simRatios[[2]]$ExprPhenotype == curExprPhenotype,]))/ 500
      
      #es.FoxP3_all <- cbind.data.frame(es.FoxP3_all, es.FoxP3)
      
      # alternative: z-score and p value
      zscore <- (realRatio_FoxP3[id,'ratio'] - mean(simRatios[[2]][simRatios[[2]]$ExprPhenotype == curExprPhenotype, 'ratio'])) / sd(simRatios[[2]][simRatios[[2]]$ExprPhenotype == curExprPhenotype, 'ratio'])
      
      # p value
      pval <- 2*pnorm(-abs(zscore))
      es.FoxP3_all <- cbind.data.frame(es.FoxP3_all, pval)
      
    }
    
    es.FoxP3_all <- es.FoxP3_all %>%
      cbind.data.frame(Core, PatientID) %>%
      `colnames<-` (c('es.FoxP3_PDL1', 'es.FoxP3_PD1', 'Core', 'PatientID'))
    
    es.FoxP3_All <- rbind.data.frame(es.FoxP3_All, es.FoxP3_all)
    
    #print(Core)
    
    
  }
  
  ##--- DO NOT DELETE THIS SECTION ------###
  if(nrow(posTumor) >= 10 & nrow(posFoxP3) >= 10 & Region == 'Central' & Group == '1'){
    
    print(Region)
    
    #Core <- 'M11'
    
    #--------- Get the tissue boundary to do the simulation-------#
    Tissue_all <- data.frame(matrix(nrow = 0, ncol = 0))
    tisse_files <- list.files(paste0('./QuPath/Tissue_Boundary/', Core))
    
    for(piece in tisse_files){
      Tissue <- fromJSON(file = paste0('./QuPath/Tissue_Boundary/', Core, '/', piece))%>%
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
      
      colnames(Tissue) <- c('x', 'y')
      
      # single polgon defined by owin
      
      Tissue_all <- rbind(Tissue_all, cbind(Tissue, piece))
    }
    
    colnames(Tissue_all) <- c('x', 'y', 'group')
    
    # convert polygon data frames to owin objects
    poly_list <- split(Tissue_all, Tissue_all$group)
    # only want lon-lats in the list, not the names
    poly_list <- lapply(poly_list, function(x) { x["group"] <- NULL; x })    
    ps <- lapply(poly_list, Polygon)
    p1 <- lapply(seq_along(ps), function(i) Polygons(list(ps[[i]]), 
                                                     ID = names(poly_list)[i]  ))
    # create SpatialPolygons object
    spatial_polys <- SpatialPolygons(p1) 
    
    owins <- maptools::as.owin.SpatialPolygons(spatial_polys)
    
    #
    #plot(owins)
    ggplot() +
      geom_polygon(aes(Tissue_all[,1], Tissue_all[,2], group = Tissue_all[,3]), fill = NA, color = 'black') +
      geom_point(data = posTumor, aes(CellXPos, CellYPos))
    
    
    ### Case 1: CD4 T cell to Macrophage
    distMatrix <- flexclust::dist2(posTumor[, c('CellXPos', 'CellYPos')], posFoxP3[, c('CellXPos', 'CellYPos')]) %>%
      as.matrix()
    
    # These CD4 T cells has at least 1 adjacent Macrophages
    adj.Tumor.id <- which(apply(distMatrix, 1, FUN = min) < 60)
    
    # Get the single-cell data for these CD4 T cells
    adj.Tumor <- posTumor[adj.Tumor.id, ]
    
    if(!(is.empty(adj.Tumor.id))){
      
      # Real ratio
      
      
      realRatio_Tumor <- table(adj.Tumor$ExprPhenotype) %>%
        data.frame() %>%
        `colnames<-` (c('ExprPhenotype', 'count'))
      
      realRatio_Tumor[realRatio_Tumor$ExprPhenotype == '4', 'count'] <- realRatio_Tumor[realRatio_Tumor$ExprPhenotype == '4', 'count']# +
      #realRatio_Tumor[realRatio_Tumor$ExprPhenotype == '68', 'count']
      
      realRatio_Tumor[realRatio_Tumor$ExprPhenotype == '64', 'count'] <- realRatio_Tumor[realRatio_Tumor$ExprPhenotype == '64', 'count'] #+
      #realRatio_Tumor[realRatio_Tumor$ExprPhenotype == '68', 'count']
      
      realRatio_Tumor <- realRatio_Tumor %>%
        mutate(ratio = count / nrow(adj.Tumor))
    }
    
    
    ### Case 2: Macrophage to CD4 T
    
    distMatrix <- flexclust::dist2(posFoxP3[, c('CellXPos', 'CellYPos')], posTumor[, c('CellXPos', 'CellYPos')]) %>%
      as.matrix()
    # These macrophages have at least 1 adjacent CD4T
    adj.FoxP3.id <- which(apply(distMatrix, 1, FUN = min) < 60)
    
    # Get the single-cell data for these Macrophages
    adj.FoxP3 <- posFoxP3[adj.FoxP3.id, ]
    
    
    
    
    if(!(is.empty(adj.FoxP3.id))){
      # Real ratio
      realRatio_FoxP3 <- table(adj.FoxP3$ExprPhenotype) %>%
        data.frame() %>%
        `colnames<-` (c('ExprPhenotype', 'count'))
      
      realRatio_FoxP3[realRatio_FoxP3$ExprPhenotype == '4', 'count'] <- realRatio_FoxP3[realRatio_FoxP3$ExprPhenotype == '4', 'count'] #+
      #realRatio_FoxP3[realRatio_FoxP3$ExprPhenotype == '68', 'count']
      
      realRatio_FoxP3[realRatio_FoxP3$ExprPhenotype == '64', 'count'] <- realRatio_FoxP3[realRatio_FoxP3$ExprPhenotype == '64', 'count'] #+
      #realRatio_FoxP3[realRatio_FoxP3$ExprPhenotype == '68', 'count']
      
      realRatio_FoxP3 <- realRatio_FoxP3 %>%
        mutate(ratio = count / nrow(adj.FoxP3))
    }
    
    
    
    
    # simulated ratios
    simRatios <- ratio.sim(posTumor, posFoxP3, spatial_polys, 500)
    
    # enrichment score for CD4 T cells clusters
    es.Tumor_all <- data.frame(matrix(nrow = 1, ncol = 0))
    for(id in c(2,3)){
      
      #id <- 2
      curExprPhenotype <- realRatio_Tumor[id, 'ExprPhenotype'] %>%
        as.character() %>%
        as.numeric()
      # enrichment score
      #es.Tumor <- (length(which((realRatio_Tumor[id,'ratio']) > simRatios[[1]][simRatios[[1]]$ExprPhenotype == curExprPhenotype,'ratio'])) 
      #             + 500 - nrow(simRatios[[1]][simRatios[[1]]$ExprPhenotype == curExprPhenotype, ]))/ 500
      
      zscore <- (realRatio_Tumor[id,'ratio'] - mean(simRatios[[1]][simRatios[[1]]$ExprPhenotype == curExprPhenotype, 'ratio'])) / sd(simRatios[[1]][simRatios[[1]]$ExprPhenotype == curExprPhenotype, 'ratio'])
      
      # p value
      pval <-  2*pnorm(-abs(zscore))
      
      
      es.Tumor_all <- cbind.data.frame(es.Tumor_all, zscore)
    }
    
    es.Tumor_all <- es.Tumor_all %>%
      cbind.data.frame(Core, PatientID) %>%
      `colnames<-` (c('es.Tumor_PDL1', 'es.Tumor_PD1', 'Core', 'PatientID'))
    
    es.Tumor_All <- rbind.data.frame(es.Tumor_All, es.Tumor_all)
    
    # enrichment score for CD4 T cells clusters
    es.FoxP3_all <- data.frame(matrix(nrow = 1, ncol = 0))
    for(id in c(2,3)){
      
      #id <- 2
      curExprPhenotype <- realRatio_FoxP3[id, 'ExprPhenotype'] %>%
        as.character() %>%
        as.numeric()
      # enrichment score
      #es.FoxP3 <- (length(which((realRatio_FoxP3[id,'ratio']) > simRatios[[2]][simRatios[[2]]$ExprPhenotype == curExprPhenotype, 'ratio'])) 
      #             + 500 - nrow(simRatios[[2]][simRatios[[2]]$ExprPhenotype == curExprPhenotype,]))/ 500
      
      # alternative: z-score and p value
      zscore <- (realRatio_FoxP3[id,'ratio'] - mean(simRatios[[2]][simRatios[[2]]$ExprPhenotype == curExprPhenotype, 'ratio'])) / sd(simRatios[[2]][simRatios[[2]]$ExprPhenotype == curExprPhenotype, 'ratio'])
      
      # p value
      pval <-  2*pnorm(-abs(zscore))
      
      es.FoxP3_all <- cbind.data.frame(es.FoxP3_all, pval)
    }
    
    es.FoxP3_all <- es.FoxP3_all %>%
      cbind.data.frame(Core, PatientID) %>%
      `colnames<-` (c('es.FoxP3_PDL1', 'es.FoxP3_PD1', 'Core', 'PatientID'))
    
    es.FoxP3_All_st <- rbind.data.frame(es.FoxP3_All_st, es.FoxP3_all)
    
    #print(Core)
    
    
  }
  
}

 
######################################



length(which(es.FoxP3_All$es.FoxP3_PDL1 < 0.01)) / nrow(es.FoxP3_All)
length(which(es.FoxP3_All$es.FoxP3_PD1 < 0.01)) / nrow(es.FoxP3_All)

length(which(es.CD163_All$es.CD163_PDL1 < 0.01)) / nrow(es.CD163_All)
length(which(es.CD163_All$es.CD163_PD1 < 0.01)) / nrow(es.CD163_All)

length(which(es.FoxP3_All_st$es.FoxP3_PDL1 < 0.01)) / nrow(es.FoxP3_All_st)
length(which(es.FoxP3_All_st$es.FoxP3_PD1 < 0.01)) / nrow(es.FoxP3_All_st)

length(which(es.Tumor_All$es.Tumor_PDL1 < 0.01)) / nrow(es.Tumor_All)
length(which(es.Tumor_All$es.Tumor_PD1 < 0.01)) / nrow(es.Tumor_All)


#----------------------Draw Bar plot ---------------#

es.FoxP3_long <- melt(es.FoxP3_All, id.vars = c('Core', 'PatientID'), variable.name = 'marker')


p<-ggplot(data=es.FoxP3_long, aes(x=marker, y= (value+0.000001), color = marker)) +
  theme_bw() +
  theme(legend.position = 'NULL',
        axis.text = element_text(size = 18),
        axis.title = element_blank()) +
  scale_y_log10() +
  geom_violin(aes(color = marker)) +
  geom_jitter(size = 2) +
  geom_hline(yintercept = 0.05, color = 'red') +
  scale_color_manual(values = c('es.FoxP3_PD1' = '#3167b6', 'es.FoxP3_PDL1' = '#f42836')) +
  ylab('p-value') +
  scale_x_discrete(labels=c("es.FoxP3_PD1" = "PD1", "es.FoxP3_PDL1" = "PD-L1"))
p
ggsave(p, file = './Figures/Simulation_results_FoxP3.jpg', width = 4, height = 6, units = "in", dpi = 300)


#------------------ CD163 ---------------------#
es.CD163_long <- melt(es.CD163_All, id.vars = c('Core', 'PatientID'), variable.name = 'marker')


p<-ggplot(data=es.CD163_long, aes(x=marker, y= (value+0.000001), color = marker)) +
  theme_bw() +
  theme(legend.position = 'NULL',
        axis.text = element_text(size = 18),
        axis.title = element_blank()) +
  scale_y_log10() +
  geom_violin(aes(color = marker)) +
  geom_jitter(size = 2) +
  geom_hline(yintercept = 0.05, color = 'red') +
  scale_color_manual(values = c('es.CD163_PD1' = '#3167b6', 'es.CD163_PDL1' = '#f42836')) +
  ylab('p-value') +
  scale_x_discrete(labels=c("es.CD163_PD1" = "PD1", "es.CD163_PDL1" = "PD-L1"))
p
ggsave(p, file = './Figures/Simulation_results_CD163.jpg', width = 4, height = 6, units = "in", dpi = 300)

#length(which(es.Tumor_All$es.Tumor_PDL1 < 0.01)) / nrow(es.Tumor_All)
#   length(which(es.Tumor_All$es.Tumor_PD1 < 0.01)) / nrow(es.Tumor_All)




#---------- CONSTRUCT RISK SCORES ---------#
#-------- Spatial Risk Score incorporates FoxP3, CD163, and Tumor

# SCORE 1: Distance from FoxP3 to CD163 and FoxP3 to Tumor
# SCORE 2: Distance from CD163 to FoxP3 and CD163 to Tumor


# How to combine SCORE 1 and SCORE 2? (Multiply first)
RiskScore1_all <- data.frame(matrix(nrow = 0, ncol = 0))
RiskScore2_all <- data.frame(matrix(nrow = 0, ncol = 0))
RiskScore3_all <- data.frame(matrix(nrow = 0, ncol = 0))


posPD1 <- NSCLCdata %>%
  #dplyr::filter(Phenotype != 'FoxP3CD8') %>%
  #dplyr::filter(Phenotype != 'Other') %>%
  dplyr::filter(Phenotype == 'FoxP3') %>%
  dplyr::filter(ExprPhenotype == 64) %>%
  merge(Core_demo, by = 'Core') %>%
  merge(patient_data_clus, by = 'PatientID') 


tertiles <- quantile(log(posPD1$MeanMembrane650), probs = c(0.33, 0.67))
posPD1$PD1_Bin <- ifelse(posPD1$MeanMembrane650 > tertiles[2], 1, 0)
posPD1$PD1_Bin <- ifelse(posPD1$MeanMembrane650 < tertiles[1], -1, posPD1$PD1_Bin)

table(posPD1$PD1_Bin)

for(Core in unique(NSCLCdata$Core)){
  

  PatientID <- Core_demo[Core_demo$Core == Core, 'PatientID']
  
  # group number
  Group <- patient_data_clus[patient_data_clus$PatientID == PatientID, 'group']
  # Region
  Region <- Core_demo[Core_demo$Core == Core, 'Region']
  
  if(is.empty(Group)){
    Group <- '0'
  }
  
  if(is.empty(PatientID)){
    PatientID <- '0'
  }
  
  if(is.empty(Region)){
    Region <- 'Zero'
  }
  
  # retrieve the core data
  core_data <- NSCLCdata[NSCLCdata$Core == Core, ]
  
  # FoxP3 data
  posFoxP3 <- core_data %>%
    dplyr::filter(Phenotype == 'FoxP3') %>%
    #dplyr::filter(PD1_Bin == '1') %>%
    dplyr::filter(ExprPhenotype == 64) %>%
    select(CellXPos, CellYPos)

  #posFoxP3 <- posFoxP3[posFoxP3$Core == Core, c('CellXPos', 'CellYPos')]
  # CD163 data
  
  posCD163 <- core_data %>%
    dplyr::filter(Phenotype == 'CD163') %>%
    dplyr::filter(ExprPhenotype != 64) %>%
    select(CellXPos, CellYPos)
  
  
  # FoxP3 data
  posTumor <- core_data %>%
    dplyr::filter(Phenotype == 'Tumor') %>%
    dplyr::filter(ExprPhenotype != 64) %>%
    select(CellXPos, CellYPos)
  
  
  #---- Compute Risk Score 1 -------#
  # Computation Method: KDTree
  # FoxP3 to CD163
  if(nrow(posTumor) >= 5 & nrow(posCD163) >= 5 & nrow(posFoxP3) >= 5 & Region == 'Central'){
    
    print(Region)
    
    nd_toFoxP3 <- RANN::nn2(posFoxP3, query = posTumor, k = 5, treetype = 'kd', searchtype = 'priority')[['nn.dists']] %>%
      rowMeans() #%>%
      #mean()
    
    # Tumor to FoxP3 *
    nd_toCD163 <- RANN::nn2(posCD163, query = posTumor, k = 5, treetype = 'kd', searchtype = 'priority')[['nn.dists']] %>%
      rowMeans() #%>%
      #mean()
    
    nd_163toFoxP3 <- RANN::nn2(posFoxP3, query = posCD163, k = 5, treetype = 'kd', searchtype = 'priority')[['nn.dists']] %>%
      rowMeans() #%>%
      #mean()
    
    nd_163toTumor <- RANN::nn2(posTumor, query = posCD163, k = 5, treetype = 'kd', searchtype = 'priority')[['nn.dists']] %>%
      rowMeans() #%>%
      #mean()
    
    nd_FoxP3toCD163 <- RANN::nn2(posCD163, query = posFoxP3, k = 5, treetype = 'kd', searchtype = 'priority')[['nn.dists']] %>%
      rowMeans() #%>%
      #mean()
    
    
    nd_FoxP3toTumor <- RANN::nn2(posTumor, query = posFoxP3, k = 5, treetype = 'kd', searchtype = 'priority')[['nn.dists']] %>%
      rowMeans()  #%>%
      #mean()
    
    
    #---------------------------#

    RiskScore1 <- sqrt(nd_FoxP3toCD163 * nd_FoxP3toTumor)

    RiskScore1_all <- rbind.data.frame(RiskScore1_all, cbind(RiskScore1, PatientID, Core, Group))
    

  }
 
  

}



RiskScore1_all$RiskScore1 <- as.numeric(RiskScore1_all$RiskScore1)
RiskScore2_all$RiskScore2 <- as.numeric(RiskScore2_all$RiskScore2)

#RiskScore1_all$RiskScore1 <- (RiskScore1_all$RiskScore1 - min(RiskScore1_all$RiskScore1)) /
#  (max(RiskScore1_all$RiskScore1) - min(RiskScore1_all$RiskScore1))


g1 <- RiskScore1_all[RiskScore1_all$Group == 1, 'RiskScore1'] %>%
  as.numeric()
mean(g1)

g2 <- RiskScore1_all[RiskScore1_all$Group == 3, 'RiskScore1'] %>%
  as.numeric()
mean(g2)

boxplot(g1, g2)
wilcox.test(g1, g2)

# statistics
#min(eff.scores[eff.scores$class == 'long', 'eff.score'])
#mean(eff.scores[eff.scores$class == 'long', 'eff.score'])
#mean(eff.scores[eff.scores$class == 'short', 'eff.score'])


RiskScore2_all <- RiskScore2_all %>%
  dplyr::filter(Group == 1 | Group == 3)

RiskScore1_all <- RiskScore1_all %>%
  dplyr::filter(Group == 1 | Group == 3)






tgc <- summarySE(RiskScore1_all, measurevar= "RiskScore1", groupvars=c("Group"))
tgc$class <- factor(tgc$Group, level = c('1', '3'))

pd <- position_dodge(0.1) # move them .05 to the left and right
#wilcox.test(eff.scores[eff.scores$class == 'short', 'eff.score'], eff.scores[eff.scores$class == 'long', 'eff.score'])
p <- ggplot(tgc, aes(x= Group, y = RiskScore1, colour=Group, group=Group)) + 
  theme_bw() +
  geom_errorbar(aes(ymin=RiskScore1 - 5*se, ymax=RiskScore1 + 5*se, color = Group), width=.05) +
  geom_point(aes(color = class), fill = 'white', position=pd, size=3, shape=21) + # 21 is filled circle
  
  scale_color_manual(values = c('3' = '#c22828', '1' = '#325698')) +
  #ylab(expression(IL-10^'+' ~ macrophage~risk~score)) +
  expand_limits(y=0) +                        # Expand y range
  ylim(50, 150) +
  geom_signif(comparisons = list(c("1", "3")), 
              map_signif_level=TRUE, size = 0.5, textsize = 10, color = 'black', 
              y_position = 130, annotations = '****'
  ) +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        legend.position = 'none') +
  ylab('Risk Scale')
p
ggsave(p, file = './Figures/Risk Scale_discovery.jpeg', width = 4, height = 5, units = "in", dpi = 300)




tgc <- summarySE(RiskScore2_all, measurevar= "RiskScore2", groupvars=c("Group"))
tgc$class <- factor(tgc$Group, level = c('1', '3'))

pd <- position_dodge(0.1) # move them .05 to the left and right
#wilcox.test(eff.scores[eff.scores$class == 'short', 'eff.score'], eff.scores[eff.scores$class == 'long', 'eff.score'])
p <- ggplot(tgc, aes(x= Group, y = RiskScore2, colour=Group, group=Group)) + 
  theme_bw() +
  geom_errorbar(aes(ymin=RiskScore2 - 20*se, ymax=RiskScore2 + 20*se, color = Group), width=.05) +
  geom_point(aes(color = class), fill = 'white', position=pd, size=3, shape=21) + # 21 is filled circle
  
  scale_color_manual(values = c('short' = '#c22828', 'long' = '#325698')) +
  #ylab(expression(IL-10^'+' ~ macrophage~risk~score)) +
  expand_limits(y=0) +                        # Expand y range
  ylim(0, 1) +
  geom_signif(comparisons = list(c("1", "3")), 
              map_signif_level=TRUE, size = 0.5, textsize = 10, color = 'black', 
              y_position = 0.9, annotations = '***'
  ) +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        legend.position = 'none') 
p






