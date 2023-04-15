#--------------- FUNCTION -------------#
# This script is used for the validation of spatial score
# Unit in Pixel
# 1 Pixel = 0.65 micron

library(spatstat); library(data.table); library(Rfast); library(RANN)
library(dplyr); library(tidyr); library(reshape2); library(tibble); library(tidyverse);
library(ggplot2); 
library(ggforce); library(rstatix)
library(tripack); library(igraph)
library(readxl); library(stringr); library(pracma); library(RANN)
library(survival); library(survminer)
# working directory
current_path <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
setwd('..')
source('./Codes/Function.R')


# read expression data (NOT the cell type)
load('./Validation Cohort/CYCIF_V4_RawZ.RData') %>%
  data.frame()


IFF.datZ_DF <- IFF.datZ %>%
  data.frame() %>%
  select(CD8A, FOXP3, CD163, KERATIN, PD1, PCNA)


# expression transformation, normalization, and standardization

for(e in 1:6){
  
  # asinh transformation
  IFF.datZ_DF[, e] <- asinh(IFF.datZ_DF[, e])
  
  # normalization, std
  IFF.datZ_DF[, e] <- IFF.datZ_DF[, e] / std(IFF.datZ_DF[, e])
  
  # Range standardization
  IFF.datZ_DF[, e] <- (IFF.datZ_DF[, e] - min(IFF.datZ_DF[, e])) / (max(IFF.datZ_DF[, e]) - min(IFF.datZ_DF[, e]))
  
  IFF.datZ_DF[, e] <- ifelse(IFF.datZ_DF[, e] >= 0.5, 1, 0)
  
  
}



# assign SID and CID
ID_lists <- strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\1 \\2', rownames(IFF.datZ_DF)), ' ') 

SID <- ID_lists %>%
  sapply("[[",1)

CID <- ID_lists %>%
  sapply("[[",2)

IFF.datZ_DF$SID <- SID
IFF.datZ_DF$CID <- CID



# read single-cell dataset
load('./Validation Cohort/CYCIF_ALL_CT.RData') 

pheA <- pheA[pheA$L2ct.T != 'UIC' & L2ct.T != 'UC', ]

# assign cell ID
CID_lists <- strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\1 \\2', rownames(pheA)), ' ') 
CID <- CID_lists %>%
  sapply("[[",2)

pheA$CID <- CID



# create a new DF to have cell coordiante and their expression profile
pheA_loc_expr <- merge(pheA, IFF.datZ_DF, by = c('SID', 'CID'))


# patient survival

clinic_protea <- read_table('./Clinic_Protea_190701.txt') %>%
  dplyr::filter(PatientID %in% unique(pheA$PID)) %>%
  select(PatientID, Survival_days, Event) # _death_1_alive_0

# get the survival threshold to determine long versus short-term survivors
thresh <- median(clinic_protea$Survival_days)




clinic_protea$Group <- ifelse(clinic_protea$Survival_days >= thresh, 1, 3)


mean(clinic_protea[clinic_protea$Group == 1,]$Survival_days)
mean(clinic_protea[clinic_protea$Group == 3,]$Survival_days)


fit <- survfit(Surv(Survival_days, Event) ~ Group, data = clinic_protea)
fit


p <- ggsurvplot(fit, data = clinic_protea, 
                palette = c('black', 'blue', 'red'),
                legend = 'none',
                legend.title = '',
                #legend.labs = c('', '),
                risk.table = TRUE,
                surv.scale = 'percent',
                font.tickslab = c(26),
                font.title = c(26),
                font.x = c(28),
                font.y = c(28),
                font.legend = c(16),
                fontsize = 10,
                tables.theme = theme(axis.text = element_text(size = 16),
                                     axis.title = element_text(size = 16),
                                     title = element_text(size = 14)),
                risk.table.y.text = FALSE,
                conf.int = FALSE,
                size = 1.5,
                censor.size = 8,
                pval = FALSE,
                break.time.by = 2000) +
  xlab('Time, (days)')
p
ggsave(file= paste("./Figures/Survival_Validation" , '.pdf'), plot = p$plot, width = 9, height = 6, dpi = 300)



res.cox <- coxph(Surv(Survival_days, Event) ~ Group, data = clinic_protea)
ggforest(res.cox, data = patient_data_clus)

surv_diff <- survdiff(Surv(Survival_days, Event) ~ Group, data = clinic_protea)
surv_diff



#---------------- Cell density comparison --------------#

for(pt in unique(vTMA$PatientID)){
  vTMApt <- vTMA %>%
    filter(PatientID == pt)
  
  # patient single-cell data
  pheApt <- pheA_loc_expr %>%
    filter(PID == pt)
  
  
  
  
  
  vTMA_circle_all <- data.frame(matrix(nrow = 0, ncol = 0))
  
  for(vCore in vTMApt$Core){
    
    
    CoreXY <- vTMApt %>%
      filter(Core == vCore) %>%
      select(X, Y)
    
    
  }
}

#-------------------------------------------------------#
#-----------------Virtual Cores ------------------------#
#-------------------------------------------------------#

# read virtual TMA cores
vTMA <- read_excel('virtualTMAs.xlsx') %>%
  drop_na()

vTMA$X <- vTMA$X / 0.65
vTMA$Y <- vTMA$Y / 0.65

# read survival data
surv <- readxl::read_xlsx('./Validation Cohort/mmc2.xlsx') %>%
  filter(`Patient ID` %in% pheA$PID)

unique(pheA$SID)
unique(pheA[, c('SID', 'PID')])

testDat <- pheA %>%
  filter(PID == 'P132666')



# --------------- Construct Virtual Cores --------------#
# P132666 - Size: 42877 x 25504 Pixels

#### Info ###
#--- Core 1: 
# X: 14162.1995 / 0.65 = 21788
# Y: 4947.7998 / 0.65 = 7612

#--- Core 2:
# X: 14663.3495 / 0.65 = 22559
# Y: 8559.1997 / 0.65 = 13168

#--- Core 3:
# X: 12251.8496 / 0.65 = 18849
# Y: 10008.0496 / 0.65 = 15397




# Create a for-loop to iterate through all patients 
radius <- 1538 / 2 # unit: pixel

# a function to draw circle on point map
circleFun <- function(center = c(),diameter = radius * 2, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

vDF_cell <- data.frame(matrix(nrow = 0, ncol = 0))

for(pt in unique(vTMA$PatientID)){
  
  pt <- 'P133537'
  # loop through each core
  vTMApt <- vTMA %>%
    filter(PatientID == pt)
  
  # patient single-cell data
  pheApt <- pheA_loc_expr %>%
    filter(PID == pt)
  
  
  
  
  
  vTMA_circle_all <- data.frame(matrix(nrow = 0, ncol = 0))
  
  for(vCore in vTMApt$Core){
    
    vCore <- '3'
    # Core X
    CoreXY <- vTMApt %>%
      filter(Core == vCore) %>%
      select(X, Y)
    
    # cells within circumscribed circle
    ccCells <- pheApt[pheApt$X >= (CoreXY$X - radius) & pheApt$X <= (CoreXY$X + radius) & pheApt$Y >= (CoreXY$Y - radius) & pheApt$Y <= (CoreXY$Y + radius),]
    
    # compute distance from all cels to center, discard cells locate outside of the circle
    dist2center <- nn2(data = ccCells[, c('X', 'Y')], query = CoreXY, k = nrow(ccCells), treetype = 'kd',
                       searchtype = 'radius', radius = radius)
    
    
    # get the ranked ids of cells that fall within the circle
    rankedID <- which(dist2center[['nn.dists']] <= radius)
    
    # use the ranked id to get the real cell id
    realID <- dist2center[['nn.idx']][,rankedID]
    
    
    # draw circle
    vTMA_circle <- circleFun(center = c(CoreXY$X, CoreXY$Y), diameter = radius * 2, npoints = 100)
    colnames(vTMA_circle) <- c('X', 'Y')
    
    #pdf('./Figures/P133537_Cores.pdf',width=6,height=6)
    #ggplot(ccCells[realID,], aes(X, Y)) +
    #  theme_void() +
    #geom_point(aes(X, Y, color = L2ct.T), size = 1) +
    #geom_path(data = vTMA_circle, aes(x,y)) +
    #  geom_voronoi(outline = vTMA_circle, aes(X, Y, fill = L2ct.T), color = '#000000', size = 0.1) +
    #ylim(max(ccCells[realID,'Y']) + 100, min(ccCells[realID,'Y']) - 100) +
    #  scale_fill_manual(values = c('CD8+T' = '#9fd5fb', 'Endothelium' = '#745d3c', 'B cell' = '#8da0cb',
    #                                'CD4+T' = '#5288cf', 'Tregs' = '#31a621', 'Tumor' = '#bd5030', 'Mast Cells' = '#66c2a5',
    #                                'Epithelium' = 'grey', 'Mac' = '#c7b0d8')) +
    #  theme(legend.position = 'NA',
    #        plot.background = element_rect(fill = "transparent",colour = NA)) 
    
    #dev.off()
    # assgin core id to circle
    
    #vTMA_circle_all <- rbind.data.frame(vTMA_circle_all, cbind.data.frame(vTMA_circle, vCore))
    
    # store in a data frame
    vDF_cell <- rbind.data.frame(vDF_cell, cbind.data.frame(ccCells[realID,], vCore))
    
    
  }
  
  # plot circles on point cloud
  
  p <- ggplot(pheApt) +
    theme_void() +
    theme(legend.position = 'none',
      rect = element_rect(fill = "transparent"),
      plot.background=element_rect(fill="white", colour=NA)) +
    geom_point(aes(X, Y, color = L2ct.T), size = 1) +
    ylim(max(pheApt$Y), 0) +
   # geom_path(data = vTMA_circle_all, aes(x,y, group = vCore), size = 4) +
    coord_fixed(ratio = 1) +
    scale_color_manual(values = c('CD8+T' = '#9fd5fb', 'Endothelium' = '#745d3c', 'B cell' = '#8da0cb',
                                  'CD4+T' = '#5288cf', 'Tregs' = '#31a621', 'Tumor' = '#bd5030', 'Mast Cells' = '#66c2a5',
                                  'Epithelium' = 'grey', 'Mac' = '#c7b0d8'))
  p
  # define plot size
  #width <- max(pheApt$X) / 1000
  #height <- max(pheApt$Y) / 1000
  
  #jpeg(paste0('./Figures/NSCLC_Havard_vTMA/', pt ,'.jpeg'), width = width, height = height, units = "in", res = 300, bg = "transparent")
  #print(p)
  #dev.off()
  
  #ggsave(p, file = paste0('./Figures/NSCLC_Havard_vTMA/', pt ,'.jpeg'), width = width, height = height, units = "in", dpi = 300)
  
}


#----------------Validation using Treg et to approximate FoxP3+ ---------------#
Q <- data.frame(matrix(nrow = 0, ncol = 0))
for(pt in unique(vDF_cell$PID)){
  
  # patient data
  patientData <- vDF_cell %>%
    dplyr::filter(PID == pt)
  
  # loop through core data
  
  for(cr in unique(patientData$vCore)){
    
    coreData <- patientData %>%
      dplyr::filter(vCore == cr)
    
    
    # 
    qCD8A <- coreData %>%
      dplyr::filter(CD8A == 1 & FOXP3 != 1 & CD163 != 1 & KERATIN != '1') %>%
      nrow()
    
    qCD8T <- coreData %>%
      dplyr::filter(L2ct.T == 'CD8+T') %>%
      nrow()  
    
    # 
    qFoxP3 <- coreData %>%
      dplyr::filter(FOXP3 == 1 & CD8A != 1 & CD163 != 1 & KERATIN != '1') %>%
      nrow()
    
    qTreg <- coreData %>%
      dplyr::filter(L2ct.T == 'Tregs') %>%
      nrow() 
    
    # 
    qCD163 <- coreData %>%
      dplyr::filter(CD163 == 1 & FOXP3 != 1 & CD8A != 1 &  KERATIN != '1') %>%
      nrow()
    
    qMacs <- coreData %>%
      dplyr::filter(L2ct.T == 'Mac') %>%
      nrow() 
    
    # 
    qKeratin <- coreData %>%
      dplyr::filter(KERATIN == 1 & FOXP3 != 1 & CD8A != 1 & CD163 != 1) %>%
      nrow()
    
    qTumor <- coreData %>%
      dplyr::filter(L2ct.T == 'Tumor') %>%
      nrow() 
    
    Q <- rbind.data.frame(Q, cbind.data.frame(qCD8A, qCD8T, qFoxP3, qTreg, qCD163, qMacs, qKeratin, qTumor))
  }
}

p <- ggplot(Q, aes(x = qCD163, y = qMacs)) +
  geom_point(size = 3) +
  theme_classic() +
  theme(axis.text = element_text(size = 22),
        axis.title = element_text(size = 24)) +
  geom_smooth(method='lm', se = FALSE, color = 'black') +
  xlab('CD8-FoxP3—CD163+Keratin-') +
  coord_fixed(ratio = 1) +
  xlim(0,1000) +
  ylab('Macrophages')

lm(Q$qMacs ~ Q$qCD163)
cor.test(Q$qMacs, Q$qCD163)

ggsave(p, file = paste0('./Figures/CD163_cor.jpeg'), width = 6, height = 6,  units = "in", dpi = 300)

p <- ggplot(Q, aes(x = qKeratin, y = qTumor)) +
  geom_point(size = 3) +
  theme_classic() +
  theme(axis.text = element_text(size = 22),
        axis.title = element_text(size = 24)) +
  geom_smooth(method='lm', se = FALSE, color = 'black') +
  xlab('CD8-FoxP3-CD163-Keratin+') +
  coord_fixed(ratio = 1) +
  xlim(0,4500) +
  ylab('Tumor cell')

lm(Q$qTumor ~ Q$qKeratin)
cor.test(Q$qKeratin, Q$qTumor)

ggsave(p, file = paste0('./Figures/Tumo_cor.jpeg'), width = 6, height = 6,  units = "in", dpi = 300)


p <- ggplot(Q, aes(x = qCD8A, y = qCD8T)) +
  geom_point(size = 3) +
  theme_classic() +
  theme(axis.text = element_text(size = 22),
        axis.title = element_text(size = 24)) +
  geom_smooth(method='lm', se = FALSE, color = 'black') +
  coord_fixed(ratio = 1) +
  xlab('CD8+FoxP3-CD163-Keratin-') +
  xlim(0,2000) +
  ylab('CD8 T')


lm(Q$qCD8T ~ Q$qCD8A)
cor.test(Q$qCD8A, Q$qCD8T)

ggsave(p, file = paste0('./Figures/CD8_cor.jpeg'), width = 6, height = 6,  units = "in", dpi = 300)

p <- ggplot(Q, aes(x = qFoxP3, y = qTreg)) +
  geom_point(size = 3) +
  theme_classic() +
  theme(axis.text = element_text(size = 22),
        axis.title = element_text(size = 24)) +
  geom_smooth(method='lm', se = FALSE, color = 'black') +
  coord_fixed(ratio = 1) +
  xlim(0,700) +
  xlab('CD8-FoxP3+CD163-Keratin-') +
  ylab('Treg')

lm(Q$qFoxP3 ~ Q$qTreg)
cor.test(Q$qFoxP3, Q$qTreg)

ggsave(p, file = paste0('./Figures/FoxP3_cor.jpeg'), width = 6, height = 6,  units = "in", dpi = 300)







densities_validation_all <- data.frame(matrix(nrow = 0, ncol = 0))

for(pt in unique(vDF_cell$PID)){
  
  # patient data
  patientData <- vDF_cell %>%
    dplyr::filter(PID == pt)
  
  # loop through core data
  
  for(cr in unique(patientData$vCore)){
    
    coreData <- patientData %>%
      dplyr::filter(vCore == cr)
    
    colnames(coreData)[7] <- 'PatientID' 
    # get the group
    
    Group <- coreData %>%
      merge(clinic_protea) %>%
      select(Group) %>%
      unique() %>%
      as.numeric()
    
    # densities
    posCD8 <- coreData %>%
      dplyr::filter(L2ct.T == 'CD8+T') 
    
    posFoxP3 <- coreData %>%
      dplyr::filter(L2ct.T == 'Tregs')
    
    # CD163 data
    posCD163 <- coreData %>%
      dplyr::filter(L2ct.T == 'Mac')
    
    # FoxP3 data
    posTumor <- coreData %>%
      dplyr::filter(L1ct.T == 'Tumor') 
    
    
    CD8 <- nrow(posCD8) / (pi*(radius/1000)^2)
    CD163 <- nrow(posCD163) / (pi*(radius/1000)^2)
    FoxP3 <- nrow(posFoxP3) / (pi*(radius/1000)^2)
    Tumor <- nrow(posTumor) / (pi*(radius/1000)^2)
    
    
    # combine all data
    densities_validation_all <- rbind.data.frame(densities_validation_all, 
                                                 cbind.data.frame(CD8, CD163, FoxP3,Tumor, cr, pt, Group))
  }
  
  
  
}


# wide table to long
for(pt in unique(densities_validation_all$pt)){
  print(densities_validation_all[densities_validation_all$pt == pt, ])
  print('CD8')
  print(mean(densities_validation_all[densities_validation_all$pt == pt, 1]))
  print(sd(densities_validation_all[densities_validation_all$pt == pt, 1]))
  
  print('CD163')
  print(mean(densities_validation_all[densities_validation_all$pt == pt, 2]))
  print(sd(densities_validation_all[densities_validation_all$pt == pt, 2]))
  
  print('Tumor')
  print(mean(densities_validation_all[densities_validation_all$pt == pt, 4]))
  print(sd(densities_validation_all[densities_validation_all$pt == pt, 4]))
  
  print('FoxP3')
  print(mean(densities_validation_all[densities_validation_all$pt == pt, 3]))
  print(sd(densities_validation_all[densities_validation_all$pt == pt, 3]))
  

}





densities_validation_all <- melt(densities_validation_all[, c('CD8', 'CD163', 'FoxP3', 'Tumor', 'pt', 'Group')], id.vars=c("pt", "Group"))
colnames(densities_validation_all)[3] <- 'ctype' 

stat.test <- densities_validation_all %>% 
  group_by(ctype) %>%
  wilcox_test(value ~ Group) %>% 
  add_significance() %>%
  add_xy_position(x = "Group")

supp.labs <- c('CD8', "Macs", "Tregs", 'Tumor')
names(supp.labs) <- c("CD8", "CD163", 'FoxP3', 'Tumor')

p <- ggplot(data = densities_validation_all, aes(x = as.factor(Group), y = value, color = as.factor(Group))) +
  stat_boxplot(geom ='errorbar', width = 0.4) +
  geom_boxplot(width = 0.6, fill = "white") +
  theme_bw() +
  stat_pvalue_manual(stat.test, size = 6) +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        #axis.text.y = element_text(angle = 90),
        strip.text = element_text(size = 20),
        legend.position = 'none',
        strip.background = element_rect(fill ='white')) +
  ylab('') +
  #ylim(0, 1.5) +
  #scale_y_continuous(breaks = seq(0, 1.5, by = 0.5)) +
  facet_grid(. ~ ctype, labeller = labeller(ctype = supp.labs)) +
  scale_color_manual(values = c('1' = '#0a377b', '3' = "#d63434")) +
  ylab(expression('Density (No.' ~ mm ^ -2 ~')')) 
p
ggsave(p, file=paste0("./Figures/Density_comparison_validation.jpeg"), width = 10, height = 5, units = "in", dpi = 300)




#------------------------------------------------#
# ----------- Network Analysis ------------------#
#------------------------------------------------#
degrees_all <- data.frame(matrix(nrow = 0, ncol = 0))
betweenness_all <- data.frame(matrix(nrow = 0, ncol = 0))
closeness_all <- data.frame(matrix(nrow = 0, ncol = 0))
transitivities_all <- data.frame(matrix(nrow = 0, ncol = 0))
assort_df <- data.frame(matrix(nrow = 0, ncol = 6))
colnames(assort_df) <- c(c('CD8+T', 'Mac', 'Tregs', 'Tumor'), 'Group', 'Core')
for(pt in unique(vDF_cell$PID)){
  
  # patient data
  patientData <- vDF_cell %>%
    dplyr::filter(PID == pt)
  
  # loop through core data
  
  for(cr in unique(patientData$vCore)){
    #cr <- 1
    coreData <- patientData %>%
      dplyr::filter(vCore == cr) %>%
      dplyr::filter(L2ct.T == 'Mac' | L2ct.T == 'CD8+T' | L2ct.T == 'Tregs' | L2ct.T == 'Tumor')
    
    # get the group
    colnames(coreData)[7] <- 'PatientID'
    Group <- coreData %>%
      merge(clinic_protea) %>%
      select(Group) %>%
      unique() %>%
      as.numeric()
    
    lists <- Network_val(coreData[, c('X', 'Y', 'L2ct.T', 'CID')], 60)
    Dual_EdgeList <- lists[[1]]
    Dual_NodeList <- lists[[2]]
    
    
    
    Dual_EdgeList$from <- Dual_NodeList[Dual_EdgeList$from, 'CellID']
    Dual_EdgeList$to <- Dual_NodeList[Dual_EdgeList$to, 'CellID']
    
    # split Dual_Edgelist to two part, part 1 contains connected nodes of same kinds
    
    # ------ part 1 -----#
    
    Dual_EdgeList_types <- data.frame(from = coreData[match(Dual_EdgeList$from, coreData$CID), 'L2ct.T'],
                                      to = coreData[match(Dual_EdgeList$to, coreData$CID), 'L2ct.T'])
    
    
    Dual_EdgeList_pt1 <- Dual_EdgeList[which(Dual_EdgeList_types$from == Dual_EdgeList_types$to), ]
    Dual_EdgeList_pt1_types <- Dual_EdgeList_types[which(Dual_EdgeList_types$from == Dual_EdgeList_types$to), ]
    
    # ------ part 2 -----#
    
    Dual_EdgeList_pt2 <- Dual_EdgeList[which(Dual_EdgeList_types$from != Dual_EdgeList_types$to), ]
    Dual_EdgeList_pt2_types <- Dual_EdgeList_types[which(Dual_EdgeList_types$from == Dual_EdgeList_types$to), ]
    
    
    # single-cell data IN the network
    #Node_singlecell <- coreData[Dual_NodeList$CellID %in% unique(),]
    Node_singlecell <- Dual_NodeList[Dual_NodeList$CellID %in% unique(c(
      as.vector(Dual_EdgeList$from), 
      as.vector(Dual_EdgeList$to))),]
    
    
    p <- ggplot() +
      #theme_bw() +
      theme_void() +
      #geom_polygon(data = mask, aes(Y, X, group = cellLabelInImage, fill = immuneGroup)) +
      geom_segment(aes(x = Dual_NodeList[match(Dual_EdgeList_pt1$from, Dual_NodeList$CellID), 2],
                       xend = Dual_NodeList[match(Dual_EdgeList_pt1$to, Dual_NodeList$CellID), 2],
                       y = Dual_NodeList[match(Dual_EdgeList_pt1$from, Dual_NodeList$CellID), 3],
                       yend = Dual_NodeList[match(Dual_EdgeList_pt1$to, Dual_NodeList$CellID), 3],
                       color = Dual_EdgeList_pt1_types$from), size = 0.4) +
      geom_segment(aes(x = Dual_NodeList[match(Dual_EdgeList_pt2$from, Dual_NodeList$CellID), 2],
                       xend = Dual_NodeList[match(Dual_EdgeList_pt2$to, Dual_NodeList$CellID), 2],
                       y = Dual_NodeList[match(Dual_EdgeList_pt2$from, Dual_NodeList$CellID), 3],
                       yend = Dual_NodeList[match(Dual_EdgeList_pt2$to, Dual_NodeList$CellID), 3],
                       color = '#b4b3b3'), size = 0.4) +
      geom_point(data = Node_singlecell, aes(x, y, color = L2ct.T), size = 1) +
      scale_color_manual(values = c('CD8+T' = '#a6cee3', 'Tregs' = '#33a02c', 'Mac' = '#cab2d6',
                                    'Tumor' = '#e31b1c')) + #f8f8f8
      theme(legend.position = 'none') +
      #ylim(670, 1398) +
      #xlim(830, 1558) +
      #ylim(700, 1428) + # 670 - 1398 (Patient 16)
      #xlim(1200, 1928) + # 830 - 1558 (Patient 16)
      coord_fixed(ratio = 1)
    p
    print(p)
    ggsave(p, file=paste0("./Figures/Spatial_Clustering_validations/Patient_", pt, '_', cr, '.jpeg'), 
           bg = 'white',width = 8, height = 8, units = "in", dpi = 300)
    
    # ------------------ #
    #      iGraph        #
    # ------------------ #
    
    
    #---------------- Compute centrality coefficients -------------#
    ig <- graph_from_data_frame(vertices = Node_singlecell[,c('CellID', 'x', 'y', 'L2ct.T')], d = Dual_EdgeList, directed = FALSE)
    
    ig <- set_vertex_attr(ig, 'cell_type', value = Node_singlecell$L2ct.T)  
    
    
    for(ctype in c('CD8+T', 'Mac', 'Tregs', 'Tumor')){
      #ctype <- 'Mac'
      degrees <- igraph::degree(ig, v = V(ig), normalized = TRUE)%>%
        as.data.frame() %>%
        cbind(Node_singlecell$L2ct.T) %>%
        `colnames<-` (c('degrees', 'Phenotype')) %>%
        dplyr::filter(Phenotype == ctype)
      
      
      closenesses <- closeness(ig, vids = V(ig), weights = NULL, normalized = TRUE)%>%
        as.data.frame() %>%
        cbind(Node_singlecell$L2ct.T) %>%
        `colnames<-` (c('closenesses', 'Phenotype')) %>%
        dplyr::filter(Phenotype == ctype)
      
      
      betweenesses <- betweenness(ig, v = V(ig), normalized = TRUE)%>%
        as.data.frame() %>%
        cbind(Node_singlecell$L2ct.T) %>%
        `colnames<-` (c('betweennesses', 'Phenotype')) %>%
        dplyr::filter(Phenotype == ctype)
      
      
      #transitivities <- transitivity(ig, type = 'localundirected', vids = V(ig)[V(ig)$`cell_type` == ctype])
      transitivities <- transitivity(ig, type = 'localundirected', vids = V(ig)) %>%
        as.data.frame() %>%
        cbind(Node_singlecell$L2ct.T) %>%
        `colnames<-` (c('transitivities', 'Phenotype')) %>%
        dplyr::filter(Phenotype == ctype)
      
      degrees_all <- rbind.data.frame(degrees_all, cbind.data.frame(degrees, cr, Group, pt))
      closeness_all <- rbind.data.frame(closeness_all, cbind.data.frame(closenesses, cr, Group, pt))
      betweenness_all <- rbind.data.frame(betweenness_all, cbind.data.frame(betweenesses, cr, Group, pt))
      transitivities_all <- rbind.data.frame(transitivities_all, cbind.data.frame(transitivities, cr, Group, pt))
      
    }
    
    
    # construct iGraph object
    ig <- graph_from_data_frame(vertices = Node_singlecell[,c('CellID', 'x', 'y', 'L2ct.T')], d = Dual_EdgeList, directed = FALSE)
    
    # assign vertex attributes
    # For this part, we set the objective cell type as type 1 and all other cell types as type 2
    assort_vec <- data.frame(matrix(nrow = 0, ncol = 4))
    colnames(assort_vec) <- c('CD8+T', 'Mac', 'Tregs', 'Tumor')
    
    for(type1 in c('CD8+T', 'Mac', 'Tregs', 'Tumor')){
      
      type1_id <- which(Node_singlecell$L2ct.T == type1)
      vertex_attr <- rep(2, nrow(Node_singlecell))
      
      vertex_attr[type1_id] <- 1
      ig <- set_vertex_attr(ig, 'cell_type', value = vertex_attr)  
      
      print(assortativity_nominal(ig, types = V(ig)$cell_type))
      
      assort_vec[1, type1] <- assortativity_nominal(ig, types = V(ig)$cell_type)
    }
    
    
    assort_vec$Group <- Group
    assort_vec$Core <- cr
    assort_df <- rbind.data.frame(assort_df, assort_vec)
    
    
  }
  
  
  
}

#--------------------------------------------------------------------------#
#-------------- Visualizations for network measurements -------------------#
#--------------------------------------------------------------------------#
library(hrbrthemes)

degrees_all_vis <- degrees_all


pd = position_dodge(width = 0.75)

signif_pairs <- degrees_all_vis %>%
  select(-c('cr', 'pt')) %>%
  melt(id.vars = c('Phenotype', 'Group')) #%>%
#select(-'core')

stat.test <- signif_pairs %>%
  group_by(Phenotype) %>%
  wilcox_test(value ~ Group) %>%
  add_significance() %>%
  add_xy_position(x = "Phenotype")

stat.test$y.position <- log(stat.test$y.position) - 2.5 # - closeness

p <- ggplot(degrees_all_vis, aes(x= Phenotype, y= log(degrees + 0.000001), color=as.factor(Group))) + 
  stat_boxplot(geom = "errorbar", width = 0.25, position=pd) + 
  scale_color_manual(values = c('1' = '#0a377b', '3' = "#d63434")) +
  #stat_pvalue_manual(stat.test, size = 5) +
  rremove("x.title") + rremove('y.title') +
  
  geom_boxplot(outlier.color = NA) +
  #ylim(-9, -4) + # degree centrality
  #ylim(2, 20) + # betweenness centrality
  #ylim(-2, 1) + # transitivity centrality
  #ylim(-5, -2) + # transitivity centrality
  theme_classic() +
  #theme(strip.text = element_text(size = 15)) +
  #ylim(800, 1300) +
  font("xy.text", size = 20) +
  font("xy.text", size = 20) +
  rremove('legend') 
p
ggsave(p, file=paste0("./Figures/closenesses_centrality.jpeg"), width = 4, height = 5, units = "in", dpi = 300)


#--------------------- Shannon Index ------------------------#

ShannonIndex_all <- data.frame(matrix(nrow = 0, ncol = 0))
for(pt in unique(vDF_cell$PID)){
  
  # patient data
  patientData <- vDF_cell %>%
    dplyr::filter(PID == pt)
  
  # loop through core data
  
  for(cr in unique(patientData$vCore)){
    
    # coordinate data
    
    coreData <- patientData %>%
      dplyr::filter(vCore == cr)
    
    colnames(coreData)[7] <- 'PatientID' 
    # get the group
    
    Group <- coreData %>%
      merge(clinic_protea) %>%
      select(Group) %>%
      unique() %>%
      as.numeric()
    # Survival group
    # Compute Shannon Index
    
    
    cell_types <- c('CD8+T', 'Tregs', 'Tumor', 'Mac')
    
    # a modified function is used since the data formats are different
    ShannonIndex <- ShannonE_val(cell_types, coreData)
    
    # 
    ShannonIndex_all <- rbind.data.frame(ShannonIndex_all, cbind.data.frame(ShannonIndex, cr, Group, pt))
  }
}


# Barplot for specific core
dataF <- cur_df %>%
  dplyr::filter(core == 'E3') %>%
  select(c('core', 'CD8', 'FoxP3', 'CD163', 'Tumor')) %>%
  melt(id.vars = c("Core"), variable.name = "ctype")



p<-ggplot(data=dataF, aes(x=ctype, y=value, fill = ctype, width = .9)) +
  theme_classic() +
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size = 32),
        axis.title.x = element_blank(),
        #axis.text.x = element_text(angle = 90, vjust = 0.6, hjust=0.6),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = 'na') +
  ylab('Nucleus type fraction') +
  geom_bar(stat="identity") +
  scale_y_continuous(breaks=c(0, 0.1, 0.2)) +
  ylim(0, 0.4) +
  scale_fill_manual(values = c('CD8' = '#a6cee3', 'FoxP3' = '#33a02c', 'CD163' = '#cab2d6','Tumor' = '#e31b1c')) 
p
ggsave(p, file="./Figures/E13_high_ShannonIndex.png", width=4, height=6)





stat.test <- ShannonIndex_all %>% 
  wilcox_test(ShannonIndex ~ Group) %>% 
  add_significance() %>%
  add_xy_position(x = "Group")

ShannonIndex_all$title <- 'Shannon Entropy '

stat.test$y.position <- stat.test$y.position + 0.5 # - closeness

p <- ggplot(ShannonIndex_all, aes(x = as.factor(Group), y = ShannonIndex, colour=as.factor(Group)), fill = NA) + 
  stat_boxplot(geom = "errorbar", width = 0.25) + 
  stat_pvalue_manual(stat.test, size = 8) + rremove("x.title") + rremove('y.title') +
  scale_color_manual(values = c('1' = '#0a377b', '3' = "#d63434")) +
  theme_bw() +
  facet_grid(.~title) +
  geom_boxplot() +
  #geom_jitter(size =3) +
  theme(axis.text = element_text(size = 22),
        axis.title = element_text(size = 24),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size = 24),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        #axis.text.x = element_text(angle = 90, vjust = 0.6, hjust=0.6),
        #axis.text.x = element_blank(),
        strip.background = element_rect(fill ='white'),
        legend.position = 'na'
  ) +
  #ylab(expression('Density (No.' ~ mm ^ -2 ~')')) +
  ylab('') +
  ylim(0, 3) 
p

ggsave(file= paste("./Figures/Shannon_Index——validation" , '.png'), plot = print(p), width = 4, height = 6, dpi = 300)



wilcox.test(ShannonIndex_all[ShannonIndex_all$group == '1', 'ShannonIndex'], ShannonIndex_all[ShannonIndex_all$group == '3', 'ShannonIndex'])


#-----------------------------------------------------------------#
#--------------- Ratio between cell types ------------------------#
#-----------------------------------------------------------------#

densities_validation_all <- melt(densities_validation_all[, c('CD8', 'CD163', 'FoxP3', 'Tumor', 'pt', 'cr', 'Group')], id.vars=c("pt", "cr", "Group"))

colnames(densities_validation_all)[1] <- 'PatientID'
colnames(densities_validation_all)[2] <- 'Core'
colnames(densities_validation_all)[4] <- 'ctype' 
colnames(densities_validation_all)[5] <- 'Density' 

cl_profile <- densities_validation_all %>%
  select(c('Group', 'Core', 'PatientID', 'ctype', 'Density')) %>%
  
  reshape(idvar = c('Group', 'PatientID', 'Core'), timevar = "ctype", direction = "wide") %>%
  `colnames<-` (c('Group', 'Core', 'PatientID', 
                  'CD8',
                  'CD163',
                  'FoxP3',
                  'Tumor'
  )) %>%
  mutate(`CD8/Tumor` = CD8/Tumor) %>%
  mutate(`CD8/FoxP3` = CD8/FoxP3) %>%
  mutate(`CD8/CD163` = CD8/CD163) %>%
  mutate(`FoxP3/CD163` = FoxP3/CD163) %>%
  mutate(`FoxP3/Tumor` = FoxP3/Tumor) %>%
  mutate(`CD163/Tumor` = CD163/Tumor) %>%
  select(c('Group', 'Core', 'PatientID', 'CD8/Tumor', 'CD8/FoxP3', 'CD8/CD163', 'FoxP3/CD163',
           'FoxP3/Tumor', 'CD163/Tumor')) %>%
  melt(id.vars = c("Core","PatientID", "Group"), variable.name = "ctype")




stat.test <- cl_profile %>% 
  group_by(ctype) %>%
  wilcox_test(value ~ Group) %>% 
  add_significance() %>%
  add_xy_position(x = "Group")

p <- ggplot(data = cl_profile, aes(x = as.factor(Group), y = value, color = as.factor(Group))) +
  stat_boxplot(geom ='errorbar', width = 0.4) +
  geom_boxplot(width = 0.6, fill = "white") +
  theme_bw() +
  stat_pvalue_manual(stat.test, size = 6) +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        #axis.text.y = element_text(angle = 90),
        strip.text = element_text(size = 20),
        legend.position = 'none',
        strip.background = element_rect(fill ='white')) +
  ylab('') +
  #ylim(0, 1.5) +
  #scale_y_continuous(breaks = seq(0, 1.5, by = 0.5)) +
  facet_grid(. ~ ctype) +
  scale_color_manual(values = c('1' = '#0a377b', '3' = "#d63434")) +
  ylab(expression('Ratio')) 
p
ggsave(p, file=paste0("./Figures/Ratio_validation.jpeg"), width = 14, height = 5, units = "in", dpi = 300)












# compute the score
RiskScore1_all <- data.frame(matrix(nrow = 0, ncol = 0))

#RiskScore2_all <- data.frame(matrix(nrow = 0, ncol = 0))

for(pt in unique(vTMA$PatientID)){
  
  #pt <- 'P132115'
  # loop through each core
  pt_cellpos <- vDF_cell %>%
    filter(PID == pt) #%>%
  #dplyr::filter(!(FOXP3 == '1' & CD8A == '1'))
  
  
  
  for(Core in unique(pt_cellpos$vCore)){
    
    colnames(pt_cellpos)[7] <- 'PatientID'
    Group <- pt_cellpos %>%
      merge(clinic_protea) %>%
      select(Group) %>%
      unique() %>%
      as.numeric()
    
    #Core <- 1
    # get the single-cell data for each cell populations
    posFoxP3 <- pt_cellpos %>%
      dplyr::filter(vCore == Core) %>%
      dplyr::filter(L2ct.T == 'Tregs') %>%
      #dplyr::filter(FOXP3 == '1') %>%
      #dplyr::filter(CD163 != '1') %>%
      dplyr::filter(PD1 == '1') %>%
      #dplyr::filter(KERATIN != '1') %>%
      #dplyr::filter(L1ct.T != 'Tumor') %>%
      select(X, Y)
    
    # CD163 data
    posCD163 <- pt_cellpos %>%
      dplyr::filter(vCore == Core) %>%
      dplyr::filter(L2ct.T == 'Mac') %>%
      #dplyr::filter(CD163 == '1') %>%
      #dplyr::filter(CD8A != '1') %>%
      #dplyr::filter(FOXP3 != '1') %>%
      dplyr::filter(PD1 == '0') %>%
      #dplyr::filter(KERATIN != '1') %>%
      
      #dplyr::filter(L1ct.T != 'Tumor') %>%
      
      select(X, Y)
    
    
    # FoxP3 data
    posTumor <- pt_cellpos %>%
      dplyr::filter(vCore == Core) %>%
      dplyr::filter(PD1 == '0') %>%
      #dplyr::filter(CD163 == '0') %>%
      #dplyr::filter(KERATIN == '1') %>%
      #dplyr::filter(CD8A != '1') %>%
      #dplyr::filter(FOXP3 != '1') %>%
      
      dplyr::filter(L1ct.T == 'Tumor') %>%
      #dplyr::filter(KERATIN == '1') %>%
      # dplyr::filter(L1ct == 'Epithelium') %>%
      select(X, Y)
    
    
    
    
    
    
    #---- Compute Risk Score 1 -------#
    # Computation Method: KDTree
    # FoxP3 to CD163
    
    tryCatch(
      expr = {
        if(nrow(posTumor) >= 5 & nrow(posCD163) >= 5 & nrow(posFoxP3) >= 5){
          
          
          # Observation 1: when FoxP3 closer to CD163, favors better survival ---- A
          #                when FoxP3 closer to Tumor, favors poor survival   ---- Bíí
          
          # to quantify A
          #nd_toCD163 <- RANN::nn2(posCD163, query = posFoxP3, k = 5, treetype = 'kd', searchtype = 'priority')[['nn.dists']] %>%
          #  rowMeans()
          # Tumor to FoxP3 *
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
          
          # Tumor to FoxP3 *
          #nd_toTumor <- RANN::nn2(posFoxP3, query = posCD163, k = 5, treetype = 'kd', searchtype = 'priority')[['nn.dists']] %>%
          #  rowMeans()
          
          #RiskScore1 <- (nd_FoxP3toCD163 / (nd_FoxP3toCD163 + nd_FoxP3toTumor)) %>%
          #  as.numeric()
          
          #RiskScore1 <- (nd_163toFoxP3 / (nd_163toFoxP3 + nd_163toTumor)) %>%
          #  as.numeric()
          
          #RiskScore1 <- (nd_toCD163 / (nd_toCD163 + nd_toFoxP3)) %>%
          #  as.numeric()
          
          RiskScore1 <- sqrt(nd_FoxP3toCD163 * nd_FoxP3toTumor) %>%
            as.numeric()
          
          #RiskScore2 <- (nd_FoxP3toCD163 / (nd_FoxP3toCD163 + nd_FoxP3toTumor)) %>%
          #  as.numeric()
          
          #RiskScore1 <- sqrt(RiskScore1 * RiskScore2)
          # Biomarker, defined as RiskScore 1 * RisckScore 2
          RiskScore1_all <- rbind.data.frame(RiskScore1_all, cbind(RiskScore1, pt, Core, Group))
          RiskScore2_all <- rbind.data.frame(RiskScore2_all, cbind(RiskScore2, pt, Core, Group))
          #RiskScore3_all <- rbind.data.frame(RiskScore3_all, cbind(RiskScore3, PatientID, Core, Group))
          
        }
      },
      error = function(e){
        
      }
    )
    
    
    
  }
  
}

colnames(RiskScore1_all) <- c('RiskScore', 'PatientID', 'Core', 'Group')
#colnames(RiskScore2_all) <- c('RiskScore', 'PatientID', 'Core')

colnames(clinic_protea) <- c('PatientID', 'Survival_days', 'Event', 'Group')

RiskScore1_all$RiskScore <- as.numeric(RiskScore1_all$RiskScore)
testScore <- merge(RiskScore1_all, clinic_protea, by = 'PatientID')

#RiskScore2_all$RiskScore <- as.numeric(RiskScore2_all$RiskScore)
#testScore <- merge(RiskScore2_all, clinic_protea, by = 'PatientID')

longData <- testScore %>%
  dplyr::filter(Group == 1) %>%
  group_by(PatientID) #%>%
#summarise_all(mean)


shortData <- testScore %>%
  dplyr::filter(Group == 3) %>%
  group_by(PatientID) #%>%
#summarise_all(mean)



tgc <- summarySE(RiskScore1_all, measurevar= "RiskScore", groupvars=c("Group"))
tgc$class <- factor(tgc$Group, level = c('1', '3'))

pd <- position_dodge(0.1) # move them .05 to the left and right
#wilcox.test(eff.scores[eff.scores$class == 'short', 'eff.score'], eff.scores[eff.scores$class == 'long', 'eff.score'])
p <- ggplot(tgc, aes(x= Group, y = RiskScore, colour=Group, group=Group)) + 
  theme_bw() +
  geom_errorbar(aes(ymin=RiskScore - 5*se, ymax=RiskScore + 5*se, color = Group), width=.05) +
  geom_point(aes(color = class), fill = 'white', position=pd, size=3, shape=21) + # 21 is filled circle
  
  scale_color_manual(values = c('3' = '#c22828', '1' = '#325698')) +
  #ylab(expression(IL-10^'+' ~ macrophage~risk~score)) +
  expand_limits(y=0) +                        # Expand y range
  ylim(100, 200) +
  geom_signif(comparisons = list(c("1", "3")), 
              map_signif_level=TRUE, size = 0.5, textsize = 10, color = 'black', 
              y_position = 180, annotations = '****'
  ) +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        legend.position = 'none') +
  ylab('Risk Scale')
p
ggsave(p, file = './Figures/Risk Scale_validation.jpeg', width = 4, height = 5, units = "in", dpi = 300)







