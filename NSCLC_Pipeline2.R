library(spatstat); library(Rfast2)
library(dplyr); library(tidyr); library(reshape2); library(tibble); library(tidyverse)
library(ggplot2); library(umap); library(ggvoronoi); library(Rtsne); library(ComplexHeatmap)
library(ggforce); library(UpSetR); library(rstatix)
library(dendextend)
library(readxl); library(stringr); 
library(survival); library(survminer)
# working directory
current_path <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))

setwd('..')

source('./Codes/Function.R')

#---------------------- DATA INPUT ----------------#
# read single-cell data
NSCLCdata <- readRDS('./NSCLCdataset.rds')
colnames(NSCLCdata)[1] <- 'Core'

# read patient data
patient_data <- read.csv('./NSCLC Clinical Data.csv')
# core demographic data
Core_demo <- read.csv('./Clinicopathology parameters.csv')
# tissue demographic data
Tissue_demo <- read.csv('./Tissue_demographics.csv')

#--------------------------------------------------#
#

# PREREQUISITE DATA: INHERITED FROM PIPELINE 1
#-------------- Cell type fraction per image with hierarchical clustering ------------#
# Create table with celltype fractions

cur_df <- NSCLCdata %>%
  # remove control cores
  dplyr::filter(Core %in% Core_demo$Core) %>%
  group_by(Core, Phenotype) %>%
  dplyr::summarise(n=n()) %>% 
  group_by(Core) %>%
  mutate(fraction = n / sum(n)) %>%
  reshape2::dcast(Core ~ Phenotype, value.var = "fraction", fill=0)


matrixrownames <- cur_df$Core 

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



#-------------- Compare the cell density for each type between two groups ------------#

# Cluster 1 data
quant_profiles_all <- readRDS('./Nucleus_Density_Core.rds')


cl_profile <- cl_tree %>%
  merge(patient_data_clus, by = 'PatientID') %>%
  merge(quant_profiles_all, by = 'Core') %>%
  select(c('group', 'Core', 'PatientID', 'ctype', 'density')) %>%
  dplyr::filter(group != 2) %>%
  reshape(idvar = c('group', 'PatientID', 'Core'), timevar = "ctype", direction = "wide") %>%
  `colnames<-` (c('Group', 'Core', 'PatientID', 
                  'CD163',
                  'FoxP3',
                  'Tumor',
                  'Other',
                  'CD8',
                  'FoxP3CD8'
  )) %>%
  mutate(`CD8/Tumor` = CD8/Tumor) %>%
  mutate(`CD8/FoxP3` = CD8/FoxP3) %>%
  mutate(`CD8/CD163` = CD8/CD163) %>%
  mutate(`FoxP3/CD163` = FoxP3/CD163) %>%
  mutate(`FoxP3/Tumor` = FoxP3/Tumor) %>%
  mutate(`CD163/Tumor` = CD163/Tumor) #%>%
  melt(id.vars = c("Core","PatientID", "Group"), variable.name = "ctype")
  



type <- 'Tumor'  
dat <- cl_profile %>%
  dplyr::filter(ctype == type) %>%
  mutate(title = type)

stat.test <- dat %>% 
  group_by(ctype) %>%
  wilcox_test(value ~ Group) %>% 
  add_significance() %>%
  add_xy_position(x = "Group")

p <- ggplot(dat, aes(x = as.factor(Group), y = value), fill = NA) + 
  theme_bw() +
  facet_grid(.~title) +
  geom_boxplot() +
  geom_jitter() +
  theme(axis.text = element_text(size = 22),
        axis.title = element_text(size = 24),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size = 24),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #axis.text.x = element_text(angle = 90, vjust = 0.6, hjust=0.6),
        #axis.text.x = element_blank(),
        legend.position = 'na'
  ) +
  #ylab(expression('Density (No.' ~ mm ^ -2 ~')')) +
  ylab('') +
  #ylim(0, 450) + # FoxP3
  #ylim(0, 700) + # CD8
  #ylim(0, 2200) + # CD163
  ylim(0, 8400) + # Tumor
  stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01, size = 6, y.position = 8200) 
p
ggsave(file= paste("./Figures/Tumor_Compare_G1_G3" , '.png'), plot = print(p), width = 4, height = 5, dpi = 300)





# --------------------------- Shannon's Index -----------------------#
# REF: PMID 20101094

ShannonIndex_all <- data.frame(matrix(nrow = 0, ncol = 0))
for(Core in cl_profile$Core){
  
  # Survival group
  group <- cl_profile[cl_profile$Core == Core, 'Group']
  
  # Survival group
  
  
  
  PatientID <- cl_profile[cl_profile$Core == Core, 'PatientID']
  
  # coordinate data
  CellProp <- cur_df %>%
    dplyr::filter(Core == Core)
  
  # Compute Shannon Index
  coredata <- NSCLCdata[NSCLCdata$Core == Core,] %>%
    dplyr::filter(Phenotype != 'Other') %>%
    dplyr::filter(Phenotype != 'FoxP3CD8')
  
  cell_types <- unique(coredata$Phenotype)
  
  ShannonIndex <- ShannonE(cell_types, coredata)
  
  # 
  ShannonIndex_all <- rbind.data.frame(ShannonIndex_all, cbind.data.frame(ShannonIndex, Core, group, PatientID))
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
  wilcox_test(ShannonIndex ~ group) %>% 
  add_significance() %>%
  add_xy_position(x = "group")

ShannonIndex_all$title <- 'Shannon Entropy '
p <- ggplot(ShannonIndex_all, aes(x = as.factor(group), y = ShannonIndex), fill = NA) + 
  theme_bw() +
  facet_grid(.~title) +
  geom_boxplot() +
  geom_jitter() +
  theme(axis.text = element_text(size = 22),
        axis.title = element_text(size = 24),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size = 24),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #axis.text.x = element_text(angle = 90, vjust = 0.6, hjust=0.6),
        #axis.text.x = element_blank(),
        legend.position = 'na'
  ) +
  #ylab(expression('Density (No.' ~ mm ^ -2 ~')')) +
  ylab('') +
  ylim(0, 2) + 
  stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01, size = 6, y.position = 1.9) 
p
ggsave(file= paste("./Figures/Shannon_Index" , '.png'), plot = print(p), width = 4, height = 6, dpi = 300)



wilcox.test(ShannonIndex_all[ShannonIndex_all$group == '1', 'ShannonIndex'], ShannonIndex_all[ShannonIndex_all$group == '3', 'ShannonIndex'])



