#--------------- FUNCTION -------------#
# This script analyze the PD-1/PD-L1

library(spatstat); library(Rfast2); library(data.table); library(Rfast)
library(dplyr); library(tidyr); library(reshape2); library(tibble); library(tidyverse)
library(ggplot2); library(umap); library(ggvoronoi); library(Rtsne);
library(ggforce); library(UpSetR); library(ComplexUpset); library(rstatix); library(dendextend)
library(tripack); library(igraph); library(ComplexHeatmap); library(circlize)
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
# read patient data
patient_data <- read.csv('./NSCLC Clinical Data.csv')
# core demographic data
Core_demo <- read.csv('./Clinicopathology parameters.csv')
# core area
tissue_area <- read.csv('./TissueArea.csv') %>%
  mutate(total = (StromaArea + TumorArea)/1000000)

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



#------------------------------------------------------------------------------#
# CONSTRUCT UpSet Diagram
# Explanations: 
#             two rows: PD1 and PDL1
#             four columns: CD8, FoxP3, Tumor, CD163

celltypes <- c('PD1', 'PDL1')

NSCLCdata_UpSet <- NSCLCdata %>%
  dplyr::filter(Phenotype != 'FoxP3CD8') %>%
  dplyr::filter(Phenotype != 'Other')

NSCLCdata_UpSet$gid <- 1:nrow(NSCLCdata_UpSet)

NSCLCdata_UpSet_R <- data.frame(gid = NSCLCdata_UpSet$gid, Phenotype = NSCLCdata_UpSet$Phenotype, PDL1 = NA, PD1 = NA)

for(id in 1:nrow(NSCLCdata_UpSet)){
  
  gid <- NSCLCdata_UpSet[id, 'gid']
  
  Expr <- NSCLCdata_UpSet[NSCLCdata_UpSet$gid == gid, 'ExprPhenotype']
  
  if(Expr == 4){
    NSCLCdata_UpSet_R[NSCLCdata_UpSet_R$gid == gid, 'PDL1'] <- 1
    NSCLCdata_UpSet_R[NSCLCdata_UpSet_R$gid == gid, 'PD1'] <- 0
  } 
  if(Expr == 64){
    NSCLCdata_UpSet_R[NSCLCdata_UpSet_R$gid == gid, 'PDL1'] <- 0
    NSCLCdata_UpSet_R[NSCLCdata_UpSet_R$gid == gid, 'PD1'] <- 1
  }
  if(Expr == 68){
    NSCLCdata_UpSet_R[NSCLCdata_UpSet_R$gid == gid, 'PDL1'] <- 1
    NSCLCdata_UpSet_R[NSCLCdata_UpSet_R$gid == gid, 'PD1'] <- 1
  }
  if(Expr == 0){
    NSCLCdata_UpSet_R[NSCLCdata_UpSet_R$gid == gid, 'PDL1'] <- 0
    NSCLCdata_UpSet_R[NSCLCdata_UpSet_R$gid == gid, 'PD1'] <- 0
  }
}

#saveRDS(NSCLCdata_UpSet_R, 'NSCLCdata_UpSet_R.rds')

NSCLCdata_UpSet_R %>%
  dplyr::filter(PD1 == 1 & Phenotype == 'CD163') %>%
  nrow()

1861/488379

NSCLCdata_UpSet_R <- readRDS('NSCLCdata_UpSet_R.rds')
pdf(file = "./Figures/Upset_PD1PDL1.pdf", width = 8, height = 8)
ComplexUpset::upset(
  NSCLCdata_UpSet_R,
  celltypes,
  base_annotations=list(
    'Intersection size'=intersection_size(
      counts=TRUE,
      mapping=aes(fill=Phenotype)
    ) +
      scale_fill_manual(values = c('CD8' = '#a6cee3', 'FoxP3' = '#33a02c', 'CD163' = '#cab2d6', 
                                   'Tumor' = '#e31b1c'))
  ),
  themes = upset_default_themes(text = element_text(size = 30)),
  set_sizes=(
    upset_set_size()
    # you can also add annotations on top of bars:
    + expand_limits(y=1100)
    + theme(axis.text.x=element_text(angle=90))
  ),
  width_ratio=0.1
)
dev.off()



nrow(NSCLCdata_UpSet_R[NSCLCdata_UpSet_R$PDL1 == 1,])
nrow(NSCLCdata_UpSet_R[NSCLCdata_UpSet_R$PD1 == 1,])

#-------------------------------------------------------------#
# Before we divide PD1/PDL1 into low, mid, and high
# we first findout positive PD-1 PD-L1 cells and compare between survival groups

###########################################
# Expr: 0: DAPI
#       4: PDL1 (Opal 520)
#      64: PD1 (High - overexhausted - Opal 650)
#      68: PDL1 PD1

NSCLCdata_eval <- NSCLCdata %>%
  dplyr::filter(Phenotype != 'FoxP3CD8') %>%
  dplyr::filter(Phenotype != 'Other')


posPDL1 <- NSCLCdata_eval %>%
  dplyr::filter(ExprPhenotype == 64 | ExprPhenotype == 68) %>%
  merge(Core_demo, by = 'Core') %>%
  merge(patient_data_clus, by = 'PatientID') %>%
  dplyr::filter(Region == 'Central')


n_posPDL1 <- posPDL1 %>%
  group_by(Core, Phenotype) %>%
  tally() %>%
  merge(Core_demo, by = 'Core') %>%
  merge(patient_data_clus, by = 'PatientID') %>%
  dplyr::filter(Region == 'Central') %>%
  dplyr::filter(Phenotype == 'Tumor') %>%
  merge(tissue_area, by = 'Core') %>%
  mutate(density = n / total)
  



stat.test <- n_posPDL1 %>% 
  wilcox_test(density ~ group) %>% 
  add_significance() %>%
  add_xy_position(x = "group")

p <- ggplot(data = n_posPDL1, aes(x = as.factor(group), y = density)) +
  stat_boxplot(geom ='errorbar', width = 0.4) +
  geom_boxplot(width = 0.6, fill = "white") +
  stat_pvalue_manual(stat.test, size = 8) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        axis.text.y = element_text(angle = 90)) +
  ylab('') +
  #scale_y_continuous(breaks=c(0, 2000, 4000, 7000)) +
  ylim(0, max(n_posPDL1$density) + 5)
p
ggsave(p, file = './Figures/Tumor_PD1.jpeg', width = 3, height = 4, units = "in", dpi = 300)


#--------------- PD-L1 expression distributions (intensity histogram) -------------#


posPDL1 <- NSCLCdata_eval %>%
  dplyr::filter(ExprPhenotype == 4 | ExprPhenotype == 68) %>%
  merge(Core_demo, by = 'Core') %>%
  merge(patient_data_clus, by = 'PatientID') %>%
  dplyr::filter(Region == 'Central')



tertiles <- log(quantile((posPDL1$MeanMembrane520), probs = c(0.33, 0.67)))

p <- ggplot(posPDL1, aes(log(MeanMembrane520))) +
  theme_bw() +
  geom_histogram(binwidth = 0.2, color = 'black', fill = '#ee1a22', alpha = 0.3) +
  geom_density(aes(y = ..density.. * nrow(posPDL1) * 0.2), color = '#ee1a22', size = 1) +
  theme(axis.text.y = element_text(angle = 90),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16)) +
  xlab('log(intensity)') +
  ylab('Nucleus count') +
  geom_vline(xintercept = tertiles, color = 'black', linetype = 'dashed', size = 1) +
  scale_y_continuous(breaks=c(0, 6250, 12500, 18750, 25000))

p
ggsave(p, file = './Figures/PD-L1_distribution.jpeg', width = 6, height = 4, units = "in", dpi = 300)

#--------------- PD-1 expression distributions -------------#
posPD1 <- NSCLCdata_eval %>%
  dplyr::filter(ExprPhenotype == 64 | ExprPhenotype == 68) %>%
  merge(Core_demo, by = 'Core') %>%
  merge(patient_data_clus, by = 'PatientID') %>%
  dplyr::filter(Region == 'Central')

tertiles <- quantile(log(posPD1$MeanMembrane650), probs = c(0.33, 0.67))


p <- ggplot(posPD1, aes(log(MeanMembrane650))) +
  theme_bw() +
  geom_histogram(binwidth = 0.2, color = 'black', fill = '#4382c4', alpha = 0.3) +
  geom_density(aes(y = ..density.. * nrow(posPD1) * 0.2), color = '#4382c4', size = 1) +
  theme(axis.text.y = element_text(angle = 90),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16)) +
  xlab('log(intensity)') +
  ylab('Nucleus count') +
  geom_vline(xintercept = tertiles, color = 'black', linetype = 'dashed', size = 1)
p
ggsave(p, file = './Figures/PD-1_distribution.jpeg', width = 6, height = 4, units = "in", dpi = 300)


# assign expression level (binned - low, medium, and high)


#------------------- Comapre the expression intensities between low, mid, and high groups

tertiles <- quantile(log(posPD1$MeanMembrane650), probs = c(0.33, 0.67))
posPD1$PD1_Bin <- ifelse(posPD1$MeanMembrane650 > tertiles[2], 1, 0)
posPD1$PD1_Bin <- ifelse(posPD1$MeanMembrane650 < tertiles[1], -1, posPD1$PD1_Bin)




stat.test <- posPD1 %>% 
  wilcox_test(MeanMembrane650 ~ PD1_Bin) %>% 
  add_significance() %>%
  add_xy_position(x = "PD1_Bin")

stat.test[1, 'y.position'] <- log(stat.test[2, 'y.position']) 
stat.test[2, 'y.position'] <- log(stat.test[2, 'y.position']) + 0.3
stat.test[3, 'y.position'] <- log(stat.test[3, 'y.position']) + 0.6

p <- ggplot(data = posPD1, aes(x = as.factor(PD1_Bin), y = log(MeanMembrane650))) +
  stat_boxplot(geom ='errorbar', width = 0.4) +
  geom_boxplot(width = 0.6, fill = "white") +
  stat_pvalue_manual(stat.test, size = 8) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 22),
        axis.title = element_text(size = 24),
        axis.text.y = element_text(angle = 90)) +
  ylab('log(Mean PD-1 membrane intensity)') +
  scale_x_discrete(breaks=c(-1, 0, 1), labels = c('Low', 'Mid', 'High')) 
p

ggsave(p, file = './Figures/PD-1_Intensity_Binned_Compare.jpeg', width = 4.5, height = 6, units = "in", dpi = 300)





tertiles <- quantile(log(posPDL1$MeanMembrane520), probs = c(0.33, 0.67))
posPDL1$PDL1_Bin <- ifelse(posPDL1$MeanMembrane520 > tertiles[2], 1, 0)
posPDL1$PDL1_Bin <- ifelse(posPDL1$MeanMembrane520 < tertiles[1], -1, posPDL1$PDL1_Bin)

stat.test <- posPDL1 %>% 
  wilcox_test(MeanMembrane520 ~ PDL1_Bin) %>% 
  add_significance() %>%
  add_xy_position(x = "PDL1_Bin")

stat.test[1, 'y.position'] <- log(stat.test[2, 'y.position']) 
stat.test[2, 'y.position'] <- log(stat.test[2, 'y.position']) + 0.3
stat.test[3, 'y.position'] <- log(stat.test[3, 'y.position']) + 0.6

p <- ggplot(data = posPDL1, aes(x = as.factor(PDL1_Bin), y = log(MeanMembrane520))) +
  stat_boxplot(geom ='errorbar', width = 0.4) +
  geom_boxplot(width = 0.6, fill = "white") +
  stat_pvalue_manual(stat.test, size = 8) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 22),
        axis.title = element_text(size = 24),
        axis.text.y = element_text(angle = 90)) +
  ylab('log(Mean PD-L1 membrane intensity)') +
  scale_x_discrete(breaks=c(-1, 0, 1), labels = c('Low', 'Mid', 'High')) 
p

ggsave(p, file = './Figures/PD-L1_Intensity_Binned_Compare.jpeg', width = 4.5, height = 6, units = "in", dpi = 300)


#-------------------- The distribution of low, mid, high PD1/PDL1 cells in patients --------#
# Prevalence map
# Present as a heatmap
# Each row will represents a patient while each column a category

posPDL1_char <- posPDL1
posPDL1_char$Category <- paste0(posPDL1_char$PDL1_Bin, 'PDL1')

posPD1_char <- posPD1
posPD1_char$Category <- paste0(posPD1_char$PD1_Bin, 'PD1')



PD1PDL1_char <- rbind.data.frame(posPD1_char[, c('Core', 'Category')], posPDL1_char[, c('Core', 'Category')])

PD1PDL1_char$Category <- as.factor(PD1PDL1_char$Category)

Prevalence <- PD1PDL1_char %>%
  group_by(Core, Category, .drop = FALSE) %>%
  tally() %>%
  spread(key = Category, value = n) %>%
  data.frame() %>%
  select(-1) %>%
  as.matrix()



Prevalence <- ifelse(Prevalence>0, 1,Prevalence) %>%
  as.data.frame() %>%
  apply(2, as.character) %>%
  `colnames<-` (c('PD-1 low', 'PD-L1 low', 'PD-1 mid', 'PD-L1 mid', 'PD-1 high', 'PD-L1 high'))


# assign color
colors_pm = structure(c('white', '#4364a9'), names = c('0', '1'))
colors_sg = structure(c('#8baad8', '#e29d4a'), names = c('1', '3'))

group <- PD1PDL1_char %>%
  group_by(Core, Category, .drop = FALSE) %>%
  tally() %>%
  spread(key = Category, value = n) %>%
  merge(Core_demo, by = 'Core') %>%
  merge(patient_data_clus, by = 'PatientID') %>%
  dplyr::filter(Region == 'Central') %>%
  select(group) %>%
  as.matrix()
  
Prevalence <- Prevalence %>%
  as.data.frame() %>%
  mutate(group = group) %>%
  apply(2, as.character)

png("./Figures/Prevalence_map.png",width=3,height=6,units="in",res=300)
ha <- ComplexHeatmap::Heatmap(Prevalence[,1:6], cluster_rows = FALSE, cluster_columns = FALSE, col = colors_pm,
                              show_heatmap_legend = FALSE) +
  Heatmap(Prevalence[, 'group'], show_heatmap_legend = FALSE, col = colors_sg)
ha

dev.off()



lgd1 <- Legend(labels = c('0', '1'), 
               title = "Presence", 
               legend_gp = gpar(fill = c('white', '#4364a9')),
               ncol = 1)

lgd2 <- Legend(labels = c('1', '3'), 
               legend_gp = gpar(fill = c('#8baad8', '#e29d4a')),
               title = "Survival group",
               ncol = 1)

pd = packLegend(lgd1, lgd2, direction = "horizontal", 
                column_gap = unit(0.5, "cm"))
draw(pd)



#----------- Chi-square test----------------#

Prevalence_DF <- Prevalence %>%
  data.frame()

x <- factor(Prevalence[,4], levels = c('0', '1'))
y <- factor(Prevalence[,5], levels = c('0', '1'))
x
y


# PD-L1 high
id <- 6
Prevalence_DF[, id] <- factor(Prevalence_DF[,id], levels = c('0', '1'))

test <- Prevalence_DF[,  c(id,7)]
table(test)

stats::chisq.test(table(test))

#----------------------------------------------------#




BinCounts <- posPDL1 %>%
  group_by(Core, PDL1_Bin) %>%
  tally() %>%
  data.frame() %>%
  merge(Core_demo, by = 'Core') %>%
  merge(patient_data_clus, by = 'PatientID') %>%
  dplyr::filter(Region == 'Central') %>%
  merge(tissue_area, by = 'Core') %>%
  mutate(density = n / total) %>%
  mutate(PDL1_Bin = recode(PDL1_Bin, 
                '-1' = 'low',
                '0' = 'mid',
                '1' = 'high'))
  
BinCounts$PDL1_Bin = factor(BinCounts$PDL1_Bin, levels=c('low','mid','high'))


stat.test <- BinCounts %>% 
  group_by(PDL1_Bin) %>%
  wilcox_test(density ~ group) %>% 
  add_significance() %>%
  add_xy_position(x = "group")



p <- ggplot(data = BinCounts, aes(x = as.factor(group), y = density)) +
  stat_boxplot(geom ='errorbar', width = 0.4) +
  geom_boxplot(width = 0.6, fill = "white") +
  stat_pvalue_manual(stat.test, size = 6) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        axis.text.y = element_text(angle = 90),
        strip.text = element_text(size = 20),
        strip.background = element_rect(fill = "white")) +
  ylab('') +
  facet_wrap(.~PDL1_Bin, scales = 'free') 
p
ggsave(p, file = './Figures/PD-L1_category_compare.jpeg', width = 8, height = 5, units = "in", dpi = 300)







BinCounts <- posPD1 %>%
  group_by(Core, PD1_Bin) %>%
  tally() %>%
  data.frame() %>%
  merge(Core_demo, by = 'Core') %>%
  merge(patient_data_clus, by = 'PatientID') %>%
  dplyr::filter(Region == 'Central') %>%
  merge(tissue_area, by = 'Core') %>%
  mutate(density = n / total) %>%
  mutate(PD1_Bin = recode(PD1_Bin, 
                           '-1' = 'low',
                           '0' = 'mid',
                           '1' = 'high'))

BinCounts$PD1_Bin = factor(BinCounts$PD1_Bin, levels=c('low','mid','high'))


stat.test <- BinCounts %>% 
  group_by(PD1_Bin) %>%
  wilcox_test(density ~ group) %>% 
  add_significance() %>%
  add_xy_position(x = "group")



p <- ggplot(data = BinCounts, aes(x = as.factor(group), y = density)) +
  stat_boxplot(geom ='errorbar', width = 0.4) +
  geom_boxplot(width = 0.6, fill = "white") +
  stat_pvalue_manual(stat.test, size = 6) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        axis.text.y = element_text(angle = 90),
        strip.text = element_text(size = 20),
        strip.background = element_rect(fill = "white")) +
  ylab('') +
  facet_wrap(.~PD1_Bin, scales = 'free') 
p
ggsave(p, file = './Figures/PD-1_category_compare.jpeg', width = 8, height = 5, units = "in", dpi = 300)




#---------------Supplementary ------------------------#

# Figure 3F further catagorized by cell type

BinCounts <- posPD1 %>%
  group_by(Core, Phenotype, PD1_Bin) %>%
  tally() %>%
  data.frame() %>%
  merge(Core_demo, by = 'Core') %>%
  merge(patient_data_clus, by = 'PatientID') %>%
  dplyr::filter(Region == 'Central') %>%
  merge(tissue_area, by = 'Core') %>%
  mutate(density = n / total) %>%
  mutate(PD1_Bin = recode(PD1_Bin, 
                          '-1' = 'low',
                          '0' = 'mid',
                          '1' = 'high')) %>%
  filter(PD1_Bin == 'low')

#BinCounts$PD1_Bin = factor(BinCounts$PD1_Bin, levels=c('low','mid','high'))


stat.test <- BinCounts %>% 
  group_by(Phenotype) %>%
  wilcox_test(density ~ group) %>% 
  add_significance() %>%
  add_xy_position(x = "group")



p <- ggplot(data = BinCounts, aes(x = as.factor(group), y = density)) +
  stat_boxplot(geom ='errorbar', width = 0.4) +
  geom_boxplot(width = 0.6, fill = "white") +
  stat_pvalue_manual(stat.test, size = 6) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        axis.text.y = element_text(angle = 90),
        strip.text = element_text(size = 20),
        strip.background = element_rect(fill = "white")) +
  ylab('') +
  facet_wrap(.~Phenotype, scales = 'free', ncol = 4) 
p
ggsave(p, file = './Figures/PD-1__category_compare.jpeg', width = 12, height = 5, units = "in", dpi = 300)






BinCounts <- posPDL1 %>%
  group_by(Core, Phenotype, PDL1_Bin) %>%
  tally() %>%
  data.frame() %>%
  merge(Core_demo, by = 'Core') %>%
  merge(patient_data_clus, by = 'PatientID') %>%
  dplyr::filter(Region == 'Central') %>%
  merge(tissue_area, by = 'Core') %>%
  mutate(density = n / total) %>%
  mutate(PDL1_Bin = recode(PDL1_Bin, 
                           '-1' = 'low',
                           '0' = 'mid',
                           '1' = 'high')) %>%
  filter(PDL1_Bin == 'high')




stat.test <- BinCounts %>% 
  group_by(PDL1_Bin) %>%
  wilcox_test(density ~ group) %>% 
  add_significance() %>%
  add_xy_position(x = "group")



p <- ggplot(data = BinCounts, aes(x = as.factor(group), y = density)) +
  stat_boxplot(geom ='errorbar', width = 0.4) +
  geom_boxplot(width = 0.6, fill = "white") +
  stat_pvalue_manual(stat.test, size = 6) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        axis.text.y = element_text(angle = 90),
        strip.text = element_text(size = 20),
        strip.background = element_rect(fill = "white")) +
  ylab('') +
  facet_wrap(.~Phenotype, scales = 'free', ncol = 4) 
p
ggsave(p, file = './Figures/PD-L1_category_compare_ct.jpeg', width = 12, height = 5, units = "in", dpi = 300)



