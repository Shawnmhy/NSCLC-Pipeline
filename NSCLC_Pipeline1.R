library(spatstat); library(Rfast2)
library(dplyr); library(tidyr); library(reshape2); library(tibble); library(tidyverse)
library(ggplot2); library(umap); library(ggvoronoi); library(Rtsne); library(ComplexHeatmap)
library(ggforce); library(UpSetR)
library(dendextend)
library(readxl); library(stringr); 
library(survival); library(survminer)
# working directory
current_path <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))

setwd('..')


#---------------------- DATA INPUT ----------------#
# read single-cell data
NSCLCdata <- readRDS('./NSCLCdataset.rds')
# read patient data
patient_data <- read.csv('./NSCLC Clinical Data.csv')
# core demographic data
Core_demo <- read.csv('./Clinicopathology parameters.csv')
# tissue demographic data
Tissue_demo <- read.csv('./Tissue_demographics.csv')

#--------------------------------------------------#



#-------------- Barplot to show overall abundance ------------#

nucleusQuant <- table(NSCLCdata$Phenotype) %>%
  data.frame()

p<-ggplot(data=nucleusQuant, aes(x=Var1, y=log(Freq), fill = Var1, width = .9)) +
  theme_classic() +
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size = 32),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.6, hjust=0.6),
        legend.position = 'na') +
  ylab('log(Total nuclei count)') +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c('CD8' = '#a6cee3', 'FoxP3' = '#33a02c', 'CD163' = '#cab2d6', 'FoxP3CD8' = '#fdbf6f',
                               'Tumor' = '#e31b1c', 'Other' = '#1f78b4')) 
p

ggsave(file= paste("./Figures/Total_Nuclei_Counts" , '.png'), plot = print(p), width = 6, height = 8, dpi = 300)


#-------------- Dimensionally reduction ------------#

# data for umap, 10% of all by sampling
NSCLCdata_UMAP <- NSCLCdata %>%
  select(-contains(c('650', '520', '480', '780')))

unique(NSCLCdata$core)

tumor_umap <- NSCLCdata_UMAP %>%
  dplyr::filter(Phenotype == 'Tumor') %>%
  sample_n(32673)

other_umap <- NSCLCdata_UMAP %>%
  dplyr::filter(Phenotype == 'Other') %>%
  sample_n(63127) 

CD163_umap <- NSCLCdata_UMAP %>%
  dplyr::filter(Phenotype == 'CD163') %>%
  sample_n(15739)

umapDat <- rbind.data.frame(tumor_umap, other_umap, CD163_umap,
                            NSCLCdata_UMAP[NSCLCdata_UMAP$Phenotype == 'CD8',],
                            NSCLCdata_UMAP[NSCLCdata_UMAP$Phenotype == 'FoxP3',],
                            NSCLCdata_UMAP[NSCLCdata_UMAP$Phenotype == 'FoxP3CD8',]) %>%
  select(c(7, 10:51)) %>%
  na.omit()


#umapRes <- umap(umapDat[, 2:43])$layout %>%
#  data.frame()

#umapRes <- data.frame(umapRes)


tsneRes <- Rtsne(umapDat[, 2:43])$Y %>%
  data.frame()

p <- ggplot() +
  theme_classic() +
  geom_point(data = tsneRes, aes(x = X1, y = X2, color = as.factor(umapDat$Phenotype))) +
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size = 32),
        axis.text.x = element_text(angle = 90, vjust = 0.6, hjust=0.6),
        legend.position = 'na') +
  xlab('tsne 1') +
  ylab('tsne 2') +
  scale_color_manual(values = c('CD8' = '#a6cee3', 'FoxP3' = '#33a02c', 'CD163' = '#cab2d6', 'FoxP3CD8' = '#fdbf6f',
                                'Tumor' = '#e31b1c', 'Other' = '#1f78b4')) 
p
ggsave(file= paste("./Figures/tSNE_subset" , '.png'), plot = print(p), width = 8, height = 8, dpi = 300)



#-------------- Comparison of cell fractions between edge and central ------------#
# 

quant_profiles_all <- readRDS('./Nucleus_Density_Core.rds')

quant_profiles_all_region <- merge(quant_profiles_all, patient_meta, by = 'Core')
# Boxplot
p <- ggplot(quant_profiles_all_region, aes(x=ctype, y = density, fill=Region)) + 
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.title.x = element_blank(),
        #axis.text.x = element_text(angle = 90, vjust = 0.6, hjust=0.6),
        #axis.text.x = element_blank(),
        #legend.position = 'na'
        ) +
  scale_fill_manual(values = c('Central' = 'black', 'Edge' = 'grey')) +
  ylab(expression('Density (No.' ~ mm ^ -2 ~')')) +
  geom_boxplot(alpha = 0.6) 

p
ggsave(file= paste("./Figures/Boxplot_Region_Compare" , '.png'), plot = print(p), width = 5, height = 4, dpi = 300)



edge_dat <- quant_profiles_all_region %>%
  dplyr::filter(Region == 'Edge') %>%
  dplyr::filter(ctype == 'Tumor') %>%
  select(density) %>%
  as.matrix()

central_dat <- quant_profiles_all_region %>%
  dplyr::filter(Region == 'Central') %>%
  dplyr::filter(ctype == 'Tumor') %>%
  select(density) %>%
  as.matrix()

wilcox.test(edge_dat, central_dat)







#-------------- Cell type fraction per image with hierarchical clustering ------------#
# Create table with celltype fractions

cur_df <- NSCLCdata %>%
  # remove control cores
  dplyr::filter(core %in% Core_demo$Core) %>%
  group_by(core, Phenotype) %>%
  summarise(n=n()) %>% 
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


row_sorted <- hc$labels

# now we generate the clinical metadata and order it in the same way as the celltype data
patient_meta <- Core_demo[, c('Core', 'Region')]

mrownames <- patient_meta$Core
patient_meta <- as.matrix(patient_meta)
rownames(patient_meta) <- mrownames
patient_meta <- data.frame(patient_meta[row_sorted,])

# generate the barplot. this is generated as the annotation for the heatmap of the Patient_ID that is generated below.
hm_dat <- hm_dat[row_sorted,]

# bring cell types in order (column order)
col_order <- c("Tumor", "CD163", "CD8", "FoxP3", "FoxP3CD8", "Other")
hm_dat <- hm_dat[, col_order]

col_vector <- structure(c("#e31b1c", "#cab2d6", '#a6cee3', '#33a02c', '#fdbf6f', '#1f78b4'), 
                        names = c("Tumor", "CD163", "CD8", "FoxP3", "FoxP3CD8", "Other"))

# rename punch locations
col_vector_region <- structure(c("black", "grey"), 
                               names = c("Central", "Edge"))

# create annotation with cell type proportions
ha <- rowAnnotation(`Cell Type Proportion` =
                      anno_barplot(hm_dat, 
                                   gp=gpar(fill=col_vector),
                                   bar_width = 1,
                                   height = unit(25,"cm"),
                                   width = unit(11,"cm"),
                                   show_row_names = FALSE),
                    col = list(`Region Location` = col_vector_region),
                    show_legend = FALSE)


dend <- as.dendrogram(hc)
dend <- color_branches(dend, k = 3, col = c("gray50", "blue", "red"), groupLabels = TRUE) # `color_branches()` returns a dendrogram object

# heatmap consisting of the patient_IDs. one color per patient

pdf(file = "./Figures/StackBarplotClustering.pdf",   # The directory you want to save the file in
    width = 15, # The width of the plot in inches
    height = 6) # The height of the plot in inches

h1 = Heatmap(patient_meta[,c("Region")], 
             col = col_vector_region, 
             width = unit(0.5, "cm"), 
             cluster_rows = dend,
             row_dend_width = unit(3, "cm"),
             height = unit(25, "cm"),
             show_heatmap_legend = FALSE, 
             heatmap_legend_param = list(title = "Core Location"),
             row_names_gp = gpar(cex=20),
             show_row_names = TRUE,
             right_annotation =  ha,
             column_labels = "Core Location")

# plot the data
ht = grid.grabExpr(draw(h1))
grid.newpage()
pushViewport(viewport(angle = 270))
grid.draw(ht)

dev.off()


lgd1 <- Legend(labels = names(col_vector), 
               title = "Nucleus Type", 
               legend_gp = gpar(fill = col_vector),
               ncol = 1)

lgd2 <- Legend(labels = names(col_vector_region), 
               legend_gp = gpar(fill = unname(col_vector_region)),
               title = "Core Location",
               ncol = 1)

pd = packLegend(lgd1, lgd2, direction = "vertical", 
                column_gap = unit(0.5, "cm"))
draw(pd)


#-------------- Plot example phenograph for each cluster ------------#

# Example: 
#         Cluster 1: C15
#         Cluster 2: A7
#         Cluster 3: B13

cellPos <- NSCLCdata %>%
  filter(Core == 'A1') %>%
  #filter(Phenotype != 'Other') %>%
  select('CellXPos', 'CellYPos', 'Phenotype', 'ExprPhenotype', 'CellID') #%>%


cellBounds <- read.csv(paste0('./polygonList/', 'TMA_1314_Core[1,1,1]_[66431,20799]_polygons.csv')) %>%
  merge(cellPos, by = 'CellID') 



p <- ggplot(cellPos) + 
  theme_void() +
  theme(legend.position = 'none',
        #panel.background = element_rect(fill = 'black')
        ) +
  #geom_polygon(aes(PosX, PosY, group = CellID, fill = as.factor(Phenotype)), color = 'black') +
  geom_point(aes(CellXPos, CellYPos, color = Phenotype)) +
  #geom_path(data=bdry_tr, aes(x=x, y=y, group=draw), color = 'white', size = 1) +  #draw is the original ring
  # group label
  #scale_fill_manual(values = c('0' = 'gray', '1' = '#d67121', '2' = '#cb181b', '3' = '#70361c',
  #                             '4' = '#2700f7', '5' = '#10006a', '6' = '#66aefa')) +
  # Cell Type
  scale_color_manual(values = c('CD8' = '#a6cee3', 'FoxP3' = '#33a02c', 'CD163' = '#cab2d6', 'FoxP3CD8' = '#fdbf6f',
                                                            'Tumor' = '#e31b1c', 'Other' = '#1f78b4')) +
  #scale_fill_manual(values = c('CD8' = '#f8f8f8', 'FoxP3' = '#f8f8f8', 'CD163' = '#cab2d6', 'Other' = '#f8f8f8',
  #                             'Tumor' = '#e31b1c', 'FoxP3CD8' = t'#f8f8f8')) + #f8f8f8  
  ylim(2008, 0)
p

ggsave(file= paste("./Figures/A1" , '_pointcloud.png'), plot = print(p), width = 8, height = 8, dpi = 300)











#------------- Survival Analysis on Clustered Tissue Architectures ------------#

# Note: 
# High tumor: Cluster 
# Low tumor: 4 + 1
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

#median(clus_dist$group)

patient_data_clus <- merge(clus_dist, patient_data, by = 'PatientID')


#patient_data_clus$label <- 1* ifelse(patient_data_clus$group >= 10, 1, 2)
#patient_data_clus$group2 <- patient_data_clus$group
#patient_data_clus$group2 <- ifelse(patient_data_clus$group == 1 | patient_data_clus$group == 4, 1, 2)

patient_data_clus <- patient_data_clus %>%
  dplyr::filter(group != 2)

fit <- survfit(Surv(Time, Event) ~ group, data = patient_data_clus)
fit


p <- ggsurvplot(fit, data = patient_data_clus, 
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
                pval = TRUE,
                break.time.by = 2000) +
  xlab('Time, (days)')
p
ggsave(file= paste("./Figures/Survival_Central_G1_3" , '.pdf'), plot = p$plot, width = 9, height = 6, dpi = 300)

res.cox <- coxph(Surv(Time, Event) ~ group, data = patient_data_clus)
ggforest(res.cox, data = patient_data_clus)

surv_diff <- survdiff(Surv(Time, Event) ~ group, data = patient_data_clus)
surv_diff





# Create dataset
data <- data.frame(
  individual=paste( "Mister ", seq(1,60), sep=""),
  group=c( rep('A', 10), rep('B', 30), rep('C', 14), rep('D', 6)) ,
  value1=sample( seq(10,100), 60, replace=T),
  value2=sample( seq(10,100), 60, replace=T),
  value3=sample( seq(10,100), 60, replace=T)
)

data <- data %>% gather(key = "Category", value="value", -c(1,2)) 

# example of list input (list of named vectors)
cl_tree <- cutree(hc, k = 3) %>%
  data.frame() %>%
  rownames_to_column() %>%
  `colnames<-` (c('Core', 'Cluster')) %>%
  merge(Core_demo, by = 'Core') %>%
  select(Cluster, PatientID, Region) %>%
  table() %>%
  data.frame() %>%
  reshape(idvar = c('Cluster', 'PatientID'), timevar = "Region", direction = "wide") %>%
  `colnames<-` (c('Cluster', 'PatientID', 'Central', 'Edge')) %>%
  mutate(Total = Central + Edge) %>%
  melt(id.vars = c("Cluster","PatientID"), variable.name = "Category")
#dplyr::filter(Region == 'Central')



# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 2
nObsType <- nlevels(as.factor(cl_tree$Cluster))
to_add <- data.frame( matrix(NA, empty_bar*nlevels(cl_tree$PatientID)*nObsType, ncol(cl_tree)) )
colnames(to_add) <- colnames(cl_tree)
to_add$PatientID <- rep(levels(cl_tree$PatientID), each=empty_bar*nObsType )
cl_tree <- rbind(cl_tree, to_add)
cl_tree <- cl_tree %>% arrange(PatientID, Category)
cl_tree$id <- rep( seq(1, nrow(cl_tree)/nObsType) , each=nObsType)

# Get the name and the y position of each label
label_data <- cl_tree %>% group_by(id, Cluster) %>% summarize(tot=sum(value))
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- cl_tree %>% 
  group_by(PatientID) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

# Make the plot


# Make the plot





#cl_tree$Cluster <- as.character(cl_tree$Cluster)
cl_tree$Category <- as.character(cl_tree$Category)

p <- ggplot(cl_tree) +      
  
  # Add the stacked bar
  geom_bar(aes(x=as.factor(id), y=value, fill=Cluster), stat="identity", alpha=0.5) +
  scale_fill_manual(values = c('1' = 'grey50', '2' = 'red', '3' = '#0f0fff')) +
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 2, xend = start, yend = 2), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 4, xend = start, yend = 4), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 6, xend = start, yend = 6), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 8, xend = start, yend = 8), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 10, xend = start, yend = 10), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  ggplot2::annotate("text", x = rep(max(cl_tree$id),6), y = c(0, 2, 4, 6, 8, 10), label = c("0", "2", "4", "6", "8", "10") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
  
  ylim(-6,max(10, na.rm=T)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,30), "cm") 
  ) +
  coord_polar() +
  
  # Add labels on top of each bar
  #geom_text(data=label_data, aes(x=id, y=tot+10, label=Cluster, hjust = hjust), color="black", fontface="bold",alpha=0.6, size=5, angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -0, xend = end, yend = 0), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title - 1, y = -1, label=PatientID), hjust=c(0,0,0,0,0,0,0,0,0,0, 0,0,0,0, 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1), colour = "black", alpha=0.8, size = 3, fontface="bold", inherit.aes = FALSE)

p

ggsave(p, file="./Figures/Circular_Stacked_Plot.png", width=6, height=6)

col_vector <- structure(c("#c9c9c9", "#ff999d", '#a1a3ff'), 
                        names = c("1", "2", "3"))

lgd1 <- Legend(labels = names(col_vector), 
               title = "Cluster type", 
               legend_gp = gpar(fill = col_vector),
               ncol = 1)


pd = packLegend(lgd1, direction = "vertical", 
                column_gap = unit(0.5, "cm"))
draw(pd)

