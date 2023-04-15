library(ggplot2)
library(RANN)
library(ggpubr)
library(rstatix)
library(dplyr); library(stringr); library(data.table)
library(readxl)

###########################################
# Expr: 0: DAPI
#       4: PDL1
#      64: PD1 (High - overexhausted - Opal 650)
#      68: PDL1 PD1
# working directory
current_path <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))

setwd('..')


allTables <- list.files('./Tables')


allTables_tr <- str_replace(allTables, '\\[', ',')
allTables_tr <- str_replace(allTables_tr, ']',',')


coreData_all <- data.frame(matrix(nrow = 0, ncol = 0))

for(id1 in seq_len(13)){
  for(id2 in seq_len(16)){
    
    
    # file id
    pattern <- paste0(',1,', id1, ',', id2, ',')
    
    core <- paste0(LETTERS[id1], id2)
    # read file data
    coreData <- read.csv(paste0('./Tables/', allTables[which(grepl(pattern,allTables_tr) == 'TRUE')]))
    
    
    # combine
    coreData_all <- rbind(coreData_all, cbind(core, coreData))
    
    
  }
}

saveRDS(coreData_all, 'NSCLCdataset.rds')


#---------------- Nucleus density per type per core

Region.area <- data.frame(read.csv('./Files_Archive/tissue_area_edge_versus_core.csv')) %>%
  select(c(Core, area_sum))

# compute density
quant_profiles_all <- data.frame(matrix(nrow = 0, ncol = 0))
for(Core in Region.area$Core){
  
  # core data
  core_data <- coreData_all %>%
    dplyr::filter(core == Core)
  
  # tissue area
  tissue_area <- Region.area[Region.area$Core == Core, 'area_sum']
  
  # get the number/density profile for the core
  quant_profiles <- table(core_data$Phenotype) %>%
    data.frame() %>%
    mutate(density = Freq / tissue_area) %>%
    mutate(Core = Core) %>%
    `colnames<-` (c('ctype', 'number', 'density', 'Core'))
  
  quant_profiles_all <- rbind.data.frame(quant_profiles_all, quant_profiles)
  
}

saveRDS(quant_profiles_all, 'Nucleus_Density_Core.rds')

quant_profiles_all <- readRDS('Nucleus_Density_Core.rds')

CD163 <- quant_profiles_all %>%
  filter(ctype == 'Other')

summary(CD163)

sd(CD163$density)



#---------------- Nucleus Polygon save all to one


for(file in list.files('./polygonList')) {
  
  #file = list.files('./boundary')[2]
  # name to Letter - Number combination
  spName <- strsplit(file, '\\]|,') # regex to extract cores
  
  core <- paste0(LETTERS[as.numeric(spName[[1]][2])], spName[[1]][3]) # pieces to cores
  
  file = fread(paste0('./polygonList/', file))
  
  # prefix
  #prefix <- strsplit(file, '_')
  #prefix <- paste0(prefix[[1]][1], '_', prefix[[1]][2], '_', prefix[[1]][3], '_', prefix[[1]][4], '_morphometrics.csv')
  
  
  # index NSCLC database to retrieve the core data
  cellPos <- NSCLCdata %>%
    filter(Core == Core) %>%
    #filter(Phenotype != 'Other') %>%
    select('CellXPos', 'CellYPos', 'Phenotype', 'ExprPhenotype', 'CellID') #%>%
  
  
  cellBounds <- read.csv(paste0('./polygonList/', 'TMA_1314_Core[1,3,1]_[66305,17745]_polygons.csv')) %>%
    merge(cellPos, by = 'CellID') 
  
  
  
}


