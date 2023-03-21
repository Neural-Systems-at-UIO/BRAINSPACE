##Script Outcomes: read in all nonlinear and rigid pilot load data, create bar/dot plots of load differences, plot changes in load between age/genotype groups
##load in packages
library(dplyr)
library(reshape2)
library(ggplot2)
library(circlize)
library(tidyr)

###function to read in load data-------------------------------------------------------------------------
getStainDat <- function(stain_name, path, folder_id, folder_name, data_id) {
  
  stain <- stain_name
  
  # remove(dataset)
  
  path <- paste0(path, stain)
  stain.files <- list.files(path, all.files = T,recursive = T, pattern = data_id)
  stain.files <- stain.files[grep(folder_id, stain.files, value = F)]
  stain.files <- stain.files[grep(paste0("/CKT-1-",stain, paste0(".*", data_id)), stain.files, value = F, perl = T)]
  for (file in stain.files) {
    
    fname <- file.path(path, file)
    brain.name <- gsub(paste0("(.*)/", folder_name, ".*"), "\\1", file)
    print(brain.name)
    
    
    if(!exists("dataset")){
      dataset = read.csv(fname, sep = ",")
      dataset$Brain.ID <- brain.name
    }
    
    if (exists("dataset")){
      temp_dataset = read.csv(fname, sep = ",")
      temp_dataset$Brain.ID <- brain.name
      dataset = rbind(dataset, temp_dataset)
      rm(temp_dataset)
    }
  }
  
  dataset <- dataset[!duplicated(dataset),]
  dataset <- dataset[colSums(is.na(dataset)) < nrow(dataset)]
  
  return(dataset)
  
}

##read in data output from QuickNII only registration (Overlay_QuickNII_Pilot) and QuickNII and VisuAlign registration (Overlay_Pilot)-------------------------------------
##load in stain data: two options since naming convention varies between stains:
##AB1-42 and NeuN have a underscore before the last ID numbers and Iba1, GFAP, and Thionine have dashes
##reading in .csv that was recreated with damaged slices removed if reading in raw data replace _PilotQNIISliceRM.csv with _Pilot_CustomRegions_All.csv

##Run for all stains or for subset or individual stains based on output needed:
stain.name <- c('NeuN', 'GFAP', 'Iba1', 'Thionine', 'AB1-42')

##adjust path = "" to the destination of parent stain folder obtained from EBRAINS portal
for (stain in stain.name) {
  if (stain != 'Iba1'){
    stain.rg <- getStainDat(stain_name = stain, 
                            path = "/CKT-1-",
                            folder_id = "/Overlay_QuickNII_Pilot/",folder_name = "Overlay_QuickNII_Pilot", data_id = "_PilotQNIISliceRM.csv")
    stain.rg$Brain.ID <- gsub('[a-zA-Z]+(\\d+\\-\\d+)?-', '', stain.rg$Brain.ID) %>% gsub('-(\\d+)$', '_\\1', .)
    stain.rg <- stain.rg[,c(1, 4:5)]
    colnames(stain.rg)[2] <- paste0('Load.',stain, '.rg')
    
    stain.nl <- getStainDat(stain_name = stain, 
                            path = "/CKT-1-",
                            folder_id = "/Overlay_Pilot/",folder_name = "Overlay_Pilot",data_id = "_PilotSliceRM.csv")
    stain.nl$Brain.ID <- gsub('[a-zA-Z]+(\\d+\\-\\d+)?-', '', stain.nl$Brain.ID) %>% gsub('-(\\d+)$', '_\\1', .)
    stain.nl <-stain.nl[,c(1, 4:5)]
    colnames(stain.nl)[2] <- paste0('Load.',stain, '.nl')
  } else {
    
    stain.rg <- getStainDat(stain_name = stain, 
                            path = "/CKT-1-",
                            folder_id = "/Overlay_QuickNII_Pilot/",folder_name = "Overlay_QuickNII_Pilot", data_id = "_PilotQNIISliceRM.csv")
    stain.rg$Brain.ID <- gsub('[a-zA-Z]+(\\d+)?-', '', stain.rg$Brain.ID) %>% gsub('-(\\d+)$', '_\\1', .)
    stain.rg <- stain.rg[,c(1, 4:5)]
    colnames(stain.rg)[2] <- paste0('Load.',stain, '.rg')
    
    stain.nl <- getStainDat(stain_name = stain, 
                            path = "/CKT-1-",
                            folder_id = "/Overlay_Pilot/",folder_name = "Overlay_Pilot",data_id = "_PilotSliceRM.csv")
    stain.nl$Brain.ID <- gsub('[a-zA-Z]+(\\d+)?-', '', stain.nl$Brain.ID) %>% gsub('-(\\d+)$', '_\\1', .)
    stain.nl <-stain.nl[,c(1, 4:5)]
    colnames(stain.nl)[2] <- paste0('Load.',stain, '.nl')
    
  }
  
  
  x <- Reduce(function(x,y) merge(x,y, by = c('Region.name', 'Brain.ID')), list(stain.rg, stain.nl))
  
  x <- x[!duplicated(x),]
  x[is.na(x)] = 0
  
  
  x3 <- dcast(x, Region.name ~ Brain.ID, value.var = c(paste0('Load.', stain, '.rg')), fill = 0)
  colnames(x3) <- c('RegionName', paste0(colnames(x3)[-1], "_RG"))
  x4 <- dcast(x, Region.name ~ Brain.ID, value.var = c(paste0('Load.', stain, '.nl')), fill = 0)
  colnames(x4) <- c('RegionName', paste0(colnames(x4)[-1], "_NL"))
  
  
  x.merged <- Reduce(function(x,y) merge(x,y, by = c('RegionName')), list(x3,x4))
  
  ##Output data in multiple formats:
  
  write.csv(x.merged, file = paste0('CustomRegistration_', stain, '_Wide.csv'), row.names = F)
  
  ########
  
  x3 <- dcast(x, Brain.ID ~ Region.name, value.var = c(paste0('Load.', stain, '.rg')), fill = 0)
  colnames(x3) <- c('Brain.ID', paste0(colnames(x3)[-1], "_RG_",stain))
  
  x4 <- dcast(x, Brain.ID ~ Region.name, value.var = c(paste0('Load.', stain, '.nl')), fill = 0)
  colnames(x4) <- c('Brain.ID', paste0(colnames(x4)[-1], "_NL_",stain))
  
  c.bind <- Reduce(function(x,y) merge(x,y, by = c('Brain.ID')), list(x3,x4))
  
  write.csv(c.bind, file = paste0('CustomRegistration_', stain, '_RM', '.csv'), row.names = F)
  
  
  #####
  
  x3 <- dcast(x, Brain.ID ~ Region.name, value.var = c(paste0('Load.', stain, '.rg')), fill = 0)
  rownames(x3) <- paste0(x3$Brain.ID, "_RG")
  
  x4 <- dcast(x, Brain.ID ~ Region.name, value.var = c(paste0('Load.', stain, '.nl')), fill = 0)
  rownames(x4) <- paste0(x4$Brain.ID, "_NL")
  
  x.bind <- rbind(x3, x4)
  
  write.csv(x.bind, file = paste0('CustomRegistration_', stain, '_Long.csv'), row.names = F)
}


##Combine load data with meta data:---------------------------------------------------------------------
###read in all individual stain files and clean up
##input path location of files created in lines 1-132:

AB <- read.csv("/CustomRegistration_AB1-42_RM.csv")
NeuN <- read.csv("/CustomRegistration_NeuN_RM.csv")
GFAP <- read.csv("/CustomRegistration_GFAP_RM.csv")
Iba1 <- read.csv("/CustomRegistration_Iba1_RM.csv")
Thionine <- read.csv("/CustomRegistration_Thionine_RM.csv")

##merge stain files to have one list with all stain/region
All_stains <- Reduce(function(x,y) merge(x,y, by = c('Brain.ID')), list(AB, NeuN, GFAP, Iba1, Thionine))

##clean up and add key to load data:
##input path location of meta data obtained from EBRAINS:
key <- read.csv("/keyxHarvestAge.csv")
names(key)[names(key) == "Harvest.Age..m."] <- "Age.Harvested"
key$Strain <- paste0("BXD", key$Strain)
key$Strain[key$Strain == "BXD2"] <- "DBA.2J"
key$Strain[key$Strain == "BXD6"] <- "C57BL.6J"

key$Genotype[key$Genotype == "5XFAD"] <- 5
key$Genotype[key$Genotype == "WT"] <- 0

names(key)[names(key) == "Gender"] <- "Sex"

key$Sex[key$Sex == "Male"] <- 1
key$Sex[key$Sex == "Female"] <- 0
key$Sex[key$Sex == "female"] <- 0
key$Sex[key$Sex == "11"] <- 1

##merge stain data and key:
df.key <- merge(key, All_stains, by = 'Brain.ID')
df.key <- df.key[,-c(2:11,14:15, 18:19)] 
df.key$Genotype <- as.numeric(df.key$Genotype)
df.key$Sex <- as.numeric(df.key$Sex)

NL_reg <- c('_NL')
NL_only <- df.key[,c(1:5,grep(paste(NL_reg, collapse = "|"),colnames(df.key)))]
colnames(NL_only)<-gsub("_NL","",colnames(NL_only))
# write.csv(NL_only, "NL_IntermediateCustomRegions_Load_Slice_RM.csv")

##reshape data frame to have long format:
long <- reshape2::melt(data = df.key, id.vars = c(1:5), measure.vars = c(6:815)) 

##new column for stain
long$Stain <- NA
long$Stain[grep('NeuN', long$variable)] = 'NeuN'
long$Stain[grep('Thionine', long$variable)] = 'Thionine'
long$Stain[grep('Iba1', long$variable)] = 'Iba1'
long$Stain[grep('GFAP', long$variable)] = 'GFAP'
long$Stain[grep('AB1.42', long$variable)] = 'AB1.42'

##new column for registration type: RG= rigid QuickNII registration, NL= VisuAlign registration
long$Registration <- NA
long$Registration[grep('RG', long$variable)] = 'RG'
long$Registration[grep('NL', long$variable)] = 'NL'

##New column for region:
long$Region<- NA
long$Region <- gsub("_RG_GFAP|_RG_NeuN|_RG_Thionine|_RG_Iba1|_RG_AB1.42|_NL_GFAP|_NL_NeuN|_NL_Thionine|_NL_Iba1|_NL_AB1.42","", long$variable)

##Rename value and convert to percentage:
names(long)[names(long) == "value"] <- "Load"
long$Load.Percentage <- (long$Load)*100

write.csv(long, 'NLvsRG_IntermediateCustomRegions_Load.csv')
##separate registration types:
rigid <- subset(long, Registration == 'RG')
names(rigid)[names(rigid) == "Load.Percentage"] <- "Load.Percentage_RG"
rigid <- rigid[,-c(6,7,9)] 

nonlinear <- subset(long, Registration == 'NL')
names(nonlinear)[names(nonlinear) == "Load.Percentage"] <- "Load.Percentage_NL"
nonlinear <- nonlinear[,-c(6,7,9)] 

##write out data:
write.csv(nonlinear, "NL_IntermediateCustomRegions_long_sliceRM_all.csv")

##combine to create a dataframe with all brains/regions and the rigid and nonlinear Load separated
registration_type <- cbind(nonlinear[,1:8],rigid[,8])
names(registration_type)[names(registration_type) == "rigid[, 8]"] <- "Load.Percentage_RG"

