library(tidyverse)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(ggfortify)

##Objective: Read in and calculate accuracy, inaccuracy, and uncertainty scores from output from QCAlign inter-rater reliability test:

##set wroking directory:
setwd("/QCAlign_InterUser_Results")

##adjust the path according to the location of files obtained from EBRAINS:
path <- ''

##List file names:
fname <- list.files(path, recursive = TRUE, pattern = "*_QControl_Thionine_stats(_SY|_BG|_AO|_HK|_IB|_MT|_SS|_TM|_TO|_US)")

##locate and read in output stats.txt files
stats_all <- lapply(fname, function(f){
  stats <- read.csv(file.path(path, f), sep = "\t")
  
  stain.id <- gsub('.*(Thionine|AB1-42|Iba1|GFAP|NeuN).*', '\\1', f)
  user.id <- gsub('.*(SY|BG|AO|HK|IB|MT|SS|TM|TO|US).*', '\\1', f)
  brain.id <- gsub('(*)\\_.*', '\\1', f)
  
  
  x <- stats %>% gather(Slice, Value, -(1:4)) %>% 
    filter(Status != 'Damaged') %>% 
    group_by(Structure, Status) %>% 
    summarise(region_total = sum(Value, na.rm = T)) %>% 
    group_by(Structure) %>% 
    
##Calculate accuracy, inaccuracy, and uncertainty scores:
    mutate(Accuracy = region_total[Status == 'Accurate']/sum(region_total[Status == 'Accurate'], region_total[Status == 'Inaccurate']),
           Inaccuracy = region_total[Status == 'Inaccurate']/sum(region_total[Status == 'Accurate'], region_total[Status == 'Inaccurate']),
           Uncertainty = region_total[Status == 'Uncertain']/sum(region_total[Status == 'Accurate'], region_total[Status == 'Inaccurate'], region_total[Status == 'Uncertain']),
           Total = sum(region_total))
  x$User = user.id
  x$Stain = stain.id
  x$Brain = brain.id
  
  stats_final <- merge(stats, x, by = c('Structure', 'Status'))
  return(x)
})

stats_all <- bind_rows(stats_all)


##---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# filter dataset to include relevant accurate, inaccurate, and uncertainty measure
# filter to remove regions that were not sampled in all the brains (Total=! 0)

stats_all <- stats_all %>% 
  filter(Status == 'Accurate' | Status == 'Uncertain' |  Status == 'Inaccurate') %>% 
  filter(Total != 0) %>% 
  filter(Structure != 'Clear Label' & Structure != 'Total' & Structure != 'root')

##Remove under-represented regions (regions that do not have load values for all brains (colSum>0, but most of rows are 0))
stats_all <- stats_all %>%
  filter(Structure != "Subparafascicular area" & Structure != "Fasciola cinerea" & Structure != "Peripeduncular nucleus" & Structure != "Cerebellar cortex" & Structure != "Medulla, behavioral state related"
         & Structure != "Medulla, sensory related"  & Structure != "Medulla" & Structure != "Induseum griseum")

##Output analysis:
write.csv(stats_all, 'QCAlign_InterRaterTest.csv')

##Summarize data by region or region and brain:--------------------------------------------------------------------------------------------
rater_data <- stats_all

##see which regions/users/brains have uncertainty scores of 100% (aka no points contributing to accuracy and inaccuracy scores and input as NAs)
uncertain_rm <- rater_data[rater_data$Uncertainty == 1, ]   
uncertain_rm <- uncertain_rm[!is.na(uncertain_rm$Uncertainty),]

##-----FOR ACCURACY AND INACCURACY ASSESSMENT ONLY: Remove rater data with only 100 inaccuracy:-----
rater_data <- rater_data[rater_data$Uncertainty != 1, ]    

##-----FOR UNCERTAINTY MEASURE ONLY: REmove any rows that have NaN for all measures to have correct error bars:----
rater_data <- rater_data[!with(rater_data,is.na(Accuracy)& is.na(Uncertainty)),]

##set which QCAlign measure you want to look at (Accurate, Inaccurate, OR Uncertain):
status <- "Accurate"

##Filter by QCAlign measure:
rater_test <- rater_data %>% 
  filter(Status == status)

##Calculate average QCAlign scores per region:-----------------------------------------------------------------------------------------------

##input which QCAlign measure you are looking at (fill in accuracy, inaccuracy, or uncertainty):
rater_test = rater_test %>%
  group_by(Structure) %>%
  summarise(mean = mean(Accuracy, na.rm = TRUE),
            sd_load = sd(Accuracy),
            n_load = n(),
            se_load = sd_load/(sqrt(n_load)))


##rename Ammon's horn to avoid syntax error:
rater_test$Structure <- gsub("Ammon's horn", "Ammons horn", rater_test$Structure)

###Rename unassigned pixels of parent regions:
rater_test$Structure[rater_test$Structure == "Midbrain"] <- "Midbrain unassigned"
rater_test$Structure[rater_test$Structure == "Hippocampal formation"] <- "Hippocampal formation unassigned"
rater_test$Structure[rater_test$Structure == "Pons"] <- "Pons unassigned"
rater_test$Structure[rater_test$Structure == "Hindbrain"] <- "Hindbrain unassigned"
rater_test$Structure[rater_test$Structure == "Thalamus"] <- "Thalamus unassigned"
rater_test$Structure[rater_test$Structure == "Cortical subplate"] <- "Cortical subplate unassigned"

