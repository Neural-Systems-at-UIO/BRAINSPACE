library(tidyverse)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(ggfortify)

##Objective: Read in and calculate accuracy, inaccuracy, and uncertainty scores from output from QCAlign inter-rater reliability test stat.txt files:

##set working directory as the location where EBRAINS data is downloaded:
setwd("/EBRAINS")

##set path where data is retrieved from:
path <- getwd()

##List file names that meet the following pattern:
fname <- list.files(path, recursive = TRUE, pattern = "*_QControl_Thionine_stats(_SY|_BG|_AO|_HK|_IB|_MT|_SS|_TM|_TO|_US)")

##locate and read in QCALign output stats.txt files per user of the inter-rater reliability test performed (QuickNII registration only (QuickNII_stats) or QuickNII and VIsuALign registration(stats)):
stats_all <- lapply(fname, function(f){
  stats <- read.csv(file.path(path, f), sep = "\t")
  
  stain.id <- gsub('.*(Thionine|AB1-42|Iba1|GFAP|NeuN).*', '\\1', f)
  user.id <- gsub('.*\\_(QuickNII_stats_SY|QuickNII_stats_BG).*', '\\1', f) %>% 
    gsub('.*(Thionine\\_)(stats_SY|stats_BG|stats_AO|stats_HK|stats_IB|stats_MT|stats_SS|stats_TM|stats_TO|stats_US).*', '\\2', .)
  brain.id <- gsub('(*)\\_.*', '\\1', f)

#remove QCAlign damage results and organize according to Structure (aka intermediate hierarchy region):    
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
#Merge data:
stats_all <- bind_rows(stats_all)

# Filter dataset to include relevant accurate, inaccurate, and uncertainty measures and remove regions that were not sampled in all the brains, including the Clear label, Total, and root terms:

stats_all <- stats_all %>% 
  filter(Status == 'Accurate' | Status == 'Uncertain' |  Status == 'Inaccurate') %>% 
  filter(Total != 0) %>% 
  filter(Structure != 'Clear Label' & Structure != 'Total' & Structure != 'root')

##remove cerebellar and hypothalamic regions (regions removed during initial dissection (hypothalamus) and sectioning (cerebellum)):
rater_data_QckNII_all <- rater_data_QckNII_all %>%
  filter(Structure != "Hypothalamus" & Structure != "Cerebellar cortex" )

##create registration column (NL: QuickNII and VisuAlign registration evaluated, RG: QuickNII only registration evaluated)
stats_all$Registration <- NA
stats_all$Registration[grep('stats_BG', stats_all$User)] = 'NL'
stats_all$Registration[grep('stats_SY', stats_all$User)] = 'NL'
stats_all$Registration[grep('stats_AO', stats_all$User)] = 'NL'
stats_all$Registration[grep('stats_HK', stats_all$User)] = 'NL'
stats_all$Registration[grep('stats_IB', stats_all$User)] = 'NL'
stats_all$Registration[grep('stats_MT', stats_all$User)] = 'NL'
stats_all$Registration[grep('stats_SS', stats_all$User)] = 'NL'
stats_all$Registration[grep('stats_TM', stats_all$User)] = 'NL'
stats_all$Registration[grep('stats_TO', stats_all$User)] = 'NL'
stats_all$Registration[grep('stats_US', stats_all$User)] = 'NL'
stats_all$Registration[grep('QuickNII_stats_BG', stats_all$User)] = 'RG'
stats_all$Registration[grep('QuickNII_stats_SY', stats_all$User)] = 'RG'

##Save dataset:
write.csv(stats_all, 'QCAlign_InterRaterTest.csv')

##Summarize data by region:--------------------------------------------------------------------------------------------------------------------------------------------------------------------
#reassign dataframe name to avoid overwriting:
rater_data <- stats_all

##see which regions/users/brains have uncertainty scores of 100% (aka no points contributing to accuracy and inaccuracy scores and input as NAs)
uncertain_rm <- rater_data[rater_data$Uncertainty == 1, ]   
uncertain_rm <- uncertain_rm[!is.na(uncertain_rm$Uncertainty),]

##set which QCAlign measure you want to look at (Accurate, Inaccurate, OR Uncertain):
status <- "Accurate"

##**run line 93 if you are evaluating accuracy or inaccuracy scores, OR run line 96 if you are evaluating uncertainty scores:
##-----FOR ACCURACY AND INACCURACY ASSESSMENT ONLY: Remove rater data with 100% uncertainty:-----
rater_data <- rater_data[rater_data$Uncertainty != 1, ]    

##-----FOR UNCERTAINTY MEASURE ONLY: Remove any rows that have NaN for all measures to have correct error bars:----
rater_data <- rater_data[!with(rater_data,is.na(Accuracy)& is.na(Uncertainty)),]

##Filter by QCAlign measure (set in line 89):
rater_test <- rater_data %>% 
  filter(Status == status)

##Calculate average QCAlign scores per region and registration type for all raters and brains evaluated in the inter-rater reliability test:-------------------------------------------------------------------

##input which QCAlign measure you are looking at (fill in accuracy, inaccuracy, or uncertainty to match the measure input in line 89):
rater_test_summary = rater_test %>%
  group_by(Structure, Registration) %>%
  summarise(mean = mean(Accuracy, na.rm = TRUE),
            sd_load = sd(Accuracy),
            n_load = n(),
            se_load = sd_load/(sqrt(n_load)))

##rename Ammon's horn to avoid syntax error:
rater_test_summary$Structure <- gsub("Ammon's horn", "Ammons horn", rater_test_summary$Structure)

###Rename unassigned pixels of parent regions:
rater_test_summary$Structure[rater_test_summary$Structure == "Midbrain"] <- "Midbrain unassigned"
rater_test_summary$Structure[rater_test_summary$Structure == "Hippocampal formation"] <- "Hippocampal formation unassigned"
rater_test_summary$Structure[rater_test_summary$Structure == "Pons"] <- "Pons unassigned"
rater_test_summary$Structure[rater_test_summary$Structure == "Hindbrain"] <- "Hindbrain unassigned"
rater_test_summary$Structure[rater_test_summary$Structure == "Thalamus"] <- "Thalamus unassigned"
rater_test_summary$Structure[rater_test_summary$Structure == "Cortical subplate"] <- "Cortical subplate unassigned"

##save output of QCALign output averaged per measure (accuracy, inaccuracy, or uncertainty) and structure 
write.csv(rater_test_summary, paste0("Averaged", status,"Score_byStructure&RegistrationType.csv"))


registration_colors <- c("green3", "navy")

p <- ggplot(rater_test_summary, aes(x = factor(Structure, level = level_order), y = mean, color = Registration, shape = Registration))
p <- p + geom_point(aes(fill = as.factor(Registration), x = factor(Structure, level = level_order)), size = 1, position=position_dodge(width=0.99) ) +
  geom_errorbar(aes(ymin=mean-se_load, ymax=mean+se_load), width = 0.3, color= "black", position = position_dodge(0.99)) +
  theme_bw() + theme(axis.line = element_line(colour = "black", linewidth = 0.5, linetype = "solid")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_text(size=6, face="bold"),
        axis.title.y = element_text(size=6, face="bold"),
        axis.text.x = element_text(angle = 90, vjust = 1, hjust=1, size=5, color = "black"),
        axis.text.y = element_text(size=5,color = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 6),
        legend.position='bottom') +  scale_fill_manual(values= registration_colors) +
  scale_color_manual(values= registration_colors) +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1.0, .05)) +
  guides(color = guide_legend("Registration"), shape = guide_legend("Registration")) +
  labs(
    title = "Inter-User Test: Regional Accuracy Score",
    x = "Region",
    y = "Accuracy Score (%)"
  )+ 
  geom_point(data = rater_test_QckNII_all, aes(fill = as.factor(Registration), x = factor(Structure, level = level_order), 
                                               y= Accuracy), size = 1, alpha = .2, position=position_dodge(width=0.99) )

p

