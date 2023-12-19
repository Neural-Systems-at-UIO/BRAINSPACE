library(ggplot2)
library(ggpmisc)
library(dplyr)
library(ggfortify)
library(sva)
library(DESeq2)

##Objective: Integrate IHC hippocampal load and gene expression data in scatterplots 

##set working directory as the directory location 'DESeq_ModelOutput_Intercept_BySample.rds' was downloaded from GitHub in 
setwd("~/MultiLevelCorr_Results_BySample")

##Set path where plots will be output:
path <- "/Multilevel_Corr_Output"

###Gene x Age and Age x Load correlation plots:--------------------------------------------------------------------------------------------------------
## Representative code applied to batch correct gene expression data: Formatting of DESeq2 intercept model used to create the input files of this script (file input on line 28). CountData consist of all gene expression count data for the 34 samples included in this analysis. Group_key consists of the hippocampal formation stain load and metadata for the 34 samples included in this analysis.
# dds_intc = DESeqDataSetFromMatrix(countData = group_count,
#                                   colData   = group_key,
#                                   design    = ~  1)
# 
# #run DESeq2:
# dds_intc = DESeq(dds_intc)
# 
# # view_results_names:
# resultsNames(dds_intc)
##that model outputs:
dds_intc <- readRDS("DESeq_ModelOutput_Intercept_BySample.rds")

##read in and organize meta file that consists of the hippocampal formation stain load and metadata for the 34 samples included in this analysis:
group_key <- read.csv("/Multilevel_Corr_Output/SampleMatch_GroupKey.csv", row.names = 1)
group_key$Age.Harvested <- factor(group_key$Age.Harvested, levels = c(6,14))
group_key$Strain <- factor(group_key$Strain)
group_key$Sex <- factor(group_key$Sex, levels = c(11,22))

##Plot:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#set plot aesthetics:
theme = theme_classic() +
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 5),
        #axis.title.x = element_blank(),
        axis.text.x=element_text(angle=0, hjust=0.5, vjust = 0.25, size = 5),
        axis.text.y=element_text(angle=0, hjust=0.5, vjust = 0.25, size = 5, color = "black"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()) +
  theme(axis.line = element_line(color = 'black')) +
  theme(legend.position = "none",
        legend.title=element_blank(),
        legend.text = element_text(size = 5)) +
  theme(plot.title=element_text(hjust=0.5))+
  theme(axis.line=element_line(size=0.25))

age_colors <- c("blue2", "goldenrod1")


##Define functions to plot DESeq2 output data as boxlplots of hippocampal stain load by age and scatterplots relating hippocampal stain load and gene expression (Figure 4d)-----------------------------------------------------------------------------

## Cor plot function:
corPlot <- function(deseq_dataset, meta_file,  goi, soi, print = c(TRUE,FALSE), assay = c('normcounts', 'vst')){
  if(assay == 'normcounts'){
    tmp <- counts(deseq_dataset, normalized = TRUE)
  } else if (assay == 'vst'){
    tmp <- DESeq2::vst(deseq_dataset, blind=FALSE)
    tmp <- assay(tmp)
  }
  
  g_vec <- tmp[goi,]
  meta_file$expr <- g_vec[match(names(g_vec), rownames(meta_file))]
  
  p1 <- ggplot(meta_file, aes(x = (100*meta_file[,soi]), y = (log2(expr+1)))) + 
    geom_point(size = 0.5) + 
    geom_smooth(method = 'lm', color = "gray40") +
    labs(y= paste(goi, "Transformed Counts"),
         x = paste(soi,"Load (%)")) + theme +
    # scale_y_continuous(limits = c(min(log2(meta_file[,'expr']+1)), max(log2(meta_file[,'expr']+1))), breaks = seq(0, max(log2(meta_file[,'expr']+1)), by = 0.1))+
    # scale_y_continuous(limits = c(5.9,9.6), breaks = seq(5.9,9.6, by = 0.1))+
    scale_color_manual(values=rev(age_colors))
  
  
  
  p2 <- ggplot(meta_file, aes(x = (100*meta_file[,soi]), y = (log2(expr+1)), color = Age.Harvested)) + 
    geom_point(size = 0.5) + 
    geom_smooth(method = 'lm') +
    labs(y= paste(goi, "Transformed Counts"),
         x = paste(soi,"Load (%)"))+ theme +
    scale_y_continuous(limits = c(min(log2(meta_file[,'expr']+1)), max(log2(meta_file[,'expr']+1))), breaks = seq(0, max(log2(meta_file[,'expr']+1)), by = 0.1))+
    # scale_x_continuous(limits = c(0, max(100*meta_file[,soi]+1)), breaks = seq(0, max(100*meta_file[,soi]+1), by = 1))+
    scale_color_manual(values=rev(age_colors))
  
  
  p_final <- cowplot::plot_grid(p1, p2, rel_widths = c(1,1))
  return(p_final)
}

###

## Age plot function:
AgePlot <- function(deseq_dataset, meta_file,  goi, soi, print = c(TRUE,FALSE), assay = c('normcounts', 'vst')){
  if(assay == 'normcounts'){
    tmp <- counts(deseq_dataset, normalized = TRUE)
  } else if (assay == 'vst'){
    tmp <- DESeq2::vst(deseq_dataset, blind=FALSE)
    tmp <- assay(tmp)
  }
  
  g_vec <- tmp[goi,]
  meta_file$expr <- g_vec[match(names(g_vec), rownames(meta_file))]
  
  meta_file$expr
  seq(0, max(log2(meta_file$expr+1)))
  
  meta_file$expr
  meta_file$expr + 1
  log2(meta_file$expr + 1)
  
  p1 <- ggplot(meta_file, aes(x = Age.Harvested, y = (100*meta_file[,soi]), color=factor(Age.Harvested))) + 
    geom_boxplot(lwd=0.25) +
    geom_point(size = 0.5) +
    theme_bw() +
    theme +
    labs(x= "Age Harvested", y=paste(soi,"Load (%)")) +
    scale_y_continuous(limits = c(min(100*meta_file[,soi]),  max(100*meta_file[,soi])), breaks = seq(0,  max(100*meta_file[,soi]), by = 2))+
    scale_color_manual(values=rev(age_colors))
  
  p2 <- ggplot(meta_file, aes(x = Age.Harvested, y = (log2(expr+1)), color=factor(Age.Harvested))) + 
    geom_boxplot(lwd=0.25) +
    geom_point(size = 0.5) +
    theme_bw() +
    theme + 
    labs(x= "Age Harvested", y= paste(goi, "Transformed Counts")) +
    scale_y_continuous(limits = c(min(log2(meta_file[,'expr']+1)), max(log2(meta_file[,'expr']+1))), breaks = seq(0, max(log2(meta_file[,'expr']+1)), by = 0.1))+
    scale_color_manual(values=rev(age_colors))
  
  p_final <- cowplot::plot_grid(p1, p2, rel_widths = c(1,1))
  return(p_final)
}


##set gene and stain of interest:---------------------------------------------------------------------------------------
##top Iba1 load gene
# #tmem39a:
g <- "ENSMUSG00000002845"

##top Iba1 age gene
# ##Galnt6
g <- "ENSMUSG00000037280"

#set stain of interest
soi <- 'Hippo_Iba1'

###########plotCorPlot----------------------
p <- corPlot(deseq_dataset = dds_intc, meta_file = group_key, goi = g, soi = soi, assay = 'normcounts')
p

ggsave(file=paste(path,soi, "Osmr_Genecorplot_tranformed.png", sep = "_"), p, width = 3, height = 1.8, units = 'in', dpi = 300, bg= NULL, device = png)

ggsave(file=paste(path,soi, "Osmr_Genecorplot_tranformed.pdf", sep = "_"), p, width = 3, height = 1.8, units = 'in', dpi = 300, bg= NULL)


###########plotAgePlot----------------------
p <- AgePlot(deseq_dataset = dds_intc, meta_file = group_key, goi = g, soi = soi, assay = 'normcounts')
p

ggsave(file=paste(path,soi, "Osmr_AgePlot_tranformed.png", sep = "_"), p, width = 3, height = 1.8, units = 'in', dpi = 300, bg= NULL, device = png)

ggsave(file=paste(path,soi, "Osmr_Ageplot_tranformed.pdf", sep = "_"), p, width = 3, height = 1.8, units = 'in', dpi = 300, bg= NULL)
