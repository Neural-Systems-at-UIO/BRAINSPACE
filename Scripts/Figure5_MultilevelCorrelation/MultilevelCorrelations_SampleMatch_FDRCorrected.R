library(dplyr)
library(ggplot2)
library(correlation)
library(SummarizedExperiment)
library("AnnotationDbi")
library("org.Mm.eg.db")

#Objective: Use output from the default R/DESeq2 pipeline to run a multilevel partial correlation adjusting for age or not adjusting for age:----------------------------  

##Set working directory to the location where the multi-level correlation input data from GitHub was downloaded:
setwd("~/MultiLevelCorr_Results/Correlation_BySample")
path <- "/Multilevel_Corr_Output"

## Output from the default R/DESeq2 pipeline. Only samples that have both IHC and RNAseq were included, n=34 when running DESeq2():
dat <- readRDS(~"/DESeq_ModelOutput_Interaction_BySample.rds")

#Separate the multiple components of the DESeq2 object:
meta <- dat$meta
counts <- dat$counts %>% as.data.frame()
vsd <- assay(dat$counts_vst) %>% as.data.frame()

##set region of interest: Hippocampal formation (Hippo) load per stain
stains <- c('Hippo_AB1.42', 'Hippo_NeuN', 'Hippo_GFAP', 'Hippo_Iba1')

#Run multi level correlation with age adjustment and without age adjustment--------------------------------------------------------------------------------------------
# multilevel: If TRUE, the factors are included as random factors. Else, if FALSE (default), they are included as fixed effects in a simple regression model
for(soi in stains){
  options(warn=-1)
  cor_res <- c()
  gene_list <- rownames(counts)
  for(goi in gene_list){
    expr_vec <- counts[goi,] %>% unlist()
    meta$expr <- expr_vec[match(meta$SampleID, names(expr_vec))]
    
    tmp <- meta[,c(soi, 'expr', 'Age.Harvested')]
    lm <- correlation(tmp) %>% as.data.frame()
    lm$gene <- goi
    lm$method <- 'nadj'
    lmm <- correlation(tmp, multilevel = TRUE) %>% as.data.frame()
    lmm$gene <- goi
    lmm$method <- 'adj'
    
    cor_res <- plyr::rbind.fill(cor_res, lm, lmm)
  }
  
  cor_res_summ.r <- reshape2::dcast(gene ~ method, value.var = 'r', data = cor_res) %>% 
    setNames(c('Gene', 'adj.r', 'nadj.r'))
  cor_res_summ.p <- reshape2::dcast(gene ~ method, value.var = 'p', data = cor_res) %>% 
    setNames(c('Gene', 'adj.p', 'nadj.p'))
  
  
  head(cor_res_summ.p)
  head(cor_res_summ.r)
  cor_res_summ <- merge(cor_res_summ.r, cor_res_summ.p, by = 'Gene')
  
  #Create a group term to classify correlations based on nominal p-value significance level
  cor_res_summ$group <- with(cor_res_summ, ifelse(adj.p < 0.05 & nadj.p > 0.05, 'load', 
                                                  ifelse(adj.p > 0.05 & nadj.p < 0.05, 'age', 
                                                         ifelse((adj.p < 0.05 & nadj.p < 0.05) & (nadj.p < adj.p), 'mostly age',
                                                                ifelse((adj.p < 0.05 & nadj.p < 0.05) & (nadj.p > adj.p), 'mostly load',
                                                                       ifelse(adj.p > 0.05 & nadj.p > 0.05, 'nonsig', 'other'))))))
  
  ##Evaluate results
  table(cor_res_summ$group)                           
  
  ##Save results
  res <- list(results_long = cor_res,
              results_wide = cor_res_summ)
  fname <- paste0(soi, '_correlations_summary_results_SampleMatch.rds')
  saveRDS(res, file.path(path, fname))
  options(warn=0)
}

##Complete FDR correction of p-values:----------------------------------------------------------------------------------------
## Per each stain: Read in the files created from multi-level corrrelation loop above, correct non-age-adjusted and age-adjusted p-values for multiple comparisons using FDR correction, and output multi-level correlation outout with added FDR-corrected p-values  

##Iba1:
Iba1 <- readRDS("/Multilevel_Corr_Output/Hippo_Iba1_correlations_summary_results_SampleMatch.rds")
Iba1 <- Iba1$results_wide
Iba1$adj.p_fdr <- p.adjust(Iba1$adj.p,  method = 'fdr') 
Iba1$nadj.p_fdr <- p.adjust(Iba1$nadj.p,  method = 'fdr') 
Iba1$Stain <- "Iba1"
# write.csv(Iba1, file.path(path, "MultilevelCorrOutput_Iba1.csv"))

#NeuN:
NeuN <- readRDS("/Multilevel_Corr_Output/Hippo_NeuN_correlations_summary_results_SampleMatch.rds")
NeuN <- NeuN$results_wide
NeuN$adj.p_fdr <- p.adjust(NeuN$adj.p,  method = 'fdr') 
NeuN$nadj.p_fdr <- p.adjust(NeuN$nadj.p,  method = 'fdr') 
NeuN$Stain <- "NeuN"
# write.csv(NeuN, file.path(path, "MultilevelCorrOutput_NeuN.csv"))

#GFAP:
GFAP <- readRDS("/Multilevel_Corr_Output/Hippo_GFAP_correlations_summary_results_SampleMatch.rds")
GFAP <- GFAP$results_wide
GFAP$adj.p_fdr <- p.adjust(GFAP$adj.p,  method = 'fdr') 
GFAP$nadj.p_fdr <- p.adjust(GFAP$nadj.p,  method = 'fdr') 
GFAP$Stain <- "GFAP"
# write.csv(GFAP, file.path(path, "MultilevelCorrOutput_GFAP.csv"))

#AB1-42
AB1.42 <- readRDS("/Multilevel_Corr_Output/Hippo_AB1.42_correlations_summary_results_SampleMatch.rds")
AB1.42 <- AB1.42$results_wide
AB1.42$adj.p_fdr <- p.adjust(AB1.42$adj.p,  method = 'fdr') 
AB1.42$nadj.p_fdr <- p.adjust(AB1.42$nadj.p,  method = 'fdr') 
AB1.42$Stain <- "AB1.42"
# write.csv(AB1.42, file.path(path, "MultilevelCorrOutput_AB1.42.csv"))

##Merge all FDR-corrected output:

###Compile the multi-level correlation resutls for all stains:-----------------------------------------------------------------------------------------------------------------------------------------------
all_stains <- rbind(Iba1, NeuN, GFAP, AB1.42)
all_stains <- all_stains[,-c(1)]

##Create comparison (adj or nadj) and direction groups:

##create a column summarizing correlations by positive or negative r values:
all_stains$Direction<- ifelse(all_stains$adj.r > 0 | all_stains$nadj.r > 0, "Positive",
                              ifelse(all_stains$adj.r < 0 | all_stains$nadj.r < 0, "Negative", "None" ))

##create a column based on nominal significance values (this "Group2" term combines correlations that were significant before and after age adjustment)
all_stains$Group2<- ifelse(all_stains$group == "mostly age" | all_stains$group == "mostly load", "Both",
                           ifelse(all_stains$group == "age" , "Age",
                                  ifelse(all_stains$group == "nonsig" , "NonSig","Load" )))

##Convert ensembl gene IDs to gene names:
all_stains$GeneName <- mapIds(org.Mm.eg.db, as.character(all_stains$Gene), 'SYMBOL', 'ENSEMBL')

#Evaluate a snippet of the data frame before plotting:
head(all_stains)

##Save compiled output:
write.csv(all_stains, "MultilevelCorrOutput_FDRAdjusted_AllStains.csv")
