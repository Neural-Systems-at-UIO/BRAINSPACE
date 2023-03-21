library(dplyr)
library(ggplot2)
library(correlation)
library(SummarizedExperiment)

##Objective: run a multilevel correlation with no age adjustment and age adjustment between hippocampal bulk RNAseq measures of gene expression and IHC measures of hippocampal cell and pathology load

setwd("~/DESeq_Test/MultiLevelCorr_Results/Correlation_BySample")
path <- "/DESeq_Test"

##dds object output from Option#2_DESeq_ModelIterations_SampleMatch.R:
  ##only samples that have both IHC and RNAseq
dat <- readRDS("/DESeq_Test/DESeq_ModelOutput_Interaction_BySample.rds")

meta <- dat$meta
counts <- dat$counts %>% as.data.frame()
vsd <- assay(dat$counts_vst) %>% as.data.frame()

stains <- c('Hippo_AB1.42', 'Hippo_NeuN', 'Hippo_GFAP', 'Hippo_Iba1')
# multilevel: If TRUE, the factors are included as random factors. Else, if FALSE (default), they are included as fixed effects in the simple regression model.
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
  cor_res_summ$group <- with(cor_res_summ, ifelse(adj.p < 0.05 & nadj.p > 0.05, 'load', 
                                                  ifelse(adj.p > 0.05 & nadj.p < 0.05, 'age', 
                                                         ifelse((adj.p < 0.05 & nadj.p < 0.05) & (nadj.p < adj.p), 'mostly age',
                                                                ifelse((adj.p < 0.05 & nadj.p < 0.05) & (nadj.p > adj.p), 'mostly load',
                                                                       ifelse(adj.p > 0.05 & nadj.p > 0.05, 'nonsig', 'other'))))))
  
  # cor_res_summ$adj.p_fdr <- p.adjust(cor_res_summ$adj.p,  method = 'fdr') 
  # cor_res_summ$nadj.p_fdr <- p.adjust(cor_res_summ$nadj.p,  method = 'fdr')
  
  
  ## evaluate results
  table(cor_res_summ$group)                           
  
  ## Save results
  res <- list(results_long = cor_res,
              results_wide = cor_res_summ)
  fname <- paste0(soi, '_correlations_summary_results_SampleMatch_TEST.rds')
  saveRDS(res, file.path(path, fname))
  options(warn=0)
}


