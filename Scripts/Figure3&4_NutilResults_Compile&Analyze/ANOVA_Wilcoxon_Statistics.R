#Script Objective: calculate statistics for age comparison and registration comparison of load values:
library(dplyr)

##set working directory as the location where EBRAINS data is downloaded:
setwd("/EBRAINS")

###ANOVA:--------------------------------------------------------------------------------------------------------------------------------------------------------

#Purpose: Evaluate differences in nonlinear stain load between 6m 5xFADs and 14m 5XFADs via ANOVA-------------------------------------------------------------------------------
##Read in long formatted data that has a stain load column and a column that specifies the registration type used ()
##This file is created in CompileNutilOutput_IntermediateRegions.R:
dat <- read.csv("IntermediateCustomRegions_Load_byRegistrationType.csv", row.names = 1)

##Subset data to be only 5XFADS (that's what the majority of the dataset is)
dat <- subset(dat, Genotype == 5)

##Remove under-represented regions (regions that do not have load values for all brains, see supplemental table 5)
dat <- dat %>%
  filter(Region != "Subparafascicular.area" & Region != "Fasciola.cinerea" & Region != "Peripeduncular.nucleus" & Region != "Cerebellar.cortex" & Region != "Medulla..behavioral.state.related"
         & Region != "Medulla..sensory.related"  & Region != "Medulla")

##set contrasts:
options(contrasts = c("contr.sum", "contr.poly"))

##list stains and regions
stains <- c("NeuN", "GFAP", "Iba1", "AB1.42", "Thionine")
regions <-  unique(dat$Region)

##Create dataframe for ANOVA output:
anova_df <- c() 

##loop to run anova for all stain and region combinations:
for (stain in stains) {
  for (region in regions) {
    tmp_df <- subset(dat, Region == region & Stain == stain)
    
    ## build a model :
    fit <- lm(Load.Percentage_NL ~ Age.Harvested + as.factor(Strain), data = tmp_df)
    
    #test normality:
    shapiro.test(residuals(fit))
    
    # fit using anova function:
    anovafit <- car::Anova(fit, type = 2)
    ## get sum squares:
    anovafitss <- anovafit$`Sum Sq`
    ## calculate percent variance explained:
    res <- cbind(anovafit, PctExp = anovafitss/sum(anovafitss) * 100)
    
    res$constrasts <- rownames(res)
    res$Region <- region
    res$Stain <- stain
    
    anova_tmp <- res
    anova_df = rbind(anova_df, anova_tmp)
    rm(anova_tmp)
    
  }
}

##save raw output from ANOVA:
write.csv(Age.Strain_ANOVA_FDR_Corrected, "LoadxAgexStrain_AnovaResults_AllRegions.csv")

##subset that dataset and run post hoc FDR correction per contrast type:

##Age:
Age_contrast <- subset(anova_df, constrasts == "Age.Harvested")
Age_contrast$FDR_adjusted_pval <- p.adjust(Age_contrast$'Pr(>F)',  method = 'fdr')

##Strain:
Strain_contrast <- subset(anova_df, constrasts == "as.factor(Strain)")
Strain_contrast$FDR_adjusted_pval <- p.adjust(Strain_contrast$'Pr(>F)',  method = 'fdr')

##Combine FDR corrected contrasts:
Age.Strain_ANOVA_FDR_Corrected <- rbind(Age_contrast, Strain_contrast)

##save FDR corrected output from ANOVA:
write.csv(Age.Strain_ANOVA_FDR_Corrected, "AllStain_NL_LoadxAgexStrain_AnovaResults_FDRCorrected.csv")

###Wilcoxon:--------------------------------------------------------------------------------------------------------------------------------------------------------

##Purpose: Complete Wilcoxon tests to evaluate the difference between stain load quantified via QuickNII registration only or quantified following QuickNII and VisuALign registration-------------------------------------------------------------------------
##"RG" is the output quantified following QuickNII registration only, "NL" is the output quantified following QuickNII and VisuALign registration

##Read in data:
##This file is created in CompileNutilOutput_IntermediateRegions.R:
reg_comparison <- read.csv("NLvsRG_IntermediateCustomRegions_Load.csv", row.names = 1)

##Subset data to be only 5XFADS and females (that's what the majority of the dataset is)
reg_comparison <- subset(reg_comparison, Genotype == 5)

##Run the Wilcoxon test per age:
reg_comparison <- subset(reg_comparison, Age.Harvested == 6)

##Remove under-represented regions (regions that do not have load values for all brains (colSum>0, but most of rows are 0))
reg_comparison <- reg_comparison %>%
  filter(Region != "Subparafascicular.area" & Region != "Fasciola.cinerea" & Region != "Peripeduncular.nucleus" & Region != "Cerebellar.cortex" & Region != "Medulla..behavioral.state.related"
         & Region != "Medulla..sensory.related"  & Region != "Medulla")

##list stains and regions
stains <- c("NeuN", "GFAP", "Iba1", "AB1.42", "Thionine")
regions <-  unique(reg_comparison$Region)

##create blank data frame to read results into:
wilcox_df <- data.frame()

for (stain in stains) {
  for (region in regions) {
    reg_temp <- subset(reg_comparison, Region == region & Stain == stain)
    N <- nrow(reg_temp)
    # wil_res <- wilcox.test(Load.Percentage ~ Registration, data = reg_temp, paired = TRUE)
    # By default (if exact is not specified), an exact p-value is computed if the samples contain less than 50 finite values and there are no ties. Otherwise, a normal approximation is used.
    wil_res <- wilcox.test(Load.Percentage ~ Registration, data = reg_temp, paired = TRUE, exact= FALSE)
    
    # Calculate the standardized Z statistic: 
    Zstat<-qnorm(wil_res$p.value/2)
    effect_size <- abs(Zstat)/sqrt(N)
    
    wilcox_tmp_df <- data.frame(pval = wil_res$p.value,
                                Region = region, 
                                Stain = stain,
                                test = wil_res$method,
                                comparison = wil_res$data.name,
                                test.stat= wil_res$statistic,
                                effect_size = effect_size)
    
    wilcox_df = rbind(wilcox_df, wilcox_tmp_df)
    rm(wilcox_tmp_df)
  }
}

##create column to identify which age was subset (change according to which age was subset in line 93):
wilcox_df$Age <- 6

##create new region_stain identifier column
wilcox_df$Region_Stain <- paste(wilcox_df$Region, wilcox_df$Stain, sep="_")

##save output of raw Wilcoxon test output:
write.csv(wilcox_df, 'RegionalRegistrationTypeComparison_Wilcoxon_AllSex_6m.csv')

##Run FDR correction of raw p-values:
wilcox_df$FDR_adjusted_pval <- p.adjust(wilcox_df$pval,  method = 'fdr')

##save FDR corrected output of Wilcoxon test:
write.csv(wilcox_df, 'RegionalRegistrationTypeComparison_Wilcoxon_AllSex_6m_FDRCorrected.csv')

##subset FDR significant results to assessment:
wilcox_df_sig <- wilcox_df %>%
  filter(FDR_adjusted_pval < 0.05)

##Create a table of regions and how many stains have significance differences at that region
Sig_regions <- as.data.frame(table(wilcox_df_sig$Region))



