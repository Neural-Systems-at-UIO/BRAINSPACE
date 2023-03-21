# Objective: calculate stats for age comparison and registration comparison of load values
library(dplyr)

setwd("Z:/Brain-wide_cell-type_maps/NeuroAssociates_20x_images_AD-BXDandNtg-BXD/Brains1-40/key and beh data/Pilot_Manuscript/StainData_IntermediateCustom/Load_Dataframes")
#-------------------------------------------------------------------------------
##Read in long data with a load column and a column that specifies the registration type used
##input location of NL-RG_RegistrationDiff_ParentRegion.csv file in line 9:

dat <- read.csv("/NL-RG_RegistrationDiff_ParentRegion.csv")
dat <- dat[, c(-1)]

##Subset data to be only 5XFADS and females (that's what the majority of the dataset is)
dat <- subset(dat, Genotype == 5)

##removing regions that do not have load values for every brain (under-represented regions in the dataset)
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
    
    # fit using anova:
    anovafit <- car::Anova(fit, type = 2)
    ## get sum squares:
    anovafitss <- anovafit$`Sum Sq`
    ## calc variance explained:
    res <- cbind(anovafit, PctExp = anovafitss/sum(anovafitss) * 100)
    
    res$constrasts <- rownames(res)
    res$Region <- region
    res$Stain <- stain
    
    anova_tmp <- res
    anova_df = rbind(anova_df, anova_tmp)
    rm(anova_tmp)
    
  }
}

# write.csv(Age.Strain_ANOVA_FDR_Corrected, "LoadxAgexStrain_AnovaResults_AllRegions.csv")

##subset dataset and run post hoc FDR correct per contrast type:

##Age:
Age_contrast <- subset(anova_df, constrasts == "Age.Harvested")
Age_contrast$FDR_adjusted_pval <- p.adjust(Age_contrast$'Pr(>F)',  method = 'fdr')
##Strain:
Strain_contrast <- subset(anova_df, constrasts == "as.factor(Strain)")
Strain_contrast$FDR_adjusted_pval <- p.adjust(Strain_contrast$'Pr(>F)',  method = 'fdr')
##Combine corrected values:
Age.Strain_ANOVA_FDR_Corrected <- rbind(Age_contrast, Strain_contrast)

write.csv(Age.Strain_ANOVA_FDR_Corrected, "AllStain_Hippo_LoadxAgexStrain_AnovaResults_FDRCorrected.csv")


##Wilcoxon test to evaluate the difference between load of 2 different registration--------------------------------------------------------------------------

##Read in data:
##input location of NLvsRG_IntermediateCustomRegions_Load.csv file in line 80 (this file is created in CompileNutilOutput_IntermediateRegions.R):

reg_comparison <- read.csv("/NLvsRG_IntermediateCustomRegions_Load.csv")
reg_comparison <- reg_comparison[, c(-1)]

##Subset data to be only 5XFADS and females (that's what the majority of the dataset is)
reg_comparison <- subset(reg_comparison, Genotype == 5 & Age.Harvested == 6)


##Remove under-represented regions (regions that do not have load values for all brains (colSum>0, but most of rows are 0))
reg_comparison <- reg_comparison %>%
  filter(Region != "Subparafascicular.area" & Region != "Fasciola.cinerea" & Region != "Peripeduncular.nucleus" & Region != "Cerebellar.cortex" & Region != "Medulla..behavioral.state.related"
         & Region != "Medulla..sensory.related"  & Region != "Medulla")


median(reg_comparison$Load.Percentage)

##set factors:
stains <- c("NeuN", "GFAP", "Iba1", "AB1.42", "Thionine")
regions <-  unique(reg_comparison$Region)

wilcox_df <- data.frame()
for (stain in stains) {
  for (region in regions) {
    reg_temp <- subset(reg_comparison, Region == region & Stain == stain)
    N <- nrow(reg_temp)
    # wil_res <- wilcox.test(Load.Percentage ~ Registration, data = reg_temp, paired = TRUE)
    wil_res <- wilcox.test(Load.Percentage ~ Registration, data = reg_temp, paired = TRUE, exact= FALSE)
    
    # By default (if exact is not specified), an exact p-value is computed if the samples contain less than 50 finite 
    # values and there are no ties. Otherwise, a normal approximation is used.
    
    # Calculate the standardised z statistic Z 
    Zstat<-qnorm(wil_res$p.value/2)
    # Calculate the effect size, 275= number of paired test (55 regions, 5 stains= 275) 
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

wilcox_df$FDR_adjusted_pval <- p.adjust(wilcox_df$pval,  method = 'fdr')
wilcox_df$Region_Stain <- paste(wilcox_df$Region, wilcox_df$Stain, sep="_")

wilcox_df_sig <- wilcox_df %>%
  filter(FDR_adjusted_pval < 0.05)

##write out according to what was subset up front (just genotype or genotype and sex)
write.csv(wilcox_df, 'RegionalRegistrationTypeComparison_Wilcoxon_AllSex_6m.csv')

##table of regions and how many stains have significance differences at that region
regions <- as.data.frame(table(wilcox_df_sig$Region))

