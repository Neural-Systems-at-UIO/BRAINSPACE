library(utils)
library(reshape2)
library(dpyr)

#Script Objective: Compile GSEA results output from WebGestalt. Results are saved according to the stain and correlation method (age-adjusted or not age-adjusted):

##set working directory as the directory location of WebGestalt output downloaded from GitHub:
setwd("/WebGestaltOutput_SampleMatch/FDR0.05")
path <- getwd()

##Unzip WebGestalt outout and list all files
fs <- list.files(path, pattern = glob2rx("*.zip"))

#create dataframe to save results into:
dat <- c()

#read in files:
for (f in fs) {
  fname <- file.path(path, f)
  Stain <- gsub("(\\w*)_.*", "\\1", f)
  Group <- gsub("\\w*_(\\w*).*", "\\1", f)
  zipped_csv_names <- grep('enrichment_results*', unzip(fname, list=TRUE)$Name, 
                           ignore.case=TRUE, value=TRUE)
  dat_tmp <- read.csv(unz(fname, zipped_csv_names), header=T, sep="\t")
#create Group column to list whether the results were from the age adjusted or non-age-adjusted multilevel correlation results input into WebGestalt:
  dat_tmp$Group <- Group
#create Stain column to list which stain the WebGestalt results are associated with:
  dat_tmp$Stain <- Stain
  dat <- rbind(dat, dat_tmp)
}

##select relevant columns ("geneSet","description", enrichmentScore","normalizedEnrichmentScore", "pValue", "FDR", "Group", "Stain")
df <- dat[,c(1,2,4:7,13:14)]

##reorganize the data frame: normalized enrichment score calculated by WebGestalt are organized per stain and correlation adjustment method (age-adjusted or non-age-adjusted):
df <- reshape2::dcast(geneSet * description ~ Group * Stain, data = dat, value.var = "normalizedEnrichmentScore")

#convert any NA values to 0:
df[is.na(df)] = 0

##read in a file that classifies child Reactome pathway descriptions into parent terms (download from Github):
rpa <- read.csv("WebGestaltOutput_SampleMatch/ReactomePA_PatenChildTable_v3.csv")

##Merge dataframes: Classify GSEA pathway results by parent term
x <- merge(df,rpa, by.x = 'description', by.y = 'child.description', all.x = T)

##save output for plotting:
write.csv(x, "WebGestaltOutput_SampleMatch_FDR005.csv")
