# load packages
{library(variancePartition)
  library(edgeR)
  library(BiocParallel)
  library(dplyr)
  library(tidyverse)
  library(readxl)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(readxl)
  library(dplyr)
  library(car)
  library(tidyverse)
  library(ggvenn)
  library(limma)
}

####
# This analysis was done on Bioconductor version 3.16 for R version 4.2.2. 
# This provides the analyses are completed with variancePartition (v1.28.9), 
# limma (3.54.2), and edgeR (3.40.2)
####

# Read in the metadata file
AGS_metadata = read_xlsx("AGS_metadata.xlsx")
AGS_metadata = as.data.frame(AGS_metadata)

# set factors
AGS_metadata$AGS_ID = as.factor(AGS_metadata$AGS_ID)
AGS_metadata$Group = as.factor(AGS_metadata$Group)
AGS_metadata$Sex = as.factor(AGS_metadata$Sex)
AGS_metadata$Age = as.factor(AGS_metadata$Age)
AGS_metadata$Sample_Date = as.factor(AGS_metadata$Sample_Date)
AGS_metadata$Individual = c(1:14)
AGS_metadata$Individual = as.factor(AGS_metadata$Individual)
rownames(AGS_metadata) = AGS_metadata$AGS_ID
rownames(AGS_metadata) == AGS_metadata$AGS_ID

# group Sex and hibernation timepoints as factors
AGS_metadata$group = factor(paste0(AGS_metadata$Group, AGS_metadata$Sex))

str(AGS_metadata) # check every variable is correctly formatted

# Read in readcount file
AGS_readcounts = read.csv("AGS_readcounts.csv", header = T)

AGS_readcounts$genes = AGS_readcounts$X
rownames(AGS_readcounts) = AGS_readcounts$genes
AGS_readcounts = AGS_readcounts[, -1] # remove first column
AGS_readcounts = AGS_readcounts[, -15] # remove last column

# make sure metadata and readcounts names match
#set column names as metadata rownames
colnames(AGS_readcounts) = rownames(AGS_metadata)
colnames(AGS_readcounts) == rownames(AGS_metadata) #double check, should be all TRUE

## Make gene expression object for dream and variancePartition
dge0 = DGEList(AGS_readcounts)
dim(dge0)
dge0$samples$group = AGS_metadata$group # add grouped factor to dge object
dge0$samples

keep = rowSums(cpm(dge0) > .5) >= 7 # filter counts in 50% of libraries
dge <- dge0[keep,]
dim(dge)
dge = calcNormFactors(dge)
dge$samples$group = AGS_metadata$group
dge$samples # double check

#dream approach
# estimate 
form <- ~ 0 + group + (1|Sample_Date) # formula for dream model

param = SnowParam(8, "SOCK", progressbar=TRUE)
vobjDream = voomWithDreamWeights(dge, form, AGS_metadata, BPPARAM=param,plot = T) # calculate gene expression values

#Extract contrast matrix for dream model
L1 = makeContrastsDream(form, AGS_metadata, contrasts = 
                          c(group = "groupLHM - groupEHM"))
L2 = makeContrastsDream(form, AGS_metadata, contrasts = 
                          c("groupLHF - groupEHF"))
L3 = makeContrastsDream(form, AGS_metadata, contrasts = 
                          c("groupLHM - groupLHF"))
L4 = makeContrastsDream(form, AGS_metadata, contrasts = 
                          c("groupEHM - groupEHF"))

L = cbind(L1, L2, L3, L4)    
plotContrasts(L) # double check all looks good

#Fit the dream model on each gene
fitmm = dream(vobjDream, form, AGS_metadata, L, ddf = "Kenward-Roger")
fitmm = eBayes(fitmm)
head(fitmm$coefficients)

# Extract contrast specific results
# 2week males late hib males
MaleTreat = topTable(fitmm, coef = 'group', sort.by = "P", n = Inf)
sum(MaleTreat$adj.P.Val < 0.05) #58
#write.csv(MaleTreat, 'dream/Male2wkTOlh.csv')

# 2 week females to late hib females
FemaleTreat = topTable(fitmm, coef = 'groupLHF - groupEHF', sort.by = "P", n = Inf)
sum(FemaleTreat$adj.P.Val < 0.05) #52
#write.csv(FemaleTreat, 'dream/Female2wkTOlh.csv')

# 2 week male to 2 week females
twoWeekSex = topTable(fitmm, coef = 'groupEHM - groupEHF', sort.by = "P", n = Inf)
sum(twoWeekSex$adj.P.Val < 0.05) #22
#write.csv(twoWeekSex, 'dream/twoWeekSex.csv')

# Late hibernation male to late hibernation female
lateHibSex = topTable(fitmm, coef = 'groupLHM - groupLHF', sort.by = "P", n = Inf)
sum(lateHibSex$adj.P.Val < 0.05) #128
#write.csv(lateHibSex, 'dream/lateHibSex.csv')

# At this point orthologs were obtained for unknown genes

# Extract and model a priori retrograde signaling genes
summaryRetro <- data.frame() #create empty data frame to populate
Sample_Date = as.factor(AGS_metadata$Sample_Date) #specify factors
group = as.factor(AGS_metadata$group) #specify factors

genes = as.data.frame(vobjDream$E) # Extract gene expression values

prioriGenes = c('Eya3', "Tshb", 'Dio2', 'Dio3', 'Kiss1', "Nes", 'Slc16a2') # specify the a priori genes to extract

for (gene_id in prioriGenes) {
  
  test_data <- vobjDream$E[gene_id,]
  test_model <- lmer(test_data ~ group + (1|(AGS_metadata$Sample_Date)))
  anova(test_model)
  output <- summary(pairs(emmeans(test_model, ~ group), adjust = "none"))
  output_line <- output$p.value
  
  output_line[7] <- gene_id
  
  summaryRetro <- rbind(summaryRetro, output_line)
}

output_colnames <- output$contrast
colnames(summaryRetro) <- output_colnames
rownames(summaryRetro) = summaryRetro[,7] # make geneID's rownames
