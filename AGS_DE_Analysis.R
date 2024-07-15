{library('variancePartition')
  library('edgeR')
  library('BiocParallel')
  library('dplyr')
  library('tidyverse')
  library('readxl')
  library(lme4)
  library(lmerTest)
  library(MuMIn)
  library(emmeans)
  library('readxl')
  library(dplyr)
  library('edgeR')
  library('emmeans')
  library('limma')
  library(DESeq2)
  library(car)
  library(tidyverse)
  library(pheatmap)
  library(network)
  library(ggnetwork)
  library(colorspace)
  library(ggvenn)
  library(patchwork)
}


# Body weight graphs
setwd("~/Library/CloudStorage/OneDrive-Colostate/CSU/Exp2 with AGS brains")
AGS_data = read_excel("CSU_AGS_UltraLowInventory.xlsx")
AGS_data = AGS_data[, c(1:10)]
AGS_data = AGS_data[c(1:25), ]

AGS_data = AGS_data[-9, ]

str(AGS_data)
AGS_data$Treatment = as.factor(AGS_data$Treatment)
AGS_data$`Body weight (g) at TC` = as.numeric(AGS_data$`Body weight (g) at TC`)
AGS_data = AGS_data %>% rename("Body weight (g) at TC" = "BodyWeight_at_TC")
AGS_data = AGS_data %>% group_by(Treatment, Sex) %>% 
  mutate(avgBW = mean(BodyWeight_at_TC))
AGS_data$Sex = as.factor(AGS_data$Sex)

AGS_data = AGS_data %>% filter(Treatment == "LH" | Treatment == "2W")
AGS_data$Treatment = gsub("2W", "EH", AGS_data$Treatment)


ggplot(AGS_data, aes(x = Treatment, y = BodyWeight_at_TC)) + 
  geom_boxplot(aes(fill = Treatment), outlier.shape = NA) + 
  geom_jitter(pch = 21, size = 4) +
  theme_bw() + 
  facet_wrap(~Sex, labeller = labeller(Sex = c(M = "Male", F = "Female"))) + 
  theme(strip.text = element_text(size = 15), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14)) + 
  scale_fill_manual(values = c("darkblue", "lightyellow"))


BWlm = aov(BodyWeight_at_TC ~ Treatment*Sex, data = AGS_data)
qqp(resid(BWlm))
anova(BWlm)

maleChange = AGS_data %>% filter(Sex == "M") # males lost on average 94.9167 g BW
femaleChange = AGS_data %>% filter(Sex == "F") # females lost on average 40.4167 g BW

## Sample Info ####
setwd("~/Library/CloudStorage/OneDrive-Colostate/CSU/Exp2 with AGS brains/Sequencing/03_output")

AGS_metadata = read_xlsx("AGS_metadata.xlsx")
AGS_metadata = as.data.frame(AGS_metadata)
AGS_metadata$AGS_ID = as.factor(AGS_metadata$AGS_ID)
AGS_metadata$Group = as.factor(AGS_metadata$Group)
AGS_metadata$Sex = as.factor(AGS_metadata$Sex)
AGS_metadata$Age = as.factor(AGS_metadata$Age)
AGS_metadata$Sample_Date = as.factor(AGS_metadata$Sample_Date)
  # Create factors for months and years
  AGS_metadata$Sample_Date2 <- as.Date(AGS_metadata$Sample_Date, format = "%Y-%m-%d")
  AGS_metadata$month <- factor(months(AGS_metadata$Sample_Date2), levels = month.name)
  AGS_metadata$year <- factor(format(AGS_metadata$Sample_Date2, "%Y"))

AGS_metadata$Extraction_batch = as.factor(AGS_metadata$Extraction_batch)
AGS_metadata$Prep_batch = as.factor(AGS_metadata$Prep_batch)
AGS_metadata$Individual = c(1:14)
AGS_metadata$Individual = as.factor(AGS_metadata$Individual)

rownames(AGS_metadata) = AGS_metadata$AGS_ID
rownames(AGS_metadata) == AGS_metadata$AGS_ID

## AGS Readcounts ####
setwd("~/Library/CloudStorage/OneDrive-Colostate/CSU/Exp2 with AGS brains/Sequencing/03_output")

AGS_readcounts = read.table("AGS_HISATwoverlapmultimap_readcounts.txt", header = TRUE)
rownames(AGS_readcounts) = AGS_readcounts$Geneid
AGS_readcounts2 = AGS_readcounts
#write.csv(AGS_readcounts2, "AGS_readcounts2.csv")
AGS_readcounts = AGS_readcounts[, -c(1:6)]

AGS_readcounts$genes = rownames(AGS_readcounts)
AGS_readcounts1 = AGS_readcounts %>% filter(genes != 'Trh')
AGS_readcounts = AGS_readcounts1 %>% select(-genes)

#AGS_readcounts = AGS_readcounts[, rownames(AGS_metadata)]
#rownames(AGS_metadata) == colnames(AGS_readcounts)
#length(colnames(AGS_readcounts))

#set column names as metadata rownames
colnames(AGS_readcounts) = rownames(AGS_metadata)
colnames(AGS_readcounts) == rownames(AGS_metadata) #double check

## Make gene expression object
dge0 = DGEList(AGS_readcounts)
dim(dge0)
dge0$samples$group = AGS_metadata$Group
dge0$samples

keep = rowSums(cpm(dge0) > .5) >= 7 
dge <- dge0[keep,]
dim(dge)
dge = calcNormFactors(dge)
dge$samples$group = AGS_metadata$Group
dge$samples

#MDS plot
AGS_metadata$group = factor(paste0(AGS_metadata$Group, AGS_metadata$Sex))

plotMDS(dge,col=as.numeric(AGS_metadata$Group))

plotMDS(dge,col=as.numeric(AGS_metadata$group))
legend('topleft', fill=c('black', 'green', 'red', 'blue'), 
       legend=c("2WF", "LHF", '2WM', 'LHM'))

write.csv(dge, "dge.csv")

#dream approach
# estimate 
form <- ~ 0 + group + (1|Sample_Date)
param = SnowParam(8, "SOCK", progressbar=TRUE)
vobjDream = voomWithDreamWeights(dge, form, AGS_metadata, BPPARAM=param,plot = T)

#Extract contrast matrix for linear mixed model
L1 = makeContrastsDream(form, AGS_metadata, contrasts = 
                          c(group = "groupLHM - groupEHM"))
L2 = makeContrastsDream(form, AGS_metadata, contrasts = 
                          c("groupLHF - groupEHF"))
L3 = makeContrastsDream(form, AGS_metadata, contrasts = 
                          c("groupLHM - groupLHF"))
L4 = makeContrastsDream(form, AGS_metadata, contrasts = 
                          c("groupEHM - groupEHF"))

L = cbind(L1, L2, L3, L4)    
plotContrasts(L) 
#Fit the dream model on each gene
fitmm = dream(vobjDream, form, AGS_metadata, L, ddf = "Kenward-Roger")
fitmm = eBayes(fitmm)
head(fitmm$coefficients)
summary(decideTests(fitmm))

##
#save workspace image
#save.image(file = "dream/AGS_dream.RData")
load("dream/AGS_dream.RData")

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

# compare DE genes
#what are the shared genes between females and males from early to late hib
MaleTreatsig = MaleTreat[which(MaleTreat$adj.P.Val < 0.05), ]
MaleTreatsig$gene_id = rownames(MaleTreatsig)
MaleTreatsigRef = as.data.frame(rownames(MaleTreatsig))
rownames(MaleTreatsigRef) = MaleTreatsigRef$`rownames(MaleTreatsig)`

intersection = intersect(rownames(MaleTreatsigRef), rownames(AGS_readcounts2))
newdf = AGS_readcounts2[intersection, ]
MaleTreatsigRef = newdf$Chr
#write.csv(MaleTreatsigRef, 'MaleTreatsigRef.txt')

FemaleTreatsig = FemaleTreat[which(FemaleTreat$adj.P.Val < 0.05), ]
FemaleTreatsig$gene_id = rownames(FemaleTreatsig)

FemaleinMale = intersect(MaleTreatsig$gene_id, FemaleTreatsig$gene_id)
FemaleinMale # 3 shared genes

#Shared genes within sexes at 2 weeks and late hib
twoWeekSexsig = twoWeekSex[which(twoWeekSex$adj.P.Val < 0.05), ]
twoWeekSexsig$gene_id = rownames(twoWeekSexsig)

lateHibSexsig = lateHibSex[which(lateHibSex$adj.P.Val < 0.05), ]
lateHibSexsig$gene_id = rownames(lateHibSexsig)

WithinSex = intersect(twoWeekSexsig$gene_id, lateHibSexsig$gene_id)
WithinSex # 10 shared genes

# Loop through retograde genes in vobj list and run mixed model on them to see
# if they are significant. 

summaryRetro <- data.frame()
all_genes <- row.names(MaleTreat)
Treatment = as.factor(AGS_metadata$Group)
Sex = as.factor(AGS_metadata$Sex)
Sample_Date = as.factor(AGS_metadata$Sample_Date)
group = as.factor(AGS_metadata$group)

prioriGenes = c('Eya3', "Tshb", 'Dio2', 'Dio3', 'Kiss1', "Nes")

for (gene_id in prioriGenes) {
  
  test_data <- vobjDream$E[gene_id,]
  test_model <- lmer(test_data ~ Treatment*Sex + (1|(AGS_metadata$Sample_Date)))
  anova(test_model)
  output <- summary(pairs(emmeans(test_model, ~ Treatment*Sex), adjust = "none"))
  output_line <- output$p.value
  output_line[7] <- gene_id
  
  summaryRetro <- rbind(summaryRetro, output_line)
}

output_colnames <- output$contrast
colnames(summaryRetro) <- output_colnames
rownames(summaryRetro) = summaryRetro[,7] # make geneID's rownames
summaryRetro = summaryRetro[,-7]
#summarynonCorrect = summary
#write.csv(summaryRetro, "summaryRetro.csv")

# Venn diagrams of shared genes
library(eulerr)

# 2W to LH between sex
MaleTreat$gene = rownames(MaleTreat)
MaleTreatsig = MaleTreat[which(MaleTreat$adj.P.Val < 0.05), ]

FemaleTreat$gene = rownames(FemaleTreat)
FemaleTreatsig = FemaleTreat[which(FemaleTreat$adj.P.Val < 0.05), ]

vennData.all <- list(`EHM vs. LHM` = MaleTreatsig$gene, 
                     `EHF vs. LHF` = FemaleTreatsig$gene)

plot(euler(vennData.all[c('EHM vs. LHM','EHF vs. LHF')]), 
     quantities = list(cex = 2),alpha=0.7,edges='black',labels=F, fill = 
       c("darkblue", "lightyellow"), 
     legend = list(side = "top", nrow = 1, ncol = 3, 
                   labels = c("EH Male \nvs. LH Male", 'EH Female \nvs. LH Female'), 
                   cex = 1.5)) #only 4 sig genes shared

# 2WK and LH within sex
twoWeekSex$gene = rownames(twoWeekSex)
twoWeekSexsig = twoWeekSex[which(twoWeekSex$adj.P.Val < 0.05), ]

lateHibSex$gene = rownames(lateHibSex)
lateHibSexsig = lateHibSex[which(lateHibSex$adj.P.Val < 0.05), ]

vennData.all <- list(`EHM vs. EHF` = twoWeekSexsig$gene, 
                     `LHM vs. LHF` = lateHibSexsig$gene)

plot(euler(vennData.all[c('EHM vs. EHF','LHM vs. LHF')]), 
     quantities = list(cex = 2),alpha=0.7,edges='black',labels=F, fill = 
       c("yellow", "forestgreen"), 
     legend = list(side = "top", nrow = 1, ncol = 2, 
                   labels = c("EH Male \nvs. EH Female", 'LH Male \nvs. LH Female'), 
                   cex = 1.5)) #only 4 sig genes shared

# Plot individual genes ####
supp.labs <- c("Male", "Female")
names(supp.labs) <- c("M", "F")
treat.labs = c("Early Hib", "Late Hib")
names(treat.labs) = c("EH", "LH")

axis_label_size <- 30
title_size <- 30

gene_List = as.data.frame(t(vobjDream$E))
gene_List$Treatment = AGS_metadata$Group
gene_List$Sex = AGS_metadata$Sex
write.csv(gene_List, 'gene_List.csv')

geneName = 'LOC113193622'

plot_gene <- function(geneName) {
  # Assuming geneList is a data frame containing the necessary data
  plot_data <- gene_List
  
  ggplot(plot_data, aes(x = Treatment, y = plot_data[[geneName]])) + 
    geom_boxplot(aes(fill = Treatment), outlier.shape = NA) +
    geom_jitter(pch = 1, size = 3) + 
    facet_wrap(~Sex, labeller = labeller(Sex = supp.labs)) + 
    scale_fill_manual(values = c("darkblue", "lightyellow")) +  # Setting manual colors
    theme_classic() + 
    labs(y = "Normalized cpm") + 
    ggtitle("Vim") +
    theme(
      plot.title = element_text(size = title_size, hjust = 0.5),  # Adjust title size
      axis.text = element_text(size = axis_label_size),  # Adjust axis label size
      axis.title = element_text(size = axis_label_size), # Adjust axis title size
      strip.text = element_text(size = axis_label_size),
      legend.title = element_text(size = axis_label_size),
      legend.text = element_text(size = axis_label_size)
    )
}
plot = plot_gene(geneName)
plot