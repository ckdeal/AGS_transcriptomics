library(ggplot2)
library(ggvenn)
library(eulerr)
library(dplyr)
library(tidyverse)
library(readxl)
# Figure Plotting code for AGS manuscript

# read in genes with orthologs
MaleTreatsig_orth = read.csv("../DE_genes/Male_LHvsEH_sig.csv")
FemaleTreatsig_orth = read.csv("../DE_genes/Female_LHvsEH_sig.csv")
twoWeekSexsig_orth = read.csv("../DE_genes/earlyHib_MvsF_sig.csv")
lateHibSexsig_orth = read.csv("../DE_genes/lateHib_MvsF_sig.csv")

# Venn diagrams of shared genes
library(eulerr)

# 2W to LH wihtin sex
MaleTreatsigDE_up = MaleTreatsig_orth %>% dplyr::filter((logFC) >= 0.5) # 25 upregulated
MaleTreatsigDE_down = MaleTreatsig_orth %>% dplyr::filter((logFC) <= -0.5) # 8 downregulated
MaleTreatsigDE = MaleTreatsig_orth %>% dplyr::filter(abs(logFC) > 0.5) #33 up and downregulated w/ abs of 0.5

FemaleTreatsigDE_up = FemaleTreatsig_orth %>% dplyr::filter((logFC) >= 0.5) # 17 upregulated
FemaleTreatsigDE_down = FemaleTreatsig_orth %>% dplyr::filter((logFC) <= -0.5) # 24 downregulated
FemaleTreatsigDE = FemaleTreatsig_orth %>% dplyr::filter(abs(logFC) > 0.5) #41 up and downregulated w/ abs of 0.5

intersect(MaleTreatsigDE$gene, FemaleTreatsigDE$gene)

vennData.all <- list(`EHM vs. LHM` = MaleTreatsigDE$gene, 
                     `EHF vs. LHF` = FemaleTreatsigDE$gene)

# plot venn diagram for within differences
plot(euler(vennData.all[c('EHM vs. LHM','EHF vs. LHF')]), 
     quantities = list(cex = 4),alpha=0.7,edges='black',labels=F, fill = 
       c("darkblue", "lightyellow"), 
     legend = list(side = "top", nrow = 1, ncol = 3, 
                   labels = c("EH Male \nvs. LH Male", 'EH Female \nvs. LH Female'), 
                   cex = 1.5)) #only 4 sig genes shared. saved 7x8 pdf

# 2WK and LH between sex
twoWeekSexsigDE_up = twoWeekSexsig_orth %>% dplyr::filter((logFC) >= 0.5) # 11 upregulated
twoWeekSexsigDE_down = twoWeekSexsig_orth %>% dplyr::filter((logFC) <= -0.5) # 7 downregulated
twoWeekSexsigDE = twoWeekSexsig_orth %>% dplyr::filter(abs(logFC) >= 0.5) # 18 DE genes

intersect(twoWeekSexsigDE_up$gene, twoWeekSexsigDE_down$gene)

lateHibSexsig_up = lateHibSexsig_orth %>% dplyr::filter((logFC) >= 0.5) # 52 upregulated
lateHibSexsig_down = lateHibSexsig_orth %>% dplyr::filter((logFC) <= -0.5) # 44 downregulated
lateHibSexsigDE = lateHibSexsig_orth %>% dplyr::filter(abs(logFC) >= 0.5) # 96 DE genes

intersect(twoWeekSexsigDE$gene, lateHibSexsigDE$gene)

vennData.all <- list(`EHM vs. EHF` = twoWeekSexsigDE$gene, 
                     `LHM vs. LHF` = lateHibSexsigDE$gene)

# plot venn diagram for between differences
plot(euler(vennData.all[c('EHM vs. EHF','LHM vs. LHF')]), 
     quantities = list(cex = 4),alpha=0.7,edges='black',labels=F, fill = 
       c("yellow", "forestgreen"), 
     legend = list(side = "top", nrow = 1, ncol = 2, 
                   labels = c("EH Male \nvs. EH Female", 'LH Male \nvs. LH Female'), 
                   cex = 1.5)) #only 4 sig genes shared

# Plot individual genes ####
supp.labs <- c("Male", "Female")
names(supp.labs) <- c("M", "F")
treat.labs = c("Early Hib", "Late Hib")
names(treat.labs) = c("EH", "LH")

axis_label_size <- 25
title_size <- 25

# get gene expression values
AGS_metadata = read_xlsx("AGS_metadata.xlsx")
AGS_metadata = as.data.frame(AGS_metadata)

# set factors
AGS_metadata$AGS_ID = as.factor(AGS_metadata$AGS_ID)
AGS_metadata$Group = as.factor(AGS_metadata$Group)
AGS_metadata$Sex = as.factor(AGS_metadata$Sex)
AGS_metadata$group = factor(paste0(AGS_metadata$Group, AGS_metadata$Sex))
rownames(AGS_metadata) = AGS_metadata$AGS_ID
rownames(AGS_metadata) == AGS_metadata$AGS_ID

# create a gene list by trabnsposing the gene expression object and saving as data frame
load("AGS_dream.RData") # load in RData, so you do not need to rerun variancePartition

gene_List = as.data.frame(t(vobjDream$E)) # create gene list
gene_List$Treatment = AGS_metadata$Group
gene_List$Sex = AGS_metadata$Sex

# specify gene name to plot
geneName = 'Nes'

# Plot specific gene
plot_gene <- function(geneName) {
  # Assuming geneList is a data frame containing the necessary data
  plot_data <- gene_List
  
  ggplot(plot_data, aes(x = Treatment, y = plot_data[[geneName]])) + 
    geom_boxplot(aes(fill = Treatment), outlier.shape = NA) +
    geom_jitter(pch = 1, size = 3, stroke =1.5, width = 0.25) + 
    facet_wrap(~Sex, labeller = labeller(Sex = supp.labs)) + 
    scale_fill_manual(values = c("darkblue", "lightyellow")) +  # Setting manual colors
    theme_minimal() + 
    #ylim(5.9,6.9) + 
    labs(y = "Normalized cpm") + 
    ggtitle(geneName) +
    theme(
      plot.title = element_text(size = title_size, hjust = 0.5),  # Adjust title size
      axis.text = element_text(size = axis_label_size),# Adjust axis label size
      axis.title = element_text(size = axis_label_size), # Adjust axis title size
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.line.x = element_line(size = 1),  # Thicker x-axis line only
      axis.line.y = element_line(size = 1),   # Thicker y-axis line only
      strip.text = element_text(size = axis_label_size),
      legend.title = element_text(size = axis_label_size),
      legend.text = element_text(size = axis_label_size)
    )
}
plot = plot_gene(geneName)
plot # 7x5 as pdf
