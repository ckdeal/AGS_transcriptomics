setwd("DEGenes_wOrthologs")

MaleTreatSig  = read.csv("Male_LHvsEH_sig.csv")
FemaleTreatsig = read.csv("Female_LHvsEH_sig.csv")
lateHibSexSig = read.csv("lateHib_MvsF_sig.csv")
earlyHibSexSig = read.csv("earlyHib_MvsF_sig.csv")


#filter out sig genes
MaleTreatsig = MaleTreatSig[which(MaleTreatSig$adj.P.Val < 0.05), ] 
sum(MaleTreatsig$adj.P.Val < 0.05 & MaleTreatsig$logFC >= 0.05) # 40 upregulated genes
sum(MaleTreatsig$adj.P.Val < 0.05 & MaleTreatsig$logFC <= 0.05) # 18 downregulated genes
rownames(MaleTreatsig) = MaleTreatsig$X

FemaleTreatsig = FemaleTreatsig[which(FemaleTreatsig$adj.P.Val < 0.05), ]
sum(FemaleTreatsig$adj.P.Val < 0.05 & FemaleTreatsig$logFC >= 0.05) # 24 upregulated genes
sum(FemaleTreatsig$adj.P.Val < 0.05 & FemaleTreatsig$logFC <= 0.05) # 28 downregulated genes
rownames(FemaleTreatsig) = FemaleTreatsig$X

earlyHibSexSigGenes = earlyHibSexSig[which(earlyHibSexSig$adj.P.Val < 0.05), ]
sum(earlyHibSexSigGenes$adj.P.Val < 0.05 & earlyHibSexSigGenes$logFC >= 0.05) # 12 upregulated genes in males
sum(earlyHibSexSigGenes$adj.P.Val < 0.05 & earlyHibSexSigGenes$logFC <= 0.05) # 10 downregulated genes in females
rownames(earlyHibSexSigGenes) = earlyHibSexSigGenes$X

lateHibSexSigGenes = lateHibSexSig[which(lateHibSexSig$adj.P.Val < 0.05), ]
sum(lateHibSexSigGenes$adj.P.Val < 0.05 & lateHibSexSigGenes$logFC >= 0.05) # 59 upregulated genes in males
sum(lateHibSexSigGenes$adj.P.Val < 0.05 & lateHibSexSigGenes$logFC <= 0.05) # 69 downregulated genes in females

rownames(lateHibSexSigGenes) = lateHibSexSigGenes$X

# Pathway Analysis ####
## Get gene orthologs ####
library(biomaRt)
# Choose an Ensembl mirror
ensembl_mirror <- useEnsembl(biomart = "ensembl", mirror = "uswest")
mmus <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org")
urop = useMart("ensembl", dataset = "uparryii_gene_ensembl", host = "https://dec2021.archive.ensembl.org")

your_genes = as.vector(lateHibSexDWNGenes) ## enter the vector of gene names here, replace with specific pairwise comparison

mouse_genes <- getLDS(attributes = c("external_gene_name"),
                      filters = "external_gene_name",
                      values = your_genes,
                      mart = urop,
                      attributesL = c("external_gene_name"),
                      martL = mmus)

setdiff(your_genes, mouse_genes$Gene.name.1)

# run enrichPathway ####
library(clusterProfiler) 
library(org.Mm.eg.db) 
library(ReactomePA)

dge = read.csv("dge.csv") # read in dream() gene expression object 
AGSbackground = as.vector(rownames(dge))
AGSbackgroundENTREZ = bitr(AGSbackground, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
AGSbackgroundENTREZ = as.vector(AGSbackgroundENTREZ$ENTREZID)
#
Genes = bitr(earlyHibSexDWNGenes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db") # Change gene list to acquire specific IDs for what pairwise comparison you are interested in 
Genes = as.vector(Genes$ENTREZID)

x <- enrichPathway(gene=Genes, pvalueCutoff = 0.05, 
                   qvalueCutoff = .2,
                   universe = AGSbackgroundENTREZ, pAdjustMethod = "BH", 
                   organism = "mouse", minGSSize = 1, maxGSSize = 500)

x = as.data.frame(x)
write.csv(x, ".csv")

## GO analysis ####
library(gprofiler2)
AGSbackground = as.vector(rownames(dge))

#GO between 2 week to LH in males
MaleTreatNEW = read.csv("Male_LHvsEH_sig.csv")
MaleTreatGO = MaleTreatNEW[which(MaleTreatNEW$adj.P.Val < 0.05), ]
rownames(MaleTreatGO) = MaleTreatGO$X
MaleTreatGO$gene = rownames(MaleTreatGO)
MaleTreatGOres = MaleTreatGO$logFC
names(MaleTreatGOres) = MaleTreatGO$gene
MaleTreatGOres = sort(MaleTreatGOres, decreasing = T)
MaleTreatGOres = as.matrix(MaleTreatGOres)
MaleTreatGOres = rownames(MaleTreatGOres)
MaleTreatGOres = as.vector(MaleTreatGOres)

MaleTreatGO = MaleTreatNEW[which(MaleTreatNEW$adj.P.Val < 0.05), ]
rownames(MaleTreatGO) = MaleTreatGO$X
MaleTreatGO$gene = rownames(MaleTreatGO)
MaleTreatGOUP = unique(MaleTreatGO[MaleTreatGO$logFC >= 0.5, ])
MaleTreatGOUPres = MaleTreatGOUP$logFC
names(MaleTreatGOUPres) = MaleTreatGOUP$gene
MaleTreatGOUP = sort(MaleTreatGOUPres, decreasing = T)
MaleTreatGOUP = as.matrix(MaleTreatGOUPres)
MaleTreatGOUP = rownames(MaleTreatGOUP)
MaleTreatGOUP = as.vector(MaleTreatGOUP)

MaleTreatGO = MaleTreatNEW[which(MaleTreatNEW$adj.P.Val < 0.05), ]
rownames(MaleTreatGO) = MaleTreatGO$X
MaleTreatGO$gene = rownames(MaleTreatGO)
MaleTreatGODWN = unique(MaleTreatGO[MaleTreatGO$logFC <= -0.5, ])
MaleTreatGODWNres = MaleTreatGODWN$logFC
names(MaleTreatGODWNres) = MaleTreatGODWN$gene
MaleTreatGODWN = sort(MaleTreatGODWNres, decreasing = F)
MaleTreatGODWN = as.matrix(MaleTreatGODWNres)
MaleTreatGODWN = rownames(MaleTreatGODWN)
MaleTreatGODWN = as.vector(MaleTreatGODWN)

# Write GO function
run_gost_and_process_results <- function(gene_list) {
  # Run gost()
  GORes <- gost(gene_list,
                organism = "uparryii",
                user_threshold = 0.05,
                custom_bg = AGSbackground,
                ordered_query = T,
                correction_method = "fdr", 
                evcodes = T)
  
  # Process results and create a dataframe
  result_df <- data.frame(Term.ID = GORes$result$term_id,
                          Term.Name = GORes$result$term_name,
                          geneid = GORes$result$intersection,
                          P.value = GORes$result$p_value,
                          Source = GORes$result$source, 
                          Term.Size = GORes$result$term_size, 
                          Precision = GORes$result$precision,
                          intersection_size = GORes$result$intersection_size,
                          query_size = GORes$result$query_size, 
                          Gene_ratio = as.numeric(GORes$result$intersection_size / GORes$result$query_size))
  
  return(result_df)
}

#run function
GOResults <- run_gost_and_process_results(MaleTreatGOres)
#write.csv(GOResults, 'dream/GOResults/MaleTreatGOResults.csv')

# Female 2W to LH
FemaleTreatNEW = read.csv("Female_LHvsEH_sig.csv")
FemaleTreatGO = FemaleTreatNEW[which(FemaleTreatNEW$adj.P.Val < 0.05), ]
rownames(FemaleTreatGO) = FemaleTreatGO$X
FemaleTreatGO$gene = rownames(FemaleTreatGO)
FemaleTreatGOres = FemaleTreatGO$logFC
names(FemaleTreatGOres) = FemaleTreatGO$gene
FemaleTreatGOres = sort(FemaleTreatGOres, decreasing = T)
FemaleTreatGOres = as.matrix(FemaleTreatGOres)
FemaleTreatGOres = rownames(FemaleTreatGOres)
FemaleTreatGOres = as.vector(FemaleTreatGOres)

FemaleTreatGO = FemaleTreatNEW[which(FemaleTreatNEW$adj.P.Val < 0.05), ]
rownames(FemaleTreatGO) = FemaleTreatGO$X
FemaleTreatGO$gene = rownames(FemaleTreatGO)
FemaleTreatGOUP = unique(FemaleTreatGO[FemaleTreatGO$logFC >= 0.5, ])
FemaleTreatGOUPres = FemaleTreatGOUP$logFC
names(FemaleTreatGOUPres) = FemaleTreatGOUP$gene
FemaleTreatGOUP = sort(FemaleTreatGOUPres, decreasing = T)
FemaleTreatGOUP = as.matrix(FemaleTreatGOUPres)
FemaleTreatGOUP = rownames(FemaleTreatGOUP)
FemaleTreatGOUP = as.vector(FemaleTreatGOUP)

FemaleTreatGO = FemaleTreatNEW[which(FemaleTreatNEW$adj.P.Val < 0.05), ]
rownames(FemaleTreatGO) = FemaleTreatGO$X
FemaleTreatGO$gene = rownames(FemaleTreatGO)
FemaleTreatGODWN = unique(FemaleTreatGO[FemaleTreatGO$logFC <= -0.5, ])
FemaleTreatGODWNres = FemaleTreatGODWN$logFC
names(FemaleTreatGODWNres) = FemaleTreatGODWN$gene
FemaleTreatGODWN = sort(FemaleTreatGODWNres, decreasing = F)
FemaleTreatGODWN = as.matrix(FemaleTreatGODWNres)
FemaleTreatGODWN = rownames(FemaleTreatGODWN)
FemaleTreatGODWN = as.vector(FemaleTreatGODWN)

# Write GO function
run_gost_and_process_results <- function(gene_list) {
  # Run gost()
  GORes <- gost(gene_list,
                organism = "uparryii",
                user_threshold = 0.05,
                custom_bg = AGSbackground,
                ordered_query = TRUE,
                correction_method = "fdr", 
                evcodes = TRUE)
  
  # Process results and create a dataframe
  result_df <- data.frame(Term.ID = GORes$result$term_id,
                          Term.Name = GORes$result$term_name,
                          geneid = GORes$result$intersection,
                          P.value = GORes$result$p_value,
                          Source = GORes$result$source, 
                          Term.Size = GORes$result$term_size, 
                          Precision = GORes$result$precision,
                          intersection_size = GORes$result$intersection_size,
                          query_size = GORes$result$query_size, 
                          Gene_ratio = as.numeric(GORes$result$intersection_size / GORes$result$query_size))
  
  return(result_df)
}

#run function
GOResults <- run_gost_and_process_results(FemaleTreatGOUP)
#write.csv(GOResults, 'dream/GOResults/FemaleTreatGOResults.csv')

# Male versus female at Late hibernation
lateHibSexNEW = read.csv("lateHib_MvsF_sig.csv")
SexLHGO = lateHibSexNEW[which(lateHibSexNEW$adj.P.Val < 0.05), ]
rownames(SexLHGO) = SexLHGO$X
SexLHGO$gene = rownames(SexLHGO) 
SexLHGOres = SexLHGO$logFC
names(SexLHGOres) = SexLHGO$gene
SexLHGOres = sort(SexLHGOres, decreasing = T)
SexLHGOres = as.matrix(SexLHGOres)
SexLHGOres = rownames(SexLHGOres)
SexLHGOres = as.vector(SexLHGOres)

SexLHGO = lateHibSexNEW[which(lateHibSexNEW$adj.P.Val < 0.05), ]
rownames(SexLHGO) = SexLHGO$X
SexLHGO$gene = rownames(SexLHGO) 
SexLHGOUP = unique(SexLHGO[SexLHGO$logFC > 0.5, ])
SexLHGOUPres = SexLHGOUP$logFC
names(SexLHGOUPres) = SexLHGOUP$gene
SexLHGOUP = sort(SexLHGOUPres, decreasing = T)
SexLHGOUP = as.matrix(SexLHGOUPres)
SexLHGOUP = rownames(SexLHGOUP)
SexLHGOUP = as.vector(SexLHGOUP)

SexLHGO = lateHibSexNEW[which(lateHibSexNEW$adj.P.Val < 0.05), ]
rownames(SexLHGO) = SexLHGO$X
SexLHGO$gene = rownames(SexLHGO) 
SexLHGODWN = unique(SexLHGO[SexLHGO$logFC <= -0.5, ])
SexLHGODWNres = SexLHGODWN$logFC
names(SexLHGODWNres) = SexLHGODWN$gene
SexLHGODWN = sort(SexLHGODWNres, decreasing = F)
SexLHGODWN = as.matrix(SexLHGODWNres)
SexLHGODWN = rownames(SexLHGODWN)
SexLHGODWN = as.vector(SexLHGODWN)

# Write GO function
run_gost_and_process_results <- function(gene_list) {
  # Run gost()
  GORes <- gost(gene_list,
                organism = "uparryii",
                user_threshold = 0.05,
                custom_bg = AGSbackground,
                ordered_query = T,
                correction_method = "fdr", 
                evcodes = TRUE)
  
  # Process results and create a dataframe
  result_df <- data.frame(Term.ID = GORes$result$term_id,
                          Term.Name = GORes$result$term_name,
                          geneid = GORes$result$intersection,
                          P.value = GORes$result$p_value,
                          Source = GORes$result$source, 
                          Term.Size = GORes$result$term_size, 
                          Precision = GORes$result$precision,
                          intersection_size = GORes$result$intersection_size,
                          query_size = GORes$result$query_size, 
                          Gene_ratio = as.numeric(GORes$result$intersection_size / GORes$result$query_size))
  
  return(result_df)
}

#run function
GOResults <- run_gost_and_process_results(SexLHGODWN)
#write.csv(GOResults, 'dream/GOResults/SexLHGODWNResults.csv')

# Male versus Female EH
earlyHibSexNEW = read.csv("earlyHib_MvsF_sig.csv")
SexEHGO = earlyHibSexNEW[which(earlyHibSexNEW$adj.P.Val < 0.05), ]
rownames(SexEHGO) = SexEHGO$X
SexEHGO$gene = rownames(SexEHGO) 
SexEHGOres = SexEHGO$logFC
names(SexEHGOres) = SexEHGO$gene
SexEHGOres = sort(SexEHGOres, decreasing = T)
SexEHGOres = as.matrix(SexEHGOres)
SexEHGOres = rownames(SexEHGOres)
SexEHGOres = as.vector(SexEHGOres)

SexEHGO = earlyHibSexNEW[which(earlyHibSexNEW$adj.P.Val < 0.05), ]
rownames(SexEHGO) = SexEHGO$X
SexEHGO$gene = rownames(SexEHGO) 
SexEHGOUP = unique(SexEHGO[SexEHGO$logFC > 0.5, ])
SexEHGOUPres = SexEHGOUP$logFC
names(SexEHGOUPres) = SexEHGOUP$gene
SexEHGOUP = sort(SexEHGOUPres, decreasing = T)
SexEHGOUP = as.matrix(SexEHGOUPres)
SexEHGOUP = rownames(SexEHGOUP)
SexEHGOUP = as.vector(SexEHGOUP)

SexEHGO = earlyHibSexNEW[which(earlyHibSexNEW$adj.P.Val < 0.05), ]
rownames(SexEHGO) = SexEHGO$X
SexEHGO$gene = rownames(SexEHGO) 
SexEHGODWN = unique(SexEHGO[SexEHGO$logFC <= -0.5, ])
SexEHGODWNres = SexEHGODWN$logFC
names(SexEHGODWNres) = SexEHGODWN$gene
SexEHGODWN = sort(SexEHGODWNres, decreasing = F)
SexEHGODWN = as.matrix(SexEHGODWNres)
SexEHGODWN = rownames(SexEHGODWN)
SexEHGODWN = as.vector(SexEHGODWN)

GOResults <- run_gost_and_process_results(SexEHGODWN)
#write.csv(GOResults, 'dream/GOResults/SexEHGODWNResults.csv')

