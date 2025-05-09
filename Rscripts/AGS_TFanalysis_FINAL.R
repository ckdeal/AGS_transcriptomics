#load packages
{library(dplyr)
  library(ggplot2)
  library(biomaRt)
  library(stringi)}


# read in genes with orthologs
MaleTreatsig_orth = read.csv("../DE_genes/Male_LHvsEH_sig.csv")
FemaleTreatsig_orth = read.csv("../DE_genes/Female_LHvsEH_sig.csv")
twoWeekSexsig_orth = read.csv("../DE_genes/earlyHib_MvsF_sig.csv")
lateHibSexsig_orth = read.csv("../DE_genes/lateHib_MvsF_sig.csv")
  
  #filter out sig genes
MaleTreatsig = MaleTreatsig_orth[which(MaleTreatsig_orth$adj.P.Val < 0.05), ] 
sum(MaleTreatsig$adj.P.Val < 0.05 & MaleTreatsig$logFC >= 0.05) # 40 upregulated genes
sum(MaleTreatsig$adj.P.Val < 0.05 & MaleTreatsig$logFC <= 0.05) # 18 downregulated genes
rownames(MaleTreatsig) = MaleTreatsig$X

FemaleTreatsig = FemaleTreatsig_orth[which(FemaleTreatsig_orth$adj.P.Val < 0.05), ]
sum(FemaleTreatsig$adj.P.Val < 0.05 & FemaleTreatsig$logFC >= 0.05) # 24 upregulated genes
sum(FemaleTreatsig$adj.P.Val < 0.05 & FemaleTreatsig$logFC <= 0.05) # 28 downregulated genes
rownames(FemaleTreatsig) = FemaleTreatsig$X

earlyHibSexSigGenes = twoWeekSexsig_orth[which(twoWeekSexsig_orth$adj.P.Val < 0.05), ]
sum(earlyHibSexSigGenes$adj.P.Val < 0.05 & earlyHibSexSigGenes$logFC >= 0.05) # 12 upregulated genes in males
sum(earlyHibSexSigGenes$adj.P.Val < 0.05 & earlyHibSexSigGenes$logFC <= 0.05) # 10 downregulated genes in females
rownames(earlyHibSexSigGenes) = earlyHibSexSigGenes$X

lateHibSexSigGenes = lateHibSexsig_orth[which(lateHibSexsig_orth$adj.P.Val < 0.05), ]
sum(lateHibSexSigGenes$adj.P.Val < 0.05 & lateHibSexSigGenes$logFC >= 0.05) # 59 upregulated genes in males
sum(lateHibSexSigGenes$adj.P.Val < 0.05 & lateHibSexSigGenes$logFC <= 0.05) # 69 downregulated genes in females
rownames(lateHibSexSigGenes) = lateHibSexSigGenes$X

# Convert gene names to ENSEMBL IDs #### ---------------------------------------
MaleTreatSigUPGenes = MaleTreatsig[which(MaleTreatsig$logFC >= 0.5), ] %>% 
  rownames(MaleTreatsig)
#write.csv(MaleTreatSigUPGenes, "DEGenes_wOrthologs/MaleTreatSigUPGenes.csv")

MaleTreatSigDWNGenes = MaleTreatsig[which(MaleTreatsig$logFC <= -0.5), ] %>% 
  rownames(MaleTreatsig)
#write.csv(MaleTreatSigDWNGenes, "DEGenes_wOrthologs/MaleTreatSigDWNGeness.csv")

MaleTreatGenes = MaleTreat[which(MaleTreat$logFC <= 0.05 & MaleTreat$logFC >= -0.05), ] %>% 
  rownames(MaleTreat) # filter genes between -0.05 and 0.05 for background

BackgroundMaleGenes = setdiff(MaleTreatGenes, MaleTreatSigGenes) # ensure DE genes are not in background

#Females
FemaleTreatSigUPGenes = FemaleTreatsig[which(FemaleTreatsig$logFC >= 0.5), ] %>% 
  rownames(FemaleTreatsig)
#write.csv(FemaleTreatSigUPGenes, "DEGenes_wOrthologs/FemaleTreatSigUPGenes.csv")

FemaleTreatSigDWNGenes = FemaleTreatsig[which(FemaleTreatsig$logFC <= -0.5), ] %>% 
  rownames(FemaleTreatsig)
#write.csv(FemaleTreatSigDWNGenes, "DEGenes_wOrthologs/FemaleTreatSigDWNGenes.csv")

FemaleTreatGenes = FemaleTreat[which(FemaleTreat$logFC <= 0.05 & FemaleTreat$logFC >= -0.05), ] %>% 
  rownames(FemaleTreat) # filter genes between -0.05 and 0.05 for background set

BackgroundGenesFemales = setdiff(FemaleTreatGenes, FemaleTreatSigGenes) # ensure DE genes are not in background

# sex differences at late hibernation
lateHibSexGenes = lateHibSexSig[which(lateHibSexSig$logFC <= 0.05 & lateHibSexSig$logFC >= -0.05), ]
rownames(lateHibSexGenes) = lateHibSexGenes$X

lateHibSexBackgroundGenes = setdiff(lateHibSexGenes, lateHibSexSigGenes)
lateHibSexBackgroundGenes = lateHibSexBackgroundGenes$X

lateHibSexUPGenes = lateHibSexSigGenes[which(lateHibSexSigGenes$logFC >= 0.5), ] %>%
  rownames(lateHibSexSigGenes)
#write.csv(lateHibSexUPGenes, "DEGenes_wOrthologs/lateHibSexUPGenes.csv")

lateHibSexDWNGenes = lateHibSexSigGenes[which(lateHibSexSigGenes$logFC <= 0.5), ] %>%
  rownames(lateHibSexSigGenes)
#write.csv(lateHibSexDWNGenes, "DEGenes_wOrthologs/lateHibSexDWNGenes.csv")

# sex differences at early hibernation
earlyHibSexGenes = earlyHibSexSig[which(earlyHibSexSig$logFC <= 0.05 & earlyHibSexSig$logFC >= -0.05), ] #background set
rownames(earlyHibSexSig) = earlyHibSexSig$X

earlyHibSexBackgroundGenes = setdiff(earlyHibSexGenes, earlyHibSexSigGenes)
earlyHibSexBackgroundGenes = earlyHibSexBackgroundGenes$X

earlyHibSexUPGenes = earlyHibSexSigGenes[which(earlyHibSexSigGenes$logFC >= 0.5), ] %>%
  rownames(earlyHibSexSigGenes)
#write.csv(lateHibSexUPGenes, "DEGenes_wOrthologs/lateHibSexUPGenes.csv")

earlyHibSexDWNGenes = earlyHibSexSigGenes[which(earlyHibSexSigGenes$logFC <= 0.5), ] %>%
  rownames(earlyHibSexSigGenes)
#write.csv(lateHibSexDWNGenes, "DEGenes_wOrthologs/lateHibSexDWNGenes.csv")

# Convert gene names to ENSEMBL IDs for CiiiDER analysis -----------------------
#
write_gene_ids_to_file <- function(values, output_file) {
  # Specify the BioMart dataset
  uro = useEnsembl('ensembl', 'uparryii_gene_ensembl', mirror = "uswest")
  # Get gene information
  results <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                   filters = "external_gene_name",
                   values = values,
                   mart = uro)
  
  # Write gene IDs to a file
  write.table(results$ensembl_gene_id,
              row.names = FALSE, col.names = FALSE,
              sep = "\t", quote = FALSE, file = output_file)
  
}

# Run ID Conversion function for CiiiDER analysis
# Replace 'YourInputValues' with the actual values you want to use,
# and 'YourOutputFile.txt' with the desired output file name.
# Here, convert 

setwd("")
write_gene_ids_to_file(values = FemaleTreatSigUPGenes, 
                       output_file = "FemaleTreatSigUPGenes.txt")

uro = useEnsembl('ensembl', 'uparryii_gene_ensembl', mirror = "uswest")
getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                 filters = "external_gene_name",
                 values = FemaleTreatSigUPGenes,
                 mart = uro)
# You can now proceed to CiiiDER enrichment analysis (https://ciiider.com/)

# Enrichment Analysis ####
# After running CiiiDER using desktop java app. Upload enrichment results as a 
# csv and filter out significant and enriched TFs. 
# ***NOTE*** You can only access the predicted targets for a priori TFs by
# opening the CiiiDER software application

### Ciiider enriched TFS Males 2W to LH ####
# Males 2W to LH Upregulated
MaleTFUP = read.csv("../CiiiDER_files/Male_LH_vs_EH/Upregulated/Enrichment_ BackgroundMaleGenes_MostSigDeficit.csv")

MaleTFUP$geneid <- stri_trans_totitle(MaleTFUP$Transcription.Factor.Name)
intersect(MaleTreat$geneid, MaleTFUP$geneid) # which overlap - does not account for TF dimers
MaleTFUP$GeneDIR = "UpRegulated" # Add column specifying these genes are UP in LH
MaleTFUP1 = MaleTFUP %>% filter(MaleTFUP$Gene.P.Value < 0.05)
MaleTFUP1 = MaleTFUP1 %>% filter(MaleTFUP1$Log2.Enrichment >= 1 | MaleTFUP1$Log2.Enrichment <= -1) #118 TFs
MaleTFUP1$GeneDIR = "UpRegulated" 

#down regulated DE genes
MaleTFDWN = read.csv("../CiiiDER_files/Male_LH_vs_EH/Downregulated/Enrichment_ BackgroundMaleGenes_MostSigDeficit.csv")

MaleTFDWN$geneid <- stri_trans_totitle(MaleTFDWN$Transcription.Factor.Name)
intersect(MaleTreat$geneid, MaleTFDWN$geneid) #445 Tfs overlap
MaleTFDWN$GeneDIR = "DownRegulated"
#MaleTFDWN1 = MaleTFDWN %>% filter(MaleTFDWN$fdrPval < 0.25)
MaleTFDWN1 = MaleTFDWN %>% filter(MaleTFDWN$Gene.P.Value < 0.05)
MaleTFDWN1 = MaleTFDWN1 %>% filter(MaleTFDWN1$Log2.Enrichment > 1 | MaleTFDWN1$Log2.Enrichment < -1) # 60 TFs
MaleTFDWN1$GeneDIR = "DownRegulated"

## Ciiider enrich Females 2W to LH ####
#Upregulated
FemaleTFUP = read.csv("../CiiiDER_files/Female_LH_vs_EH/Upregulated/Up_Filtered_Enrichment_BackgroundFemalesGenes_MostSigDeficit.csv")

FemaleTFUP$geneid <- stri_trans_totitle(FemaleTFUP$Transcription.Factor.Name)
intersect(FemaleTreat$geneid, FemaleTFUP$geneid) #445
FemaleTFUP$GeneDIR = "UpRegulated"
FemaleTFUP1 = FemaleTFUP %>% filter(FemaleTFUP$Gene.P.Value < 0.05)
FemaleTFUP1 = FemaleTFUP1 %>% filter(FemaleTFUP1$Log2.Enrichment > 1 | FemaleTFUP1$Log2.Enrichment < -1) #76
FemaleTFUP1$GeneDIR = "UpRegulated"

# downregulated
FemaleTFDWN = read.csv("../CiiiDER_files/Female_LH_vs_EH/Downregulated/DownFiltered_Enrichment_BackgroundFemalesGenes_MostSigDeficit.csv")

FemaleTFDWN$geneid <- stri_trans_totitle(FemaleTFDWN$Transcription.Factor.Name)
intersect(FemaleTreat$gene_id, FemaleTFDWN$geneid)
FemaleTFDWN$GeneDIR = "DownRegulated"
FemaleTFDWN1 = FemaleTFDWN %>% filter(FemaleTFDWN$Gene.P.Value < 0.05)
FemaleTFDWN1 = FemaleTFDWN1 %>% filter(FemaleTFDWN1$Log2.Enrichment > 1 | FemaleTFDWN1$Log2.Enrichment < -1) #97
FemaleTFDWN1$GeneDIR = "DownRegulated"

## Ciiider enrich sexes at late hibernation ####
# upregulated DE genes in males at LH
lateHibTFUP = read.csv("../CiiiDER_files/LateHibernation/Upregulated/Enrichment_lateHibSexBackgroundGenes_MostSigDeficit.csv")

lateHibTFUP$geneid <- stri_trans_totitle(lateHibTFUP$Transcription.Factor.Name)
intersect(lateHibSex$gene_id, lateHibTFUP$geneid)
lateHibTFUP$GeneDIR = "UpRegulated"
lateHibTFUP1 = lateHibTFUP %>% filter(lateHibTFUP$Gene.P.Value < 0.05)
lateHibTFUP1 = lateHibTFUP1 %>% filter(lateHibTFUP1$Log2.Enrichment >= 1 | lateHibTFUP1$Log2.Enrichment <= -1) #63
lateHibTFUP1$GeneDIR = "UpRegulated"

#down regulated DE genes (i.e., upregulated in females)
lateHibTFDWN = read.csv("../CiiiDER_files/LateHibernation/Downregulated/Enrichment_lateHibSexBackgroundGenes_MostSigDeficit.csv")

lateHibTFDWN$geneid <- stri_trans_totitle(lateHibTFDWN$Transcription.Factor.Name)
intersect(lateHibSex$gene_id, lateHibTFDWN$geneid)
lateHibTFDWN$GeneDIR = "DownRegulated"
lateHibTFDWN1 = lateHibTFDWN %>% filter(lateHibTFDWN$Gene.P.Value < 0.05)
lateHibTFDWN1 = lateHibTFDWN1 %>% filter(lateHibTFDWN1$Log2.Enrichment > 1 | lateHibTFDWN1$Log2.Enrichment < -1) #90 TFs
lateHibTFDWN1$GeneDIR = "DownRegulated"

## Overlap of TFs ####

# male versus female across hibernation
MalevsFemaleUPEHtoLH = intersect(MaleTFUP1$geneid, FemaleTFUP1$geneid)

MalevsFemaleDWNEHtoLH = intersect(MaleTFDWN1$geneid, FemaleTFDWN1$geneid)

LHSexes = intersect(lateHibTFUP1$geneid, lateHibTFDWN1$geneid)

MaleUPvsDWNEHtoLH = intersect(MaleTFUP1$geneid, MaleTFDWN1$geneid)

FemaleUPvsDWNEHtoLH = intersect(FemaleTFUP1$geneid, FemaleTFDWN1$geneid)

# UP In all comparisons only
UPEHLHbothSexes = intersect(MaleTFUP1$geneid, FemaleTFUP1$geneid)
LHSexcombined = c(lateHibTFUP1$geneid, lateHibTFDWN1$geneid)
duplicates = duplicated(LHSexcombined)
LHSexcombined = LHSexcombined[!duplicates]
intersect(UPEHLHbothSexes, LHSexcombined)

# DOWN in all comparisons
DWNEHLHbothSexes = intersect(MaleTFDWN1$geneid, FemaleTFDWN1$geneid)
FEMALELHvsEHLHbothSexes = intersect(DWNEHLHbothSexes, LHSexcombined)

# see where there is overlap across groups where enrichment is > 1 in order
# to look at DE targets. 
#UP across hibernation
MaleTFUP2 = MaleTFUP1 %>% 
  filter(MaleTFUP1$Gene.P.Value <= 0.05 & MaleTFUP1$Log2.Enrichment > 1)

intersect(MaleTFUP1$geneid, lateHibTFUP1$geneid) # 14 TFs
intersect(MaleTFUP2$geneid, lateHibTFUP2$geneid) # 1 TFs Srebf1

FemaleTFUP2 = FemaleTFUP1 %>% 
  filter(FemaleTFUP1$Gene.P.Value <= 0.05 & FemaleTFUP1$Log2.Enrichment > 1)

intersect(FemaleTFUP1$geneid, lateHibTFDWN1$geneid)
intersect(FemaleTFUP2$geneid, lateHibTFDWN2$geneid) # 3 TFs: Gcm1, Maf, Mef2b

intersect(MaleTFUP2$geneid, FemaleTFUP2$geneid) # 2 TFs Dmbx1, Gcm2

lateHibTFUP2 = lateHibTFUP1 %>% 
  filter(lateHibTFUP1$Gene.P.Value <= 0.05 & lateHibTFUP1$Log2.Enrichment > 1)

lateHibTFDWN2 = lateHibTFDWN1 %>% 
  filter(lateHibTFDWN1$Gene.P.Value <= 0.05 & lateHibTFDWN1$Log2.Enrichment > 1)

LHSexcombined = c(lateHibTFUP2$geneid, lateHibTFDWN2$geneid)
duplicates = duplicated(LHSexcombined)
LHSexcombined = LHSexcombined[!duplicates] #97

# Down
MaleTFDWN2 = MaleTFDWN1 %>% 
  filter(MaleTFDWN1$Gene.P.Value < 0.05 & MaleTFDWN1$Log2.Enrichment > 1)

FemaleTFDWN2 = FemaleTFDWN1 %>% 
  filter(FemaleTFDWN1$Gene.P.Value < 0.05 & FemaleTFDWN1$Log2.Enrichment > 1)

intersect(MaleTFDWN2$geneid, FemaleTFDWN2$geneid) # 7 TFs

intersect(maleFemEHLHDNW, LHSexcombined) # 1 Hoxc11
