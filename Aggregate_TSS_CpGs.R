# Astrid M Manuel
# Last Update: 09/07/2021
# Version 2 of OnlyPromoters Script
# Obtaining methylation gene-level p-values with consideration of only promoter regions (as annotated by Illumina)
# MS PBMCs methylation

# Set your working directory in the location of the files: annotation file of Human Methylation 450 BeadChip & the Methylaiton data file
setwd("C:/Users/astri/Documents/ZhaoLab/MS_PBMCs/Methylation")

# GPL13534 annotations file for HumanMethylation450 BeadChip 
GPL <- read.delim("GPL13534-11288.txt", as.is = T)

# To store the indeces of GPL probe records containing an annotation for transcription start sites (TSS) of genes
idx.TSS <- c()

# Spliting the "UCSC_RefGene_Group" column of GPL which contains several different annotations for a particular probe
region_annotated <- strsplit(as.vector(GPL$UCSC_RefGene_Group), split = ";")

# Getting the idx.TSS
for(i in seq(1:length(region_annotated))){
  
  if("TSS1500" %in% region_annotated[[i]] | "TSS200" %in% region_annotated[[i]]){
    idx.TSS <- c(idx.TSS, i)
  }
  
}

# Getting all TSS probes
GPL_TSS <- GPL[idx.TSS,]

#-------------------------------------------------------------------------------
# Read in the differentially methylated CpGs, as calculated by Limma R package
Methy <- read.table("MS_PBMCs_DifferentialMethylation.txt", as.is = T)
# Extracting TSS probes from Methy data
# MethyTSS will store the CpG gene-level annotations (as annotated by Illumina manifest file: GPL13534)
MethyTSS <- Methy[(which(rownames(Methy) %in% GPL_TSS$ID)),]
idx.MethyTSS <- match(rownames(MethyTSS), GPL_TSS$ID)
MethyTSS$Gene_Name <- GPL_TSS$UCSC_RefGene_Name[idx.MethyTSS]
MethyTSS$Group <- GPL_TSS$UCSC_RefGene_Group[idx.MethyTSS]

# Spliting the groups from TSS regions
TSS_regions <- strsplit(as.vector(MethyTSS$Group), split = ";")
TSS_genes <- strsplit(as.vector(MethyTSS$Gene), split = ";")

TSS1500 <- c()
TSS1500_genes <- c()
TSS200 <- c()
TSS200_genes <- c()

for(i in seq(1:length(TSS_regions))){
  for(j in seq(1:length(TSS_regions[[i]]))){
    if("TSS1500" == TSS_regions[[i]][j]){
      TSS1500 <- c(TSS1500, i, j)
      TSS1500_genes <- c(TSS1500_genes, TSS_genes[[i]][j])
    }else if("TSS200" == TSS_regions[[i]][j]){
      TSS200 <- c(TSS200, i, j)
      TSS200_genes <- c(TSS200_genes, TSS_genes[[i]][j])
    }
  }
}


# Getting indeces for probes from TSS1500
# Odd numbers store the list index for obtaining the corresponding probe in TSS
TSS1500_probe.idx <- c()

for(i in seq(0,length(TSS1500_genes)-1)){
  
  TSS1500_probe.idx <- c(TSS1500_probe.idx, TSS1500[2*i + 1])
  
}

# Storing the Methy values from TSS1500
TSS1500_probes <- rownames(MethyTSS)[TSS1500_probe.idx]
TSS1500_raw.pvalue <- MethyTSS$P.Value[TSS1500_probe.idx]
TSS1500_logFC <- MethyTSS$logFC[TSS1500_probe.idx]
TSS1500_group <- rep("TSS1500", length(TSS1500_genes))

TSS1500_MethyValues <- data.frame(cbind(TSS1500_genes, TSS1500_probes, TSS1500_raw.pvalue, TSS1500_logFC,  TSS1500_group))

# Writing file with the TSS1500 probes.
# ****Name it according to the type of methylation data that is being  studied:
write.table(TSS1500_MethyValues,"TSS1500_MS_PBMCs.txt", quote = F, sep = "\t", row.names = F, col.names = F)


#Getting indeces for probes from TSS200
TSS200_probe.idx <- c()

for(i in seq(0,length(TSS200_genes)-1)){
  
  TSS200_probe.idx <- c(TSS200_probe.idx, TSS1500[2*i + 1])
  
}

TSS200_probes <- rownames(MethyTSS)[TSS200_probe.idx]
TSS200_raw.pvalue <- MethyTSS$P.Value[TSS200_probe.idx]
TSS200_logFC <- MethyTSS$logFC[TSS200_probe.idx]
TSS200_group <- rep("TSS200", length(TSS200_genes))
TSS200_MethyValues <- data.frame(cbind(TSS200_genes, TSS200_probes, TSS200_raw.pvalue, TSS200_logFC, TSS200_group))

# Writing file with the TSS200 probes.
# ****Name it according to the type of methylation data that is being  studied:
write.table(TSS200_MethyValues,"TSS200_MS_PBMCs.txt", quote = F, sep = "\t", row.names = F, col.names = F)

#-------------------------------------------------------------------------------
# Getting all the methylation probes annotated and calculating gene level methylation scores from unique annotations

colnames(TSS1500_MethyValues) <- c("Genes", "Probes", "Raw_P_value", "logFC", "TSS_Group")
colnames(TSS200_MethyValues) <- c("Genes", "Probes", "Raw_P_value", "logFC", "TSS_Group")
TSS_MethyValues <- rbind(TSS1500_MethyValues, TSS200_MethyValues)
TSS_MethyValues <- TSS_MethyValues[which(!duplicated(TSS_MethyValues)), ]
TSS_MethyValues <- TSS_MethyValues[-which(is.na(TSS_MethyValues$Genes)), ]


# Getting the gene-level methylation scores
# All unique genes from TSS annotations:
genes <- unique(TSS_MethyValues$Genes)
Z.scores.genelevel <- c()

for(i in seq(1, length(genes))){
 
  records <- TSS_MethyValues[which(TSS_MethyValues$Genes == genes[i]),]
  
  # Obtaining z-scores for each probe records within the corresponding gene
  p.values <- as.numeric(records$Raw_P_value)
  Z.scores <- qnorm(1-(p.values/2))
  
  # Getting direction of methylation from logFC
  FC <- as.numeric(records$logFC)
  for(j in seq(1, length(FC))){
    if(FC[j] < 0){
      Z.scores[j] <- Z.scores[j]*(-1)
    }
  }
  
  # Combining Z-scores
  genelevel.z <- sum(Z.scores)/ sqrt(length(Z.scores))
  Z.scores.genelevel <- c(Z.scores.genelevel, genelevel.z)
  
}

forNodes <- abs(Z.scores.genelevel)
MethyNodeWeights <- cbind(genes, forNodes)
MethyDirection <- cbind(genes, Z.scores.genelevel)

# ****Name it according to the type of methylation data that is being  studied:
write.table(MethyNodeWeights, "MethyNodeWeights_MS_PBMCs.txt", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(MethyDirection, "MethyDirection_MS_PBMCs.txt", row.names = F, col.names = F, quote = F, sep = "\t")

par(mfrow=c(1,2))
hist(Methy$P.Value, main = "Histogram of CpG Probe-Level Raw P Values", xlab = "CpG Probe-Level P-values")
hist(as.numeric(MethyNodeWeights[,2]), main= "Histogram for Methylation Gene-Level Z Scores",  xlab = "Gene-Level Z scores")