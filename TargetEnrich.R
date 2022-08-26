# Astrid M Manuel
# Last Update: 02/22/2022
# Handling the raw data from Therapeutic Target Database (TTD): http://db.idrblab.net/ttd/full-data-download
# Preparing TTD Drug Targets for drug signature enrichment analyses

# Reading the TTD drug targets
rawData <- read.delim("C:/Users/astri/Documents/ZhaoLab/TTD_DrugTargets/TTD_DrugTargets.txt", header = F)

# Preparing the drug targets
targetID_idx <- which(rawData$V2 == "TARGETID")
TargetIds <- rawData$V3[targetID_idx]
geneNames_idx <- which(rawData$V2 == "GENENAME")
geneNames <- rawData$V3[geneNames_idx]

MS_GRN_genes <- read.table("C:/Users/astri/Documents/ZhaoLab/MS_PBMCs/MS_methy_dmGWAS/GRN_NW_top50.txt", header = T)
MS_GRN_genes <- MS_GRN_genes$GeneName

targets_enriched <-  MS_GRN_genes[which(MS_GRN_genes %in% geneNames)]

targetIDs <- c()
genenames <- c()
targnames <- c()
targtypes <- c()
bioclasses <- c()
druginfos <- c()
clinicalstatuses <- c()
idx_approved <- c()
idx_phase3 <- c()

for(i in seq(1, length(targets_enriched))){
  target_idx <- which(rawData$V3 == targets_enriched[i])
  targetID <- rawData$V1[target_idx]
  for(x in targetID){
    targetIDs <- c(x, targetIDs)
    genenames <- c(targets_enriched[i], genenames)
    targname <- rawData$V3[which(rawData$V1==x & rawData$V2=="TARGNAME")]
    targnames <- c(targname, targnames)
    targtype <- rawData$V3[which(rawData$V1==x & rawData$V2=="TARGTYPE")]
    targtypes <- c(targtype, targtypes)
    bioclass <- rawData$V3[which(rawData$V1==x & rawData$V2=="BIOCLASS")]
    bioclasses <- c(bioclass, bioclasses)
    druginfo <- list(rawData$V4[which(rawData$V1 == x & rawData$V2 == "DRUGINFO")])
    druginfos <- c(druginfo, druginfos)
    clinicalstatus <- list(rawData$V5[which(rawData$V1 == x & rawData$V2 == "DRUGINFO")])
    clinicalstatuses <- c(clinicalstatus, clinicalstatuses)
    uniqueCS <- (unique(clinicalstatus[[1]]))
    approved_drugs <- c()
    phase3_drugs <- c()
    if("Approved" %in% uniqueCS){
      idx_A <- which(clinicalstatus[[1]] == "Approved")
      print(targets_enriched[i])
      print(clinicalstatus[[1]][idx_A])
      print(druginfo[[1]][idx_A])
      idx_approved <- append(idx_approved, i)
    } 
    approved_drugs <- toString(approved_drugs)
    if("Phase 3" %in% uniqueCS){
      idx_P3 <- which(clinicalstatus[[1]] == "Phase 3")
      print(targets_enriched[i])
      print(clinicalstatus[[1]][idx_P3])
      print(druginfo[[1]][idx_P3])
      idx_phase3 <- append(idx_phase3, i) 
    }
    
    
  }
}

enrichmentDetails <- cbind(targetIDs, genenames, targnames, targtypes, bioclasses)
colnames(enrichmentDetails) <- c("TTD ID", "Gene Name", "Target Name", "Target Type", "Bio Class")
write.table(enrichmentDetails, "enrichmentDetails.txt", row.names = F, quote = F, sep = "\t")






