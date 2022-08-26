#Astrid M Manuel
#07/02/2019

#Install and load Limma R package
#R version 3.5 is needed

setwd("C:/Users/amanuel1/Dev/ZhaoLab/MS/Methylation")

na_val <- c()

for(i in seq(1,nrow(Methy_Beta))){if(length(which(is.na(Methy_Beta[i,]))) > 0){na_val <- c(na_val,i)}}


# Methylation matrix with M-values after removal of probes with NA values
Methy <- read.table("MS_Methy_Mvalues_rmNA.txt")
Methy_Beta <- read.table("MS_Methy_BetaValues_rmNA.txt", as.is=T)

Sample_ID <- colnames(Methy_Beta)

groups <- c(rep("MS", 28), rep("Control", 19))

age <- c(63,74,64,50,42,54,35,49,51,75,50,49,51,60,61,57,59,61,55,45,59,78,46,60,54,55,44,53,54,76,71,72,79,59,73,64,58,64,54,66,68,58,64,81,71,53,71)

sex <- c(rep("F", 17), rep("M", 11), rep("F", 7), rep("M", 12))

sample_info <- data.frame(cbind(Sample_ID,groups,age,sex))
sample_info$age <- as.numeric(age)

f <- factor(sample_info$groups)

design <- model.matrix(~0 + f + age + sex, sample_info)

fit <- lmFit(Methy, design)

cont.matrix_new <- makeContrasts(P1="fMS-fControl", levels=design)

fit <- contrasts.fit(fit, cont.matrix_new)
fit <- eBayes(fit, trend=T)

res_new <- topTable(fit, adjust="BH", number=nrow(Methy))


#-------------------------------
#MS hippocampus methylation

MS_HC <- read.delim("MS_hippocampus_methy.txt")

na_val <- c()

for(i in seq(1,nrow(MS_HC))){if(length(which(is.na(MS_HC[i,]))) > 0){na_val <- c(na_val,i)}}


rownames(MS_HC) <- MS_HC$ID_REF
MS_HC <- MS_HC[,c(seq(2,16))]

Sample_ID <- colnames(MS_HC)
groups <- c(rep("Myelinated", 8), rep("Demyelinated", 7))
sample_info <- data.frame(cbind(Sample_ID,groups))

f <- factor(sample_info$groups)

design <- model.matrix(~0 + f , sample_info)

fit <- lmFit(Methy, design)

cont.matrix_new <- makeContrasts(P1="fMyelinated-fDemyelinated", levels=design)

fit <- contrasts.fit(fit, cont.matrix_new)
fit <- eBayes(fit)

res_new <- topTable(fit, adjust="BH", number=nrow(MS_HC))

MS_HC_DMCpGs <- read.table("MS_hippocampus_DM_CpGs.txt", as.is = T)