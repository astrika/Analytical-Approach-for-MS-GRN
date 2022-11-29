# Astrid M Manuel
# Created: 11/29/2022
# Aggregating genetic and epigenetic weights from GWAS summary statistics and disease-specific DNA methylation datasets
# Genetic weights: gene-level GWAS-based weights obtained from MAGMA
# Epigenetic weights: gene-level methylation-based weights obtained by applying the 'Aggregate_TSS_CpGs.R' scripts
# Here using example datasets from Alzheimer's Disease

# Read in GWAS-based weights
genetic_Weights = read.table('C:/Users/ALFON/Documents/Astrid/AD_GWAS_KeithBoyan/NodeWeights_ADGWAS.txt')
genetic_Weights$V2 = abs(genetic_Weights$V2)

# Read in the methylation-based weights
epigenetic_weights = read.table('C:/Users/ALFON/Documents/Astrid/AD_GWAS_KeithBoyan/2022_dmGWAS/MethyNodeWeights_AD_PostmortemBrain.txt')
epigenetic_weights_directions = read.table('C:/Users/ALFON/Documents/Astrid/AD_GWAS_KeithBoyan/2022_dmGWAS/MethyDirection_AD_PostmortemBrain.txt')

# Matching the genes
match_idx = match(epigenetic_weights$V1, genetic_Weights$V1)
idx_na = which(is.na(match_idx))
match_idx = match_idx[-idx_na]

epigenetic_weights = epigenetic_weights[-idx_na, ]
epigenetic_weights_directions = epigenetic_weights_directions[-idx_na, ]
genetic_Weights = genetic_Weights[match_idx, ]

# Calculate the lambda = variance ratio, scaling factor to minimize variance during integration
var_genetic = var(genetic_Weights$V2)
var_epigenetic = var(epigenetic_weights$V2)
lambda_varianceRatio = var_genetic/var_epigenetic

# Caluculating the new genetic and epigenetic weights for dmGWAS run
methy_varianceRatio = lambda_varianceRatio * epigenetic_weights$V2
nodeWeights = genetic_Weights$V2 + lambda_varianceRatio * epigenetic_weights$V2

# create tables for summary file and node weight file:
summary = cbind(genetic_Weights$V1, genetic_Weights$V2, epigenetic_weights_directions$V2, epigenetic_weights$V2,  methy_varianceRatio, nodeWeights)
colnames(summary) = c('Gene Symbol', 'GWAS weight', 'Methy Direction Weight', 'Abs Methy Weight', 'Abs Methy Weight Scaled', 'Node Weight')
setwd('C:/Users/ALFON/Documents/Astrid/AD_GWAS_KeithBoyan/2022_dmGWAS')
write.table(summary, 'summary_AD_GWAS_Methy_ROSMAP.txt', quote = F, row.names = F, sep = '\t')

node_weight_file = cbind(genetic_Weights$V1, nodeWeights)
write.table(node_weight_file, 'NodeWeights_AD_geneticNepigenetic.txt', quote = F, row.names = F, col.names = F, sep = '\t')
