# Analytical-Approach-for-MS-GRN
Reference: Astrid M Manuel, Yulin Dai, Peilin Jia, Leorah A Freeman, Zhongming Zhao, A gene regulatory network approach harmonizes genetic and epigenetic signals and reveals repurposable drug candidates for multiple sclerosis, Human Molecular Genetics, 2022;, ddac265, https://doi.org/10.1093/hmg/ddac265


The analytical approach for identification and interpretation of a disease-specific gene regulatory network (GRN):
![Figure1_AnalyticalApproach](https://user-images.githubusercontent.com/10716674/186982697-94b12168-0421-4bc1-89fe-63dd012a207e.png)

Custom scripts of the  pipeline are contained within this repository, including:
* Aggregate_TSS_CpGs.R
* Aggregate_geneticNepigeneticWeights.R
* TargetEnrich.R

Also contained in this repository:
* Script for running Limma R package to optain differentially methylated CpG probes
* Output from Limma
* The dmGWAS java script version we used in this pipeline
* Result files for all genetic and epigenetic weights used in running dmGWAS


For downloading data from the Target Therapeutic Database please follow this link: http://db.idrblab.net/ttd/full-data-download

For running WebCSEA please follow this link:https://bioinfo.uth.edu/webcsea/

For questions you may e-mail me: Astrid.M.Manuel@uth.tmc.edu
