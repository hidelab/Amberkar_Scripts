########################################################################################################################
# The following script performs a WGCNA analysis on the Mt. Sinai RNAseq dataset (https://www.synapse.org/%23!Synapse:syn3219494).
# Briefly, the 448 samples across 3 brain regions were segregated into "low" and "high" pathology samples. WGCNA was
# run on the low and high samples for each brain regions. The modules generated from both pathologies will then be
# compared to determine coexpression dynamics. 
# This script uses the WGCNA functions that are wrappers to common tasks in the WGCNA protocols and are sourced from (https://github.com/SamBuckberry/RUN-WGCNA)
########################################################################################################################
library(limma)
library(Biobase)
library(metaArray)
library(doMC)
library(pheatmap)
library(clusterProfiler)
library(WGCNA)

library(org.Hs.eg.db)
library(edgeR)
source("/shared/hidelab2/user/md4zsa/Work/Scripts/RUN-WGCNA/WGCNA_functions.R")
#Set active directory and read data
setwd("/shared/hidelab2/user/md4zsa/Work/Data/MSMM_RNAseq/MSMM_RNAseq_WGCNA")
msmm_rnaseq_data=read.delim2("msmm_rnaseq_B1_raw_counts.txt",header = T,row.names = 1,as.is = T)
ensembl_entrezMap=select(x = org.Hs.eg.db,keys=rownames(msmm_rnaseq_data),keytype = "ENSEMBL",columns = "ENTREZID")
ensembl_entrezMap=ensembl_entrezMap[-which(is.na(ensembl_entrezMap$ENTREZID)==T),]
msmm_rnaseq_data=msmm_rnaseq_data[which(rownames(msmm_rnaseq_data)%in%ensembl_entrezMap$ENSEMBL),]

#Read in covariates and regroup data according to brain region
msmm_rnaseq.covariates=read.delim("msmm_rnaseq_B1_covariates.tsv",header = T,as.is = T)
msmm_rnaseq=vector(mode="list",length=length(unique(msmm_rnaseq.covariates$TISSUE)))
names(msmm_rnaseq)=c("BM_10","BM_22","BM_36")
msmm_rnaseq[[1]]=msmm_rnaseq_data[,grep("BM_10",colnames(msmm_rnaseq_data))]
msmm_rnaseq[[2]]=msmm_rnaseq_data[,grep("BM_22",colnames(msmm_rnaseq_data))]
msmm_rnaseq[[3]]=msmm_rnaseq_data[,grep("BM_36",colnames(msmm_rnaseq_data))]

#Set AOD as 90 for patients with AOD=89+
msmm_rnaseq.covariates$AOD[which(msmm_rnaseq.covariates$AOD=="89+")]="90"
msmm_rnaseq.covariates$AOD=as.numeric(msmm_rnaseq.covariates$AOD)
msmm_rnaseq.covariates=msmm_rnaseq.covariates[which(is.na(msmm_rnaseq.covariates$bbscore)==F),]

#Define low and high pathology samples and regroup into individual lists
low_Tangle_Samples=high_Tangle_Samples=low_Plaque_Samples=high_Plaque_Samples=vector(mode="list",length = 3)
names(low_Tangle_Samples)=names(high_Tangle_Samples)=names(low_Plaque_Samples)=names(high_Plaque_Samples)=names(msmm_rnaseq)
low_Tangle_Samples[[1]]=grep(pattern = "BM_10",msmm_rnaseq.covariates$Sample.ID[which(msmm_rnaseq.covariates$bbscore<=2)],value = T)
low_Tangle_Samples[[2]]=grep(pattern = "BM_22",msmm_rnaseq.covariates$Sample.ID[which(msmm_rnaseq.covariates$bbscore<=2)],value = T)
low_Tangle_Samples[[3]]=grep(pattern = "BM_36",msmm_rnaseq.covariates$Sample.ID[which(msmm_rnaseq.covariates$bbscore<=2)],value = T)

low_Plaque_Samples[[1]]=grep(pattern = "BM_10",msmm_rnaseq.covariates$Sample.ID[which(msmm_rnaseq.covariates$PlaqueMean<=2)],value = T)
low_Plaque_Samples[[2]]=grep(pattern = "BM_22",msmm_rnaseq.covariates$Sample.ID[which(msmm_rnaseq.covariates$PlaqueMean<=2)],value = T)
low_Plaque_Samples[[3]]=grep(pattern = "BM_36",msmm_rnaseq.covariates$Sample.ID[which(msmm_rnaseq.covariates$PlaqueMean<=2)],value = T)

high_Tangle_Samples[[1]]=grep(pattern = "BM_10",msmm_rnaseq.covariates$Sample.ID[which(msmm_rnaseq.covariates$bbscore>=5)],value = T)
high_Tangle_Samples[[2]]=grep(pattern = "BM_22",msmm_rnaseq.covariates$Sample.ID[which(msmm_rnaseq.covariates$bbscore>=5)],value = T)
high_Tangle_Samples[[3]]=grep(pattern = "BM_36",msmm_rnaseq.covariates$Sample.ID[which(msmm_rnaseq.covariates$bbscore>=5)],value = T)

high_Plaque_Samples[[1]]=grep(pattern = "BM_10",msmm_rnaseq.covariates$Sample.ID[which(msmm_rnaseq.covariates$PlaqueMean>=20)],value = T)
high_Plaque_Samples[[2]]=grep(pattern = "BM_22",msmm_rnaseq.covariates$Sample.ID[which(msmm_rnaseq.covariates$PlaqueMean>=20)],value = T)
high_Plaque_Samples[[3]]=grep(pattern = "BM_36",msmm_rnaseq.covariates$Sample.ID[which(msmm_rnaseq.covariates$PlaqueMean>=20)],value = T)

#Normalise counts for each pathology/brain region
msmm_rnaseq.y_Tangle=msmm_rnaseq.y_Plaques=vector(mode = "list",length = 3)
names(msmm_rnaseq.y_Tangle)=names(msmm_rnaseq.y_Plaques)=names(msmm_rnaseq)
for (m in 1:length(names(msmm_rnaseq.y_Tangle))){
  msmm_rnaseq.y_Tangle[[m]]=DGEList(counts = msmm_rnaseq[[m]][,c(low_Tangle_Samples[[m]],high_Tangle_Samples[[m]])],group = factor(c(rep(1,length(low_Tangle_Samples[[m]])),rep(2,length(high_Tangle_Samples[[m]])))))
  msmm_rnaseq.y_Plaques[[m]]=DGEList(counts = msmm_rnaseq[[m]][,c(low_Plaque_Samples[[m]],high_Plaque_Samples[[m]])],group = factor(c(rep(1,length(low_Plaque_Samples[[m]])),rep(2,length(high_Plaque_Samples[[m]])))))
}
msmm_rnaseq.y_Tangle=lapply(msmm_rnaseq.y_Tangle,calcNormFactors,method="TMM")
msmm_rnaseq.y_Plaques=lapply(msmm_rnaseq.y_Plaques,calcNormFactors,method="TMM")
msmm_rnaseq.Tangle_counts=lapply(msmm_rnaseq.y_Tangle,cpm,normalized.lib.sizes = T)
msmm_rnaseq.Plaque_counts=lapply(msmm_rnaseq.y_Plaques,cpm,normalized.lib.sizes = T)

#Rearrange counts data to suit WGCNA functions
msmm_rnaseq.TangleExprs2=vector(mode = "list",length = 3)
msmm_rnaseq.TangleExprs2[[1]]=msmm_rnaseq.TangleExprs2[[2]]=msmm_rnaseq.TangleExprs2[[3]]=vector(mode = "list",length = 2)
names(msmm_rnaseq.TangleExprs2)=names(msmm_rnaseq)
names(msmm_rnaseq.TangleExprs2[[1]])=names(msmm_rnaseq.TangleExprs2[[2]])=names(msmm_rnaseq.TangleExprs2[[3]])=c("low","high")
msmm_rnaseq.PlaqueExprs2=vector(mode = "list",length = 3)
names(msmm_rnaseq.PlaqueExprs2)=names(msmm_rnaseq)
msmm_rnaseq.PlaqueExprs2[[1]]=msmm_rnaseq.PlaqueExprs2[[2]]=msmm_rnaseq.PlaqueExprs2[[3]]=vector(mode="list",length=2)
names(msmm_rnaseq.PlaqueExprs2[[1]])=names(msmm_rnaseq.PlaqueExprs2[[2]])=names(msmm_rnaseq.PlaqueExprs2[[3]])=c("low","high")

for (j in 1:length(msmm_rnaseq.TangleExprs2)){
  msmm_rnaseq.TangleExprs2[[j]][[1]]=msmm_rnaseq.Tangle_counts[[j]][,low_Tangle_Samples[[j]]]
  msmm_rnaseq.TangleExprs2[[j]][[2]]=msmm_rnaseq.Tangle_counts[[j]][,high_Tangle_Samples[[j]]]
  write.table(msmm_rnaseq.TangleExprs2[[j]][[1]],paste(names(msmm_rnaseq.TangleExprs2)[j],"_TangleExprs_low.csv",sep=""),sep = ",",col.names = T,row.names = T,quote=F)
  write.table(msmm_rnaseq.TangleExprs2[[j]][[2]],paste(names(msmm_rnaseq.TangleExprs2)[j],"_TangleExprs_high.csv",sep=""),sep = ",",col.names = T,row.names = T,quote=F)
}
for (j in 1:length(names(msmm_rnaseq.PlaqueExprs2))){
  msmm_rnaseq.PlaqueExprs2[[j]][[1]]=msmm_rnaseq.Plaque_counts[[j]][,low_Plaque_Samples[[j]]]
  msmm_rnaseq.PlaqueExprs2[[j]][[2]]=msmm_rnaseq.Plaque_counts[[j]][,high_Plaque_Samples[[j]]]
  write.table(msmm_rnaseq.PlaqueExprs2[[j]][[1]],paste(names(msmm_rnaseq.PlaqueExprs2)[j],"_PlaqueExprs_low.csv",sep=""),sep = ",",col.names = T,row.names = T,quote=F)
  write.table(msmm_rnaseq.PlaqueExprs2[[j]][[2]],paste(names(msmm_rnaseq.PlaqueExprs2)[j],"_PlaqueExprs_high.csv",sep=""),sep = ",",col.names = T,row.names = T,quote=F)
}
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# msmm_rnaseq.TangleExprs2_low.sft=lapply(lapply(lapply(msmm_rnaseq.TangleExprs2,`[[`,1),t),pickSoftThreshold,dataIsExpr = TRUE,powerVector = powers,corFnc = cor,corOptions = list(use = 'p'),networkType = "unsigned")
# msmm_rnaseq.TangleExprs2_high.sft=lapply(lapply(lapply(msmm_rnaseq.TangleExprs2,`[[`,2),t),pickSoftThreshold,dataIsExpr = TRUE,powerVector = powers,corFnc = cor,corOptions = list(use = 'p'),networkType = "unsigned")
# msmm_rnaseq.PlaqueExprs2_low.sft=lapply(lapply(lapply(msmm_rnaseq.PlaqueExprs2,`[[`,1),t),pickSoftThreshold,dataIsExpr = TRUE,powerVector = powers,corFnc = cor,corOptions = list(use = 'p'),networkType = "unsigned")
# msmm_rnaseq.PlaqueExprs2_high.sft=lapply(lapply(lapply(msmm_rnaseq.PlaqueExprs2,`[[`,2),t),pickSoftThreshold,dataIsExpr = TRUE,powerVector = powers,corFnc = cor,corOptions = list(use = 'p'),networkType = "unsigned")

msmm_rnaseq.TangleExprs.sft=vector(mode = "list",length = 3)
msmm_rnaseq.TangleExprs.sft[[1]]=msmm_rnaseq.TangleExprs.sft[[2]]=msmm_rnaseq.TangleExprs.sft[[3]]=vector(mode = "list",length = 2)
names(msmm_rnaseq.TangleExprs.sft[[1]])=names(msmm_rnaseq.TangleExprs.sft[[2]])=names(msmm_rnaseq.TangleExprs.sft[[3]])=c("low","high")
names(msmm_rnaseq.TangleExprs.sft)=names(msmm_rnaseq)

msmm_rnaseq.PlaqueExprs.sft=vector(mode = "list",length = 3)
msmm_rnaseq.PlaqueExprs.sft[[1]]=msmm_rnaseq.PlaqueExprs.sft[[2]]=msmm_rnaseq.PlaqueExprs.sft[[3]]=vector(mode = "list",length = 2)
names(msmm_rnaseq.PlaqueExprs.sft[[1]])=names(msmm_rnaseq.PlaqueExprs.sft[[2]])=names(msmm_rnaseq.PlaqueExprs.sft[[3]])=c("low","high")
names(msmm_rnaseq.PlaqueExprs.sft)=names(msmm_rnaseq)

determineSoftPowerWGCNA(data = msmm_rnaseq.TangleExprs2$BM_10$low,"msmm_rnaseq_Tangle_FP_low.png",0.75)
determineSoftPowerWGCNA(data = msmm_rnaseq.TangleExprs2$BM_10$high,"msmm_rnaseq_Tangle_FP_high.png",0.75)
determineSoftPowerWGCNA(data = msmm_rnaseq.TangleExprs2$BM_22$low,"msmm_rnaseq_Tangle_STG_low.png",0.75)
determineSoftPowerWGCNA(data = msmm_rnaseq.TangleExprs2$BM_22$high,"msmm_rnaseq_Tangle_STG_high.png",0.75)
determineSoftPowerWGCNA(data = msmm_rnaseq.TangleExprs2$BM_36$low,"msmm_rnaseq_Tangle_PHG_low.png",0.75)
determineSoftPowerWGCNA(data = msmm_rnaseq.TangleExprs2$BM_36$high,"msmm_rnaseq_Tangle_PHG_high.png",0.75)

# determineSoftPowerWGCNA(data = msmm_rnaseq.PlaqueExprs2$BM_10$low,"msmm_rnaseq_PLQ_FP_low",1)
# determineSoftPowerWGCNA(data = msmm_rnaseq.PlaqueExprs2$BM_10$high,"msmm_rnaseq_PLQ_FP_high",1) Too few genes
determineSoftPowerWGCNA(data = msmm_rnaseq.PlaqueExprs2$BM_22$low,"msmm_rnaseq_PLQ_STG_low.png",0.75)
determineSoftPowerWGCNA(data = msmm_rnaseq.PlaqueExprs2$BM_22$high,"msmm_rnaseq_PLQ_STG_high.png",0.75)
determineSoftPowerWGCNA(data = msmm_rnaseq.PlaqueExprs2$BM_36$low,"msmm_rnaseq_PLQ_PHG_low.png",0.75)
determineSoftPowerWGCNA(data = msmm_rnaseq.PlaqueExprs2$BM_36$high,"msmm_rnaseq_PLQ_PHG_high.png",0.75)

msmm_rnaseq.TangleExprs.sft$BM_10$low=6 #8 for data including non-Entrez
msmm_rnaseq.TangleExprs.sft$BM_10$high=8 #16
msmm_rnaseq.TangleExprs.sft$BM_22$low=4 #5
msmm_rnaseq.TangleExprs.sft$BM_22$high=6 #14
msmm_rnaseq.TangleExprs.sft$BM_36$low=9 #8
msmm_rnaseq.TangleExprs.sft$BM_36$high=5 #12

msmm_rnaseq.PlaqueExprs.sft$BM_22$low=4 #3
msmm_rnaseq.PlaqueExprs.sft$BM_22$high=20 #14
msmm_rnaseq.PlaqueExprs.sft$BM_36$low=5 #5
msmm_rnaseq.PlaqueExprs.sft$BM_36$high=18 #10

msmm_rnaseq.TangleExprs.wgcna=vector(mode = "list",length = 3)
msmm_rnaseq.TangleExprs.wgcna[[1]]=msmm_rnaseq.TangleExprs.wgcna[[2]]=msmm_rnaseq.TangleExprs.wgcna[[3]]=vector(mode = "list",length = 2)
names(msmm_rnaseq.TangleExprs.wgcna[[1]])=names(msmm_rnaseq.TangleExprs.wgcna[[2]])=names(msmm_rnaseq.TangleExprs.wgcna[[3]])=c("low","high")
names(msmm_rnaseq.TangleExprs.wgcna)=names(msmm_rnaseq)

msmm_rnaseq.PlaqueExprs.wgcna=vector(mode = "list",length = 3)
msmm_rnaseq.PlaqueExprs.wgcna[[1]]=msmm_rnaseq.PlaqueExprs.wgcna[[2]]=msmm_rnaseq.PlaqueExprs.wgcna[[3]]=vector(mode = "list",length = 2)
names(msmm_rnaseq.PlaqueExprs.wgcna[[1]])=names(msmm_rnaseq.PlaqueExprs.wgcna[[2]])=names(msmm_rnaseq.PlaqueExprs.wgcna[[3]])=c("low","high")
names(msmm_rnaseq.PlaqueExprs.wgcna)=names(msmm_rnaseq)

msmm_rnaseq.TangleExprs.wgcna$BM_10$low=runWGCNA(data1 = msmm_rnaseq.TangleExprs2$BM_10$low,propGenes = 0.75,softPower = msmm_rnaseq.TangleExprs.sft$BM_10$low,signedNetwork = T)
msmm_rnaseq.TangleExprs.wgcna$BM_10$high=runWGCNA(data1 = msmm_rnaseq.TangleExprs2$BM_10$high,propGenes = 0.75,softPower = msmm_rnaseq.TangleExprs.sft$BM_10$high,signedNetwork = T)
msmm_rnaseq.TangleExprs.wgcna$BM_22$low=runWGCNA(data1 = msmm_rnaseq.TangleExprs2$BM_22$low,propGenes = 0.75,softPower = msmm_rnaseq.TangleExprs.sft$BM_22$low,signedNetwork = T)
msmm_rnaseq.TangleExprs.wgcna$BM_22$high=runWGCNA(data1 = msmm_rnaseq.TangleExprs2$BM_22$high,propGenes = 0.75,softPower = msmm_rnaseq.TangleExprs.sft$BM_22$high,signedNetwork = T)
msmm_rnaseq.TangleExprs.wgcna$BM_36$low=runWGCNA(data1 = msmm_rnaseq.TangleExprs2$BM_36$low,propGenes = 0.75,softPower = msmm_rnaseq.TangleExprs.sft$BM_36$low,signedNetwork = T)
msmm_rnaseq.TangleExprs.wgcna$BM_36$high=runWGCNA(data1 = msmm_rnaseq.TangleExprs2$BM_36$high,propGenes = 0.75,softPower = msmm_rnaseq.TangleExprs.sft$BM_36$high,signedNetwork = T)

msmm_rnaseq.PlaqueExprs.wgcna$BM_22$low=runWGCNA(data1 = msmm_rnaseq.PlaqueExprs2$BM_22$low,propGenes = 0.75,softPower = msmm_rnaseq.PlaqueExprs.sft$BM_22$low,signedNetwork = T)
msmm_rnaseq.PlaqueExprs.wgcna$BM_22$high=runWGCNA(data1 = msmm_rnaseq.PlaqueExprs2$BM_22$high,propGenes = 0.75,softPower = msmm_rnaseq.PlaqueExprs.sft$BM_22$high,signedNetwork = T)
msmm_rnaseq.PlaqueExprs.wgcna$BM_36$low=runWGCNA(data1 = msmm_rnaseq.PlaqueExprs2$BM_36$low,propGenes = 0.75,softPower = msmm_rnaseq.PlaqueExprs.sft$BM_36$low,signedNetwork = T)
msmm_rnaseq.PlaqueExprs.wgcna$BM_36$high=runWGCNA(data1 = msmm_rnaseq.PlaqueExprs2$BM_36$high,propGenes = 0.75,softPower = msmm_rnaseq.PlaqueExprs.sft$BM_36$high,signedNetwork = T)

msmm_rnaseq.TangleExprs.EigenGenes=vector(mode = "list",length = 3)
msmm_rnaseq.TangleExprs.EigenGenes[[1]]=msmm_rnaseq.TangleExprs.EigenGenes[[2]]=msmm_rnaseq.TangleExprs.EigenGenes[[3]]=vector(mode = "list",length = 2)
names(msmm_rnaseq.TangleExprs.EigenGenes)=names(msmm_rnaseq)
names(msmm_rnaseq.TangleExprs.EigenGenes[[1]])=names(msmm_rnaseq.TangleExprs.EigenGenes[[2]])=names(msmm_rnaseq.TangleExprs.EigenGenes[[3]])=c("low","high")
for(k in 1:length(names(msmm_rnaseq.TangleExprs.EigenGenes))){
  msmm_rnaseq.TangleExprs.EigenGenes[[k]][[1]]=calculateModuleEigengenes(referenceDataset=msmm_rnaseq.TangleExprs.wgcna[[k]][[1]],split = 2,minClusterSize = 20)
  msmm_rnaseq.TangleExprs.EigenGenes[[k]][[2]]=calculateModuleEigengenes(referenceDataset=msmm_rnaseq.TangleExprs.wgcna[[k]][[2]],split = 2,minClusterSize = 20)
}

msmm_rnaseq.PlaqueExprs.EigenGenes=vector(mode = "list",length = 2)
msmm_rnaseq.PlaqueExprs.EigenGenes[[1]]=msmm_rnaseq.PlaqueExprs.EigenGenes[[2]]=vector(mode = "list",length = 2)
names(msmm_rnaseq.PlaqueExprs.EigenGenes)=names(msmm_rnaseq)[2:3]
names(msmm_rnaseq.PlaqueExprs.EigenGenes[[1]])=names(msmm_rnaseq.PlaqueExprs.EigenGenes[[2]])=c("low","high")
for(k in 1:length(names(msmm_rnaseq.PlaqueExprs.EigenGenes))){
  msmm_rnaseq.PlaqueExprs.EigenGenes[[k]][[1]]=calculateModuleEigengenes(referenceDataset=msmm_rnaseq.PlaqueExprs.wgcna[[k]][[1]],split = 2,minClusterSize = 20)
  msmm_rnaseq.PlaqueExprs.EigenGenes[[k]][[2]]=calculateModuleEigengenes(referenceDataset=msmm_rnaseq.PlaqueExprs.wgcna[[k]][[2]],split = 2,minClusterSize = 20)
}

msmm_rnaseq.TangleExprs.moduleMembership=vector(mode = "list",length = 3)
msmm_rnaseq.TangleExprs.moduleMembership[[1]]=msmm_rnaseq.TangleExprs.moduleMembership[[2]]=msmm_rnaseq.TangleExprs.moduleMembership[[3]]=vector(mode = "list",length = 2)
names(msmm_rnaseq.TangleExprs.moduleMembership)=names(msmm_rnaseq)
names(msmm_rnaseq.TangleExprs.moduleMembership[[1]])=names(msmm_rnaseq.TangleExprs.moduleMembership[[2]])=names(msmm_rnaseq.TangleExprs.moduleMembership[[3]])=c("low","high")
for (k in 1:length(names(msmm_rnaseq.TangleExprs.moduleMembership))){
  msmm_rnaseq.TangleExprs.moduleMembership[[k]][[1]]=moduleMembership(referenceDataset = msmm_rnaseq.TangleExprs.wgcna[[k]][[1]],MEs = msmm_rnaseq.TangleExprs.EigenGenes[[k]][[1]],split = 2,kME_threshold = 0.75)
  msmm_rnaseq.TangleExprs.moduleMembership[[k]][[2]]=moduleMembership(referenceDataset = msmm_rnaseq.TangleExprs.wgcna[[k]][[2]],MEs = msmm_rnaseq.TangleExprs.EigenGenes[[k]][[2]],split = 2,kME_threshold = 0.75)
}

msmm_rnaseq.PlaqueExprs.moduleMembership=vector(mode = "list",length = 2)
msmm_rnaseq.PlaqueExprs.moduleMembership[[1]]=msmm_rnaseq.PlaqueExprs.moduleMembership[[2]]=vector(mode = "list",length = 2)
names(msmm_rnaseq.PlaqueExprs.moduleMembership)=names(msmm_rnaseq)[2:3]
names(msmm_rnaseq.PlaqueExprs.moduleMembership[[1]])=names(msmm_rnaseq.PlaqueExprs.moduleMembership[[2]])=c("low","high")
for (k in 1:length(names(msmm_rnaseq.PlaqueExprs.moduleMembership))){
  msmm_rnaseq.PlaqueExprs.moduleMembership[[k]][[1]]=moduleMembership(referenceDataset = msmm_rnaseq.PlaqueExprs.wgcna[[k]][[1]],MEs = msmm_rnaseq.PlaqueExprs.EigenGenes[[k]][[1]],split = 2,kME_threshold = 0.75)
  msmm_rnaseq.PlaqueExprs.moduleMembership[[k]][[2]]=moduleMembership(referenceDataset = msmm_rnaseq.PlaqueExprs.wgcna[[k]][[2]],MEs = msmm_rnaseq.PlaqueExprs.EigenGenes[[k]][[2]],split = 2,kME_threshold = 0.75)
}

save.image(file = "msmm_rnaseq_wgcna.RData")
