library(limma)
library(Biobase)
library(jetset)
library(pathprint)
library(metaArray)
library(doMC)
library(pheatmap)
library(entropy)
library(Rarity)
library(org.Hs.eg.db)
source("/Users/sandeepamberkar/Work/Software/RUN-WGCNA/WGCNA_functions.R")
setwd("/Users/sandeepamberkar/Work/Data/MSMM-MSBB-HBTRC-Data/msbb_array19/Normalised_Data/")
load(file="msbb_corr_PLQ_TrackGenes005.RData") 
load(file="msbb_corr_NTr_TrackGenes005.RData") 
msbb_array19.NTr_corrData=msbb_array19.PLQ_corrData=vector(mode = "list",length = 3)
names(msbb_array19.NTr_corrData)=names(msbb_array19.PLQ_corrData)=c("FP","PHG","STG")
msbb_array19.NTr_corrData[[1]]=read.table("msbb_array19.corr_TrackGenes_NTr_Genes005_FP.txt",header = T,row.names = 1,as.is = T)
msbb_array19.NTr_corrData[[2]]=read.table("msbb_array19.corr_TrackGenes_NTr_Genes005_PHG.txt",header = T,row.names = 1,as.is = T)
msbb_array19.NTr_corrData[[3]]=read.table("msbb_array19.corr_TrackGenes_NTr_Genes005_STG.txt",header = T,row.names = 1,as.is = T)

msbb_array19.PLQ_corrData[[1]]=read.table("msbb_array19.corr_TrackGenes_PLQ_Genes005_FP.txt",header = T,row.names = 1,as.is = T)
msbb_array19.PLQ_corrData[[2]]=read.table("msbb_array19.corr_TrackGenes_PLQ_Genes005_PHG.txt",header = T,row.names = 1,as.is = T)
msbb_array19.PLQ_corrData[[3]]=read.table("msbb_array19.corr_TrackGenes_PLQ_Genes005_STG.txt",header = T,row.names = 1,as.is = T)

for (n in 1:length(names(msbb_array19.NTr_corrData))){
   rownames(msbb_array19.NTr_corrData[[n]])=paste(rep("ENTREZ",dim(msbb_array19.NTr_corrData[[n]])[1]),rownames(msbb_array19.NTr_corrData[[n]]),sep = "_")
   rownames(msbb_array19.PLQ_corrData[[n]])=paste(rep("ENTREZ",dim(msbb_array19.PLQ_corrData[[n]])[1]),rownames(msbb_array19.PLQ_corrData[[n]]),sep = "_")
}

msbb_array19.files=list.files(pattern = "*.tsv",full.names = T)
msbb_array19.files=msbb_array19.files[-c(18:19)]
msbb_array19=mapply(FUN = read.table,msbb_array19.files,MoreArgs = list(header=T,sep="\t",as.is=T))
names(msbb_array19)=c("AC","CN","DLPFC","FP","HP","IFG","ITG","MTG","OVC","PHG","PCC","PCG","PFC","PTMN","SPL","STG","TP")
#Select regions that show known regions to be affected!!!
#msbb_array19=msbb_array19[c("PHG","HP","AC","STG","SPL","DLPFC","OVC")]
#msbb_array19.SampleperBrainRegion=list()
msbb_array19.SampleperBrainRegion=list()
# for (a in 1:length(msbb_array19.covariates$BrainBank)){
#   msbb_array19.SampleperBrainRegion[[a]]=paste(names(unlist(lapply(lapply(lapply(lapply(msbb_array19,`[`,-c(1:4)),colnames),sort),function(x)which(msbb_array19.covariates$BrainBank[a]%in%x)))),collapse = ",")
# }
msbb_array19.covariates=read.delim2("../AMP-AD_MSBB_MSSM_covariates_mRNA_AffymetrixU133AB.tsv",header = T,as.is = T)
msbb_array19.covariates$Age[which(msbb_array19.covariates$Age=="89+")]=90 #For simplicity, all samples with AOD as 89+ were labeled 90
msbb_array19.covariates$Age=as.integer(msbb_array19.covariates$Age)
msbb_array19.covariates$pH=as.numeric(msbb_array19.covariates$pH)
msbb_array19.covariates$CDR=as.numeric(msbb_array19.covariates$CDR)
msbb_array19.covariates$PLQ_Mn=as.numeric(msbb_array19.covariates$PLQ_Mn)
msbb_array19.covariates$BrainBank=paste("X",msbb_array19.covariates$BrainBank,sep="")
for (a in 1:length(msbb_array19.covariates$BrainBank)){
  msbb_array19.SampleperBrainRegion[[a]]=paste(names(unlist(lapply(lapply(lapply(lapply(msbb_array19,`[`,-c(1:4)),colnames),sort),function(x)which(msbb_array19.covariates$BrainBank[a]%in%x)))),collapse = ",")
  names(msbb_array19.SampleperBrainRegion)[a]=msbb_array19.covariates$BrainBank[a]
}

noRegion=names(unlist(lapply(msbb_array19.SampleperBrainRegion,function(x)x[which(x=="")])))
msbb_array19.SampleperBrainRegion=msbb_array19.SampleperBrainRegion[-which(names(msbb_array19.SampleperBrainRegion)%in%noRegion)]
msbbFiltered_indices=lapply(lapply(msbb_array19,colnames),function(x)which(x%in%names(msbb_array19.SampleperBrainRegion)))
msbb_array19.covariates2=msbb_array19.covariates[-which(msbb_array19.covariates$BrainBank%in%noRegion),]
msbb_array19.Samples_FP=names(which(unlist(lapply(lapply(msbb_array19.SampleperBrainRegion,grep,pattern="FP",value=T),length))>0))
msbb_array19.Samples_PHG=names(which(unlist(lapply(lapply(msbb_array19.SampleperBrainRegion,grep,pattern="PHG",value=T),length))>0))
msbb_array19.Samples_STG=names(which(unlist(lapply(lapply(msbb_array19.SampleperBrainRegion,grep,pattern="STG",value=T),length))>0))

low_Tangle_Samples=high_Tangle_Samples=low_Plaque_Samples=high_Plaque_Samples=vector(mode="list",length = 3)
names(low_Tangle_Samples)=names(high_Tangle_Samples)=names(low_Plaque_Samples)=names(high_Plaque_Samples)=names(msbb_array19)[c(4,10,16)]
low_Tangle_Samples[[1]]=msbb_array19.covariates2$BrainBank[which((msbb_array19.covariates2$BrainBank%in%msbb_array19.Samples_FP)&(msbb_array19.covariates2$Braak<=2))]
low_Tangle_Samples[[2]]=msbb_array19.covariates2$BrainBank[which((msbb_array19.covariates2$BrainBank%in%msbb_array19.Samples_PHG)&(msbb_array19.covariates2$Braak<=2))]
low_Tangle_Samples[[3]]=msbb_array19.covariates2$BrainBank[which((msbb_array19.covariates2$BrainBank%in%msbb_array19.Samples_STG)&(msbb_array19.covariates2$Braak<=2))]

high_Tangle_Samples[[1]]=msbb_array19.covariates2$BrainBank[which((msbb_array19.covariates2$BrainBank%in%msbb_array19.Samples_FP)&(msbb_array19.covariates2$Braak>=5))]
high_Tangle_Samples[[2]]=msbb_array19.covariates2$BrainBank[which((msbb_array19.covariates2$BrainBank%in%msbb_array19.Samples_PHG)&(msbb_array19.covariates2$Braak>=5))]
high_Tangle_Samples[[3]]=msbb_array19.covariates2$BrainBank[which((msbb_array19.covariates2$BrainBank%in%msbb_array19.Samples_STG)&(msbb_array19.covariates2$Braak>=5))]

low_Plaque_Samples[[1]]=msbb_array19.covariates2$BrainBank[which((msbb_array19.covariates2$BrainBank%in%msbb_array19.Samples_FP)&(msbb_array19.covariates2$PLQ_Mn<=2))]
low_Plaque_Samples[[2]]=msbb_array19.covariates2$BrainBank[which((msbb_array19.covariates2$BrainBank%in%msbb_array19.Samples_PHG)&(msbb_array19.covariates2$PLQ_Mn<=2))]
low_Plaque_Samples[[3]]=msbb_array19.covariates2$BrainBank[which((msbb_array19.covariates2$BrainBank%in%msbb_array19.Samples_STG)&(msbb_array19.covariates2$PLQ_Mn<=2))]

high_Plaque_Samples[[1]]=msbb_array19.covariates2$BrainBank[which((msbb_array19.covariates2$BrainBank%in%msbb_array19.Samples_FP)&(msbb_array19.covariates2$PLQ_Mn>=20))]
high_Plaque_Samples[[2]]=msbb_array19.covariates2$BrainBank[which((msbb_array19.covariates2$BrainBank%in%msbb_array19.Samples_PHG)&(msbb_array19.covariates2$PLQ_Mn>=20))]
high_Plaque_Samples[[3]]=msbb_array19.covariates2$BrainBank[which((msbb_array19.covariates2$BrainBank%in%msbb_array19.Samples_STG)&(msbb_array19.covariates2$PLQ_Mn>=20))]

msbb_array.TangleExprs=msbb_array.PlaqueExprs=vector(mode = "list",length = 3)
names(msbb_array.TangleExprs)=names(msbb_array.PlaqueExprs)=names(msbb_array19.NTr_corrData)
msbb_array.TangleExprs[[1]]=msbb_array.TangleExprs[[2]]=msbb_array.TangleExprs[[3]]=msbb_array.PlaqueExprs[[1]]=msbb_array.PlaqueExprs[[2]]=msbb_array.PlaqueExprs[[3]]=vector(mode = "list",length = 2)
names(msbb_array.TangleExprs[[1]])=names(msbb_array.TangleExprs[[2]])=names(msbb_array.TangleExprs[[3]])=names(msbb_array.PlaqueExprs[[1]])=names(msbb_array.PlaqueExprs[[2]])=names(msbb_array.PlaqueExprs[[3]])=c("low","high")
for (i in 1:length(names(msbb_array.TangleExprs))){
  msbb_array.TangleExprs[[i]][[1]]=msbb_array19.NTr_corrData[[i]][,low_Tangle_Samples[[i]]]
  msbb_array.TangleExprs[[i]][[2]]=msbb_array19.NTr_corrData[[i]][,high_Tangle_Samples[[i]]]
  msbb_array.PlaqueExprs[[i]][[1]]=msbb_array19.PLQ_corrData[[i]][,low_Plaque_Samples[[i]]]
  msbb_array.PlaqueExprs[[i]][[2]]=msbb_array19.PLQ_corrData[[i]][,high_Plaque_Samples[[i]]]
}
# 
# determineSoftPowerWGCNA(data = msbb_array.TangleExprs$FP$low,"msbb_array_Tangle_FP_low.png",1)
# determineSoftPowerWGCNA(data = msbb_array.TangleExprs$FP$high,"msbb_array_Tangle_FP_high.png",1)
# determineSoftPowerWGCNA(data = msbb_array.TangleExprs$STG$low,"msbb_array_Tangle_STG_low.png",1)
# determineSoftPowerWGCNA(data = msbb_array.TangleExprs$STG$high,"msbb_array_Tangle_STG_high.png",1)
# determineSoftPowerWGCNA(data = msbb_array.TangleExprs$PHG$low,"msbb_array_Tangle_PHG_low.png",1)
# determineSoftPowerWGCNA(data = msbb_array.TangleExprs$PHG$high,"msbb_array_Tangle_PHG_high.png",1)
# 
# determineSoftPowerWGCNA(data = msbb_array.PlaqueExprs$FP$low,"msbb_array_Plaque_FP_low.png",1)
# determineSoftPowerWGCNA(data = msbb_array.PlaqueExprs$FP$high,"msbb_array_Plaque_FP_high.png",1)
# determineSoftPowerWGCNA(data = msbb_array.PlaqueExprs$STG$low,"msbb_array_Plaque_STG_low.png",1)
# determineSoftPowerWGCNA(data = msbb_array.PlaqueExprs$STG$high,"msbb_array_Plaque_STG_high.png",1)
# determineSoftPowerWGCNA(data = msbb_array.PlaqueExprs$PHG$low,"msbb_array_Plaque_PHG_low.png",1)
# determineSoftPowerWGCNA(data = msbb_array.PlaqueExprs$PHG$high,"msbb_array_Plaque_PHG_high.png",1)
# 
# msbb_array.TangleExprs.sft=vector(mode = "list",length = 3)
# msbb_array.TangleExprs.sft[[1]]=msbb_array.TangleExprs.sft[[2]]=msbb_array.TangleExprs.sft[[3]]=vector(mode = "list",length = 2)
# names(msbb_array.TangleExprs.sft[[1]])=names(msbb_array.TangleExprs.sft[[2]])=names(msbb_array.TangleExprs.sft[[3]])=c("low","high")
# names(msbb_array.TangleExprs.sft)=names(msbb_array19.NTr_corrData)
# 
# msbb_array.PlaqueExprs.sft=vector(mode = "list",length = 3)
# msbb_array.PlaqueExprs.sft[[1]]=msbb_array.PlaqueExprs.sft[[2]]=msbb_array.PlaqueExprs.sft[[3]]=vector(mode = "list",length = 2)
# names(msbb_array.PlaqueExprs.sft[[1]])=names(msbb_array.PlaqueExprs.sft[[2]])=names(msbb_array.PlaqueExprs.sft[[3]])=c("low","high")
# names(msbb_array.PlaqueExprs.sft)=names(msbb_array19.PLQ_corrData)
# 
# msbb_array.TangleExprs.sft$FP$low=7 #8 for data including non-Entrez
# msbb_array.TangleExprs.sft$FP$high=9 #16
# msbb_array.TangleExprs.sft$STG$low=12 #5
# msbb_array.TangleExprs.sft$STG$high=14 #14
# msbb_array.TangleExprs.sft$PHG$low=5 #8
# msbb_array.TangleExprs.sft$PHG$high=12 #12
# 
# msbb_array.PlaqueExprs.sft$FP$low=7#3
# msbb_array.PlaqueExprs.sft$FP$high=7
# msbb_array.PlaqueExprs.sft$STG$low=7 #3
# msbb_array.PlaqueExprs.sft$STG$high=14 #14
# msbb_array.PlaqueExprs.sft$PHG$low=8 #5
# msbb_array.PlaqueExprs.sft$PHG$high=16 #10
# 
# msbb_array.TangleExprs.wgcna=vector(mode = "list",length = 3)
# msbb_array.TangleExprs.wgcna[[1]]=msbb_array.TangleExprs.wgcna[[2]]=msbb_array.TangleExprs.wgcna[[3]]=vector(mode = "list",length = 2)
# names(msbb_array.TangleExprs.wgcna[[1]])=names(msbb_array.TangleExprs.wgcna[[2]])=names(msbb_array.TangleExprs.wgcna[[3]])=c("low","high")
# names(msbb_array.TangleExprs.wgcna)=names(msbb_array19.NTr_corrData)
# 
# msbb_array.PlaqueExprs.wgcna=vector(mode = "list",length = 3)
# msbb_array.PlaqueExprs.wgcna[[1]]=msbb_array.PlaqueExprs.wgcna[[2]]=msbb_array.PlaqueExprs.wgcna[[3]]=vector(mode = "list",length = 2)
# names(msbb_array.PlaqueExprs.wgcna[[1]])=names(msbb_array.PlaqueExprs.wgcna[[2]])=names(msbb_array.PlaqueExprs.wgcna[[3]])=c("low","high")
# names(msbb_array.PlaqueExprs.wgcna)=names(msbb_array19.PLQ_corrData)
# 
# msbb_array.TangleExprs.wgcna$FP$low=runWGCNA(data1 = msbb_array.TangleExprs$FP$low,propGenes = 1,softPower = msbb_array.TangleExprs.sft$FP$low,signedNetwork = T)
# msbb_array.TangleExprs.wgcna$FP$high=runWGCNA(data1 = msbb_array.TangleExprs$FP$high,propGenes = 1,softPower = msbb_array.TangleExprs.sft$FP$high,signedNetwork = T)
# msbb_array.TangleExprs.wgcna$STG$low=runWGCNA(data1 = msbb_array.TangleExprs$STG$low,propGenes = 1,softPower = msbb_array.TangleExprs.sft$STG$low,signedNetwork = T)
# msbb_array.TangleExprs.wgcna$STG$high=runWGCNA(data1 = msbb_array.TangleExprs$STG$high,propGenes = 1,softPower = msbb_array.TangleExprs.sft$STG$high,signedNetwork = T)
# msbb_array.TangleExprs.wgcna$PHG$low=runWGCNA(data1 = msbb_array.TangleExprs$PHG$low,propGenes = 1,softPower = msbb_array.TangleExprs.sft$PHG$low,signedNetwork = T)
# msbb_array.TangleExprs.wgcna$PHG$high=runWGCNA(data1 = msbb_array.TangleExprs$PHG$high,propGenes = 1,softPower = msbb_array.TangleExprs.sft$PHG$high,signedNetwork = T)
# 
# msbb_array.PlaqueExprs.wgcna$FP$low=runWGCNA(data1 = msbb_array.PlaqueExprs$FP$low,propGenes = 1,softPower = msbb_array.PlaqueExprs.sft$FP$low,signedNetwork = T)
# msbb_array.PlaqueExprs.wgcna$FP$high=runWGCNA(data1 = msbb_array.PlaqueExprs$FP$high,propGenes = 1,softPower = msbb_array.PlaqueExprs.sft$FP$high,signedNetwork = T)
# msbb_array.PlaqueExprs.wgcna$STG$low=runWGCNA(data1 = msbb_array.PlaqueExprs$STG$low,propGenes = 1,softPower = msbb_array.PlaqueExprs.sft$STG$low,signedNetwork = T)
# msbb_array.PlaqueExprs.wgcna$STG$high=runWGCNA(data1 = msbb_array.PlaqueExprs$STG$high,propGenes = 1,softPower = msbb_array.PlaqueExprs.sft$STG$high,signedNetwork = T)
# msbb_array.PlaqueExprs.wgcna$PHG$low=runWGCNA(data1 = msbb_array.PlaqueExprs$PHG$low,propGenes = 1,softPower = msbb_array.PlaqueExprs.sft$PHG$low,signedNetwork = T)
# msbb_array.PlaqueExprs.wgcna$PHG$high=runWGCNA(data1 = msbb_array.PlaqueExprs$PHG$high,propGenes = 1,softPower = msbb_array.PlaqueExprs.sft$PHG$high,signedNetwork = T)
# 
# msbb_array.TangleExprs.EigenGenes=vector(mode = "list",length = 3)
# msbb_array.TangleExprs.EigenGenes[[1]]=msbb_array.TangleExprs.EigenGenes[[2]]=msbb_array.TangleExprs.EigenGenes[[3]]=vector(mode = "list",length = 2)
# names(msbb_array.TangleExprs.EigenGenes)=names(msbb_array19.NTr_corrData)
# names(msbb_array.TangleExprs.EigenGenes[[1]])=names(msbb_array.TangleExprs.EigenGenes[[2]])=names(msbb_array.TangleExprs.EigenGenes[[3]])=c("low","high")
# for(k in 1:length(names(msbb_array.TangleExprs.EigenGenes))){
#   msbb_array.TangleExprs.EigenGenes[[k]][[1]]=calculateModuleEigengenes(referenceDataset=msbb_array.TangleExprs.wgcna[[k]][[1]],split = 2,minClusterSize = 50)
#   msbb_array.TangleExprs.EigenGenes[[k]][[2]]=calculateModuleEigengenes(referenceDataset=msbb_array.TangleExprs.wgcna[[k]][[2]],split = 2,minClusterSize = 50)
# }
# 
# msbb_array.PlaqueExprs.EigenGenes=vector(mode = "list",length = 2)
# msbb_array.PlaqueExprs.EigenGenes[[1]]=msbb_array.PlaqueExprs.EigenGenes[[2]]=msbb_array.PlaqueExprs.EigenGenes[[3]]=vector(mode = "list",length = 2)
# names(msbb_array.PlaqueExprs.EigenGenes)=names(msbb_array19.PLQ_corrData)
# names(msbb_array.PlaqueExprs.EigenGenes[[1]])=names(msbb_array.PlaqueExprs.EigenGenes[[2]])=names(msbb_array.PlaqueExprs.EigenGenes[[3]])=c("low","high")
# for(k in 1:length(names(msbb_array.PlaqueExprs.EigenGenes))){
#   msbb_array.PlaqueExprs.EigenGenes[[k]][[1]]=calculateModuleEigengenes(referenceDataset=msbb_array.PlaqueExprs.wgcna[[k]][[1]],split = 2,minClusterSize = 50)
#   msbb_array.PlaqueExprs.EigenGenes[[k]][[2]]=calculateModuleEigengenes(referenceDataset=msbb_array.PlaqueExprs.wgcna[[k]][[2]],split = 2,minClusterSize = 50)
# }
# 
# msbb_array.TangleExprs.moduleMembership=vector(mode = "list",length = 3)
# msbb_array.TangleExprs.moduleMembership[[1]]=msbb_array.TangleExprs.moduleMembership[[2]]=msbb_array.TangleExprs.moduleMembership[[3]]=vector(mode = "list",length = 2)
# names(msbb_array.TangleExprs.moduleMembership)=names(msbb_array19.NTr_corrData)
# names(msbb_array.TangleExprs.moduleMembership[[1]])=names(msbb_array.TangleExprs.moduleMembership[[2]])=names(msbb_array.TangleExprs.moduleMembership[[3]])=c("low","high")
# for (k in 1:length(names(msbb_array.TangleExprs.moduleMembership))){
#   msbb_array.TangleExprs.moduleMembership[[k]][[1]]=moduleMembership(referenceDataset = msbb_array.TangleExprs.wgcna[[k]][[1]],MEs = msbb_array.TangleExprs.EigenGenes[[k]][[1]],split = 2,kME_threshold = 0.8)
#   msbb_array.TangleExprs.moduleMembership[[k]][[2]]=moduleMembership(referenceDataset = msbb_array.TangleExprs.wgcna[[k]][[2]],MEs = msbb_array.TangleExprs.EigenGenes[[k]][[2]],split = 2,kME_threshold = 0.8)
# }
# 
# msbb_array.PlaqueExprs.moduleMembership=vector(mode = "list",length = 2)
# msbb_array.PlaqueExprs.moduleMembership[[1]]=msbb_array.PlaqueExprs.moduleMembership[[2]]=msbb_array.PlaqueExprs.moduleMembership[[3]]=vector(mode = "list",length = 2)
# names(msbb_array.PlaqueExprs.moduleMembership)=names(msbb_array19.PLQ_corrData)
# names(msbb_array.PlaqueExprs.moduleMembership[[1]])=names(msbb_array.PlaqueExprs.moduleMembership[[2]])=names(msbb_array.PlaqueExprs.moduleMembership[[3]])=c("low","high")
# for (k in 1:length(names(msbb_array.PlaqueExprs.moduleMembership))){
#   msbb_array.PlaqueExprs.moduleMembership[[k]][[1]]=moduleMembership(referenceDataset = msbb_array.PlaqueExprs.wgcna[[k]][[1]],MEs = msbb_array.PlaqueExprs.EigenGenes[[k]][[1]],split = 2,kME_threshold = 0.8)
#   msbb_array.PlaqueExprs.moduleMembership[[k]][[2]]=moduleMembership(referenceDataset = msbb_array.PlaqueExprs.wgcna[[k]][[2]],MEs = msbb_array.PlaqueExprs.EigenGenes[[k]][[2]],split = 2,kME_threshold = 0.8)
# }
# 
# msbb_array.PlaqueExprs.moduleMembership.compareDF=vector(mode = "list",length = 3)
# names(msbb_array.PlaqueExprs.moduleMembership.compareDF)=names(msbb_array.PlaqueExprs)
# for (n in 1:length(names(msbb_array.PlaqueExprs.moduleMembership))){
#   module_colours=intersect(names(msbb_array.PlaqueExprs.moduleMembership[[n]][[1]]),names(msbb_array.PlaqueExprs.moduleMembership[[n]][[2]]))
#   msbb_array.PlaqueExprs.moduleMembership.compareDF[[n]]=matrix(NA,nrow=length(module_colours),ncol = length(module_colours))
#   rownames(msbb_array.PlaqueExprs.moduleMembership.compareDF[[n]])=colnames(msbb_array.PlaqueExprs.moduleMembership.compareDF[[n]])=module_colours
#   for (i in 1:length(module_colours)){
#       msbb_array.PlaqueExprs.moduleMembership.compareDF[[n]][i,]=unlist(unname(lapply(lapply(msbb_array.PlaqueExprs.moduleMembership[[n]][[1]],intersect,msbb_array.PlaqueExprs.moduleMembership[[n]][[2]][[module_colours[i]]]),length)))
#       msbb_array.PlaqueExprs.moduleMembership.compareDF[[n]][,i]=unlist(unname(lapply(lapply(msbb_array.PlaqueExprs.moduleMembership[[n]][[1]],intersect,msbb_array.PlaqueExprs.moduleMembership[[n]][[2]][[module_colours[i]]]),length)))
#     }
# }
# 
# msbb_array.TangleExprs.moduleMembership.compareDF=vector(mode = "list",length = 3)
# names(msbb_array.TangleExprs.moduleMembership.compareDF)=names(msbb_array.TangleExprs)
# for (n in 1:length(names(msbb_array.TangleExprs.moduleMembership))){
#   module_colours=intersect(names(msbb_array.TangleExprs.moduleMembership[[n]][[1]]),names(msbb_array.TangleExprs.moduleMembership[[n]][[2]]))
#   msbb_array.TangleExprs.moduleMembership.compareDF[[n]]=matrix(NA,nrow=length(module_colours),ncol = length(module_colours))
#   rownames(msbb_array.TangleExprs.moduleMembership.compareDF[[n]])=colnames(msbb_array.TangleExprs.moduleMembership.compareDF[[n]])=module_colours
#   for (i in 1:length(module_colours)){
#     msbb_array.TangleExprs.moduleMembership.compareDF[[n]][i,]=unlist(unname(lapply(lapply(msbb_array.TangleExprs.moduleMembership[[n]][[1]],intersect,msbb_array.TangleExprs.moduleMembership[[n]][[2]][[module_colours[i]]]),length)))
#     msbb_array.TangleExprs.moduleMembership.compareDF[[n]][,i]=unlist(unname(lapply(lapply(msbb_array.TangleExprs.moduleMembership[[n]][[1]],intersect,msbb_array.TangleExprs.moduleMembership[[n]][[2]][[module_colours[i]]]),length)))
#   }
# }

save.image(file = "msbb_array_corrGenes_wgcna.RData")

str(lapply(lapply(lapply(lapply(lapply(msbb_array.TangleExprs.moduleMembership,`[[`,2),`[[`,3),gsub,pattern = "ENTREZ_",replacement=""),enrichKEGG,organism="hsa",pvalueCutoff=0.05),summary))