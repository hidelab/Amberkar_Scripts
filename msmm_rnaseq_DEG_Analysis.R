library(limma)
library(Biobase)
library(metaArray)
library(doMC)
library(pheatmap)
library(clusterProfiler)
library(WGCNA)
library(DESeq2)
library(org.Hs.eg.db)
library(edgeR)
source("/shared/hidelab2/user/md4zsa/Work/Scripts/RUN-WGCNA/WGCNA_functions.R")
#Set active directory and read data
setwd("/shared/hidelab2/user/md4zsa/Work/Data/MSMM_RNAseq/MSMM_RNAseq_DEG_Analysis")
msmm_rnaseq_data=read.delim2("msmm_rnaseq_B1_raw_counts.txt",header = T,row.names = 1,as.is = T)

#Map Ensembl gene identifiers from raw data to known Entrez IDs. Mapping results in 23158 Entrez IDs from 
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
########################################################################################################################
# DEG Analysis
########################################################################################################################
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

#Build design matrix for DEG analysis to be used as input to DESeq2
msmm_rnaseq.colData1=msmm_rnaseq.DEG1=msmm_rnaseq.colData2=msmm_rnaseq.DEG2=vector(mode = "list",length=3)
names(msmm_rnaseq.colData1)=names(msmm_rnaseq.DEG1)=names(msmm_rnaseq.colData2)=names(msmm_rnaseq.DEG2)=c("BM_10","BM_22","BM_36")
for (i in 1:length(names(msmm_rnaseq.colData1))){
  x=matrix(NA,nrow=sum(length(low_Tangle_Samples[[i]]),length(high_Tangle_Samples[[i]])),ncol=2)
  rownames(x)=c(low_Tangle_Samples[[i]],high_Tangle_Samples[[i]])
  x[,1]=c(rep("low",length(low_Tangle_Samples[[i]])),rep("high",length(high_Tangle_Samples[[i]])))
  x[,2]=c(rep("paired-end",sum(length(low_Tangle_Samples[[i]]),length(high_Tangle_Samples[[i]]))))
  msmm_rnaseq.colData1[[i]]=data.frame(x,stringsAsFactors = F)
  #rownames(msmm_rnaseq.colData[[i]])=sort(c(low_Tangle_Samples[grep(pattern = names(msmm_rnaseq.colData)[i],low_Tangle_Samples)],high_Tangle_Samples[grep(pattern = names(msmm_rnaseq.colData)[i],high_Tangle_Samples)]))
  colnames(msmm_rnaseq.colData1[[i]])=c("condition","type")
  rm(x)
}
for(j in 1:length(names(msmm_rnaseq))){
  dds=DESeqDataSetFromMatrix(countData = msmm_rnaseq[[j]][,rownames(msmm_rnaseq.colData1[[j]])],colData = msmm_rnaseq.colData1[[j]],design=~condition)
  msmm_rnaseq.DEG1[[j]]=dds
}

msmm_rnaseq.colData2=msmm_rnaseq.DEG=vector(mode = "list",length=3)
names(msmm_rnaseq.colData2)=names(msmm_rnaseq.DEG2)=c("BM_10","BM_22","BM_36")
for (i in 1:length(names(msmm_rnaseq.colData2))){
  x=matrix(NA,nrow=sum(length(low_Plaque_Samples[[i]]),length(high_Plaque_Samples[[i]])),ncol=2)
  rownames(x)=c(low_Plaque_Samples[[i]],high_Plaque_Samples[[i]])
  x[,1]=c(rep("low",length(low_Plaque_Samples[[i]])),rep("high",length(high_Plaque_Samples[[i]])))
  x[,2]=c(rep("paired-end",sum(length(low_Plaque_Samples[[i]]),length(high_Plaque_Samples[[i]]))))
  msmm_rnaseq.colData2[[i]]=data.frame(x,stringsAsFactors = F)
  #rownames(msmm_rnaseq.colData[[i]])=sort(c(low_Tangle_Samples[grep(pattern = names(msmm_rnaseq.colData)[i],low_Tangle_Samples)],high_Tangle_Samples[grep(pattern = names(msmm_rnaseq.colData)[i],high_Tangle_Samples)]))
  colnames(msmm_rnaseq.colData2[[i]])=c("condition","type")
  rm(x)
}
for(j in 1:length(names(msmm_rnaseq))){
  dds=DESeqDataSetFromMatrix(countData = msmm_rnaseq[[j]][,rownames(msmm_rnaseq.colData2[[j]])],colData = msmm_rnaseq.colData2[[j]],design=~condition)
  msmm_rnaseq.DEG2[[j]]=dds
}
msmm_rnaseq.DEG1=lapply(msmm_rnaseq.DEG1,DESeq)
msmm_rnaseq.DEG2=lapply(msmm_rnaseq.DEG2,DESeq)
msmm_rnaseq.results1=lapply(msmm_rnaseq.DEG1,results)
msmm_rnaseq.results2=lapply(msmm_rnaseq.DEG2,results)

#Find cpm for each region/pathology using edgeR as DESeq2 
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

msmm_rnaseq.DEG_Tangle_Genes=lapply(lapply(lapply(msmm_rnaseq.results1,as.data.frame),function(x)x[which((x$padj<0.1)&(abs(x$log2FoldChange)>log2(1.3))),]),rownames)
msmm_rnaseq.DEG_PLQ_Genes=lapply(lapply(lapply(msmm_rnaseq.results2,as.data.frame),function(x)x[which((x$padj<0.1)&(abs(x$log2FoldChange)>log2(1.3))),]),rownames)
msmm_rnaseq.DEG_Tangle_GenesMap=msmm_rnaseq.DEG_Tangle_GenesMap2=msmm_rnaseq.DEG_PLQ_GenesMap=msmm_rnaseq.DEG_PLQ_GenesMap2=vector(mode = "list",length = length(msmm_rnaseq.DEG_PLQ_Genes))
names(msmm_rnaseq.DEG_Tangle_GenesMap)=names(msmm_rnaseq.DEG_Tangle_GenesMap2)=names(msmm_rnaseq.DEG_PLQ_GenesMap)=names(msmm_rnaseq.DEG_PLQ_GenesMap2)=names(msmm_rnaseq.results1)
for (k in 1:length(msmm_rnaseq.DEG_Tangle_Genes)){
  msmm_rnaseq.DEG_Tangle_GenesMap[[k]]=select(x = org.Hs.eg.db,keys = msmm_rnaseq.DEG_Tangle_Genes[[k]],keytype = "ENSEMBL",columns = c("ENTREZID","SYMBOL"))
}
indices1=lapply(lapply(msmm_rnaseq.DEG_Tangle_GenesMap,`[`,3),function(x)which(is.na(x)==F))
for (k in 1:length(msmm_rnaseq.DEG_Tangle_GenesMap)){
  msmm_rnaseq.DEG_Tangle_GenesMap2[[k]]=msmm_rnaseq.DEG_Tangle_GenesMap[[k]][indices1[[k]],]
}
for (l in 1:length(msmm_rnaseq.DEG_PLQ_Genes)){
  msmm_rnaseq.DEG_PLQ_GenesMap[[l]]=select(x = org.Hs.eg.db,keys = msmm_rnaseq.DEG_PLQ_Genes[[l]],keytype = "ENSEMBL",columns = c("ENTREZID","SYMBOL"))
}
indices2=lapply(lapply(msmm_rnaseq.DEG_PLQ_GenesMap,`[`,3),function(x)which(is.na(x)==F))
for (l in 1:length(msmm_rnaseq.DEG_PLQ_GenesMap)){
  msmm_rnaseq.DEG_PLQ_GenesMap2[[l]]=msmm_rnaseq.DEG_PLQ_GenesMap[[l]][indices2[[l]],]
}
#Enrichment analysis using KEGG Pathways
msmm_rnaseq.DEG_Tangle_KEGG=lapply(lapply(lapply(msmm_rnaseq.DEG_Tangle_GenesMap2,`[[`,2),enrichKEGG,organism = "hsa",pvalueCutoff = 0.05,pAdjustMethod = "BH"),summary)
msmm_rnaseq.DEG_PLQ_KEGG=lapply(lapply(lapply(msmm_rnaseq.DEG_PLQ_GenesMap2[2:3],`[[`,2),enrichKEGG,organism = "hsa",pvalueCutoff = 0.05,pAdjustMethod = "BH"),summary)
