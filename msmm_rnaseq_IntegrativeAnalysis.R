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
source("~/Work/Software/RUN-WGCNA/WGCNA_functions.R")
#Set active directory and read data
setwd("/Users/sandeepamberkar/Work/Data/MSMM-MSBB-HBTRC-Data/MSMM_RNAseq/")
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

for(f in 1:length(names(msmm_rnaseq.Plaque_counts))){
  nf=select(x = org.Hs.eg.db,keys=rownames(msmm_rnaseq.Tangle_counts[[f]]),keytype = "ENSEMBL",columns = "ENTREZID")
  write.table(data.frame(Name=nf[!duplicated(nf$ENSEMBL),2],Description=rep("NA",dim(msmm_rnaseq.Plaque_counts[[f]])[1]),msmm_rnaseq.Tangle_counts[[f]]),paste(names(msmm_rnaseq.Tangle_counts)[f],"_TangleCounts_entrez.txt",sep=""),sep = "\t",col.names = T,row.names = T,quote=F)
  #write.table(data.frame(Name=select(x = org.Hs.eg.db,keys=rownames(msmm_rnaseq.DEG_TangleExprs[[f]]),keytype = "ENSEMBL",columns = "SYMBOL")[,2],Description=rep("NA",dim(msmm_rnaseq.DEG_TangleExprs[[f]])[1]),msmm_rnaseq.DEG_TangleExprs[[f]]),paste(names(msmm_rnaseq.DEG_TangleExprs)[f],"DEG_TangleCounts_symbol.txt",sep=""),sep = "\t",col.names = T,row.names = T,quote=F)
}
for(f in 1:length(names(msmm_rnaseq.Plaque_counts))){
  nf=select(x = org.Hs.eg.db,keys=rownames(msmm_rnaseq.Plaque_counts[[f]]),keytype = "ENSEMBL",columns = "ENTREZID")
  write.table(data.frame(Name=nf[!duplicated(nf$ENSEMBL),2],Description=rep("NA",dim(msmm_rnaseq.Plaque_counts[[f]])[1]),msmm_rnaseq.Plaque_counts[[f]]),paste(names(msmm_rnaseq.Plaque_counts)[f],"_PlaqueCounts_entrez.txt",sep=""),sep = "\t",col.names = T,row.names = T,quote=F)
  #write.table(data.frame(Name=select(x = org.Hs.eg.db,keys=rownames(msmm_rnaseq.DEG_PlaqueExprs[[f]]),keytype = "ENSEMBL",columns = "SYMBOL",multiVals="asNA")[,2],Description=rep("NA",dim(msmm_rnaseq.DEG_PlaqueExprs[[f]])[1]),msmm_rnaseq.DEG_PlaqueExprs[[f]]),paste(names(msmm_rnaseq.DEG_PlaqueExprs)[f],"DEG_PlaqueCounts_symbol.txt",sep=""),sep = "\t",col.names = T,row.names = T,quote=F)
}
annotation_df=vector(mode = "list",length = 3)
annotation_df[[1]]=annotation_df[[2]]=annotation_df[[3]]=vector(mode = "list",length = 2)
names(annotation_df)=names(msmm_rnaseq)
names(annotation_df[[1]])=names(annotation_df[[2]])=names(annotation_df[[3]])=c("Tangles","Plaques")
for(i in 1:length(names(annotation_df))){
  
    annotation_df[[i]][[1]]=data.frame(Pathology=factor(c(rep("low",length(low_Tangle_Samples[[i]])),rep("high",length(high_Tangle_Samples[[i]]))),levels = c("low","high")))
    annotation_df[[i]][[2]]=data.frame(Pathology=factor(c(rep("low",length(low_Plaque_Samples[[i]])),rep("high",length(high_Plaque_Samples[[i]]))),levels = c("low","high")))
    rownames(annotation_df[[i]][[1]])=c(low_Tangle_Samples[[i]],high_Tangle_Samples[[i]])
    rownames(annotation_df[[i]][[2]])=c(low_Plaque_Samples[[i]],high_Plaque_Samples[[i]])
}


#For CoGa analysis - Generate required files
for (f in 1:length(names(msmm_rnaseq.DEG_TangleExprs))){
  write(paste(sum(length(low_Tangle_Samples[[f]]),length(high_Tangle_Samples[[f]])),2,1,sep=" "),file = paste("msmm_rnaseq_",names(msmm_rnaseq.DEG_TangleExprs)[f],"_Tangles_phenoLabels.cls",sep=""),append = T)
  write(paste("#","low","high",sep=" "),file = paste("msmm_rnaseq_",names(msmm_rnaseq.DEG_TangleExprs)[f],"_Tangles_phenoLabels.cls",sep=""),append = T)
  write(paste(c(rep(0,length(low_Tangle_Samples[[f]])),(rep(1,length(high_Tangle_Samples[[f]])))),collapse = " "),file = paste("msmm_rnaseq_",names(msmm_rnaseq.DEG_TangleExprs)[f],"_Tangles_phenoLabels.cls",sep=""),append = T)
}
for (f in 1:length(names(msmm_rnaseq.DEG_PlaqueExprs))){
  write(paste(sum(length(low_Plaque_Samples[[f]]),length(high_Plaque_Samples[[f]])),2,1,sep=" "),file = paste("msmm_rnaseq_",names(msmm_rnaseq.DEG_PlaqueExprs)[f],"_Plaques_phenoLabels.cls",sep=""),append = T)
  write(paste("#","low","high",sep=" "),file = paste("msmm_rnaseq_",names(msmm_rnaseq.DEG_PlaqueExprs)[f],"_Plaques_phenoLabels.cls",sep=""),append = T)
  write(paste(c(rep(0,length(low_Plaque_Samples[[f]])),(rep(1,length(high_Plaque_Samples[[f]])))),collapse = " "),file = paste("msmm_rnaseq_",names(msmm_rnaseq.DEG_PlaqueExprs)[f],"_Plaques_phenoLabels.cls",sep=""),append = T)
}
#Cogena analysis
# annoGMT="c2.cp.kegg.v5.0.symbols.gmt.xz"
# annofile=system.file("extdata", annoGMT, package="cogena")
# nClust=10
# sampleLabels_Tangles=vector(mode = "list",length = 3)
# sampleLabels_Plaques=vector(mode = "list",length = 3)
# names(sampleLabels_Tangles)=names(sampleLabels_Plaques)=names(msmm_rnaseq)
# for (s in 1:length(msmm_rnaseq)){
#   sampleLabels_Tangles[[s]]=c(rep("low",length(low_Tangle_Samples[[s]])),rep("high",length(high_Tangle_Samples[[s]])))
#   names(sampleLabels_Tangles[[s]])=c(low_Tangle_Samples[[s]],high_Tangle_Samples[[s]])
#   sampleLabels_Plaques[[s]]=c(rep("low",length(low_Plaque_Samples[[s]])),rep("high",length(high_Plaque_Samples[[s]])))
#   names(sampleLabels_Plaques[[s]])=c(low_Plaque_Samples[[s]],high_Plaque_Samples[[s]])
#   sampleLabels_Tangles[[s]]=factor(sampleLabels_Tangles[[s]], levels=c("low","high"))
#   sampleLabels_Plaques[[s]]=factor(sampleLabels_Plaques[[s]], levels=c("low","high"))
#   
# }
# ncore=2
# clMethods=c("hierarchical","pam")
# metric="correlation"
# method="complete"
# 
 msmm_rnaseq.DEG_TangleExprs=msmm_rnaseq.DEG_PlaqueExprs=msmm_rnaseq.DEG_TangleExprs2=msmm_rnaseq.DEG_PlaqueExprs2=vector(mode = "list",length = 3)
# 
 for (i in 1:length(names(msmm_rnaseq))){
  msmm_rnaseq.DEG_TangleExprs[[i]]=msmm_rnaseq.Tangle_counts[[i]][which(rownames(msmm_rnaseq.Tangle_counts[[i]])%in%msmm_rnaseq.DEG_Tangle_GenesMap2[[i]]$ENSEMBL),]
  msmm_rnaseq.DEG_PlaqueExprs[[i]]=msmm_rnaseq.Plaque_counts[[i]][which(rownames(msmm_rnaseq.Plaque_counts[[i]])%in%msmm_rnaseq.DEG_PLQ_GenesMap2[[i]]$ENSEMBL),]
}
# genecl_result_DEG_Tangle=genecl_result_DEG_Plaques=vector(mode = "list",length = 3)
# genecl_result_DEG_Plaques=lapply(msmm_rnaseq.DEG_PlaqueExprs[2:3],coExp,nClust,clMethods,metric,method,ncore)
# genecl_result_DEG_Tangle=lapply(msmm_rnaseq.DEG_TangleExprs,coExp,nClust,clMethods,metric,method,ncore)
# 
# clen_res_DEG_Tangles=vector(mode = "list",length = 3)
# clen_res_DEG_Plaques=vector(mode = "list",length = 2)
# names(clen_res_DEG_Plaques)=names(msmm_rnaseq)[2:3]
# names(clen_res_DEG_Tangles)=names(msmm_rnaseq)
# for(j in 1:2){
#   clen_res_DEG_Plaques[[j]]=clEnrich(genecl_result_DEG_Plaques[[j]],annofile = annofile,sampleLabel = sampleLabels_Plaques[[j]])
# }
# for(j in 1:length(clen_res_DEG_Tangles)){
#   clen_res_DEG_Tangles[[j]]=clEnrich(genecl_result_DEG_Tangle[[j]],annofile = annofile,sampleLabel = sampleLabels_Tangles[[j]])
# }
# clen_res_DEG_Plaques=lapply(genecl_result_DEG_Plaques[2:3],clEnrich,annofile,lapply(sampleLabels_Plaques[2:3],function(x)x))
# clen_res_DEG_Tangles=lapply(genecl_result_DEG_Tangle,clEnrich,annofile,lapply(sampleLabels_Tangles,function(x)x))

#WGCNA analysis
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

determineSoftPowerWGCNA(data = msmm_rnaseq.TangleExprs2$BM_10$low,"msmm_rnaseq_Tangle_FP_low",0.75)
determineSoftPowerWGCNA(data = msmm_rnaseq.TangleExprs2$BM_10$high,"msmm_rnaseq_Tangle_FP_high",0.75)
determineSoftPowerWGCNA(data = msmm_rnaseq.TangleExprs2$BM_22$low,"msmm_rnaseq_Tangle_STG_low",0.75)
determineSoftPowerWGCNA(data = msmm_rnaseq.TangleExprs2$BM_22$high,"msmm_rnaseq_Tangle_STG_high",0.75)
determineSoftPowerWGCNA(data = msmm_rnaseq.TangleExprs2$BM_36$low,"msmm_rnaseq_Tangle_PHG_low",0.75)
determineSoftPowerWGCNA(data = msmm_rnaseq.TangleExprs2$BM_36$high,"msmm_rnaseq_Tangle_PHG_high",0.75)

# determineSoftPowerWGCNA(data = msmm_rnaseq.PlaqueExprs2$BM_10$low,"msmm_rnaseq_PLQ_FP_low",1)
# determineSoftPowerWGCNA(data = msmm_rnaseq.PlaqueExprs2$BM_10$high,"msmm_rnaseq_PLQ_FP_high",1) Too few genes
determineSoftPowerWGCNA(data = msmm_rnaseq.PlaqueExprs2$BM_22$low,"msmm_rnaseq_PLQ_STG_low",0.75)
determineSoftPowerWGCNA(data = msmm_rnaseq.PlaqueExprs2$BM_22$high,"msmm_rnaseq_PLQ_STG_high",0.75)
determineSoftPowerWGCNA(data = msmm_rnaseq.PlaqueExprs2$BM_36$low,"msmm_rnaseq_PLQ_PHG_low",0.75)
determineSoftPowerWGCNA(data = msmm_rnaseq.PlaqueExprs2$BM_36$high,"msmm_rnaseq_PLQ_PHG_high",0.75)

msmm_rnaseq.TangleExprs.sft$BM_10$low=6 #8 for data including non-Entrez
msmm_rnaseq.TangleExprs.sft$BM_10$high=10 #16
msmm_rnaseq.TangleExprs.sft$BM_22$low=4 #5
msmm_rnaseq.TangleExprs.sft$BM_22$high=9 #14
msmm_rnaseq.TangleExprs.sft$BM_36$low=6 #8
msmm_rnaseq.TangleExprs.sft$BM_36$high=5 #12

msmm_rnaseq.PlaqueExprs.sft$BM_22$low=4 #3
msmm_rnaseq.PlaqueExprs.sft$BM_22$high=9 #14
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
