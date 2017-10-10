library(limma)
library(Biobase)
library(jetset)
library(metaArray)
library(doMC)
library(pheatmap)
library(entropy)
library(clusterProfiler)
library(WGCNA)
library(DESeq2)

#Read data
setwd("/Users/sandeepamberkar/Work/Data/MSMM-MSBB-HBTRC-Data/MSMM_RNAseq/")
msbb_array19.files=list.files(pattern="*.tsv",full.names=T)
msbb_array19.files=msbb_array19.files[-c(18:19)]
msbb_array19=mapply(FUN=read.table,msbb_array19.files,MoreArgs=list(header=T,sep="\t",as.is=T))
names(msbb_array19)=c("AC","CN","DLPFC","FP","HP","IFG","ITG","MTG","OVC","PHG","PCC","PCG","PFC","PTMN","SPL","STG","TP")
#Resolve genes with multiple IDs and assign single,unique identifier,usually the first
multiID.u133a=msbb_array19[[1]]$ENTREZ_GENE_ID[grep(pattern="///",msbb_array19[[1]]$ENTREZ_GENE_ID)]
u133a_universe=c(msbb_array19[[1]]$ENTREZ_GENE_ID[-grep(pattern="///",msbb_array19[[1]]$ENTREZ_GENE_ID)],gsub(pattern=" ",replacement="",lapply(strsplit(x=multiID.u133a,split="///"),`[`,1)))
msbb_array19=lapply(msbb_array19,function(x){x$ENTREZ_GENE_ID<-u133a_universe;x})
msbb_array19.agg=lapply(lapply(msbb_array19,`[`,-c(1:4)),aggregate,by=list(EntrezID=msbb_array19[[1]]$ENTREZ_GENE_ID),mean)
msbb_array19.agg2=lapply(lapply(msbb_array19,`[`,-c(1:4)),aggregate,by=list(Gene.Symbol=msbb_array19[[1]]$Gene.Symbol),mean)
for (i in 1:length(names(msbb_array19.agg))){
  msbb_array19.agg[[i]]=msbb_array19.agg[[i]][-1,]
  rownames(msbb_array19.agg[[i]])=msbb_array19.agg[[i]][,1]
  msbb_array19.agg[[i]]=msbb_array19.agg[[i]][,-1]
}

# for (i in 1:length(names(msbb_array19.agg2))){
#   msbb_array19.agg2[[i]]=msbb_array19.agg2[[i]][-1,]
#   rownames(msbb_array19.agg2[[i]])=msbb_array19.agg2[[i]][,1]
#   msbb_array19.agg2[[i]]=msbb_array19.agg2[[i]][,-1]
# }
msbb_array19.covariates=read.delim2("../AMP-AD_MSBB_MSSM_covariates_mRNA_AffymetrixU133AB.tsv",header=T,as.is=T)
msbb_array19.covariates$Age[which(msbb_array19.covariates$Age=="89+")]=90 #For simplicity,all samples with AOD as 89+ were labeled 90
msbb_array19.covariates$Age=as.integer(msbb_array19.covariates$Age)
msbb_array19.covariates$pH=as.numeric(msbb_array19.covariates$pH)
msbb_array19.covariates$CDR=as.numeric(msbb_array19.covariates$CDR)
msbb_array19.covariates$PLQ_Mn=as.numeric(msbb_array19.covariates$PLQ_Mn)
msbb_array19.covariates$BrainBank=paste("X",msbb_array19.covariates$BrainBank,sep="")

#msbb_array19.2=lapply(msbb_array19,function(x){rownames(x)<-x$IDx})
#msbb_array19.2=lapply(msbb_array19,`[`,-c(1:4))
msbb_array19.NTr_PLQ_phenoData=msbb_array19.covariates[,c(1,6,10,12)]
msbb_array19.phenoDataIndices=lapply(lapply(msbb_array19.agg2,colnames),function(x)which(msbb_array19.NTr_PLQ_phenoData$BrainBank%in%x))
msbb_array19.NTr_PLQ_phenoData2=lapply(msbb_array19.phenoDataIndices,function(x)msbb_array19.NTr_PLQ_phenoData[x,])
msbb_array19.NTr_PLQ_phenoData2=lapply(msbb_array19.NTr_PLQ_phenoData2,function(x){rownames(x)<-x$BrainBank;x})
msbb_array19.NTr_PLQ_phenoData2=lapply(msbb_array19.NTr_PLQ_phenoData2,`[`,c(-1))

msbb_array19.eset=msbb_array19.DEG_PLQ=msbb_array19.DEG_NTR=list()
control.plq=lapply(msbb_array19.NTr_PLQ_phenoData2,function(x)which(x$Age>65 & x$PLQ_Mn<5))
disease.plq=lapply(msbb_array19.NTr_PLQ_phenoData2,function(x)which(x$Age>65 & x$PLQ_Mn>10))
control.ntr=lapply(msbb_array19.NTr_PLQ_phenoData2,function(x)which(x$Age>65 & x$NTrSum<5))
disease.ntr=lapply(msbb_array19.NTr_PLQ_phenoData2,function(x)which(x$Age>65 & x$NTrSum>10))
for (i in 1:length(names(msbb_array19.NTr_PLQ_phenoData2))){
  design.df=matrix(NA,nrow=sum(length(control.plq[[i]]),length(disease.plq[[i]])),ncol=2)
  rownames(design.df)=c(rownames(msbb_array19.NTr_PLQ_phenoData2[[i]][control.plq[[i]],]),rownames(msbb_array19.NTr_PLQ_phenoData2[[i]][disease.plq[[i]],]))
  colnames(design.df)=c("control.plq","disease.plq")
  design.df[which(rownames(design.df)%in%rownames(msbb_array19.NTr_PLQ_phenoData2[[i]][control.plq[[i]],])),1]=1
  design.df[-which(rownames(design.df)%in%rownames(msbb_array19.NTr_PLQ_phenoData2[[i]][control.plq[[i]],])),1]=0
  design.df[which(rownames(design.df)%in%rownames(msbb_array19.NTr_PLQ_phenoData2[[i]][disease.plq[[i]],])),2]=1
  design.df[-which(rownames(design.df)%in%rownames(msbb_array19.NTr_PLQ_phenoData2[[i]][disease.plq[[i]],])),2]=0
  phenoData=AnnotatedDataFrame(data=msbb_array19.NTr_PLQ_phenoData2[[i]][which(rownames(msbb_array19.NTr_PLQ_phenoData2[[i]])%in%rownames(design.df)),])
  msbb_array19.eset[[i]]=ExpressionSet(assayData=as.matrix(msbb_array19.agg[[i]][,which(colnames(msbb_array19.agg[[i]])%in%rownames(design.df))]),phenoData=phenoData)  
  fit=lmFit(msbb_array19.eset[[i]],design=design.df)
  fit=eBayes(fit)
  contMatrix=makeContrasts(CtrlvsDisease=disease.plq-control.plq,levels=design.df)
  fit2=contrasts.fit(fit,contMatrix)
  fit2=eBayes(fit2)
  msbb_array19.DEG_PLQ[[i]]=topTableF(fit2,adjust.method="BH",number=Inf)
}

for (j in 1:length(names(msbb_array19.NTr_PLQ_phenoData2))){
  design.df=matrix(NA,nrow=sum(length(control.ntr[[j]]),length(disease.ntr[[j]])),ncol=2)
  rownames(design.df)=c(rownames(msbb_array19.NTr_PLQ_phenoData2[[j]][control.ntr[[j]],]),rownames(msbb_array19.NTr_PLQ_phenoData2[[j]][disease.ntr[[j]],]))
  colnames(design.df)=c("control.ntr","disease.ntr")
  design.df[which(rownames(design.df)%in%rownames(msbb_array19.NTr_PLQ_phenoData2[[j]][control.ntr[[j]],])),1]=1
  design.df[-which(rownames(design.df)%in%rownames(msbb_array19.NTr_PLQ_phenoData2[[j]][control.ntr[[j]],])),1]=0
  design.df[which(rownames(design.df)%in%rownames(msbb_array19.NTr_PLQ_phenoData2[[j]][disease.ntr[[j]],])),2]=1
  design.df[-which(rownames(design.df)%in%rownames(msbb_array19.NTr_PLQ_phenoData2[[j]][disease.ntr[[j]],])),2]=0
  phenoData=AnnotatedDataFrame(data=msbb_array19.NTr_PLQ_phenoData2[[j]][which(rownames(msbb_array19.NTr_PLQ_phenoData2[[j]])%in%rownames(design.df)),])
  msbb_array19.eset[[j]]=ExpressionSet(assayData=as.matrix(msbb_array19.agg[[j]][,which(colnames(msbb_array19.agg[[j]])%in%rownames(design.df))]),phenoData=phenoData)  
  fit=lmFit(msbb_array19.eset[[j]],design=design.df)
  fit=eBayes(fit)
  contMatrix=makeContrasts(CtrlvsDisease=disease.ntr-control.ntr,levels=design.df)
  fit2=contrasts.fit(fit,contMatrix)
  fit2=eBayes(fit2)
  msbb_array19.DEG_NTR[[j]]=topTableF(fit2,adjust.method="BH",number=Inf)
}
names(msbb_array19.eset)=names(msbb_array19.DEG_PLQ)=names(msbb_array19)

#Check regions that have highest DEGs at FDR<=0.2,for plaques
unlist(lapply(lapply(msbb_array19.DEG_PLQ,function(x)which(x$adj.P.Val<=0.25)),length))
#Check regions that have highest DEGs at FDR<=0.2,for tangles
unlist(lapply(lapply(msbb_array19.DEG_NTR,function(x)which(x$adj.P.Val<=0.25)),length))

#Since Hippocampus (HP) has highest DEGs for both plaques and tangles,we further process it to determine functional enrichment
#We would also use HP for performing WGCNA in both control and samples
msbb_array19.DEG_PLQ_HPGenes=rownames(msbb_array19.DEG_PLQ[[5]][which(msbb_array19.DEG_PLQ[[5]]$adj.P.Val<=0.25),])
msbb_array19.DEG_NTR_HPGenes=rownames(msbb_array19.DEG_NTR[[5]][which(msbb_array19.DEG_NTR[[5]]$adj.P.Val<=0.25),])


multiID.DEG_PLQ_HPGenes=msbb_array19.DEG_PLQ_HPGenes[grep(pattern="///",msbb_array19.DEG_PLQ_HPGenes)]
multiID.DEG_NTR_HPGenes=msbb_array19.DEG_NTR_HPGenes[grep(pattern="///",msbb_array19.DEG_NTR_HPGenes)]
msbb_array19.DEG_NTR_HPGenes[grep(pattern="///",msbb_array19.DEG_NTR_HPGenes)]=gsub(pattern=" ",replacement="",lapply(strsplit(x=multiID.DEG_NTR_HPGenes,split="///"),`[`,1))
msbb_array19.DEG_PLQ_HPGenes[grep(pattern="///",msbb_array19.DEG_PLQ_HPGenes)]=gsub(pattern=" ",replacement="",lapply(strsplit(x=multiID.DEG_PLQ_HPGenes,split="///"),`[`,1))
#Repeat it for all the EntrezIDs of the original dataset

DEG_PLQ_HPGenes.KEGG=summary(enrichKEGG(gene=as.numeric(msbb_array19.DEG_PLQ_HPGenes),organism="hsa",pvalueCutoff=0.05))
DEG_NTR_HPGenes.KEGG=summary(enrichKEGG(gene=msbb_array19.DEG_NTR_HPGenes,organism="hsa",pvalueCutoff=0.05))
###########################################################################################################################
#Run WGCNA on HP,using control and disease samples
nSets=2
setLabels=c("Ctrl_NTr","Dis_NTr")
multiExpr.NTr=vector(mode="list",length=nSets)
multiExpr.NTr[[1]]=list(data=as.data.frame(t(msbb_array19.agg[[5]][,control.ntr[[5]]])))
multiExpr.NTr[[2]]=list(data=as.data.frame(t(msbb_array19.agg[[5]][,disease.ntr[[5]]])))
#rownames(multiExpr.NTr[[1]]$data)
# #Remove outliers using the data removal procedure highlighted in WCGNA workflow.
# #Accordingly,one sample from control.NTR (X616) and one sample from disease.NTR (X375) were removed
sampleTrees=list()
for (set in 1:nSets){sampleTrees[[set]]=hclust(dist(multiExpr.NTr[[set]]$data),method="average")}
par(mfrow=c(2,1))
par(mar=c(0,4,2,0))
for (set in 1:nSets){
  plot(sampleTrees[[set]],main=paste("Sample clustering on all genes in",setLabels[set]),xlab="",sub="",cex=0.7)
}
#   
multiExpr.NTr[[1]]$data=multiExpr.NTr[[1]]$data[-which(rownames(multiExpr.NTr[[1]]$data)=="X616"),]
multiExpr.NTr[[2]]$data=multiExpr.NTr[[2]]$data[-which(rownames(multiExpr.NTr[[2]]$data)=="X375"),]
powers = c(seq(4,10,by=1), seq(12,20, by=2))
powerTables = vector(mode = "list", length = nSets)
for (set in 1:nSets){
  powerTables[[set]]=list(data=pickSoftThreshold(multiExpr.NTr[[set]]$data[which(rownames(multiExpr.NTr[[set]]$data)%in%rownames(msbb_array19.DEG_PLQ_HP),],powerVector=powers,verbose=2,corFnc=corFast)[[2]])
}
collectGarbage()
colNames=c("Scale Free Topology Model Fit","Mean connectivity","Median connectivity","Max connectivity")
plotCols = c(2,5,6,7)
colors = c("black", "red")
ylim=matrix(NA,nrow=2,ncol=4)
for (set in 1:nSets){
  for (col in 1:length(plotCols)){
    ylim[1,col]=min(ylim[1,col],powerTables[[set]]$data[,plotCols[col]],na.rm=TRUE)
    ylim[2,col]=max(ylim[2,col],powerTables[[set]]$data[,plotCols[col]],na.rm=TRUE)
    }
}

sizeGrWindow(8,6)
par(mfcol=c(2,2))
par(mar=c(4.2,4.2 ,2.2,0.5))
cex1=0.7
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
  if (set==1)
  {
    plot(powerTables[[set]]$data[,1],-sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         xlab="Soft Threshold (power)",ylab=colNames[col],type="n",ylim=ylim[,col],
         main=colNames[col])
    addGrid()
  }
  if (col==1)
  {
    text(powerTables[[set]]$data[,1],-sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         labels=powers,cex=cex1,col=colors[set])
  } else
    text(powerTables[[set]]$data[,1],powerTables[[set]]$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[set])
  if (col==1)
  {
    legend("bottomleft",legend=setLabels,col=colors,pch=20)
  } else
    legend("topleft",legend=setLabels,col=colors,pch=20)
}
NTR_HP.bnet=blockwiseConsensusModules(multiExpr.NTr,maxBlockSize=2000,power=8,minModuleSize=30,deepSplit=2,pamRespectsDendro=FALSE,mergeCutHeight=0.25,numericLabels=TRUE,minKMEtoStay=0,saveTOMs=TRUE,verbose=5)
determineSoftPowerWGCNA(msbb_array19.agg[[5]][which(rownames(msbb_array19.agg[[5]])%in%rownames(msbb_array19.DEG_PLQ_HP)),disease.ntr[[5]]],"msbb_array19_NTR_DEG_HP_SoftPowers.png",propGenes=1)
NTR_HP_coexpp=runWGCNA(msbb_array19.agg[[5]][,disease.ntr[[5]]],1,softPower=7,signedNetwork=TRUE)
NTR_HP_MEs=calculateModuleEigengenes(NTR_HP_coexpp,1,30)
