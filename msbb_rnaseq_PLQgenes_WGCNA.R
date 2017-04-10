library(igraph)
library(org.Hs.eg.db)
source("https://raw.githubusercontent.com/SamBuckberry/RUN-WGCNA/master/WGCNA_functions.R")

msmm_rnaseq_PLQGenes_exprs=msbb_rnaseq_covariates.merged_final=lowPlaque_samples=highPlaque_samples=vector(mode='list',length=4)
names(msmm_rnaseq_PLQGenes_exprs)=names(msbb_rnaseq_covariates.merged_final)=names(lowPlaque_samples)=names(highPlaque_samples)=c("FP","IFG","PHG","STG") 
msbb_rnaseq_covariates.merged_final$FP=read.table("MSBB_RNAseq2016_FP_covariates.txt",sep = "\t",header = T,as.is = T)
msbb_rnaseq_covariates.merged_final$IFG=read.table("MSBB_RNAseq2016_IFG_covariates.txt",sep = "\t",header = T,as.is = T)
msbb_rnaseq_covariates.merged_final$PHG=read.table("MSBB_RNAseq2016_PHG_covariates.txt",sep = "\t",header = T,as.is = T)
msbb_rnaseq_covariates.merged_final$STG=read.table("MSBB_RNAseq2016_STG_covariates.txt",sep = "\t",header = T,as.is = T)

msmm_rnaseq_PLQGenes_exprs$FP=read.table('msmm_rnaseq_PLQGenes_FP_exprs.txt',sep='\t',header=T,row.names=1,as.is=T)
msmm_rnaseq_PLQGenes_exprs$IFG=read.table('msmm_rnaseq_PLQGenes_IFG_exprs.txt',sep='\t',header=T,row.names=1,as.is=T)
msmm_rnaseq_PLQGenes_exprs$PHG=read.table('msmm_rnaseq_PLQGenes_PHG_exprs.txt',sep='\t',header=T,row.names=1,as.is=T)
msmm_rnaseq_PLQGenes_exprs$STG=read.table('msmm_rnaseq_PLQGenes_STG_exprs.txt',sep='\t',header=T,row.names=1,as.is=T)

lowPlaque_samples=highPlaque_samples=msbb_rnaseq.colData=vector(mode = "list",length = 4)
names(lowPlaque_samples)=names(highPlaque_samples)=names(msbb_rnaseq.colData)=names(msmm_rnaseq_PLQGenes_exprs)
lowPlaque_samples=lapply(msbb_rnaseq_covariates.merged_final,function(x)x$sampleIdentifier[which(x$PlaqueMean<=1)])
highPlaque_samples=lapply(msbb_rnaseq_covariates.merged_final,function(x)x$sampleIdentifier[which(x$PlaqueMean>=15)])

exprs_rank=vector(mode = "list",length = 4)
names(exprs_rank)=names(msmm_rnaseq_PLQGenes_exprs)

for (j in 1:4){
  colnames(msmm_rnaseq_PLQGenes_exprs[[j]])=gsub(pattern='X',replacement='',colnames(msmm_rnaseq_PLQGenes_exprs[[j]]))
  exprs_rank[[j]]=msmm_rnaseq_PLQGenes_exprs[[j]]
  c_exprs_rank=exprs_rank[[j]][,which(colnames(exprs_rank[[j]])%in%unlist(lapply(strsplit(lowPlaque_samples[[j]],split="_"),`[[`,3)))]
  t_exprs_rank=exprs_rank[[j]][,which(colnames(exprs_rank[[j]])%in%unlist(lapply(strsplit(highPlaque_samples[[j]],split="_"),`[[`,3)))]
  cat(paste('Computing soft power in Low-plaque samples for brain region', names(msmm_rnaseq_PLQGenes_exprs)[j],'...\n',sep=' '))
  determineSoftPowerWGCNA(data1=c_exprs_rank, outFile=paste('PLQGenes_Low',names(msmm_rnaseq_PLQGenes_exprs)[j],'powerPlots.png',sep='_'), propGenes=1)
  cat(paste('Computing soft power in High-plaque samples for brain region', names(msmm_rnaseq_PLQGenes_exprs)[j],'...\n',sep=' '))
  determineSoftPowerWGCNA(data1=t_exprs_rank, outFile=paste('PLQGenes_High',names(msmm_rnaseq_PLQGenes_exprs)[j],'powerPlots.png',sep='_'), propGenes=1)
}
msmm_rnaseq_softPowers_LowPLQ=c(4,4,10,4)
msmm_rnaseq_softPowers_HighPLQ=c(7,5,7,5)
msmm_rnaseq_PLQgenes.CoexpNet_WGCNA=vector(mode = 'list',length = 4)
names(msmm_rnaseq_PLQgenes.CoexpNet_WGCNA)=names(msmm_rnaseq_PLQGenes_exprs)
msmm_rnaseq_PLQgenes.CoexpNet_WGCNA$FP=msmm_rnaseq_PLQgenes.CoexpNet_WGCNA$IFG=msmm_rnaseq_PLQgenes.CoexpNet_WGCNA$PHG=msmm_rnaseq_PLQgenes.CoexpNet_WGCNA$STG=vector(mode = 'list',length = 2)
names(msmm_rnaseq_PLQgenes.CoexpNet_WGCNA$FP)=names(msmm_rnaseq_PLQgenes.CoexpNet_WGCNA$IFG)=names(msmm_rnaseq_PLQgenes.CoexpNet_WGCNA$PHG)=names(msmm_rnaseq_PLQgenes.CoexpNet_WGCNA$STG)=c('Low','High')


for (j in 1:4){
  
  msmm_rnaseq_PLQgenes.CoexpNet_WGCNA[[j]][[1]]=runWGCNA(data1=msmm_rnaseq_PLQGenes_exprs[[j]], propGenes=1, softPower=msmm_rnaseq_softPowers_LowPLQ[j], signedNetwork=TRUE)
  msmm_rnaseq_PLQgenes.CoexpNet_WGCNA[[j]][[2]]=runWGCNA(data1=msmm_rnaseq_PLQGenes_exprs[[j]], propGenes=1, softPower=msmm_rnaseq_softPowers_HighPLQ[j], signedNetwork=TRUE)
  
}


#WGCNA analysis
msmm_rnaseq_PLQgenes_CoexpNet.WGCNA=readRDS('msmm_rnaseq_PLQgenes_CoexpNet_WGCNA.RDS')
msmm_rnaseq_LowPLQgenes_CoexpNet_WGCNA.graph=lapply(msmm_rnaseq_PLQgenes_CoexpNet.WGCNA,function(x)graph_from_adjacency_matrix(adjmatrix = x$Low$dissTOMA1,mode = 'undirected',weighted = T))
msmm_rnaseq_HighPLQgenes_CoexpNet_WGCNA.graph=lapply(msmm_rnaseq_PLQgenes_CoexpNet.WGCNA,function(x)graph_from_adjacency_matrix(adjmatrix = x$High$dissTOMA1,mode = 'undirected',weighted = T))

msmm_rnaseq_LowPLQgenes_CoexpNet_WGCNA.Betweenness=lapply(msmm_rnaseq_LowPLQgenes_CoexpNet_WGCNA.graph,function(x)betweenness(graph=x,weights=E(x)$weight))
