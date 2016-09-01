#In this script, we're trying to build a correlation matrix of genes from RNAseq datasets for 3 brain regions.
#The script computes pearson correlation using the corFast function from the WGCNA package that uses multithreading to compute the Pearson correlation
#coefficient of gene-gene pairs.

library(WGCNA)
library(NetComp)
library(doMC)
library(foreach)
registerDoMC(cores = 16)
#Initiate lists to store data
msmm_rnaseq=msmm_rnaseq.corr.threshold=list()
#Initialise multilevel lists
msmm_rnaseq.corr.threshold[[1]]=msmm_rnaseq.corr.threshold[[2]]=msmm_rnaseq.corr.threshold[[3]]=list()
setwd("/shared/hidelab2/user/md4zsa/Work/Data/MSMM_RNAseq/")
#Read data
msmm_rnaseq[[1]]=read.table("MSMM_RNAseq_FP_B1.tsv",sep="\t",header=T,as.is=T,row.names=1)
msmm_rnaseq[[2]]=read.table("MSMM_RNAseq_PHG_B1.tsv",sep="\t",header=T,as.is=T,row.names=1)
msmm_rnaseq[[3]]=read.table("MSMM_RNAseq_STG_B1.tsv",sep="\t",header=T,as.is=T,row.names=1)
msmm_rnaseq=lapply(msmm_rnaseq,t)
msmm_rnaseq.corr=lapply(msmm_rnaseq,corFast)
msmm_rnaseq.corr=lapply(msmm_rnaseq.corr,abs)
corr.threshold=c(0.5,0.9)
for (m in 2:3)
  msmm_rnaseq.corr.threshold[[m]] <- foreach(i=corr.threshold) %dopar% matrix_threshold(as.matrix(msmm_rnaseq.corr[[m]]),threshold=i,abs = T)
save.image(file="msmm_rnaseq_corr.RData")