library(cocor)
library(data.table)
require(parallel)
library(org.Hs.eg.db)
library(foreach)

#Read data and covariates
setwd("/shared/hidelab2/user/md4zsa/Work/Data/MSMM_RNAseq/MSMM_RNAseq_FinalRelease2/")
msbb_rnaseq2016_data=read.table("AMP-AD_MSBB_MSSM_IlluminaHiSeq2500_normalized_counts_September_2016.txt",sep="\t",header = T,as.is = T)
msbb_rnaseq_covariates=read.csv("MSBB_RNAseq_covariates.csv",header = T,as.is = T)
msbb_rnaseq_clinical_covariates=read.csv("MSBB_clinical.csv",header = T,as.is = T)
colnames(msbb_rnaseq2016_data)=unlist(lapply(strsplit(x = colnames(msbb_rnaseq2016_data),split = "X"),`[[`,2))
msbb_ensembl_symbol=data.frame(Ensembl=names(mapIds(x = org.Hs.eg.db,keys=rownames(msbb_rnaseq2016_data),column = "SYMBOL",keytype = "ENSEMBL",multiVals = "first")),Symbol=mapIds(x = org.Hs.eg.db,keys=rownames(msbb_rnaseq2016_data),column = "SYMBOL",keytype = "ENSEMBL",multiVals = "first"),stringsAsFactors = F)
msbb_rnaseq_covariates.merged=merge(x = msbb_rnaseq_clinical_covariates,y=msbb_rnaseq_covariates,by=c("individualIdentifier","individualIdentifier"))
msbb_rnaseq_covariates.merged2=msbb_rnaseq_covariates.merged[grep(pattern = "unmapped|resequenced",x = msbb_rnaseq_covariates.merged$fileName,invert = T),]#include resequenced samples

#Initialise datastructre to store preprocessed data
msbb_rnaseq2016_byRegion=msbb_rnaseq_covariates.merged_final=msbb_rnaseq2016_PLQGenes=vector(mode = "list",length = 4)
names(msbb_rnaseq2016_byRegion)=names(msbb_rnaseq_covariates.merged_final)=names(msbb_rnaseq2016_PLQGenes)=c("FP","IFG","PHG","STG") 
msbb_rnaseq2016_data2=msbb_rnaseq2016_data[-which(rownames(msbb_rnaseq2016_data)%in%msbb_ensembl_symbol$Ensembl[which(is.na(msbb_ensembl_symbol$Symbol)==T)]),]
msbb_rnaseq2016_data2=data.frame(cbind(Symbol=msbb_ensembl_symbol$Symbol[-which(is.na(msbb_ensembl_symbol$Symbol)==T)],msbb_rnaseq2016_data2),stringsAsFactors = F)
#Compute mean expression for  gene symbol
msbb_rnaseq2016_data2.agg=aggregate(x=msbb_rnaseq2016_data2[,-1],by=list(geneSymbol=msbb_rnaseq2016_data2$Symbol),mean)
rownames(msbb_rnaseq2016_data2.agg)=msbb_rnaseq2016_data2.agg$geneSymbol
msbb_rnaseq2016_data2.agg=msbb_rnaseq2016_data2.agg[,-1]
colnames(msbb_rnaseq2016_data2.agg)=gsub(pattern = "X",replacement = "",x = colnames(msbb_rnaseq2016_data2.agg))

msbb_rnaseq_covariates.merged_final$FP=read.table("MSBB_RNAseq2016_FP_covariates.txt",sep = "\t",header = T,as.is = T)
msbb_rnaseq_covariates.merged_final$IFG=read.table("MSBB_RNAseq2016_IFG_covariates.txt",sep = "\t",header = T,as.is = T)
msbb_rnaseq_covariates.merged_final$PHG=read.table("MSBB_RNAseq2016_PHG_covariates.txt",sep = "\t",header = T,as.is = T)
msbb_rnaseq_covariates.merged_final$STG=read.table("MSBB_RNAseq2016_STG_covariates.txt",sep = "\t",header = T,as.is = T)

#Segregate preprocessed data into datastructure
msbb_rnaseq2016_byRegion$FP=msbb_rnaseq2016_data2.agg[,which(colnames(msbb_rnaseq2016_data2.agg)%in%unlist(lapply(strsplit(x = msbb_rnaseq_covariates.merged_final$FP$sampleIdentifier,split = "_"),`[[`,3)))]
msbb_rnaseq2016_byRegion$IFG=msbb_rnaseq2016_data2.agg[,which(colnames(msbb_rnaseq2016_data2.agg)%in%unlist(lapply(strsplit(x = msbb_rnaseq_covariates.merged_final$IFG$sampleIdentifier,split = "_"),`[[`,3)))]
msbb_rnaseq2016_byRegion$PHG=msbb_rnaseq2016_data2.agg[,which(colnames(msbb_rnaseq2016_data2.agg)%in%unlist(lapply(strsplit(x = msbb_rnaseq_covariates.merged_final$PHG$sampleIdentifier,split = "_"),`[[`,3)))]
msbb_rnaseq2016_byRegion$STG=msbb_rnaseq2016_data2.agg[,which(colnames(msbb_rnaseq2016_data2.agg)%in%unlist(lapply(strsplit(x = msbb_rnaseq_covariates.merged_final$STG$sampleIdentifier,split = "_"),`[[`,3)))]

#Separate samples into Low (Control) and High (AD) plaque samples
braak_12=lapply(msbb_rnaseq_covariates.merged_final,function(x)x$sampleIdentifier[which(x$bbscore>0&x$bbscore<=2&x$bbscore!='NA')])
braak_34=lapply(msbb_rnaseq_covariates.merged_final,function(x)x$sampleIdentifier[which(x$bbscore>2&x$bbscore<=4&x$bbscore!='NA')])
braak_56=lapply(msbb_rnaseq_covariates.merged_final,function(x)x$sampleIdentifier[which(x$bbscore>4&x$bbscore<=6&x$bbscore!='NA')])

nc = 15 # number of cores
blocksize=100000 #Gene pair blocks to compute in parallel

#Define function to compute diffcorr
ProcessElement <- function(ic){
  A = ceiling((sqrt(8*(ic+1)-7)+1)/2)
  B = ic-choose(floor(1/2+sqrt(2*ic)),2)
  
  c_A = c_exprs_rank[A,]
  c_B = c_exprs_rank[B,]
  
  t_A = t_exprs_rank[A,]
  t_B = t_exprs_rank[B,]
  
  # get correlation between the summaries for the unique genes
  tmp = data.frame(Gene.A=gene.names[A],Gene.B=gene.names[B])
  c_cortest<-cor.test(unname(unlist(c_A)), unname(unlist(c_B)), method="spearman")
  t_cortest<-cor.test(unname(unlist(t_A)), unname(unlist(t_B)), method="spearman")
  rc<-unname(cor.test(unlist(c_A),unlist(c_B))$estimate)
  rt<-unname(cor.test(unlist(t_A),unlist(t_B))$estimate)
  # diffcor<-cocor.indep.groups(rc, rt, n.c, n.t)
  tmp$r.c<-rc
  tmp$p.c<-c_cortest$p.value
  tmp$n.c<-n.c
  tmp$r.t<-rt
  tmp$p.t<-t_cortest$p.value
  tmp$n.t<-n.t
  tmp$p.cocor<-NA
  tmp$abs.corr.change<-abs(rt-rc)
  if ( (!is.na(tmp$r.c)) && (!is.na(tmp$r.t)) )
  {
    diffcor<-cocor.indep.groups(tmp$r.c, tmp$r.t, tmp$n.c, tmp$n.t)
    tmp$p.cocor<-diffcor@fisher1925$p.value
  }
  setTxtProgressBar(pb,ic %% blocksize)
  return(tmp)
}

#Keep genes with counts>0 and occuring in at least 1/3 samples
msbb_rnaseq2016_byRegion.final_keep=lapply(msbb_rnaseq2016_byRegion,function(x)rowSums(x>0)>=ncol(x)/3)

exprs_rank=vector(mode = "list",length = 4)
names(exprs_rank)=names(msbb_rnaseq2016_byRegion)
#Outer for loop to loop over tissue type
for(j in 1:4){
  
  exprs_rank[[j]]=msbb_rnaseq2016_byRegion[[j]][msbb_rnaseq2016_byRegion.final_keep[[j]],][1:50,]
  number_of_combinations<-choose(nrow(exprs_rank[[j]][1:50,]),2)
  c_exprs_rank=exprs_rank[[j]][1:50,which(colnames(exprs_rank[[j]])%in%unlist(lapply(strsplit(braak_12[[j]],split="_"),`[[`,3)))]
  t_exprs_rank=exprs_rank[[j]][1:50,which(colnames(exprs_rank[[j]])%in%unlist(lapply(strsplit(braak_34[[j]],split="_"),`[[`,3)))]
  # t2_exprs_rank=exprs_rank[[j]][1:50,which(colnames(exprs_rank[[j]])%in%unlist(lapply(strsplit(braak_56[[j]],split="_"),`[[`,3)))]
  n.c<-ncol(c_exprs_rank)
  n.t<-ncol(t_exprs_rank)
  # n.t2<-ncol(t2_exprs_rank)
  gene.names<-rownames(exprs_rank[[j]][1:50,])
  #Parallel diffcorr computation
  i<-0
  start<-i*blocksize+1
  end<-min((i+1)*blocksize, number_of_combinations)
  while( start < number_of_combinations)
  {
    input<-start:end
    pb = txtProgressBar(min=0,max=length(input),style=3,initial=0)
    cat("\n")
    res = mclapply(input,ProcessElement,mc.cores=nc)
    close(pb)	
    result <- rbindlist(res)
    result <- as.data.frame(result)
    result <- data.frame(result,stringsAsFactors = F)
    write.table(result, file=paste("/shared/hidelab2/user/md4zsa/Work/Data/MSMM_RNAseq/MSMM_RNAseq_FinalRelease2/",names(msbb_rnaseq2016_byRegion)[j],"_R2_",i,".txt",sep=""),sep="\t",col.names = T,row.names=FALSE, quote = FALSE)
    i<-i+1
    start<-i*blocksize+1
    end<-min((i+1)*blocksize, number_of_combinations)
  }
  
}
cat(paste("Done!\n"))
system.time()




