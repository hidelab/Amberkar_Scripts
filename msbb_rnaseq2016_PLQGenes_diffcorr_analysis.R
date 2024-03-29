library(cocor)
library(data.table)
require(parallel)
library(org.Hs.eg.db)

setwd("/shared/hidelab2/user/md4zsa/Work/Data/MSMM_RNAseq/MSMM_RNAseq_FinalRelease2/")
#Read data and preprocess - remove unmapped/reseq samples
msbb_rnaseq2016_data=read.table("AMP-AD_MSBB_MSSM_IlluminaHiSeq2500_normalized_counts_September_2016.txt",sep="\t",header = T,as.is = T)
msbb_rnaseq_covariates=read.csv("MSBB_RNAseq_covariates.csv",header = T,as.is = T)
msbb_rnaseq_clinical_covariates=read.csv("MSBB_clinical.csv",header = T,as.is = T)
colnames(msbb_rnaseq2016_data)=unlist(lapply(strsplit(x = colnames(msbb_rnaseq2016_data),split = "X"),`[[`,2))
msbb_ensembl_symbol=data.frame(Ensembl=names(mapIds(x = org.Hs.eg.db,keys=rownames(msbb_rnaseq2016_data),column = "SYMBOL",keytype = "ENSEMBL",multiVals = "first")),Symbol=mapIds(x = org.Hs.eg.db,keys=rownames(msbb_rnaseq2016_data),column = "SYMBOL",keytype = "ENSEMBL",multiVals = "first"),stringsAsFactors = F)
msbb_rnaseq_covariates.merged=merge(x = msbb_rnaseq_clinical_covariates,y=msbb_rnaseq_covariates,by=c("individualIdentifier","individualIdentifier"))
msbb_rnaseq_covariates.merged2=msbb_rnaseq_covariates.merged[grep(pattern = "unmapped|resequenced",x = msbb_rnaseq_covariates.merged$fileName,invert = T),]

msbb_rnaseq2016_byRegion=msbb_rnaseq_covariates.merged_final=msbb_rnaseq2016_PLQGenes=vector(mode = "list",length = 4)
names(msbb_rnaseq2016_byRegion)=names(msbb_rnaseq_covariates.merged_final)=names(msbb_rnaseq2016_PLQGenes)=c("FP","IFG","PHG","STG") 
msbb_rnaseq2016_data2=msbb_rnaseq2016_data[-which(rownames(msbb_rnaseq2016_data)%in%msbb_ensembl_symbol$Ensembl[which(is.na(msbb_ensembl_symbol$Symbol)==T)]),]
msbb_rnaseq2016_data2=data.frame(cbind(Symbol=msbb_ensembl_symbol$Symbol[-which(is.na(msbb_ensembl_symbol$Symbol)==T)],msbb_rnaseq2016_data2),stringsAsFactors = F)
msbb_rnaseq2016_data2.agg=aggregate(x=msbb_rnaseq2016_data2[,-1],by=list(geneSymbol=msbb_rnaseq2016_data2$Symbol),mean)
#rm_genes=grep(pattern = "_|\\.|^RP|-",msbb_rnaseq2016_data2.agg$geneSymbol,invert = F)
rownames(msbb_rnaseq2016_data2.agg)=msbb_rnaseq2016_data2.agg$geneSymbol
msbb_rnaseq2016_data2.agg=msbb_rnaseq2016_data2.agg[,-1]
colnames(msbb_rnaseq2016_data2.agg)=gsub(pattern = "X",replacement = "",x = colnames(msbb_rnaseq2016_data2.agg))

#Read region-wise covariates independently from text files to avoid whitespace issues
msbb_rnaseq_covariates.merged_final$FP=read.table("MSBB_RNAseq2016_FP_covariates.txt",sep = "\t",header = T,as.is = T)
msbb_rnaseq_covariates.merged_final$IFG=read.table("MSBB_RNAseq2016_IFG_covariates.txt",sep = "\t",header = T,as.is = T)
msbb_rnaseq_covariates.merged_final$PHG=read.table("MSBB_RNAseq2016_PHG_covariates.txt",sep = "\t",header = T,as.is = T)
msbb_rnaseq_covariates.merged_final$STG=read.table("MSBB_RNAseq2016_STG_covariates.txt",sep = "\t",header = T,as.is = T)

#Read preprocessed gene counts by region
msbb_rnaseq2016_byRegion$FP=msbb_rnaseq2016_data2.agg[,which(colnames(msbb_rnaseq2016_data2.agg)%in%unlist(lapply(strsplit(x = msbb_rnaseq_covariates.merged2$sampleIdentifier[which(msbb_rnaseq_covariates.merged2$BrodmannArea=="BM10")],split = "_"),`[[`,3)))]
msbb_rnaseq2016_byRegion$IFG=msbb_rnaseq2016_data2.agg[,which(colnames(msbb_rnaseq2016_data2.agg)%in%unlist(lapply(strsplit(x = msbb_rnaseq_covariates.merged2$sampleIdentifier[which(msbb_rnaseq_covariates.merged2$BrodmannArea=="BM44")],split = "_"),`[[`,3)))]
msbb_rnaseq2016_byRegion$PHG=msbb_rnaseq2016_data2.agg[,which(colnames(msbb_rnaseq2016_data2.agg)%in%unlist(lapply(strsplit(x = msbb_rnaseq_covariates.merged2$sampleIdentifier[which(msbb_rnaseq_covariates.merged2$BrodmannArea=="BM36")],split = "_"),`[[`,3)))]
msbb_rnaseq2016_byRegion$STG=msbb_rnaseq2016_data2.agg[,which(colnames(msbb_rnaseq2016_data2.agg)%in%unlist(lapply(strsplit(x = msbb_rnaseq_covariates.merged2$sampleIdentifier[which(msbb_rnaseq_covariates.merged2$BrodmannArea=="BM22")],split = "_"),`[[`,3)))]

#Define Low and High PLQ samples i.e. Control and Disease samples
lowPlaque_samples=highPlaque_samples=msbb_rnaseq.colData=vector(mode = "list",length = 4)
names(lowPlaque_samples)=names(highPlaque_samples)=names(msbb_rnaseq.colData)=names(msbb_rnaseq2016_byRegion)
lowPlaque_samples=lapply(msbb_rnaseq_covariates.merged_final,function(x)x$sampleIdentifier[which(x$PlaqueMean<=1)])
highPlaque_samples=lapply(msbb_rnaseq_covariates.merged_final,function(x)x$sampleIdentifier[which(x$PlaqueMean>=15)])

#Read PLQ associated genes
msbb_rnaseq2016_PLQGenes=vector(mode = "list",length = 4)
names(msbb_rnaseq2016_PLQGenes)=names(msbb_rnaseq2016_byRegion)
msbb_rnaseq2016_PLQGenes$FP=read.table("PLQ_Assoc_Genes/FP_allPLQ_AssocGenes.txt",sep = "\t",header = T,as.is = T)
msbb_rnaseq2016_PLQGenes$IFG=read.table("PLQ_Assoc_Genes/IFG_allPLQ_AssocGenes.txt",sep = "\t",header = T,as.is = T)
msbb_rnaseq2016_PLQGenes$PHG=read.table("PLQ_Assoc_Genes/PHG_allPLQ_AssocGenes.txt",sep = "\t",header = T,as.is = T)
msbb_rnaseq2016_PLQGenes$STG=read.table("PLQ_Assoc_Genes/STG_allPLQ_AssocGenes.txt",sep = "\t",header = T,as.is = T)
#Compute FDR for all regions
msbb_rnaseq2016_PLQGenes_Rho.padj=lapply(msbb_rnaseq2016_PLQGenes,function(x)p.adjust(x$Rho.p,method="fdr"))
msbb_rnaseq2016_PLQGenes$FP$Rho.padj=msbb_rnaseq2016_PLQGenes_Rho.padj$FP
msbb_rnaseq2016_PLQGenes$IFG$Rho.padj=msbb_rnaseq2016_PLQGenes_Rho.padj$IFG
msbb_rnaseq2016_PLQGenes$PHG$Rho.padj=msbb_rnaseq2016_PLQGenes_Rho.padj$PHG
msbb_rnaseq2016_PLQGenes$STG$Rho.padj=msbb_rnaseq2016_PLQGenes_Rho.padj$STG

#Apply FDR filter of 5%
msbb_rnaseq2016_PLQGenes2=lapply(msbb_rnaseq2016_PLQGenes,function(x)x[which(x$Rho.padj<=0.05),])

#Setup parallel cocor analysis
nc = 8
blocksize=100000
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

msbb_rnaseq2016_byRegion.final_keep=lapply(msbb_rnaseq2016_byRegion,function(x)x[rowSums(x>0)>=ncol(x)/3,])

exprs_rank=vector(mode = "list",length = 4)
names(exprs_rank)=names(msbb_rnaseq2016_byRegion)
for (j in 1:4){
  
  exprs_rank[[j]]=msbb_rnaseq2016_byRegion.final_keep[[j]]
  exprs_rank[[j]]=exprs_rank[[j]][which(rownames(exprs_rank[[j]])%in%msbb_rnaseq2016_PLQGenes2[[j]]$Genes),]
  number_of_combinations<-choose(nrow(exprs_rank[[j]]),2)
  c_exprs_rank=exprs_rank[[j]][,which(colnames(exprs_rank[[j]])%in%unlist(lapply(strsplit(lowPlaque_samples[[j]],split="_"),`[[`,3)))]
  t_exprs_rank=exprs_rank[[j]][,which(colnames(exprs_rank[[j]])%in%unlist(lapply(strsplit(highPlaque_samples[[j]],split="_"),`[[`,3)))]
  n.c<-ncol(c_exprs_rank)
  n.t<-ncol(t_exprs_rank)
  gene.names<-rownames(exprs_rank[[j]])
  
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
    write.table(result, file=paste("/shared/hidelab2/user/md4zsa/Work/Data/MSMM_RNAseq/MSMM_RNAseq_FinalRelease2/",names(msbb_rnaseq2016_byRegion)[j],"_PLQGenes_DiffCorr_",i,".txt",sep=""),sep="\t",col.names = T,row.names=FALSE, quote = FALSE)
    i<-i+1
    start<-i*blocksize+1
    end<-min((i+1)*blocksize, number_of_combinations)
  }
  
}
cat(paste("Done!\n"))

