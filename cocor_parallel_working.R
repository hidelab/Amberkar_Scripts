library(cocor)
library(doMC)
library(foreach)
library(data.table)
#Set number of cores
nc=16

#This RData object has to be preloaded else the code will crash
load("msmm_rnaseq_FinalDataset_contPLQ.RData")
#Set number of cores
nc=16

#Set up the work directory
msmm_rnaseq_data=read.table("AMP-AD_MSBB_MSSM_IlluminaHiSeq2500_final.TMM_normalized.race_age_RIN_PMI_batch_corrected.tsv",sep = "\t",header = T,as.is = T)
msmm_rnaseq_raw_data=read.table("AMP-AD_MSBB_MSSM_IlluminaHiSeq2500_raw_counts_final.tsv",sep = "\t",header = T,as.is = T)
msmm_rnaseq_clinical=read.table("MSBB_clinical.csv",header = T,sep = ",",as.is = T)
msmm_rnaseq_metaSample=read.table("AMP-AD_MSBB_MSSM_IlluminaHiSeq2500.sample.meta.tsv",header = T,sep = "\t",as.is = T)
msmm_rnaseq_covariates=read.table("AMP-AD_MSBB_MSSM_meta.traits.tsv",header = T,sep = "\t",as.is = T)
msmm_rnaseq_covariates2=msmm_rnaseq_covariates[which(msmm_rnaseq_covariates$Sample.ID%in%msmm_rnaseq_metaSample$LibID[grep(pattern = "^BM_",msmm_rnaseq_metaSample$LibID)]),]
msmm_rnaseq_covariates2$AOD[which(msmm_rnaseq_covariates2$AOD=="90+")]=90
msmm_rnaseq_covariates2=msmm_rnaseq_covariates2[-which(is.na(msmm_rnaseq_covariates2$bbscore)==T),]

msmm_rnaseq=msmm_rnaseq_alog=msmm_rnaseq_raw=vector(mode = "list",length = 3)
names(msmm_rnaseq)=names(msmm_rnaseq_alog)=names(msmm_rnaseq_raw)=c("BM_10","BM_22","BM_36")
for (l in 1:length(names(msmm_rnaseq))){
  msmm_rnaseq[[l]]=msmm_rnaseq_data2.agg[,c(1,grep(pattern = names(msmm_rnaseq)[l],colnames(msmm_rnaseq_data2.agg)))]  
  rownames(msmm_rnaseq[[l]])=msmm_rnaseq[[l]]$geneSymbol
  msmm_rnaseq[[l]]=msmm_rnaseq[[l]][,-1]
}

cat(paste("Grouping samples by brain region ...\n"))
low_Plaque_Samples=high_Plaque_Samples=vector(mode="list",length = 3)
names(low_Plaque_Samples)=names(high_Plaque_Samples)=names(msmm_rnaseq)

low_Plaque_Samples[[1]]=grep(pattern = "^BM_10",msmm_rnaseq_covariates2$Sample.ID[which(msmm_rnaseq_covariates2$PlaqueMean<=1)],value = T)
low_Plaque_Samples[[2]]=grep(pattern = "^BM_22",msmm_rnaseq_covariates2$Sample.ID[which(msmm_rnaseq_covariates2$PlaqueMean<=1)],value = T)
low_Plaque_Samples[[3]]=grep(pattern = "^BM_36",msmm_rnaseq_covariates2$Sample.ID[which(msmm_rnaseq_covariates2$PlaqueMean<=1)],value = T)

high_Plaque_Samples[[1]]=grep(pattern = "^BM_10",msmm_rnaseq_covariates2$Sample.ID[which(msmm_rnaseq_covariates2$PlaqueMean>=15)],value = T)
high_Plaque_Samples[[2]]=grep(pattern = "^BM_22",msmm_rnaseq_covariates2$Sample.ID[which(msmm_rnaseq_covariates2$PlaqueMean>=15)],value = T)
high_Plaque_Samples[[3]]=grep(pattern = "^BM_36",msmm_rnaseq_covariates2$Sample.ID[which(msmm_rnaseq_covariates2$PlaqueMean>=15)],value = T)

msmm_rnaseq.cocor=vector(mode = "list",length = 3)
msmm_rnaseq.cocor[[1]]=msmm_rnaseq.cocor[[2]]=msmm_rnaseq.cocor[[3]]=vector(mode = "list",length = 2)
names(msmm_rnaseq.cocor)=names(msmm_rnaseq)
names(msmm_rnaseq.cocor[[1]])=names(msmm_rnaseq.cocor[[2]])=names(msmm_rnaseq.cocor[[3]])

#Function to compute diff correlation between gene pairs
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
  rc<-c_cortest$estimate
  rt<-t_cortest$estimate
  diffcor<-cocor.indep.groups(rc, rt, n.c, n.t)
  tmp$r.c<-rc
  tmp$p.c<-c_cortest$p.value
  tmp$n.c<-n.c
  tmp$r.t<-rt
  tmp$p.t<-t_cortest$p.value
  tmp$n.t<-n.t
  tmp$p.cocor<-diffcor@fisher1925$p.value
  tmp$abs.corr.change<-abs(rt-rc)
  
  setTxtProgressBar(pb,ic)
  return(tmp)
}

#This is the differing point with the test script wherein a specific subset of genes/brain region is used for diff. correlation computation that follows
msmm_rnaseq_Plaque.Corr005.Genes=lapply(msmm_rnaseq_Plaque.Corr,function(x)x[which(x$Pval<=0.05),1])
msmm_rnaseq_Plaque.Corr001.Genes=lapply(msmm_rnaseq_Plaque.Corr,function(x)x[which(x$Pval<=0.01),1])
for (t in 1:length(names(msmm_rnaseq))){
  
  exprs_rank<-msmm_rnaseq[[t]][which(rownames(msmm_rnaseq[[t]])%in%msmm_rnaseq_Plaque.Corr005.Genes[[t]]),]
  number_of_combinations<-choose(nrow(exprs_rank),2)
  c_exprs_rank=exprs_rank[,which(colnames(exprs_rank)%in%low_Plaque_Samples[[t]])]
  t_exprs_rank=exprs_rank[,which(colnames(exprs_rank)%in%high_Plaque_Samples[[t]])]
  n.c<-ncol(c_exprs_rank)
  n.t<-ncol(t_exprs_rank)
  gene.names<-rownames(exprs_rank)
  
  input = 1:number_of_combinations
  pb = txtProgressBar(min=0,max=number_of_combinations,style=3,initial=0)
  cat("\n")
  res = mclapply(input,ProcessElement,mc.cores=nc)
  close(pb)
  
  result <- rbindlist(res)
  result <- data.frame(result,stringsAsFactors = F)
  result$FDR<-p.adjust(result$p.cocor, method="fdr")
  result[,c(3:8,10:11)]=round(result[,c(3:8,10:11)],digits = 3)
  msmm_rnaseq.cocor[[t]]=result
}
