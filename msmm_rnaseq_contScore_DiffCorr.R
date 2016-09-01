library(DiffCorr)

setwd("/shared/hidelab2/user/md4zsa/Work/Data/MSMM_RNAseq/MSMM_RNAseq_FinalRelease/")
cat(paste("Reading data ...\n"))
load("msmm_rnaseq_FinalDataset_contPLQ.RData")
cat(paste("Applying significance threshold of 0.05 ...\n"))

msmm_rnaseq.diffcorr=vector(mode = "list",length = 3)
names(msmm_rnaseq.diffcorr)=names(msmm_rnaseq)


for (t in 1:length(names(msmm_rnaseq))){
  
  iqr<-apply(msmm_rnaseq[[t]], 1, IQR,8)
  qt<-quantile(iqr, probs=.995)
  exprs_rank<-msmm_rnaseq[[t]][iqr>qt,]
  c_exprs_rank=exprs_rank[,which(colnames(msmm_rnaseq[[t]])%in%low_Plaque_Samples[[t]])]
  t_exprs_rank=exprs_rank[,which(colnames(msmm_rnaseq[[t]])%in%high_Plaque_Samples[[t]])]
  comp.2.cc.fdr(output.file = paste("res",names(msmm_rnaseq.diffcorr)[t],"diffcorr.txt",sep = "_"),data1 = c_exprs_rank,data2 = t_exprs_rank,method = "spearman",p.adjust.methods = "fdr",threshold = 0.05)
}
save.image("msmm_rnaseq_contScore_DiffCorr.RData")