library(org.Hs.eg.db)
library(cocor)
library(data.table)
library(parallel)

setwd("/shared/hidelab2/user/md4zsa/Work/Data/MSMM_RNAseq/MSMM_RNAseq_FinalRelease/")
cat(paste("Reading data ...\n"))

msmm_rnaseq_data=read.table("AMP-AD_MSBB_MSSM_IlluminaHiSeq2500_final.TMM_normalized.race_age_RIN_PMI_batch_corrected.tsv",sep = "\t",header = T,as.is = T)
msmm_rnaseq_raw_data=read.table("AMP-AD_MSBB_MSSM_IlluminaHiSeq2500_raw_counts_final.tsv",sep = "\t",header = T,as.is = T)
msmm_rnaseq_clinical=read.table("MSBB_clinical.csv",header = T,sep = ",",as.is = T)
msmm_rnaseq_metaSample=read.table("AMP-AD_MSBB_MSSM_IlluminaHiSeq2500.sample.meta.tsv",header = T,sep = "\t",as.is = T)
msmm_rnaseq_covariates=read.table("AMP-AD_MSBB_MSSM_meta.traits.tsv",header = T,sep = "\t",as.is = T)
msmm_rnaseq_covariates2=msmm_rnaseq_covariates[which(msmm_rnaseq_covariates$Sample.ID%in%msmm_rnaseq_metaSample$LibID[grep(pattern = "^BM_",msmm_rnaseq_metaSample$LibID)]),]
msmm_rnaseq_covariates2$AOD[which(msmm_rnaseq_covariates2$AOD=="90+")]=90
msmm_rnaseq_covariates2=msmm_rnaseq_covariates2[-which(is.na(msmm_rnaseq_covariates2$bbscore)==T),]

cat(paste("Averaging mean expression for multiple Ensembl probeIDs to Gene Symbols ...\n"))
msmm_rnaseq_data2=msmm_rnaseq_data[,c(1,7,which(colnames(msmm_rnaseq_data)%in%msmm_rnaseq_metaSample$NewBarcode[grep(pattern = "^BM_",msmm_rnaseq_metaSample$LibID)]))]
msmm_rnaseq_raw_data2=msmm_rnaseq_raw_data[,c(1,7,which(colnames(msmm_rnaseq_data)%in%msmm_rnaseq_metaSample$NewBarcode[grep(pattern = "^BM_",msmm_rnaseq_metaSample$LibID)]))]
colnames(msmm_rnaseq_data2)[-c(1:2)]=msmm_rnaseq_metaSample$LibID[1:469]
#colnames(msmm_rnaseq_raw_data2)[-c(1:2)]=msmm_rnaseq_metaSample$LibID[1:469]
rm_genes=grep(pattern = "_|\\.|^RP|-",msmm_rnaseq_data2$geneSymbol)
msmm_rnaseq_data2.agg=aggregate(x = msmm_rnaseq_data2[-rm_genes,-c(1:2)],by=list(geneSymbol=msmm_rnaseq_data2$geneSymbol[-rm_genes]),mean)
#msmm_rnaseq_raw_data2.agg=aggregate(x = msmm_rnaseq_raw_data2[-rm_genes,-c(1:2)],by=list(geneSymbol=msmm_rnaseq_raw_data2$geneSymbol),mean)

msmm_rnaseq=vector(mode = "list",length = 3)
names(msmm_rnaseq)=c("BM_10","BM_22","BM_36")
for (l in 1:length(names(msmm_rnaseq))){
  msmm_rnaseq[[l]]=msmm_rnaseq_data2.agg[,c(1,grep(pattern = names(msmm_rnaseq)[l],colnames(msmm_rnaseq_data2.agg)))]  
  rownames(msmm_rnaseq[[l]])=msmm_rnaseq[[l]]$geneSymbol
  msmm_rnaseq[[l]]=msmm_rnaseq[[l]][,-1]
}
nc = detectCores()
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
  c_cortest<-cor.test(unname(unlist(c_A)), unname(unlist(c_B)), method="pearson")
  t_cortest<-cor.test(unname(unlist(t_A)), unname(unlist(t_B)), method="pearson")
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
low_Plaque_Samples=high_Plaque_Samples=vector(mode="list",length = 3)
names(low_Plaque_Samples)=names(high_Plaque_Samples)=names(msmm_rnaseq)

low_Plaque_Samples[[1]]=grep(pattern = "^BM_10",msmm_rnaseq_covariates2$Sample.ID[which(msmm_rnaseq_covariates2$PlaqueMean<=1)],value = T)
low_Plaque_Samples[[2]]=grep(pattern = "^BM_22",msmm_rnaseq_covariates2$Sample.ID[which(msmm_rnaseq_covariates2$PlaqueMean<=1)],value = T)
low_Plaque_Samples[[3]]=grep(pattern = "^BM_36",msmm_rnaseq_covariates2$Sample.ID[which(msmm_rnaseq_covariates2$PlaqueMean<=1)],value = T)

high_Plaque_Samples[[1]]=grep(pattern = "^BM_10",msmm_rnaseq_covariates2$Sample.ID[which(msmm_rnaseq_covariates2$PlaqueMean>=15)],value = T)
high_Plaque_Samples[[2]]=grep(pattern = "^BM_22",msmm_rnaseq_covariates2$Sample.ID[which(msmm_rnaseq_covariates2$PlaqueMean>=15)],value = T)
high_Plaque_Samples[[3]]=grep(pattern = "^BM_36",msmm_rnaseq_covariates2$Sample.ID[which(msmm_rnaseq_covariates2$PlaqueMean>=15)],value = T)

msmm_rnaseq.final_keep=lapply(msmm_rnaseq,function(x)rowSums(x>0)>=ncol(x)/3)
exprs_rank=vector(mode = "list",length = 3)
names(exprs_rank)=c("FP","STG","PHG")
for (t in 1:3){
  exprs_rank[[t]]<-msmm_rnaseq[[t]][msmm_rnaseq.final_keep[[t]],]
  number_of_combinations<-choose(nrow(exprs_rank[[t]]),2)
  c_exprs_rank=exprs_rank[[t]][,low_Plaque_Samples[[t]]]
  t_exprs_rank=exprs_rank[[t]][,high_Plaque_Samples[[t]]]
  n.c<-ncol(c_exprs_rank)
  n.t<-ncol(t_exprs_rank)
  gene.names<-rownames(exprs_rank[[t]])
  
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
    write.table(result, file=paste0("/shared/hidelab2/user/md4zsa/Work/Data/MSMM_RNAseq/MSMM_RNAseq_FinalRelease/FullDataset_DiffCorr/",names(exprs_rank)[[t]],"/", i, ".txt"), sep="\t", row.names=FALSE, quote = FALSE)
    i<-i+1
    start<-i*blocksize+1
    end<-min((i+1)*blocksize, number_of_combinations)
  }
  cat(paste("Done!"))
}
