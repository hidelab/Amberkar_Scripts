library(data.table)
library(parallel)
library(cocor)
library(org.Hs.eg.db)
library(igraph)
library(gtools)

setwd("/shared/hidelab2/user/md4zsa/Work/Data/AMP-AD_RNAseq_ReSeq/Normalised_covariate_corrected_NoResiduals")
mayo_reseq_data=fread("MAYO/MAYO_CBE_TCX_netResidualExpression.tsv",sep="\t",header=T,data.table=F)
rownames(mayo_reseq_data)=mayo_reseq_data$ensembl_gene_id
mayo_covariates_CER=fread("MAYO/MayoRNAseq_RNAseq_CBE_covariates.csv",sep="\t",header=T,data.table=F)
mayo_covariates_TCX=fread("MAYO/MayoRNAseq_RNAseq_TCX_covariates.csv",sep="\t",header=T,data.table=F)
mayo_reseq_dataList=vector(mode="list",length=2)
names(mayo_reseq_dataList)=c("CER","TCX")
mayo_reseq_dataList$CER=mayo_reseq_dataList$TCX=vector(mode="list",length=2)
names(mayo_reseq_dataList$CER)=names(mayo_reseq_dataList$TCX)=c("Control","AD")
mayo_reseq_dataList$CER$Control=mayo_reseq_data[,grep(pattern = paste(mayo_covariates_CER$SampleID[grep(pattern = "Control",x = mayo_covariates_CER$Diagnosis)],collapse = "|"),x = colnames(mayo_reseq_data))]
mayo_reseq_dataList$CER$AD=mayo_reseq_data[,grep(pattern = paste(mayo_covariates_CER$SampleID[grep(pattern = "AD",x = mayo_covariates_CER$Diagnosis)],collapse = "|"),x = colnames(mayo_reseq_data))]
mayo_reseq_dataList$TCX$Control=mayo_reseq_data[,grep(pattern = paste(mayo_covariates_TCX$ID[grep(pattern = "Control",x = mayo_covariates_TCX$Diagnosis)],collapse = "|"),x = colnames(mayo_reseq_data))]
mayo_reseq_dataList$TCX$AD=mayo_reseq_data[,grep(pattern = paste(mayo_covariates_TCX$ID[grep(pattern = "AD",x = mayo_covariates_TCX$Diagnosis)],collapse = "|"),x = colnames(mayo_reseq_data))]

mayo_reseq_data.final_keep=vector(mode = "list",length = 2)
names(mayo_reseq_data.final_keep)=c("CER","TCX")
d1=do.call("cbind",mayo_reseq_dataList$CER)
d2=do.call("cbind",mayo_reseq_dataList$TCX)
rownames(d1)=rownames(mayo_reseq_dataList$CER$Control)
rownames(d2)=rownames(mayo_reseq_dataList$TCX$Control)
mayo_reseq_data.final_keep$CER=d1[which((rowSums(d1>0)>=ncol(d1)/3)==T),]
mayo_reseq_data.final_keep$TCX=d2[which((rowSums(d2>0)>=ncol(d2)/3)==T),]

ncore = detectCores()
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
for(t in 1:2){
  exprs_rank=mayo_reseq_data.final_keep[[t]]
  number_of_combinations<-choose(nrow(exprs_rank),2)
  Control_samples=grep("Control",colnames(mayo_reseq_data.final_keep[[t]]))
  AD_samples=grep("AD",colnames(mayo_reseq_data.final_keep[[t]]))
  c_exprs_rank=exprs_rank[,Control_samples]
  t_exprs_rank=exprs_rank[,AD_samples]
  n.c<-ncol(c_exprs_rank)
  n.t<-ncol(t_exprs_rank)
  gene.names<-rownames(exprs_rank)

  i<-0
  start<-i*blocksize+1
  end<-min((i+1)*blocksize, number_of_combinations)

  while(start < number_of_combinations){
  input<-start:end
  pb = txtProgressBar(min=0,max=length(input),style=3,initial=0)
  cat("\n")
  res = mclapply(input,ProcessElement,mc.cores=ncore)
  close(pb)	
  # save results
  #saveRDS(res,paste0("/shared/hidelab2/user/md4zsa/Work/Data/ALS/fcx_c9orf72_als_", i, ".RDS"))
  #rbindlist is faster than rbind.fill
  result <- rbindlist(res)
  result <- as.data.frame(result)
  result <- data.frame(result,stringsAsFactors = F)
  #result$FDR<-p.adjust(result$p.cocor, method="fdr")
  #result[,c(3:8,10:11)]=round(result[,c(3:8,10:11)],digits = 3)
  write.table(result, file=paste0("MAYO/results/mayo_reseq_",names(mayo_reseq_data.final_keep)[t],"_cocor_tmp", i, ".txt"), sep="\t", row.names=FALSE, quote = FALSE)
  i<-i+1
  start<-i*blocksize+1
  end<-min((i+1)*blocksize, number_of_combinations)
}
  cat(paste("Done!"))
}

setwd("/shared/hidelab2/user/md4zsa/Work/Data/AMP-AD_RNAseq_ReSeq/Normalised_covariate_corrected_NoResiduals/MAYO/results")
system("awk 'FNR==1 && NR!=1 { while (/^<header>/) getline; }    1 {print}' mayo_reseq_CER*.txt >MAYO_ReSeq_CER_DiffCorr_Results.txt")
system("awk 'FNR==1 && NR!=1 { while (/^<header>/) getline; }    1 {print}' mayo_reseq_TCX*.txt >MAYO_ReSeq_TCX_DiffCorr_Results.txt")
allResults_CER=fread("MAYO_ReSeq_CER_DiffCorr_Results.txt",sep = "\t",header = T,data.table = T,showProgress = T)
allResults_TCX=fread("MAYO_ReSeq_TCX_DiffCorr_Results.txt",sep = "\t",header = T,data.table = T,showProgress = T)
allResults_CER$FDR=p.adjust(p = allResults_CER$p.cocor,method = "fdr")
allResults_CER$FDR.c=p.adjust(p = allResults_CER$p.c,method = "fdr")
allResults_CER$FDR.t=p.adjust(p = allResults_CER$p.t,method = "fdr")
allResults_TCX$FDR=p.adjust(p = allResults_TCX$p.cocor,method = "fdr")
allResults_TCX$FDR.c=p.adjust(p = allResults_TCX$p.c,method = "fdr")
allResults_TCX$FDR.t=p.adjust(p = allResults_TCX$p.t,method = "fdr")


fwrite(allResults_CER,"MAYO_ReSeq_CER_DiffCorr_Results_FDR.txt",sep="\t",col.names = T,row.names = F,nThread = 12,buffMB = 100,showProgress = T)
fwrite(allResults_TCX,"MAYO_ReSeq_TCX_DiffCorr_Results_FDR.txt",sep="\t",col.names = T,row.names = F,nThread = 12,buffMB = 100,showProgress = T)

