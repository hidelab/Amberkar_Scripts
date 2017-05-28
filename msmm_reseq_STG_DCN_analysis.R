library(cocor)
library(org.Hs.eg.db)
library(parallel)
library(data.table)

ncore = detectCores()
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
setwd("/shared/hidelab2/user/md4zsa/Work/Data/AMP-AD_RNAseq_ReSeq/Normalised_covariate_corrected_NoResiduals")
msmm_data=fread("./MSMM/MSSM_FP_STG_PHG_IFG_netResidualExpression.tsv",sep="\t",header=T,data.table=F,showProgress=T)
rownames(msmm_data)=msmm_data$ensembl_gene_id
#Remove genes with at least 20 NA counts
msmm_data2=msmm_data[rowSums(is.na(msmm_data))<20,]
msmm_reseq_design_matrix=read.table("./MSMM/MSSM_FP_STG_PHG_IFG_Design.tsv",sep = "\t",header = T,as.is = T)
number_of_combinations=choose(nrow(msmm_data),2)

msmm_reseq_design_matrix=read.table("./MSMM/MSSM_FP_STG_PHG_IFG_Design.tsv",sep = "\t",header = T,as.is = T)
number_of_combinations=choose(nrow(msmm_data),2)
c_exprs_rank=msmm_data2[,which(colnames(msmm_data)%in%msmm_reseq_design_matrix$SampleID[msmm_reseq_design_matrix$BrainRegion.DiagnosisSTG.CONTROL==1])]
t_exprs_rank=msmm_data2[,which(colnames(msmm_data)%in%msmm_reseq_design_matrix$SampleID[msmm_reseq_design_matrix$BrainRegion.DiagnosisSTG.AD==1])]
n.c<-ncol(c_exprs_rank)
n.t<-ncol(t_exprs_rank)
gene.names<-rownames(msmm_data)
dir.create("./MSMM/results_STG",showWarnings = T,mode = "0777")
i<-442
blocksize=100000
start<-i*blocksize+1
end<-min((i+1)*blocksize, number_of_combinations)
#setwd("./MSMM/results_STG")
#str(Reduce(intersect,lapply(lapply(exprs_rank,function(x)x[which((rowSums(x>0)>=ncol(x)/3)==T),]),rownames)))
while(start < number_of_combinations){
  input<-start:end
  pb = txtProgressBar(min=0,max=length(input),style=3,initial=0)
  cat("\n")
  res = mclapply(input,ProcessElement,mc.cores=ncore)
  close(pb)	
  #save results
  #rbindlist is faster than rbind.fill
  result <- rbindlist(res)
  result <- as.data.frame(result)
  result <- data.frame(result,stringsAsFactors = F)
  #Write results, use fwrite instead of write.table for faster processing
  fwrite(result, file=paste("./MSMM/results_STG/mssm_reseq_STG_cocor_tmp", i, ".txt",sep = ""), sep="\t",col.names = T,row.names = F,buffMB = 100,nThread = 16,quote = F)
  i<-i+1
  start<-i*blocksize+1
  end<-min((i+1)*blocksize, number_of_combinations)
}

#Remove headers from tmp files and combine in a single one
system(paste("find . -name '*.txt'|grep 'tmp'|xargs -n 1 tail -n +2 ",paste(">MSMM_ReSeq_STG",sep = "_",collapse = "_"),"_allResults_DiffCorr.txt",sep = ""))
allResults_FDR=fread(input = list.files(pattern = "*allResults_DiffCorr.txt"),sep = "\t",header = F,showProgress = T,data.table = F,strip.white = T,stringsAsFactors = F)
#Assign single header for the combined file
colnames(allResults_FDR)=colnames(result)
#Compute FDR values for all p-values
allResults_FDR$FDR.cocor=p.adjust(p = as.numeric(allResults_FDR$p.cocor),method = "fdr")
allResults_FDR$FDR.c=p.adjust(p = as.numeric(allResults_FDR$p.c),method = "fdr")
allResults_FDR$FDR.t=p.adjust(p = as.numeric(allResults_FDR$p.t),method = "fdr")
fwrite(allResults_FDR,"msmm_reseq_STG_allResults_DiffCorr_FDR.txt",sep = "\t",col.names = T,row.names = F,buffMB = 100,nThread = 12,showProgress = T)
#Generate DCNs adn write to R objects
control_DCN=AD_DCN=graph.data.frame(d = allResults_FDR[which(allResults_FDR$FDR.cocor<=0.1),c(1:2)],directed = F)
E(control_DCN)$weight=allResults_FDR[which(allResults_FDR$FDR.cocor<=0.1),3]+1
E(AD_DCN)$weight=allResults_FDR[which(allResults_FDR$FDR.cocor<=0.1),6]+1
saveRDS(control_DCN,file = paste("mssm_reseq_STG_Control.RDS",sep = ""))
saveRDS(AD_DCN,file = paste("mssm_reseq_STG_AD.RDS",sep = ""))
