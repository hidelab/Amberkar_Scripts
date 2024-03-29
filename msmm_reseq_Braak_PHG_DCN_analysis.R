library(cocor)
library(org.Hs.eg.db)
library(parallel)
library(data.table)
library(centiserve)

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
setwd("/shared/hidelab2/user/md4zsa/Work/Data/AMP-AD_RNAseq_ReSeq/Normalised_covariate_corrected_NoResiduals/")
msmm_data=fread("./MSMM/MSSM_FP_STG_PHG_IFG_netResidualExpression.tsv",sep="\t",header=T,data.table=F,showProgress=T)
msmm_data2=msmm_data[,-grep(pattern="resequenced",x=colnames(msmm_data))]
rownames(msmm_data2)=msmm_data$ensembl_gene_id

#Merge clinical and pathological covariates
msmm_rnaseq_covariates=read.csv("./MSMM/MSBB_RNAseq_covariates.csv",header = T,as.is = T)
msmm_rnaseq_clinical_covariates=read.csv("./MSMM/MSBB_clinical.csv",header = T,as.is = T)
msmm_rnaseq_covariates.merged=merge(x = msmm_rnaseq_clinical_covariates,y=msmm_rnaseq_covariates,by=c("individualIdentifier","individualIdentifier"))
msmm_rnaseq_covariates.merged2=msmm_rnaseq_covariates.merged[grep(pattern = "unmapped|resequenced",x = msmm_rnaseq_covariates.merged$fileName,invert = T),]
msmm_rnaseq_covariates.merged_allBraak=msmm_rnaseq_covariates.merged2[grep(pattern="0|9",x=msmm_rnaseq_covariates.merged2$bbscore,invert=T),]

msmm_rnaseq_sampleList.Braak=vector(mode = "list",length = 4)
names(msmm_rnaseq_sampleList.Braak)=c("FP","IFG","STG","PHG")
msmm_rnaseq_sampleList.Braak$FP=msmm_rnaseq_sampleList.Braak$IFG=msmm_rnaseq_sampleList.Braak$STG=msmm_rnaseq_sampleList.Braak$PHG=vector(mode = "list",length = 6)
names(msmm_rnaseq_sampleList.Braak$FP)=names(msmm_rnaseq_sampleList.Braak$IFG)=names(msmm_rnaseq_sampleList.Braak$STG)=names(msmm_rnaseq_sampleList.Braak$PHG)=paste("B",c(1:6),sep = "")

braak_scores=c(1:6)
brain_regions=c("BM10","BM44","BM22","BM36")

for(r in 1:length(brain_regions)){
  for(b in 1:length(braak_scores)){
    msmm_rnaseq_sampleList.Braak[[r]][[b]]=msmm_rnaseq_covariates.merged_allBraak$sampleIdentifier[which((msmm_rnaseq_covariates.merged_allBraak$BrodmannArea==brain_regions[r])&(msmm_rnaseq_covariates.merged_allBraak$bbscore==braak_scores[b]))]
  }
}


exprs_rank=lapply(msmm_rnaseq_sampleList.Braak$PHG,function(x)msmm_data2[,which(colnames(msmm_data2)%in%x)])
common_genes=Reduce(intersect,lapply(lapply(exprs_rank,function(x)x[which((rowSums(x>0)>=ncol(x)/3)==T),]),rownames))
exprs_rank2=lapply(lapply(exprs_rank,function(x)x[which((rowSums(x>0)>=ncol(x)/3)==T),]),function(y)y[rownames(y)%in%common_genes,])
number_of_combinations=choose(length(common_genes),2)
braak_scores_combn=matrix(NA,nrow=5,ncol=2)
for (i in 1:5){
  braak_scores_combn[i,]=c(i,i+1)  
}
dir.create("results")
for(c in 1:5){
  c_exprs_rank=exprs_rank2[braak_scores_combn[c,1]][[1]]
  t_exprs_rank=exprs_rank2[braak_scores_combn[c,2]][[1]]
  n.c<-ncol(c_exprs_rank)
  n.t<-ncol(t_exprs_rank)
  gene.names<-rownames(exprs_rank2$B1)
  dir.create(paste("results/Braak",braak_scores_combn[c,1],"_","Braak",braak_scores_combn[c,2],sep = ""),showWarnings = T,mode = "0777")
  setwd(paste("Braak",braak_scores_combn[c,1],"_","Braak",braak_scores_combn[c,2],"/",sep = ""))
  i<-0
  blocksize=50
  start<-i*blocksize+1
  end<-min((i+1)*blocksize, number_of_combinations)
  
  #str(Reduce(intersect,lapply(lapply(exprs_rank,function(x)x[which((rowSums(x>0)>=ncol(x)/3)==T),]),rownames)))
  while(start < 100){
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
    fwrite(result, file=paste("mssm_reseq_PHG_",paste("Braak",braak_scores_combn[c,],sep = "_",collapse = "_"),"_cocor_tmp", i, ".txt",sep = ""), sep="\t",col.names = T,row.names = F,buffMB = 10,nThread = 4,quote = F)
    i<-i+1
    start<-i*blocksize+1
    end<-min((i+1)*blocksize, number_of_combinations)
  }
  #Remove headers from tmp files and combine in a single one
  system(paste("find . -name '*.txt'|grep 'tmp'|xargs -n 1 tail -n +2",paste(">MSMM_ReSeq_PHG_",paste("Braak",braak_scores_combn[c,],sep = "_",collapse = "_"),"_allResults_DiffCorr.txt",sep = "")))
  allResults_FDR=fread(input = list.files(pattern = "*allResults_DiffCorr.txt"),sep = "\t",header = F,showProgress = T,data.table = F,strip.white = T,stringsAsFactors = F)
  #Assign single header for the combined file
  colnames(allResults_FDR)=colnames(result)
  #Compute FDR values for all p-values
  allResults_FDR$FDR.cocor=p.adjust(p = as.numeric(allResults_FDR$p.cocor),method = "fdr")
  allResults_FDR$FDR.c=p.adjust(p = as.numeric(allResults_FDR$p.c),method = "fdr")
  allResults_FDR$FDR.t=p.adjust(p = as.numeric(allResults_FDR$p.t),method = "fdr")
  fwrite(allResults_FDR,file = paste("mssm_reseq_PHG_",paste("Braak",braak_scores_combn[c,],sep = "_",collapse = "_"),".txt",sep="",collapse = "_"),sep = "\t",col.names = T,row.names = F,buffMB = 100,nThread = 12)
  #Generate DCNs adn write to R objects
  control_DCN=disease_DCN=graph.data.frame(d = allResults_FDR[which(allResults_FDR$FDR.cocor<=0.1),c(1:2)],directed = F)
  E(control_DCN)$weight=allResults_FDR[which(allResults_FDR$FDR.cocor<=0.1),3]+1
  E(disease_DCN)$weight=allResults_FDR[which(allResults_FDR$FDR.cocor<=0.1),6]+1
  saveRDS(control_DCN,file = paste("mssm_reseq_PHG_",paste("Braak",braak_scores_combn[c,],sep = "_",collapse = "_"),"_Control.RDS",sep = ""))
  saveRDS(disease_DCN,file = paste("mssm_reseq_PHG_",paste("Braak",braak_scores_combn[c,],sep = "_",collapse = "_"),"_Disease.RDS",sep = ""))
  setwd("../")
}
