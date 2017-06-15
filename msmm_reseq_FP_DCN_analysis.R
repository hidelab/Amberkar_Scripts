library(cocor)
library(org.Hs.eg.db)
library(parallel)
library(data.table)
library(igraph)

ncore = detectCores()
ProcessElement <- function(ic){
  A = ceiling((sqrt(8*(ic+1)-7)+1)/2)
  B = ic-choose(floor(1/2+sqrt(2*ic)),2)
  
  c_A = as.numeric(c_counts[A,])
  c_B = as.numeric(c_counts[B,])
  
  t_A = as.numeric(t_counts[A,])
  t_B = as.numeric(t_counts[B,])
  
  # get correlation between the summaries for the unique genes
  tmp = data.frame(IC= ic, Gene.A=gene.names[A], Gene.B=gene.names[B], r.c=NA, p.c=NA, n.c=NA, r.t=NA, p.t=NA, n.t=NA, p.cocor=NA)
  #if( (var(c_A) * var(c_B))!=0)
  tmp$n.c<-sum(!is.na(c_A + c_B))
  if (tmp$n.c >=10)
  {
    c_cortest<-cor.test(c_A, c_B, method="spearman")
    tmp$r.c<-c_cortest$estimate
    tmp$p.c<-c_cortest$p.value
  }
  
  #if( (var(t_A) * var(t_B))!=0)
  tmp$n.t<-sum(!is.na(t_A + t_B))
  if(tmp$n.t >=10)
  {
    t_cortest<-cor.test(t_A, t_B, method="spearman")
    tmp$r.t<-t_cortest$estimate
    tmp$p.t<-t_cortest$p.value
  }
  
  if ( (!is.na(tmp$r.c)) && (!is.na(tmp$r.t)) )
  {
    diffcor<-cocor.indep.groups(tmp$r.c, tmp$r.t, tmp$n.c, tmp$n.t)
    tmp$p.cocor<-diffcor@fisher1925$p.value
  }
  
  #setTxtProgressBar(pb,ic %% n_part)
  setTxtProgressBar(pb,ic %% blocksize)
  return(tmp)
}
mapIds2<-function(IDs,IDFrom,IDTo){
  idmap=mapIds(x = org.Hs.eg.db,keys = IDs,column = IDTo,keytype = IDFrom,multiVals = "first")
  na_vec=names(idmap[is.na(idmap)==T])
  idmap=idmap[is.na(idmap)==F]
  idmap_df=data.frame("From"=names(idmap),"To"=unlist(unname(idmap)),stringsAsFactors = F)
  return(list(map=idmap_df,noMap=na_vec))
}
setwd("/shared/hidelab2/user/md4zsa/Work/Data/AMP-AD_RNAseq_ReSeq/Normalised_covariate_corrected_NoResiduals")
msmm_data=fread("./MSMM/MSSM_FP_STG_PHG_IFG_netResidualExpression.tsv",sep="\t",header=T,data.table=F,showProgress=T)
rownames(msmm_data)=msmm_data$ensembl_gene_id
ensembl_geneSymbol_map=mapIds2(IDs = msmm_data$ensembl_gene_id,IDFrom = "ENSEMBL",IDTo = "SYMBOL")
msmm_data2=msmm_data[-which(rownames(msmm_data)%in%mapIds2(IDs = msmm_data$ensembl_gene_id,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[2]]),]
msmm_data2$gene_symbol=mapIds2(IDs = rownames(msmm_data2),IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[1]][,2]
msmm_data2.agg=aggregate(x=msmm_data2[,-c(1,957)],by=list(Symbol=msmm_data2$gene_symbol),mean)
rownames(msmm_data2.agg)=msmm_data2.agg$Symbol
msmm_data2.agg=msmm_data2.agg[,-1]
msmm_reseq_design_matrix=read.table("./MSMM/MSSM_FP_STG_PHG_IFG_Design.tsv",sep = "\t",header = T,as.is = T)
number_of_combinations=choose(nrow(msmm_data2.agg),2)

c_exprs_rank=msmm_data2.agg[,which(colnames(msmm_data2.agg)%in%msmm_reseq_design_matrix$SampleID[msmm_reseq_design_matrix$BrainRegion.DiagnosisFP.CONTROL==1])]#16035
t_exprs_rank=msmm_data2.agg[,which(colnames(msmm_data2.agg)%in%msmm_reseq_design_matrix$SampleID[msmm_reseq_design_matrix$BrainRegion.DiagnosisFP.AD==1])]#12000
n.c<-ncol(c_exprs_rank)
n.t<-ncol(t_exprs_rank)
gene.names<-rownames(msmm_data2.agg)
dir.create("./MSMM/results_FP",showWarnings = T,mode = "0777")
#setwd("./MSMM/results_FP")
i<-0
blocksize=100000
start<-i*blocksize+1
end<-min((i+1)*blocksize, number_of_combinations)

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
  fwrite(result, file=paste("./MSMM/results_FP/mssm_reseq_FP_cocor_corrPval_tmp", i, ".txt",sep = ""), sep="\t",col.names = T,row.names = F,buffMB = 100,nThread = 16,quote = F)
  i<-i+1
  start<-i*blocksize+1
  end<-min((i+1)*blocksize, number_of_combinations)
}
setwd("./MSMM/results_FP")
#Remove headers from tmp files and combine in a single one
system(paste("find . -name '*.txt'|grep 'tmp'|xargs -n 1 tail -n +2 ",paste(">MSMM_ReSeq_FP_corrPval",sep = "_",collapse = "_"),"_allResults_DiffCorr.txt",sep = ""))
allResults_FDR=fread(input = list.files(pattern = "*allResults_DiffCorr.txt"),sep = "\t",header = F,showProgress = T,data.table = F,strip.white = T,stringsAsFactors = F)
#Assign single header for the combined file
colnames(allResults_FDR)=colnames(result)
#Compute FDR values for all p-values
allResults_FDR$FDR.cocor=p.adjust(p = as.numeric(allResults_FDR$p.cocor),method = "fdr")
allResults_FDR$FDR.c=p.adjust(p = as.numeric(allResults_FDR$p.c),method = "fdr")
allResults_FDR$FDR.t=p.adjust(p = as.numeric(allResults_FDR$p.t),method = "fdr")
write.table(allResults_FDR,"msmm_reseq_FP_corrPval_allResults_DiffCorr_FDR.txt",sep = "\t",col.names = T,row.names = F)
#Generate DCNs adn write to R objects
control_DCN=AD_DCN=graph.data.frame(d = allResults_FDR[which(allResults_FDR$FDR.cocor<=0.1),c(1:2)],directed = F)
E(control_DCN)$weight=allResults_FDR[which(allResults_FDR$FDR.cocor<=0.1),3]+1
E(AD_DCN)$weight=allResults_FDR[which(allResults_FDR$FDR.cocor<=0.1),6]+1
saveRDS(control_DCN,file = paste("mssm_reseq_FP_corrPval_Control.RDS",sep = ""))
saveRDS(AD_DCN,file = paste("mssm_reseq_FP_corrPval_AD.RDS",sep = ""))
