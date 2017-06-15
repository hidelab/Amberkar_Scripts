library(data.table)
library(parallel)
library(cocor)
library(org.Hs.eg.db)
library(igraph)
library(gtools)
ncore = detectCores()
blocksize=100000
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
mayo_reseq_data=fread("MAYO/MAYO_CBE_TCX_netResidualExpression.tsv",sep="\t",header=T,data.table=F)
rownames(mayo_reseq_data)=mayo_reseq_data$ensembl_gene_id
ensembl_geneSymbol_map=mapIds2(IDs = mayo_reseq_data$ensembl_gene_id,IDFrom = "ENSEMBL",IDTo = "SYMBOL")
mayo_reseq_data2=mayo_reseq_data[-which(mayo_reseq_data$ensembl_gene_id%in%mapIds2(IDs = mayo_reseq_data$ensembl_gene_id,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[2]]),]
mayo_reseq_data2$gene_symbol=mapIds2(IDs = mayo_reseq_data$ensembl_gene_id,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[1]][,2]
mayo_reseq_data2.agg=aggregate(x=mayo_reseq_data2[,-c(1,529)],by=list(Symbol=mayo_reseq_data2$gene_symbol),mean)
rownames(mayo_reseq_data2.agg)=mayo_reseq_data2.agg$Symbol
mayo_reseq_data2.agg=mayo_reseq_data2.agg[,-1]
mayo_covariates_TCX=fread("MAYO/MayoRNAseq_RNAseq_TCX_covariates.csv",sep="\t",header=T,data.table=F)

exprs_rank=mayo_reseq_data2.agg
c_counts=mayo_reseq_data2.agg[,grep(pattern = paste(mayo_covariates_TCX$ID[grep(pattern = "Control",x = mayo_covariates_TCX$Diagnosis)],collapse = "|"),x = colnames(mayo_reseq_data2.agg))]
t_counts=mayo_reseq_data2.agg[,grep(pattern = paste(mayo_covariates_TCX$ID[grep(pattern = "AD",x = mayo_covariates_TCX$Diagnosis)],collapse = "|"),x = colnames(mayo_reseq_data2.agg))]
number_of_combinations<-choose(nrow(exprs_rank),2)
n.c<-ncol(c_counts)
n.t<-ncol(t_counts)
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
  write.table(result, file=paste0("MAYO/results/mayo_reseq_TCX_cocor_tmp", i, ".txt"), sep="\t", row.names=FALSE, quote = FALSE)
  i<-i+1
  start<-i*blocksize+1
  end<-min((i+1)*blocksize, number_of_combinations)
}

setwd("/shared/hidelab2/user/md4zsa/Work/Data/AMP-AD_RNAseq_ReSeq/Normalised_covariate_corrected_NoResiduals/MAYO/results")
system("awk 'FNR==1 && NR!=1 { while (/^<header>/) getline; }    1 {print}' mayo_reseq_TCX*.txt >MAYO_ReSeq_TCX_DiffCorr_Results.txt")

allResults_TCX=fread("MAYO_ReSeq_TCX_DiffCorr_Results.txt",sep = "\t",header = T,data.table = T,showProgress = T)
allResults_TCX$FDR=p.adjust(p = allResults_TCX$p.cocor,method = "fdr")
allResults_TCX$FDR.c=p.adjust(p = allResults_TCX$p.c,method = "fdr")
allResults_TCX$FDR.t=p.adjust(p = allResults_TCX$p.t,method = "fdr")
fwrite(allResults_TCX,"MAYO_ReSeq_TCX_DiffCorr_Results_FDR.txt",sep="\t",col.names = T,row.names = F,nThread = 12,buffMB = 100,showProgress = T)


