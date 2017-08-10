library(cocor)
library(org.Hs.eg.db)
library(data.table)
library(parallel)

mapIds2<-function(IDs,IDFrom,IDTo){
  require(org.Hs.eg.db)
  idmap=mapIds(x = org.Hs.eg.db,keys = IDs,column = IDTo,keytype = IDFrom,multiVals = "first")
  na_vec=names(idmap[is.na(idmap)==T])
  idmap=idmap[is.na(idmap)==F]
  idmap_df=data.frame("From"=names(idmap),"To"=unlist(unname(idmap)),stringsAsFactors = F)
  return(list(map=idmap_df,noMap=na_vec))
}
jaccard=function(A,B){
  jc=set_cardinality(intersect(A,B))/set_cardinality(union(A,B))
  return(jc)
}
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
ncore=12
setwd("/shared/hidelab2/user/md4zsa/Work/Data/IPAH")

ipah_counts.normalised=read.table("IPAH_Pilot_NormalisedCounts.txt",header = T,sep = "\t",row.names = 1)
ipah_metadata=readRDS("lawrie_sample_group.RDS")
ipah_metadata$External.ID[ipah_metadata$group=="HV"]=gsub(pattern = "_v1",replacement = "",x = ipah_metadata$External.ID[ipah_metadata$group=="HV"])
ipah_counts.filtered=ipah_counts.normalised[rowSums(ipah_counts.normalised>0)>=ncol(ipah_counts.normalised)/3,]

c_counts=ipah_counts.filtered[,grep(pattern = paste(ipah_metadata$External.ID[ipah_metadata$group=="HV"],collapse = "|"),x = colnames(ipah_counts.filtered))]
t_counts=ipah_counts.filtered[,grep(pattern = paste(ipah_metadata$External.ID[ipah_metadata$group=="IPAH"],collapse = "|"),x = colnames(ipah_counts.filtered))]
n.c=ncol(c_counts)
n.t=ncol(t_counts)
gene.names=rownames(ipah_counts.filtered)
gene.names2=gene.names[which(gene.names%in%mapIds2(IDs = gene.names,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[1]][,1])]
number_of_combinations=choose(length(gene.names2),2)
dir.create("cocor_results",showWarnings = T,mode = "0777")

i=601
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
  fwrite(result, file=paste("./cocor_results/ipah_cocor_tmp", i, ".txt",sep = ""), sep="\t",col.names = T,row.names = F,buffMB = 100,nThread = 16,quote = F)
  i<-i+1
  start<-i*blocksize+1
  end<-min((i+1)*blocksize, number_of_combinations)
}
setwd("./cocor_results")
system(paste("find . -name '*.txt'|grep 'tmp'|xargs -n 1 tail -n +2 ",paste(">IPAH",sep = "_",collapse = "_"),"_allResults_DiffCorr.txt",sep = ""))
allResults_FDR=fread(input = list.files(pattern = "*allResults_DiffCorr.txt"),sep = "\t",header = F,showProgress = T,data.table = F,strip.white = T,stringsAsFactors = F)
colnames(allResults_FDR)=colnames(result)
#Compute FDR values for all p-values
allResults_FDR$FDR.cocor=p.adjust(p = as.numeric(allResults_FDR$p.cocor),method = "fdr")
allResults_FDR$FDR.c=p.adjust(p = as.numeric(allResults_FDR$p.c),method = "fdr")
allResults_FDR$FDR.t=p.adjust(p = as.numeric(allResults_FDR$p.t),method = "fdr")
write.table(allResults_FDR,"IPAH_corrPval_allResults_DiffCorr_FDR.txt",sep = "\t",col.names = T,row.names = F)

