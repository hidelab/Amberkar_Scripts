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
  if (tmp$n.c >=5)
  {
    c_cortest<-cor.test(c_A, c_B, method="spearman")
    tmp$r.c<-c_cortest$estimate
    tmp$p.c<-c_cortest$p.value
  }
  
  #if( (var(t_A) * var(t_B))!=0)
  tmp$n.t<-sum(!is.na(t_A + t_B))
  if(tmp$n.t >=5)
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

ipah_150bp_counts.normalised=data.frame(readRDS("IPAH_150bp_normCounts.RDS"),stringsAsFactors = F)
rownames(ipah_150bp_counts.normalised)=unlist(lapply(strsplit(x = rownames(ipah_150bp_counts.normalised),split = "\\."),`[[`,1))
ipah_metadata=readRDS("lawrie_sample_group.RDS")
ipah_metadata$External.ID[ipah_metadata$group=="HV"]=gsub(pattern = "_v1",replacement = "",x = ipah_metadata$External.ID[ipah_metadata$group=="HV"])
ipah_counts.filtered1=ipah_150bp_counts.normalised[-which(rownames(ipah_150bp_counts.normalised)%in%mapIds2(IDs = rownames(ipah_150bp_counts.normalised),IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[2]]),]
ipah_counts.filtered2=ipah_counts.filtered1[rowSums(ipah_counts.filtered1>0)>=ncol(ipah_counts.filtered1)/3,]

c_counts=ipah_counts.filtered2[,grep(pattern = paste(ipah_metadata$External.ID[ipah_metadata$group=="HV"],collapse = "|"),x = colnames(ipah_counts.filtered2))]
t_counts=ipah_counts.filtered2[,grep(pattern = paste(ipah_metadata$External.ID[ipah_metadata$group=="IPAH"],collapse = "|"),x = colnames(ipah_counts.filtered2))]
n.c=ncol(c_counts)
n.t=ncol(t_counts)
gene.names=rownames(ipah_counts.filtered2)
number_of_combinations=choose(length(gene.names),2)
dir.create("cocor_results_150bp",showWarnings = T,mode = "0777")

i=84
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
  fwrite(result, file=paste("./cocor_results_150bp/ipah_cocor_tmp_150bp_", i, ".txt",sep = ""), sep="\t",col.names = T,row.names = F,buffMB = 100,nThread = 16,quote = F)
  i<-i+1
  start<-i*blocksize+1
  end<-min((i+1)*blocksize, number_of_combinations)
}
setwd("./cocor_results_150bp")
system(paste("find . -name '*.txt'|grep 'tmp'|xargs -n 1 tail -n +2 ",paste(">IPAH",sep = "_",collapse = "_"),"_allResults_DiffCorr.txt",sep = ""))
allResults_FDR=fread(input = list.files(pattern = "*allResults_DiffCorr.txt"),sep = "\t",header = F,showProgress = T,data.table = F,strip.white = T,stringsAsFactors = F)
colnames(allResults_FDR)=colnames(result)
#Compute FDR values for all p-values
allResults_FDR$FDR.cocor=p.adjust(p = as.numeric(allResults_FDR$p.cocor),method = "fdr")
allResults_FDR$FDR.c=p.adjust(p = as.numeric(allResults_FDR$p.c),method = "fdr")
allResults_FDR$FDR.t=p.adjust(p = as.numeric(allResults_FDR$p.t),method = "fdr")
fwrite(allResults_FDR,"IPAH_150bp_allResults_DiffCorr_FDR.txt",sep = "\t",col.names = T,row.names = F,buffMB = 200,nThread = 18)
HV_DCN=IPAH_DCN=graph.data.frame(d = allResults_FDR[which(allResults_FDR$FDR.cocor<=0.1),c(1:2)],directed = F)
E(HV_DCN)$weight=allResults_FDR[which(allResults_FDR$FDR.cocor<=0.05),4]+1
E(IPAH_DCN)$weight=allResults_FDR[which(allResults_FDR$FDR.cocor<=0.05),7]+1
saveRDS(HV_DCN,file = paste("HV_150bp_DCN05.RDS",sep = ""))
saveRDS(IPAH_DCN,file = paste("IPAH_150bp_DCN05.RDS",sep = ""))

