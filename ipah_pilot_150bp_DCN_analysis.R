library(org.Hs.eg.db)
library(data.table)
library(parallel)
library(diffcoexp)
library(gProfileR)
library(dplyr)
library(magrittr)
library(DCGL)

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

setwd("/shared/hidelab2/user/md4zsa/Work/Data/IPAH")

ipah_150bp_counts.normalised=data.frame(readRDS("IPAH_150bp_normCounts.RDS"),stringsAsFactors = F)
rownames(ipah_150bp_counts.normalised)=unlist(lapply(strsplit(x = rownames(ipah_150bp_counts.normalised),split = "\\."),`[[`,1))
ipah_metadata=readRDS("lawrie_sample_group.RDS")
ipah_metadata$External.ID[ipah_metadata$group=="HV"]=gsub(pattern = "_v1",replacement = "",x = ipah_metadata$External.ID[ipah_metadata$group=="HV"])
ipah_counts.filtered1=ipah_150bp_counts.normalised[-which(rownames(ipah_150bp_counts.normalised)%in%mapIds2(IDs = rownames(ipah_150bp_counts.normalised),IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[2]]),]
ipah_counts.filtered2=ipah_counts.filtered1[rowSums(ipah_counts.filtered1>0)>=ncol(ipah_counts.filtered1)/3,]
ipah_counts.filtered2$symbol=mapIds2(IDs=rownames(ipah_counts.filtered2),IDFrom="ENSEMBL",IDTo="SYMBOL")[[1]][,2]
ipah_counts.filtered2$entrez=mapIds2(IDs=rownames(ipah_counts.filtered2),IDFrom="ENSEMBL",IDTo="ENTREZID")[[1]][,2]
ipah_counts_filtered2.agg=aggregate(x=ipah_counts.filtered2[,-43],by=list(Symbol=ipah_counts.filtered2$symbol),mean)
rownames(ipah_counts_filtered2.agg)=ipah_counts_filtered2.agg$Symbol
ipah_counts_filtered2.agg=ipah_counts_filtered2.agg[,-1]

c_counts=ipah_counts.filtered2[,grep(pattern = paste(ipah_metadata$External.ID[ipah_metadata$group=="HV"],collapse = "|"),x = colnames(ipah_counts_filtered2.agg))]
t_counts=ipah_counts.filtered2[,grep(pattern = paste(ipah_metadata$External.ID[ipah_metadata$group=="IPAH"],collapse = "|"),x = colnames(ipah_counts_filtered2.agg))]
n.c=ncol(c_counts)
n.t=ncol(t_counts)
gene.names=rownames(ipah_counts.filtered2)

#Read TF-Target interactions, preprocess data
regnet_tf2target=fread("/shared/hidelab2/user/md4zsa/Work/Data/TF_Databases/RegNetwork_human_regulators2.txt",header = T,sep = "\t",showProgress = T,data.table = F)%>%filter(evidence=="Experimental")%>%select(c(regulator_symbol,target_symbol))

ipah_150_diffcoexp=diffcoexp(exprs.1 = c_counts,exprs.2 = t_counts,rth=0.6, qth=0.1, r.diffth=0.1, q.diffth=0.1)
ipah_150_diffcoexp_DRsort=DRsort(DCGs = ipah_150_diffcoexp$DCGs,DCLs = ipah_150_diffcoexp$DCLs,tf2target = regnet_tf2target,expGenes = rownames(ipah_counts_filtered2.agg))
saveRDS(ipah_150_diffcoexp,"IPAH_150_diffcoexp_noLFC_bugfix.res.RDS")
saveRDS(ipah_150_diffcoexp_DRsort,"IPAH_150_diffcoexp_DRsort_noLFC_bugfix.res.RDS")
