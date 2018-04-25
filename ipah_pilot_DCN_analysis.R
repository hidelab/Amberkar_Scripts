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

setwd("/shared/hidelab2/user/md4zsa/Work/Data/IPAH")

ipah_counts.normalised=read.table("IPAH_Pilot_NormalisedCounts.txt",header = T,sep = "\t",row.names = 1)
rownames(ipah_counts.normalised)=unlist(lapply(strsplit(x = rownames(ipah_counts.normalised),split = "\\."),`[[`,1))
ipah_metadata=readRDS("lawrie_sample_group.RDS")
ipah_metadata$External.ID[ipah_metadata$group=="HV"]=gsub(pattern = "_v1",replacement = "",x = ipah_metadata$External.ID[ipah_metadata$group=="HV"])
ipah_counts.filtered1=ipah_counts.normalised[-which(rownames(ipah_counts.normalised)%in%mapIds2(IDs = rownames(ipah_counts.normalised),IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[2]]),]
ipah_counts.filtered2=ipah_counts.filtered1[rowSums(ipah_counts.filtered1>0)>=ncol(ipah_counts.filtered1)/3,]

c_counts=ipah_counts.filtered2[,grep(pattern = paste(ipah_metadata$External.ID[ipah_metadata$group=="HV"],collapse = "|"),x = colnames(ipah_counts.filtered2))]
t_counts=ipah_counts.filtered2[,grep(pattern = paste(ipah_metadata$External.ID[ipah_metadata$group=="IPAH"],collapse = "|"),x = colnames(ipah_counts.filtered2))]
n.c=ncol(c_counts)
n.t=ncol(t_counts)
gene.names=rownames(ipah_counts.filtered2)
number_of_combinations=choose(length(gene.names),2)
dir.create("cocor_results",showWarnings = T,mode = "0777")

