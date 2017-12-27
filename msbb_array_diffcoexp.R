library(org.Hs.eg.db)
library(dplyr)
library(magrittr)
library(data.table)
library(parallel)
#library(clusterProfiler)
library(doParallel)
library(diffcoexp)

mapIds2<-function(IDs,IDFrom,IDTo){
  require(org.Hs.eg.db)
  idmap=mapIds(x = org.Hs.eg.db,keys = IDs,column = IDTo,keytype = IDFrom,multiVals = "first")
  na_vec=names(idmap[is.na(idmap)==T])
  idmap=idmap[is.na(idmap)==F]
  idmap_df=data.frame("From"=names(idmap),"To"=unlist(unname(idmap)),stringsAsFactors = F)
  return(list(map=idmap_df,noMap=na_vec))
}

setwd("/shared/hidelab2/user/md4zsa/Work/Data/MSMM_Array/MSBB_Array19/Normalised_Data")
msbb_array19.files=list.files(pattern = "*.tsv",full.names = T)
msbb_array19.files=msbb_array19.files[-c(18:19)]
msbb_array19=mapply(FUN = fread,msbb_array19.files,MoreArgs = list(header=T,sep="\t",data.table=F,showProgress=T))
names(msbb_array19)=c("AC","CN","DLPFC","FP","HP","IFG","ITG","MTG","OVC","PHG","PCC","PCG","PFC","PTMN","SPL","STG","TP")
msbb_array19.covariates=read.delim2("AMP-AD_MSBB_MSSM_covariates_mRNA_AffymetrixU133AB.txt",header = T,as.is = T)%>%
  mutate(Age=replace(Age,Age=="89+",90))%>%
  mutate(pH=as.numeric(pH))%>%
  mutate(CDR=as.numeric(CDR))%>%
  mutate(PLQ_Mn=as.numeric(PLQ_Mn))%>%
  mutate(BrainBank=paste("X",BrainBank,sep=""))
msbb_array19.2=lapply(msbb_array19[-18],function(x){rownames(x)<-x$ID;x})
#msbb_array19.3=lapply(msbb_array19.2,function(x)x[,-c(1:4)])
multiID.u133a=msbb_array19.2[[1]]$ENTREZ_GENE_ID[grep(pattern="///",msbb_array19.2[[1]]$ENTREZ_GENE_ID)]
msbb_array19.2.agg=lapply(msbb_array19.2[-18],function(y)aggregate(x=y[,-c(1:4)],by=list(EntrezID=y$ENTREZ_GENE_ID),mean))
msbb_array19.2.agg2=lapply(msbb_array19.2.agg,function(x){rownames(x)<- x$EntrezID;x <- x[-grep(pattern = "///",x = rownames(x)),];x})
msbb_array19.covariates$SampleType="OTHER"

msbb_array19.covariates$SampleType[(msbb_array19.covariates$CDR<=0.5&msbb_array19.covariates$NP1<=1&msbb_array19.covariates$Braak<=3)]="CONTROL"
msbb_array19.covariates$SampleType[(msbb_array19.covariates$CDR>=1&msbb_array19.covariates$NP1>=2&msbb_array19.covariates$Braak>=4)]="AD"

AD_sample.vector=paste(gsub(pattern = "X",replacement = "",msbb_array19.covariates$BrainBank[msbb_array19.covariates$SampleType=="AD"]),collapse = "|")
Control_sample.vector=paste(gsub(pattern = "X",replacement = "",msbb_array19.covariates$BrainBank[msbb_array19.covariates$SampleType=="CONTROL"]),collapse = "|")

#Regroup brain regions by lobes
# lobe_bm_area.map=matrix(nrow=17,ncol=2)
# lobe_bm_area.map[1,]=c("Frontal_Lobe","FP")
# lobe_bm_area.map[2,]=c("Frontal_Lobe","AC")
# lobe_bm_area.map[3,]=c("Frontal_Lobe","PFC")
# lobe_bm_area.map[4,]=c("Occipetal_Lobe","OVC")
# lobe_bm_area.map[5,]=c("Temporal_Lobe","ITG")
# lobe_bm_area.map[6,]=c("Temporal_Lobe","MTG")
# lobe_bm_area.map[7,]=c("Temporal_Lobe","STG")
# lobe_bm_area.map[8,]=c("Parietal_Lobe","PCC")
# lobe_bm_area.map[9,]=c("Temporal_Lobe","PHG")
# lobe_bm_area.map[10,]=c("Temporal_Lobe","TP")
# lobe_bm_area.map[11,]=c("Frontal_Lobe","PCG")
# lobe_bm_area.map[12,]=c("Frontal_Lobe","IFG")
# lobe_bm_area.map[13,]=c("Frontal_Lobe","DLPFC")
# lobe_bm_area.map[14,]=c("Parietal_Lobe","SPL")
# lobe_bm_area.map[16,]=c("Temporal_Lobe","HP")
# lobe_bm_area.map[15,]=c("Dorsal_striatum","PTMN")
# lobe_bm_area.map[17,]=c("Dorsal_striatum","CN")
# lobe_bm_area.map=data.frame(Lobe=lobe_bm_area.map[,1],BrainRegion=lobe_bm_area.map[,2],stringsAsFactors = F)
# 
# 
# msbb_array.byLobe=msbb_array.byLobe2=vector(mode = "list",length = length(unique(lobe_bm_area.map$Lobe)))
# names(msbb_array.byLobe)=names(msbb_array.byLobe2)=unique(lobe_bm_area.map$Lobe)
# msbb_array.byLobe=foreach(i=1:length(names(msbb_array.byLobe)))%dopar%{
#   array_byLobe=do.call("cbind",msbb_array19.2.agg2[which(names(msbb_array19.2.agg2)%in%lobe_bm_area.map$BrainRegion[lobe_bm_area.map$Lobe==names(msbb_array.byLobe)[i]])])  
# } 

msbb_array_region.diffcoexp=vector(mode = "list",length = 17)
names(msbb_array_region.diffcoexp)=names(msbb_array19)
#msbb_array_lobe.diffcoexp=vector(mode = "list",length = 5)
for(i in 6:17){
  msbb_array_region.diffcoexp[[i]]=diffcoexp(exprs.1 = as.data.frame(msbb_array19.2.agg2[[i]])[,grep(pattern = Control_sample.vector,x = colnames(as.data.frame(msbb_array19.2.agg2[[i]])))],
                 exprs.2 = as.data.frame(msbb_array19.2.agg2[[i]])[,grep(pattern = AD_sample.vector,x = colnames(as.data.frame(msbb_array19.2.agg2[[i]])))],
                 rth=0.6, qth=0.2, r.diffth=0.1, q.diffth=0.1)
  saveRDS(msbb_array_region.diffcoexp[[i]],paste(names(msbb_array_region.diffcoexp)[i],"diffcoexp.RDS",sep="_"))
}
names(diffcoexp.res)=names(msbb_array19)
msbb_array.DCG=lapply(diffcoexp.res,function(y)y$DCGs)[lapply(diffcoexp.res,function(y)dim(y$DCGs)[1])>0]
msbb_array.DCL=lapply(diffcoexp.res,function(y)y$DCLs)[lapply(diffcoexp.res,function(y)dim(y$DCLs)[1])>0]
saveRDS(msbb_array.DCG,"msbb_array_DCG.RDS")
saveRDS(msbb_array.DCL,"msbb_array_DCL.RDS")



names(msbb_array.diffcoexp)=names(msbb_array.byLobe)
msbb_array.DRsort=foreach(i=1:17,.export = "DCGL")%dopar%{
  res2=DRsort(DCGs = as.data.frame(msbb_array.DCG[i]),DCLs = as.data.frame(msbb_array.DCL[i]),tf2target = regnet_tf2target,expGenes = rownames(msbb_array19.2.agg2$AC))
}
msbb_array.DRsort=mapply(FUN = function(msbb_array.DCG), function(x)DRsort(DCGs = x$DCGs,DCLs = x$DCLs,tf2target = regnet_tf2target,expGenes = rownames(msbb_array19.2.agg2$AC)))
#Read TF-Target interactions, preprocess data
regnet_tf2target=fread("/shared/hidelab2/user/md4zsa/Work/Data/TF_Databases/RegNetwork_human_regulators2.txt",header = T,sep = "\t",showProgress = T,data.table = F)%>%filter(evidence=="Experimental")%>%select(c(regulator_symbol,target_symbol))



