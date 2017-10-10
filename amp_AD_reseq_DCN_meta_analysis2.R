library(org.Hs.eg.db)
library(igraph)
library(data.table)
library(centiserve)
library(sets)
library(clusterProfiler)
library(dplyr)
library(gdata)


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
setwd("/Users/sandeepamberkar/Work/Data/AMP-AD_RNAseq_ReSeq/Normalised_covariate_corrected_NoResiduals/")

ConsensusDEGs=vector(mode = "list",length = 7)
ConsensusDEGs_df=vector(mode = "list",length = 7)
names(ConsensusDEGs)=names(ConsensusDEGs_df)=c("DLPFC","FP","IFG","PHG","STG","CER","TCX")
ConsensusDEGs_df$DLPFC=fread("../DEG_Analyses/ROSMAP_DLPFC_DiffExpression.tsv",header = T,stringsAsFactors = F,data.table = F) %>%filter(Comparison=="cogdx1-cogdx4")
ConsensusDEGs_df$CER=fread("../DEG_Analyses/MAYO_CBE_TCX_DiffExpression.tsv",header = T,stringsAsFactors = F,data.table = F) %>%filter(BrainRegion=="CER"&Comparison=="CONTROL-AD")
ConsensusDEGs_df$TCX=fread("../DEG_Analyses/MAYO_CBE_TCX_DiffExpression.tsv",header = T,stringsAsFactors = F,data.table = F) %>%filter(BrainRegion=="TCX"&Comparison=="CONTROL-AD")
ConsensusDEGs_df$FP=fread("../DEG_Analyses/MSSM_FP_STG_PHG_IFG_DiffExpression.tsv",header = T,stringsAsFactors = F,data.table = F)%>%filter(BrainRegion=="FP"&Comparison=="CONTROL-AD")
ConsensusDEGs_df$IFG=fread("../DEG_Analyses/MSSM_FP_STG_PHG_IFG_DiffExpression.tsv",header = T,stringsAsFactors = F,data.table = F)%>%filter(BrainRegion=="IFG"&Comparison=="CONTROL-AD")
ConsensusDEGs_df$PHG=fread("../DEG_Analyses/MSSM_FP_STG_PHG_IFG_DiffExpression.tsv",header = T,stringsAsFactors = F,data.table = F)%>%filter(BrainRegion=="PHG"&Comparison=="CONTROL-AD")
ConsensusDEGs_df$STG=fread("../DEG_Analyses/MSSM_FP_STG_PHG_IFG_DiffExpression.tsv",header = T,stringsAsFactors = F,data.table = F)%>%filter(BrainRegion=="STG"&Comparison=="CONTROL-AD")
ConsensusDEGs=lapply(ConsensusDEGs_df,function(x)x%>%filter(Direction=="UP"|Direction=="DOWN")%>%filter(hgnc_symbol!="")%>%select(hgnc_symbol))

amp_ad.MetaDCN=vector(mode = "list",length = 6)
names(amp_ad.MetaDCN)=c("DLPFC","TCX","FP","IFG","PHG","STG")
amp_ad.MetaDCN$DLPFC=V(rosmap_DLPFC_DCN$AD2)$symbol
amp_ad.MetaDCN$FP=V(msmm_DCN$FP$AD2)$name
amp_ad.MetaDCN$IFG=V(msmm_DCN$IFG$AD2)$name
amp_ad.MetaDCN$PHG=V(msmm_DCN$PHG$AD2)$name
amp_ad.MetaDCN$STG=V(msmm_DCN$STG$AD2)$name
amp_ad.MetaDCN$TCX=V(mayo_TCX_DCN$AD2)$name

metaAD_DCLs=metaAD_DCGs=metaAD_DRLs=metaAD_DRGs=metaAD_DRsort_results=vector(mode = "list",length=7)
names(metaAD_DCLs)=names(metaAD_DCGs)=names(metaAD_DRLs)=names(metaAD_DRGs)=names(metaAD_DRsort_results)=c("Mayo_TCX","Mayo_CER","ROSMAP_DLPFC","MSBB_FP","MSBB_IFG","MSBB_PHG","MSBB_STG")
metaAD_DRsort_results$Mayo_TCX=readRDS("MAYO/DCGL_results/TCX_DRsort_Results.RDS")
metaAD_DRsort_results$Mayo_CER=readRDS("MAYO/DCGL_results/CER_DRsort_Results.RDS")
metaAD_DRsort_results$ROSMAP_DLPFC=readRDS("ROSMAP/DCGL_results/ROSMAP_DCe_DRsort_q005.res.RDS")
metaAD_DRsort_results$MSBB_FP=readRDS("MSMM/DCGL_results/MSBB_FP_DCe_DRsort.res.RDS")
metaAD_DRsort_results$MSBB_IFG=readRDS("MSMM/DCGL_results/MSBB_IFG_DCe_DRsort.res.RDS")
metaAD_DRsort_results$MSBB_PHG=readRDS("MSMM/DCGL_results/MSBB_PHG_DCe_DRsort.res.RDS")
metaAD_DRsort_results$MSBB_STG=readRDS("MSMM/DCGL_results/MSBB_STG_DCe_DRsort.res.RDS")

metaAD_DCLs=lapply(metaAD_DRsort_results,function(x)x$DCLs)
metaAD_DCGs=lapply(metaAD_DRsort_results,function(x)x$DCGs)
metaAD_DRLs=lapply(metaAD_DRsort_results,function(x)x$DRLs)
metaAD_DRGs=lapply(metaAD_DRsort_results,function(x)x$DRGs)
metaAD_TF_bridged_DCLs=lapply(metaAD_DRsort_results,function(x)x$TF_bridged_DCL)
metaAD_TFDCLs
function(y)mapIds2(IDs = levels(y),IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2]),function(z)data.frame(enrichKEGG(gene = z,organism="hsa",pvalueCutoff=0.05,pAdjustMethod="BH"))[,2])
