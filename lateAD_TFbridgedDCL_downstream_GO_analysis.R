library(clusterProfiler)
library(dplyr)
library(org.Hs.eg.db)



setwd("/shared/hidelab2/user/md4zsa/Work/Data/MSBB_Array19/GSE84422/lateAD_diffcoexp/")
common_lateAD_sigTF_regulated_bridgedDCLs.df=read.table("common_lateAD_sigTF_regulated_bridgedDCLs_df.txt",header = T,sep = "\t",as.is = T)
lateAD_TF_bridged_DCL.TF_downstream=lateAD_TF_bridged_DCL.TF_downstream_subgraph=lateAD_TF_bridged_DCL.TF_downstream_subgraph_entrez=lateAD_TF_bridged_DCL.TF_downstream_subgraph_GO=lateAD_TF_bridged_DCL.TF_downstream_subgraph_GO_simplify=lateAD_TF_bridged_DCL.TF_downstream_GOsummaries=vector(mode = "list",length = dim(common_lateAD_sigTF_regulated_bridgedDCLs.df)[2])
names(lateAD_TF_bridged_DCL.TF_downstream)=names(lateAD_TF_bridged_DCL.TF_downstream_subgraph)=names(lateAD_TF_bridged_DCL.TF_downstream_subgraph_entrez)=names(lateAD_TF_bridged_DCL.TF_downstream_subgraph_GO)=names(lateAD_TF_bridged_DCL.TF_downstream_subgraph_GO_simplify)=names(lateAD_TF_bridged_DCL.TF_downstream_GOsummaries)=colnames(common_lateAD_sigTF_regulated_bridgedDCLs.df)

lateAD_TF_bridged_DCL.TF_downstream=readRDS("lateAD_TF_bridged_DCL_TF_downstream.RDS")
lateAD_TF_bridged_DCL.TF_downstream_subgraph=readRDS("lateAD_TF_bridged_DCL_TF_downstream_subgraph.RDS")
lateAD_TF_bridged_DCL.TF_downstream_subgraph_entrez=readRDS("lateAD_TF_bridged_DCL_TF_downstream_subgraph_entrez.RDS")

for(l in 1:length(lateAD_TF_bridged_DCL.TF_downstream_subgraph_GO)){
  cat(paste("Processing TF",names(lateAD_TF_bridged_DCL.TF_downstream_subgraph_GO)[l],"...\n",sep = " "))
  lateAD_TF_bridged_DCL.TF_downstream_subgraph_GO[[l]]=lapply(lateAD_TF_bridged_DCL.TF_downstream_subgraph_entrez[[l]],function(x)enrichGO(gene = x,OrgDb = org.Hs.eg.db,ont = "BP",pvalueCutoff = 0.05,pAdjustMethod = "BH"))
  lateAD_TF_bridged_DCL.TF_downstream_subgraph_GO_simplify[[l]]=lapply(lateAD_TF_bridged_DCL.TF_downstream_subgraph_GO[[l]],function(x)data.frame(clusterProfiler::simplify(x,cutoff=0.8,by="p.adjust",select_fun = min, measure = "Wang")))
  
}
saveRDS(lateAD_TF_bridged_DCL.TF_downstream_subgraph_GO,"lateAD_TF_bridged_DCL_TF_downstream_subgraph_GO.RDS")
saveRDS(lateAD_TF_bridged_DCL.TF_downstream_subgraph_GO_simplify,"lateAD_TF_bridged_DCL_TF_downstream_subgraph_GO_simplify.RDS")
proc.time()
