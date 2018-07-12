library(clusterProfiler)
library(dplyr)
library(org.Hs.eg.db)


earlyAD_TF_bridged_DCL.TF_downstream=earlyAD_TF_bridged_DCL.TF_downstream_subgraph=earlyAD_TF_bridged_DCL.TF_downstream_subgraph_entrez=earlyAD_TF_bridged_DCL.TF_downstream_subgraph_GO=earlyAD_TF_bridged_DCL.TF_downstream_subgraph_GO_simplify=earlyAD_TF_bridged_DCL.TF_downstream_GOsummaries=vector(mode = "list",length = dim(common_earlyAD_sigTF_regulated_bridgedDCLs.df)[2])
names(earlyAD_TF_bridged_DCL.TF_downstream)=names(earlyAD_TF_bridged_DCL.TF_downstream_subgraph)=names(earlyAD_TF_bridged_DCL.TF_downstream_subgraph_entrez)=names(earlyAD_TF_bridged_DCL.TF_downstream_subgraph_GO)=names(earlyAD_TF_bridged_DCL.TF_downstream_subgraph_GO_simplify)=names(earlyAD_TF_bridged_DCL.TF_downstream_GOsummaries)=colnames(common_earlyAD_sigTF_regulated_bridgedDCLs.df)

setwd("/shared/hidelab2/user/md4zsa/Work/Data/MSBB_Array19/GSE84422/EarlyAD_diffcoexp/")
common_earlyAD_sigTF_regulated_bridgedDCLs.df=read.table("common_earlyAD_sigTF_regulated_bridgedDCLs_df.txt",header = T,sep = "\t",as.is = T)
earlyAD_TF_bridged_DCL.TF_downstream=readRDS("earlyAD_TF_bridged_DCL_TF_downstream.RDS")
earlyAD_TF_bridged_DCL.TF_downstream_subgraph=readRDS("earlyAD_TF_bridged_DCL_TF_downstream_subgraph.RDS")
earlyAD_TF_bridged_DCL.TF_downstream_subgraph_entrez=readRDS("earlyAD_TF_bridged_DCL_TF_downstream_subgraph_entrez.RDS")

  for(l in 1:length(earlyAD_TF_bridged_DCL.TF_downstream_subgraph_GO)){
    cat(paste("Processing TF",names(earlyAD_TF_bridged_DCL.TF_downstream_subgraph_GO)[l],"...\n",sep = " "))
    earlyAD_TF_bridged_DCL.TF_downstream_subgraph_GO[[l]]=lapply(earlyAD_TF_bridged_DCL.TF_downstream_subgraph_entrez[[l]],function(x)enrichGO(gene = x,OrgDb = org.Hs.eg.db,ont = "BP",pvalueCutoff = 0.01,pAdjustMethod = "BH"))
    earlyAD_TF_bridged_DCL.TF_downstream_subgraph_GO_simplify[[l]]=lapply(earlyAD_TF_bridged_DCL.TF_downstream_subgraph_GO[[l]],function(x)data.frame(clusterProfiler::simplify(x,cutoff=0.8,by="p.adjust",select_fun = min, measure = "Wang")))
    
  }
saveRDS(earlyAD_TF_bridged_DCL.TF_downstream_subgraph_GO,"earlyAD_TF_bridged_DCL_TF_downstream_subgraph_GO.RDS")
saveRDS(earlyAD_TF_bridged_DCL.TF_downstream_subgraph_GO_simplify,"earlyAD_TF_bridged_DCL_TF_downstream_subgraph_GO_simplify.RDS")
proc.time()
  