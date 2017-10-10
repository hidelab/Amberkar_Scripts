library(org.Hs.eg.db)
library(igraph)
library(data.table)
library(Pi)
library(foreach)
library(doParallel)
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
tanzi_selectGenes=scan("Tanzi_select44.txt",what = "char",sep = "\n")
names(tanzi_selectGenes)=mapIds2(IDs = tanzi_selectGenes,IDFrom = "SYMBOL",IDTo = "ENSEMBL")[[1]][,2]
mayo_CER_DCN.AD=readRDS("MAYO/MAYO_ReSeq_CER_DCN_AD.RDS")
mayo_CER_DCN.AD2=delete_vertices(graph = mayo_CER_DCN.AD,v = grep(pattern = paste(mapIds2(IDs = V(mayo_CER_DCN.AD)$name,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[2]],collapse = "|"),x = V(mayo_CER_DCN.AD)$name))
V(mayo_CER_DCN.AD2)$symbol=mapIds2(IDs = V(mayo_CER_DCN.AD2)$name,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[1]][,2]
V(mayo_CER_DCN.AD2)$entrez=mapIds2(IDs = V(mayo_CER_DCN.AD2)$name,IDFrom = "ENSEMBL",IDTo = "ENTREZID")[[1]][,2]
mayo_CER_DCN_AD2.SeedData=data.frame(V(mayo_CER_DCN.AD2)$name,rep(0,length(V(mayo_CER_DCN.AD2)$name)),stringsAsFactors = F)
mayo_CER_DCN_AD2.SeedData[grep(pattern = paste(tanzi_selectGenes,collapse = "|"),x = V(mayo_CER_DCN.AD2)$symbol),2]=1
colnames(mayo_CER_DCN_AD2.SeedData)=c()
mayo_CER_DCN_AD2.pData=xPierGenes(data=mayo_CER_DCN_AD2.SeedData,weighted=T,restart=0.75,network.customised=mayo_CER_DCN.AD2,parallel = F)
mayo_CER_DCN_AD2.pSubnet=xPierSubnet(pNode = mayo_CER_DCN_AD2.pData,network.customised = mayo_CER_DCN.AD2,subnet.significance = 0.05,verbose = T,priority.quantile = NA)
V(mayo_CER_DCN_AD2.pSubnet)$name=mapIds2(IDs = V(mayo_CER_DCN_AD2.pSubnet)$name,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[1]][,2]
V(mayo_CER_DCN_AD2.pSubnet)$significance=as.numeric(V(mayo_CER_DCN_AD2.pSubnet)$significance)
V(mayo_CER_DCN_AD2.pSubnet)$score=as.numeric(V(mayo_CER_DCN_AD2.pSubnet)$score)
V(mayo_CER_DCN_AD2.pSubnet)$priority=as.numeric(V(mayo_CER_DCN_AD2.pSubnet)$priority)
V(mayo_CER_DCN_AD2.pSubnet)$colour[grep(pattern = paste(tanzi_selectGenes,collapse = "|"),x = V(mayo_CER_DCN_AD2.pSubnet)$name)]="#EE4000"
V(mayo_CER_DCN_AD2.pSubnet)$colour[-grep(pattern = paste(tanzi_selectGenes,collapse = "|"),x = V(mayo_CER_DCN_AD2.pSubnet)$name)]="#8B8989"
write.graph(graph = mayo_CER_DCN_AD2.pSubnet,"MAYO_ReSeq_CER_DCN_AD_pSubnet.gml",format = "gml")

mayo_TCX_diffcorr=fread("MAYO/MAYO_ReSeq_TCX_corrPval_DiffCorr_Results_FDR01.txt",sep = "\t",header = T,data.table = F,showProgress = T)
mayo_TCX_DCN.AD=graph.data.frame(d = mayo_TCX_diffcorr[,c(2:3)],directed = F)
mayo_TCX_DCN.AD2=delete_vertices(graph = mayo_TCX_DCN.AD,v = grep(pattern = paste(mapIds2(IDs = V(mayo_TCX_DCN.AD)$name,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[2]],collapse = "|"),x = V(mayo_TCX_DCN.AD)$name))
# V(mayo_TCX_DCN.AD2)$symbol=mapIds2(IDs = V(mayo_TCX_DCN.AD2)$name,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[1]][,2]
V(mayo_TCX_DCN.AD2)$entrez=mapIds2(IDs = V(mayo_TCX_DCN.AD2)$name,IDFrom = "ENSEMBL",IDTo = "ENTREZID")[[1]][,2]
mayo_TCX_DCN_AD2.SeedData=data.frame(V(mayo_TCX_DCN.AD2)$name,rep(0,length(V(mayo_TCX_DCN.AD2)$name)),stringsAsFactors = F)
mayo_TCX_DCN_AD2.SeedData[grep(pattern = paste(tanzi_selectGenes,collapse = "|"),x = V(mayo_TCX_DCN.AD2)$symbol),2]=1
colnames(mayo_TCX_DCN_AD2.SeedData)=c()
mayo_TCX_DCN_AD2.pData=xPierGenes(data=mayo_TCX_DCN_AD2.SeedData,weighted=T,restart=0.75,network.customised=mayo_TCX_DCN.AD2,parallel = F)
mayo_TCX_DCN_AD2.pSubnet=xPierSubnet(pNode = mayo_TCX_DCN_AD2.pData,network.customised = mayo_TCX_DCN.AD2,subnet.significance = 0.05,verbose = T,priority.quantile = NA)
V(mayo_TCX_DCN_AD2.pSubnet)$name=mapIds2(IDs = V(mayo_TCX_DCN_AD2.pSubnet)$name,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[1]][,2]
V(mayo_TCX_DCN_AD2.pSubnet)$significance=as.numeric(V(mayo_TCX_DCN_AD2.pSubnet)$significance)
V(mayo_TCX_DCN_AD2.pSubnet)$score=as.numeric(V(mayo_TCX_DCN_AD2.pSubnet)$score)
V(mayo_TCX_DCN_AD2.pSubnet)$priority=as.numeric(V(mayo_TCX_DCN_AD2.pSubnet)$priority)
V(mayo_TCX_DCN_AD2.pSubnet)$colour[grep(pattern = paste(tanzi_selectGenes,collapse = "|"),x = V(mayo_TCX_DCN_AD2.pSubnet)$name)]="#EE4000"
V(mayo_TCX_DCN_AD2.pSubnet)$colour[-grep(pattern = paste(tanzi_selectGenes,collapse = "|"),x = V(mayo_TCX_DCN_AD2.pSubnet)$name)]="#8B8989"
write.graph(graph = mayo_TCX_DCN_AD2.pSubnet,"MAYO_ReSeq_TCX_DCN_AD_pSubnet.gml",format = "gml")

mayo_CER_DCN.Control=readRDS("MAYO/MAYO_ReSeq_CER_DCN_Control.RDS")
mayo_CER_DCN.Control2=delete_vertices(graph = mayo_CER_DCN.Control,v = grep(pattern = paste(mapIds2(IDs = V(mayo_CER_DCN.Control)$name,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[2]],collapse = "|"),x = V(mayo_CER_DCN.Control)$name))
V(mayo_CER_DCN.Control2)$symbol=mapIds2(IDs = V(mayo_CER_DCN.Control2)$name,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[1]][,2]
V(mayo_CER_DCN.Control2)$entrez=mapIds2(IDs = V(mayo_CER_DCN.Control2)$name,IDFrom = "ENSEMBL",IDTo = "ENTREZID")[[1]][,2]

mayo_TCX_DCN.Control=readRDS("MAYO/MAYO_ReSeq_TCX_DCN_Control.RDS")
mayo_TCX_DCN.Control2=delete_vertices(graph = mayo_TCX_DCN.Control,v = grep(pattern = paste(mapIds2(IDs = V(mayo_TCX_DCN.Control)$name,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[2]],collapse = "|"),x = V(mayo_TCX_DCN.Control)$name))
V(mayo_TCX_DCN.Control2)$symbol=mapIds2(IDs = V(mayo_TCX_DCN.Control2)$name,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[1]][,2]
V(mayo_TCX_DCN.Control2)$entrez=mapIds2(IDs = V(mayo_TCX_DCN.Control2)$name,IDFrom = "ENSEMBL",IDTo = "ENTREZID")[[1]][,2]
####################################################################################################
rosmap_diffcorr=fread("ROSMAP/ROSMAP_ReSeq_DLPFC_corrPval_diffcorr_FDR01.txt",sep = "\t",header = T,showProgress = T,data.table = F)
rosmap_DCN.Control=readRDS("ROSMAP/rosmap_DCN_Control_pValCorr.RDS")
rosmap_DCN.AD=readRDS("ROSMAP/rosmap_DCN_AD_pValCorr.RDS")
# rosmap_DCN_pCorr.Control=readRDS("ROSMAP/rosmap_DCN_Control_pValCorr.RDS")
rosmap_covariates=read.table("ROSMAP/ROSMAP_DLPFC_Covariates.tsv",sep = "\t",header = T,as.is = T)
rosmap_DCN.Control2=delete_vertices(graph = rosmap_DCN.Control,v = grep(pattern = paste(mapIds2(IDs = V(rosmap_DCN.Control)$name,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[2]],collapse = "|"),x = V(rosmap_DCN.Control)$name))
V(rosmap_DCN.Control2)$symbol=mapIds2(IDs = V(rosmap_DCN.Control2)$name,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[1]][,2]
V(rosmap_DCN.Control2)$entrez=mapIds2(IDs = V(rosmap_DCN.Control2)$name,IDFrom = "ENSEMBL",IDTo = "ENTREZID")[[1]][,2]

rosmap_DCN.AD=readRDS("ROSMAP/rosmap_DCN_AD_pValCorr.RDS")
# rosmap_DCN_pCorr.AD=readRDS("ROSMAP/rosmap_DCN_AD_pValCorr.RDS")
rosmap_DCN.AD2=delete_vertices(graph = rosmap_DCN.AD,v = grep(pattern = paste(mapIds2(IDs = V(rosmap_DCN.AD)$name,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[2]],collapse = "|"),x = V(rosmap_DCN.AD)$name))
V(rosmap_DCN.AD2)$symbol=mapIds2(IDs = V(rosmap_DCN.AD2)$name,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[1]][,2]
V(rosmap_DCN.AD2)$entrez=mapIds2(IDs = V(rosmap_DCN.AD2)$name,IDFrom = "ENSEMBL",IDTo = "ENTREZID")[[1]][,2]
V(rosmap_DCN.AD2)$lobby=lobby(graph = rosmap_DCN.AD2,loops = F)
#names(V(rosmap_DCN.AD2)$lobby)=c(rank(-V(rosmap_DCN.AD2)$lobby,ties.method = "first"))
V(rosmap_DCN.AD2)$degree=degree(graph = rosmap_DCN.AD2,loops = F)
#names(V(rosmap_DCN.AD2)$degree)=c(rank(-V(rosmap_DCN.AD2)$degree,ties.method = "first"))

rosmap_DCN_AD2.SeedData=data.frame(V(rosmap_DCN.AD2)$name,rep(0,length(V(rosmap_DCN.AD2)$name)),stringsAsFactors = F)
rosmap_DCN_AD2.SeedData[grep(pattern = paste(tanzi_selectGenes,collapse = "|"),x = V(rosmap_DCN.AD2)$symbol),2]=1
colnames(rosmap_DCN_AD2.SeedData)=c()
rosmap_DCN_AD2.pData=xPierGenes(data=rosmap_DCN_AD2.SeedData,weighted=T,restart=0.75,network.customised=rosmap_DCN.AD2,parallel = T,multicores = 12)
rosmap_DCN_AD2.pData$priority$name=mapIds2(IDs = rosmap_DCN_AD2.pData$priority$name,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[1]][,2]
rosmap_DCN_AD2.pData$priority$description=mapIds2(IDs = rosmap_DCN_AD2.pData$priority$name,IDFrom = "SYMBOL",IDTo = "GENENAME")[[1]][,2]
rosmap_DCN_AD2.pSubnet=xPierSubnet(pNode = rosmap_DCN_AD2.pData,network.customised = rosmap_DCN.AD2,subnet.significance = 0.01,verbose = T)
V(rosmap_DCN_AD2.pSubnet)$name=mapIds2(IDs = V(rosmap_DCN_AD2.pSubnet)$name,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[1]][,2]
V(rosmap_DCN_AD2.pSubnet)$significance=as.numeric(V(rosmap_DCN_AD2.pSubnet)$significance)
V(rosmap_DCN_AD2.pSubnet)$score=as.numeric(V(rosmap_DCN_AD2.pSubnet)$score)
V(rosmap_DCN_AD2.pSubnet)$priority=as.numeric(V(rosmap_DCN_AD2.pSubnet)$priority)
#V(rosmap_DCN_AD2.pSubnet)$colour=c("#8B8989")
V(rosmap_DCN_AD2.pSubnet)$colour[grep(pattern = paste(tanzi_selectGenes,collapse = "|"),x = V(rosmap_DCN_AD2.pSubnet)$name)]="#EE4000"
V(rosmap_DCN_AD2.pSubnet)$colour[-grep(pattern = paste(tanzi_selectGenes,collapse = "|"),x = V(rosmap_DCN_AD2.pSubnet)$name)]="#8B8989"
V(rosmap_DCN_AD2.pSubnet)$lobby=lobby(rosmap_DCN_AD2.pSubnet)
V(rosmap_DCN_AD2.pSubnet)$degree=degree(rosmap_DCN_AD2.pSubnet)
write.graph(graph = rosmap_DCN_AD2.pSubnet,"rosmap_ReSeq_DCN_AD_pSubnet.gml",format = "gml")
rosmap_dcn_degree_df=data.frame(DCN.Degree=degree(graph = rosmap_DCN.AD2,v = which(V(rosmap_DCN.AD2)$symbol%in%V(rosmap_DCN_AD2.pSubnet)$name)),DCN_Tanzi.Degree=degree(graph = rosmap_DCN_AD2.pSubnet,v = which(V(rosmap_DCN_AD2.pSubnet)$name%in%V(rosmap_DCN.AD2)$symbol)),stringsAsFactors = F)
rosmap_dcn_degree_df=data.frame(rosmap_dcn_degree_df[,1],rank(-rosmap_dcn_degree_df[,1],ties.method = "first"),rosmap_dcn_degree_df[,2],rank(-rosmap_dcn_degree_df[,2],ties.method = "first"),stringsAsFactors = F)
colnames(rosmap_dcn_degree_df)=c("DCN.Degree","Rank_DCN.Degree","DCN_Tanzi.Degree","Rank_DCN_Tanzi.Degree")
rosmap_dcn_lobby_df=data.frame(DCN.Lobby=lobby(graph = rosmap_DCN.AD2,vids = which(V(rosmap_DCN.AD2)$symbol%in%V(rosmap_DCN_AD2.pSubnet)$name)),DCN_Tanzi.Lobby=lobby(graph = rosmap_DCN_AD2.pSubnet,vids = which(V(rosmap_DCN_AD2.pSubnet)$name%in%V(rosmap_DCN.AD2)$symbol)),stringsAsFactors = F)
rosmap_dcn_lobby_df=data.frame(rosmap_dcn_lobby_df[,1],rank(-rosmap_dcn_lobby_df[,1],ties.method = "first"),rosmap_dcn_lobby_df[,2],rank(-rosmap_dcn_lobby_df[,2],ties.method = "first"),stringsAsFactors = F)
colnames(rosmap_dcn_lobby_df)=c("DCN.Lobby","Rank_DCN.Lobby","DCN_Tanzi.Lobby","Rank_DCN_Tanzi.Lobby")
#rosmap_dcn_bwn_df=cbind(DCN.Betweenness=betweenness(graph = rosmap_DCN.AD2,v = which(V(rosmap_DCN.AD2)$symbol%in%V(rosmap_DCN_AD2.pSubnet)$name),normalized = T),DCN_Tanzi.Betweenness=betweenness(graph = rosmap_DCN_AD2.pSubnet,v = which(V(rosmap_DCN_AD2.pSubnet)$name%in%V(rosmap_DCN.AD2)$symbol),normalized = T))
rownames(rosmap_dcn_degree_df)=V(rosmap_DCN_AD2.pSubnet)$name
rownames(rosmap_dcn_lobby_df)=V(rosmap_DCN_AD2.pSubnet)$name
#rownames(rosmap_dcn_bwn_df)=V(rosmap_DCN_AD2.pSubnet)$name
write.table(rosmap_dcn_degree_df,"./ROSMAP/ROSMAP_DCN_Degree_df.txt",sep = "\t",col.names = T,row.names = T,quote = F)
write.table(rosmap_dcn_lobby_df,"./ROSMAP/ROSMAP_DCN_Lobby_df.txt",sep = "\t",col.names = T,row.names = T,quote = F)
#write.table(rosmap_dcn_bwn_df,"./ROSMAP/ROSMAP_DCN_Betweenness_df.txt",sep = "\t",col.names = T,row.names = T,quote = F)
rosmap_dcn_topology_df=data.frame(rosmap_dcn_degree_df,rosmap_dcn_lobby_df,stringsAsFactors = F)
rosmap_dcn_topology_df$TanziVariant="N"
rosmap_dcn_topology_df$TanziVariant[rownames(rosmap_dcn_topology_df)%in%tanzi_selectGenes]="Y"
write.table(rosmap_dcn_topology_df,"./ROSMAP/ROSMAP_DCN_pSubnet_Topology_df.txt",sep = "\t",col.names = T,row.names = T,quote = F)
rosmap_DCN_AD2.Topology_df=data.frame(Genes=V(rosmap_DCN.AD2)$symbol,Lobby=V(rosmap_DCN.AD2)$lobby,Lobby.Rank=rank(-V(rosmap_DCN.AD2)$lobby,ties.method = "first"),Degree=V(rosmap_DCN.AD2)$degree,Degree.Rank=rank(-V(rosmap_DCN.AD2)$degree,ties.method = "first"),stringsAsFactors = F)
write.table(rosmap_DCN_AD2.Topology_df,"./ROSMAP/ROSMAP_DCN_Topology_df.txt",sep = "\t",col.names = T,row.names = F,quote = F)
####################################################################################################
msmm_DCN=msmm_DCN.SeedData=msmm_DCN_AD.pData=msmm_DCN_AD.pSubnet=vector(mode = "list",length=4)
names(msmm_DCN)=names(msmm_DCN.SeedData)=names(msmm_DCN_AD.pData)=names(msmm_DCN_AD.pSubnet)=c("FP","IFG","PHG","STG")
msmm_DCN$FP=msmm_DCN$IFG=msmm_DCN$PHG=msmm_DCN$STG=vector(mode = "list",length = 2)
names(msmm_DCN$FP)=names(msmm_DCN$IFG)=names(msmm_DCN$PHG)=names(msmm_DCN$STG)=c("Control","AD")
msmm_fp_diffcorr=read.table("MSMM/MSMM_ReSeq_FP_corrPval_allResults_DiffCorr_FDR01.txt",sep = "\t",header = T,as.is = T)
msmm_ifg_diffcorr=read.table("MSMM/MSMM_ReSeq_IFG_corrPval_allResults_DiffCorr_FDR01.txt",sep = "\t",header = T,as.is = T)
msmm_phg_diffcorr=read.table("MSMM/MSMM_ReSeq_PHG_corrPval_allResults_DiffCorr_FDR01.txt",sep = "\t",header = T,as.is = T)
msmm_stg_diffcorr=read.table("MSMM/MSMM_ReSeq_STG_corrPval_allResults_DiffCorr_FDR01.txt",sep = "\t",header = T,as.is = T)
msmm_DCN$FP$Control=msmm_DCN$FP$AD=graph.data.frame(d = msmm_fp_diffcorr[msmm_fp_diffcorr$FDR<=0.1,c(2:3)],directed = F)
msmm_DCN$IFG$Control=msmm_DCN$IFG$AD=graph.data.frame(d = msmm_ifg_diffcorr[msmm_ifg_diffcorr$FDR<=0.1,c(2:3)],directed = F)
msmm_DCN$PHG$Control=msmm_DCN$PHG$AD=graph.data.frame(d = msmm_phg_diffcorr[msmm_phg_diffcorr$FDR<=0.1,c(2:3)],directed = F)
msmm_DCN$STG$Control=msmm_DCN$STG$AD=graph.data.frame(d = msmm_stg_diffcorr[msmm_stg_diffcorr$FDR<=0.1,c(2:3)],directed = F)
E(msmm_DCN$FP$Control)$weight=msmm_fp_diffcorr[msmm_fp_diffcorr$FDR<=0.1,4]+1
E(msmm_DCN$FP$AD)$weight=msmm_fp_diffcorr[msmm_fp_diffcorr$FDR<=0.1,7]+1
E(msmm_DCN$IFG$Control)$weight=msmm_ifg_diffcorr[msmm_ifg_diffcorr$FDR<=0.1,4]+1
E(msmm_DCN$IFG$AD)$weight=msmm_ifg_diffcorr[msmm_ifg_diffcorr$FDR<=0.1,7]+1
E(msmm_DCN$PHG$Control)$weight=msmm_phg_diffcorr[msmm_phg_diffcorr$FDR<=0.1,4]+1
E(msmm_DCN$PHG$AD)$weight=msmm_phg_diffcorr[msmm_phg_diffcorr$FDR<=0.1,7]+1
E(msmm_DCN$STG$Control)$weight=msmm_stg_diffcorr[msmm_stg_diffcorr$FDR<=0.1,4]+1
E(msmm_DCN$STG$AD)$weight=msmm_stg_diffcorr[msmm_stg_diffcorr$FDR<=0.1,7]+1
for(i in 1:4){
  for(t in 1:2){
    #msmm_DCN[[i]][[t]]=delete_vertices(graph = msmm_DCN[[i]][[t]],v = grep(pattern = paste(mapIds2(IDs = V(msmm_DCN[[i]][[t]])$name,IDFrom = "SYMBOL",IDTo = "SYMBOL")[[2]],collapse = "|"),x = V(msmm_DCN[[i]][[t]])$name))  
    V(msmm_DCN[[i]][[t]])$symbol=mapIds2(IDs = V(msmm_DCN[[i]][[t]])$name,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[1]][,2]  
    #V(msmm_DCN[[i]][[t]])$entrez=mapIds2(IDs = V(msmm_DCN[[i]][[t]])$name,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2]
    V(msmm_DCN[[i]][[t]])$lobby=lobby(graph = msmm_DCN[[i]][[t]],loops = F)
    V(msmm_DCN[[i]][[t]])$degree=degree(graph = msmm_DCN[[i]][[t]],loops = F)
    
  }
  topology.df=data.frame(Genes=V(msmm_DCN[[i]]$AD)$name,
                         Degree=V(msmm_DCN[[i]]$AD)$degree,Lobby=V(msmm_DCN[[i]]$AD)$lobby,
                         Degree.Rank=rank(-V(msmm_DCN[[i]]$AD)$degree,ties.method = "first"),
                         Lobby.Rank=rank(-V(msmm_DCN[[i]]$AD)$lobby,ties.method = "first"),
                         stringsAsFactors = F)
  write.table(topology.df,paste(names(msmm_DCN)[i],"topology_df.txt",sep = "_"),sep = "\t",col.names = T,row.names = F,quote = F)
}
for(i in 1:4){
  msmm_DCN.SeedData[[i]]=data.frame(V(msmm_DCN[[i]][[2]])$name,rep(0,length(V(msmm_DCN[[i]][[2]])$name)),stringsAsFactors = F)
  msmm_DCN.SeedData[[i]][grep(pattern = paste(tanzi_selectGenes,collapse = "|"),x = V(msmm_DCN[[i]][[2]])$symbol),2]=1
  colnames(msmm_DCN.SeedData[[i]])=c()
  msmm_DCN_AD.pData[[i]]=xPierGenes(data=msmm_DCN.SeedData[[i]],weighted=T,restart=0.75,network.customised=msmm_DCN[[i]][[2]],parallel = T,multicores = 12)
  msmm_DCN_AD.pSubnet[[i]]=xPierSubnet(pNode = msmm_DCN_AD.pData[[i]],network.customised = msmm_DCN[[i]][[2]],subnet.significance = 0.05,verbose = T)
  V(msmm_DCN_AD.pSubnet[[i]])$name=mapIds2(IDs = V(msmm_DCN_AD.pSubnet[[i]])$name,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[1]][,2]
  V(msmm_DCN_AD.pSubnet[[i]])$colour="#8B8989"
  V(msmm_DCN_AD.pSubnet[[i]])$colour[grep(pattern = paste(tanzi_selectGenes,collapse = "|"),x = V(msmm_DCN_AD.pSubnet[[i]])$name)]="#EE4000"
  V(msmm_DCN_AD.pSubnet[[i]])$lobby=lobby(msmm_DCN_AD.pSubnet[[i]])
  V(msmm_DCN_AD.pSubnet[[i]])$degree=degree(msmm_DCN_AD.pSubnet[[i]])
  write.graph(graph = msmm_DCN_AD.pSubnet[[i]],paste("msmm",names(msmm_DCN_AD.pSubnet)[i],"ReSeq_DCN_AD_pSubnet.gml",sep = "_"),format = "gml")
  
}

V(msmm_DCN$FP$Control)$entrez=V(msmm_DCN$FP$AD)$entrez=mapIds2(IDs = V(msmm_DCN$FP$Control)$name,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2]
V(msmm_DCN$IFG$Control)$entrez=V(msmm_DCN$IFG$AD)$entrez=mapIds2(IDs = V(msmm_DCN$IFG$Control)$name,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2]
V(msmm_DCN$PHG$Control)$entrez=V(msmm_DCN$PHG$AD)$entrez=mapIds2(IDs = V(msmm_DCN$PHG$Control)$name,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2]
V(msmm_DCN$STG$Control)$entrez=V(msmm_DCN$STG$AD)$entrez=mapIds2(IDs = V(msmm_DCN$STG$Control)$name,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2]




amp_ad.MetaDCN=vector(mode = "list",length = 7)
names(amp_ad.MetaDCN)=c("MAYO_CER","MAYO_TCX","MSMM_FP","MSMM_IFG","MSMM_PHG","MSMM_STG","ROSMAP_DLPFC")
amp_ad.MetaDCN$MAYO_CER=V(mayo_CER_DCN.AD2)$entrez
amp_ad.MetaDCN$MAYO_TCX=V(mayo_TCX_DCN.AD2)$entrez
amp_ad.MetaDCN$MSMM_FP=V(msmm_DCN$FP$AD)$entrez
amp_ad.MetaDCN$MSMM_IFG=V(msmm_DCN$IFG$AD)$entrez
amp_ad.MetaDCN$MSMM_PHG=V(msmm_DCN$PHG$AD)$entrez
amp_ad.MetaDCN$MSMM_STG=V(msmm_DCN$STG$AD)$entrez
amp_ad.MetaDCN$ROSMAP_DLPFC=V(rosmap_DCN.AD2)$entrez


#Pathway DCN for ROSMAP
rosmap_pathway_diffcorr=readRDS("./ROSMAP/ROSMAP_AD_vs_Control_cp_mean_rpkm_corr_change.RDS")
rosmap_pathway_diffcorr2=rbindlist(rosmap_pathway_diffcorr)
rosmap_pathway_diffcorr2$FDR=p.adjust(p = rosmap_pathway_diffcorr2$p.cocor,method = "fdr")
rosmap_pathway_diffcorr2$abs.corr.change=abs(rosmap_pathway_diffcorr2$r.c-rosmap_pathway_diffcorr2$r.t)
rosmap_pathway.DCN=graph.data.frame(d = rosmap_pathway_diffcorr2[which(rosmap_pathway_diffcorr2$FDR<0.1),c(2:3)],directed = F)
rosmap_pathway_top10.DCN=graph.data.frame(d = rosmap_pathway_diffcorr2[which(rosmap_pathway_diffcorr2$FDR<0.1),][order(rosmap_pathway_diffcorr2[which(rosmap_pathway_diffcorr2$FDR<0.1),13],decreasing = T)[1:10],c(2:3)],directed = F)
canonical_pathways=read.table("../c2.cp.v6.0.symbols.gmt",header = F,sep = "\t",fill = T,stringsAsFactors = F)
for(i in 1:length(V(rosmap_pathway.DCN)$name)){
  
  cnp_gene=unname(unlist(canonical_pathways[which(canonical_pathways$V1==V(rosmap_pathway.DCN)$name[i]),-c(1:2)]))
  V(rosmap_pathway.DCN)$olp[i]=paste(intersect(cnp_gene[cnp_gene!=""],V(rosmap_DCN.AD2)$symbol),collapse = ",")
  V(rosmap_pathway.DCN)$jc[i]=round(jaccard(A = cnp_gene[cnp_gene!=""],B = V(rosmap_DCN.AD2)$symbol),digits = 4)
  
}
rosmap_pathway_DCN.topPathways=V(rosmap_pathway.DCN)$name[which(V(rosmap_pathway.DCN)$jc%in%sort(V(rosmap_pathway.DCN)$jc,decreasing = T)[1:10])]
rosmap_pathway_DCN.SeedData=data.frame(V(rosmap_pathway.DCN)$name,rep(0,length(V(rosmap_pathway.DCN)$name)),stringsAsFactors = F)
rosmap_pathway_DCN.SeedData[grep(pattern = paste(rosmap_pathway_DCN.topPathways,collapse = "|"),x = V(rosmap_pathway.DCN)$name),2]=1
colnames(rosmap_pathway_DCN.SeedData)=c()
rosmap_pathway_DCN.SeedData=xPierGenes(data=rosmap_pathway_DCN.SeedData,weighted=T,restart=0.75,network.customised=rosmap_pathway.DCN,parallel = T,multicores = 2)
rosmap_pathway_DCN.pSubnet=xPierSubnet(pNode = rosmap_pathway_DCN.SeedData,network.customised = rosmap_pathway.DCN,subnet.significance = 0.01,verbose = T)
write.graph(rosmap_pathway_DCN.pSubnet,"./ROSMAP/rosmap_ReSeq_pathway_DCN_pSubnet.gml",format = "gml")
write.graph(rosmap_pathway_top10.DCN,"./ROSMAP/rosmap_ReSeq_pathway_DCN_topPathways.gml",format = "gml")
for (j in 1:length(V(rosmap_pathway_DCN.pSubnet))){
  write.graph(graph = induced_subgraph(graph = rosmap_DCN.AD2,vids = which(V(rosmap_DCN.AD2)$symbol%in%unlist(strsplit(x = V(rosmap_pathway.DCN)$olp[which(V(rosmap_pathway.DCN)$name%in%V(rosmap_pathway_DCN.pSubnet)$name[j])],split = ",")))),file = paste("rosmap_pathway_DCN.pSubnet_P",j,"gene_DCN.gml",sep = "_"),format = "gml")
}


############################################################################################################################
#Determine enrichment of variants within 3 hops of top 10 Superhubs
rosmap_DCN_AD.top10_Superhubs=V(rosmap_DCN.AD2)$symbol[which(rank(-V(rosmap_DCN.AD2)$lobby)%in%sort(rank(-V(rosmap_DCN.AD2)$lobby))[1:10])]
rosmap_DCN_AD.top10_Hubs=V(rosmap_DCN.AD2)$symbol[which(rank(-V(rosmap_DCN.AD2)$degree)%in%sort(rank(-V(rosmap_DCN.AD2)$degree))[1:10])]
rosmap_DCN_AD_top10_Superhubs.Hops=matrix(NA,nrow=10,ncol = 6)
rosmap_DCN_AD_top10_Superhubs.Hops2=matrix(NA,nrow=10,ncol = 6)
tanzi_selectGenes.Damaging=scan("../../../Collaborations/Tanzi_WGS/Damaging_exonic_GeneSymbol.txt",what = "character",sep = "\n")
tanzi_selectGenes.Protective=scan("../../../Collaborations/Tanzi_WGS/Protective_exonic_GeneSymbol.txt",what = "character",sep = "\n")
colnames(rosmap_DCN_AD_top10_Superhubs.Hops)=colnames(rosmap_DCN_AD_top10_Superhubs.Hops2)=c(1:6)
rownames(rosmap_DCN_AD_top10_Superhubs.Hops)=rownames(rosmap_DCN_AD_top10_Superhubs.Hops2)=rosmap_DCN_AD.top10_Superhubs
for(j in 1:6){
  for(i in 1:10){
    rosmap_DCN_AD_top10_Superhubs.Hops[i,j]=length(intersect(mapIds2(IDs = names(neighborhood(graph = rosmap_DCN.AD2,order = j,nodes = which(V(rosmap_DCN.AD2)$symbol%in%rosmap_DCN_AD.top10_Superhubs[i]))[[1]]),IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[1]][,2],tanzi_selectGenes))
    rosmap_DCN_AD_top10_Superhubs.Hops2[i,j]=paste("Protective=",length(intersect(mapIds2(IDs = names(neighborhood(graph = rosmap_DCN.AD2,order = j,nodes = which(V(rosmap_DCN.AD2)$symbol%in%rosmap_DCN_AD.top10_Superhubs[i]))[[1]]),IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[1]][,2],tanzi_selectGenes.Protective)),
                                                   ";Damaging=",length(intersect(mapIds2(IDs = names(neighborhood(graph = rosmap_DCN.AD2,order = j,nodes = which(V(rosmap_DCN.AD2)$symbol%in%rosmap_DCN_AD.top10_Superhubs[i]))[[1]]),IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[1]][,2],tanzi_selectGenes.Damaging)),sep = "")
  }
}
write.table(rosmap_DCN_AD_top10_Superhubs.Hops2,"ROSMAP/ROSMAP_ReSeq_DCN_AD_top10_Superhubs_Hops_DamagingProtective.txt",sep = "\t",col.names = T,row.names = T,quote = F)

random_genes.Hops=random_genes_damaging.Hops=random_genes_protective.Hops=vector(mode = "list",length = 6)
names(random_genes.Hops)=names(random_genes_damaging.Hops)=names(random_genes_protective.Hops)=c(1:6)
random_genes.Hops[[1]]=random_genes.Hops[[2]]=random_genes.Hops[[3]]=random_genes.Hops[[4]]=random_genes.Hops[[5]]=random_genes.Hops[[6]]=rep(0,100)
random_genes_damaging.Hops[[1]]=random_genes_damaging.Hops[[2]]=random_genes_damaging.Hops[[3]]=random_genes_damaging.Hops[[4]]=random_genes_damaging.Hops[[5]]=random_genes_damaging.Hops[[6]]=rep(0,100)
random_genes_protective.Hops[[1]]=random_genes_protective.Hops[[2]]=random_genes_protective.Hops[[3]]=random_genes_protective.Hops[[4]]=random_genes_protective.Hops[[5]]=random_genes_protective.Hops[[6]]=rep(0,100)
hs_gene_table=toTable(org.Hs.egSYMBOL)
rand_Superhub.damaging=rand_Superhub.protective=vector(mode = "list",length = 10)
names(rand_Superhub.damaging)=names(rand_Superhub.protective)=rosmap_DCN_AD.top10_Superhubs
rand_Superhub.damaging[[1]]=rand_Superhub.damaging[[2]]=rand_Superhub.damaging[[3]]=rand_Superhub.damaging[[4]]=rand_Superhub.damaging[[5]]=rand_Superhub.damaging[[6]]=rand_Superhub.damaging[[7]]=rand_Superhub.damaging[[8]]=rand_Superhub.damaging[[9]]=rand_Superhub.damaging[[10]]=vector(mode = "list",length = 4)
rand_Superhub.protective[[1]]=rand_Superhub.protective[[2]]=rand_Superhub.protective[[3]]=rand_Superhub.protective[[4]]=rand_Superhub.protective[[5]]=rand_Superhub.protective[[6]]=rand_Superhub.protective[[7]]=rand_Superhub.protective[[8]]=rand_Superhub.protective[[9]]=rand_Superhub.protective[[10]]=vector(mode = "list",length = 4)
for(sb in 1:10){
  for(j in 1:4){
    tmp=tmp.d=tmp.p=c()
    for(n in 1:500){
      
      rand_damaging_genes=sample(x = hs_gene_table$symbol,size = length(tanzi_selectGenes.Damaging),replace = F)
      rand_protective_genes=sample(x = hs_gene_table$symbol,size = length(tanzi_selectGenes.Protective),replace = F)
      #tmp[i,j]=length(intersect(mapIds2(IDs = names(neighborhood(graph = rosmap_DCN.AD2,order = j,nodes = which(V(rosmap_DCN.AD2)$symbol%in%rand_Superhub))[[i]]),IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[1]][,2],tanzi_selectGenes))
      tmp.damaging=length(intersect(mapIds2(IDs = names(neighborhood(graph = rosmap_DCN.AD2,order = j,nodes = which(V(rosmap_DCN.AD2)$symbol%in%rosmap_DCN_AD.top10_Hubs[sb]))[[1]]),IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[1]][,2],rand_damaging_genes))
      tmp.protective=length(intersect(mapIds2(IDs = names(neighborhood(graph = rosmap_DCN.AD2,order = j,nodes = which(V(rosmap_DCN.AD2)$symbol%in%rosmap_DCN_AD.top10_Hubs[sb]))[[1]]),IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[1]][,2],rand_protective_genes))
      tmp.d=append(x = tmp.d,values=tmp.damaging)
      tmp.p=append(x = tmp.p,values=tmp.protective)
    }
    rand_Superhub.damaging[[sb]][[j]]=tmp.d
    rand_Superhub.protective[[sb]][[j]]=tmp.p
  }
  # rand_Superhub.damaging[[sb]][[j]]=random_genes_damaging.Hops[[j]]
  # rand_Superhub.protective[[sb]][[j]]=random_genes_protective.Hops[[j]]
}

dmg1=data.frame(olp=rand_Superhub.damaging$`THAP9-AS1`[[1]])
dmg2=data.frame(olp=rand_Superhub.damaging$`THAP9-AS1`[[2]])
dmg3=data.frame(olp=rand_Superhub.damaging$`THAP9-AS1`[[3]])
dmg4=data.frame(olp=rand_Superhub.damaging$`THAP9-AS1`[[4]])
dmg1$hop="Hop1"
dmg2$hop="Hop2"
dmg3$hop="Hop3"
dmg4$hop="Hop4"
prt1=data.frame(olp=rand_Superhub.protective$`THAP9-AS1`[[1]])
prt2=data.frame(olp=rand_Superhub.protective$`THAP9-AS1`[[2]])
prt3=data.frame(olp=rand_Superhub.protective$`THAP9-AS1`[[3]])
prt4=data.frame(olp=rand_Superhub.protective$`THAP9-AS1`[[4]])
prt1$hop="Hop1"
prt2$hop="Hop2"
prt3$hop="Hop3"
prt4$hop="Hop4"

par(mfrow=c(2,2))
mtext(text ="THAP9-AS1, Overlap distribution w/ Damaging variants",outer = T,side = 3,cex = 1.5)
hist(rand_Superhub.damaging$`THAP9-AS1`[[1]],col="blue1",xlab = "#Overlaps",main = "Hop1",xlim = c(0,10),cex.lab=1.25)
abline(v=6,col="red")
text(x = 7,y = 20,labels = "pval=0.01")
hist(rand_Superhub.damaging$`THAP9-AS1`[[2]],col="blue1",xlab = "#Overlaps",main = "Hop2",xlim = c(0,40),cex.lab=1.25)
abline(v=34,col="red")
text(x = 30,y = 15,labels = "pval=5.67e-08")
hist(rand_Superhub.damaging$`THAP9-AS1`[[3]],col="blue1",xlab = "#Overlaps",main = "Hop3",xlim=c(0,120),cex.lab=1.25)
abline(v=103,col="red")
text(x = 90,y = 15,labels = "pval=6.04e-05")
hist(rand_Superhub.damaging$`THAP9-AS1`[[4]],col="blue1",xlab = "#Overlaps",main = "Hop4",xlim=c(0,140),cex.lab=1.25)
abline(v=131,col="red")
text(x = 120,y = 15,labels = "pval=0.013")

par(mfrow=c(2,2))
plot(x = 1:500,y = rand_Superhub.damaging$`THAP9-AS1`[[1]],ylim = c(0,140),ylab="#Overlaps",xlab="#Iteractions",main="THAP9-AS1,Hop1")
lines(lowess(x = 1:500,y = rand_Superhub.damaging$`THAP9-AS1`[[1]]),col="blue")
points(400,6,col="red",pch=19)
text(400,4,labels = "Damaging variants",cex = 0.75)
points(400,5.5,col="green",pch=19)
text(390,3.5,labels = "Protective variants",cex = 0.75)

plot(x = 1:500,y = rand_Superhub.damaging$`THAP9-AS1`[[2]],ylim = c(0,140),ylab="#Overlaps",xlab="#Iteractions",main="THAP9-AS1,Hop2")
lines(lowess(x = 1:500,y = rand_Superhub.damaging$`THAP9-AS1`[[2]]),col="blue")
points(400,34,col="red",pch=19,)
#text(400,30,labels = "Damaging variants",cex = 0.75)
points(400,43,col="green",pch=19)
text(400,47,labels = "Protective variants",cex = 0.75)

plot(x = 1:500,y = rand_Superhub.damaging$`THAP9-AS1`[[3]],ylim = c(0,130),ylab="#Overlaps",xlab="#Iterations",cex=1.5)
lines(lowess(x = 1:500,y = rand_Superhub.damaging$`THAP9-AS1`[[3]]),col="blue")
points(400,103,col="red",pch=19,cex=1.5)
text(400,108,labels = "Damaging variants",cex = 0.75)
points(400,88,col="green",pch=19,cex=1.5)
text(400,84,labels = "Protective variants",cex = 0.75)
title(main = "Overlap distribution of Damaging/Protective variants",cex=1.5)
mtext(text = "in neighbourhood at 3rd Hop of THAP9-AS1",outer = F)

plot(x = 1:500,y = rand_Superhub.damaging$`THAP9-AS1`[[4]],ylim = c(0,140),ylab="#Overlaps",xlab="#Iteractions",main="THAP9-AS1,Hop4")
lines(lowess(x = 1:500,y = rand_Superhub.damaging$`THAP9-AS1`[[4]]),col="blue")
points(400,131,col="red",pch=19)
text(400,135,labels = "Damaging variants",cex = 0.75)
points(400,103,col="green",pch=19)
text(400,99,labels = "Protective variants",cex = 0.75)
title(main = "THAP9-AS1, neighbourhood enrichment \n with Damaging and Protective variants",line = -3,outer = T)

plot(x = 1:2000,y = unique(unname(unlist(rand_Superhub.damaging$`THAP9-AS1`))),ylim=c(0,140),ylab="#Overlaps",xlab="#Iteractions")
title("THAP9-AS1, neighbourhood enrichment",cex=1.5)
mtext("with Damaging and Protective variants")
points(400,6,col="red",pch=19)
text(400,3,labels = "Hop1,pval=0.01",cex = 0.75)
points(400,5.5,col="green",pch=19)
points(400,34,col="red",pch=19)
text(400,38,labels = "Hop2,pval=5.67e-08",cex = 0.75)
points(400,43,col="green",pch=19)
points(400,103,col="red",pch=19)
text(400,95,labels = "Hop3,pval=6.04e-05",cex = 0.75)
points(400,88,col="green",pch=19)
points(400,131,col="red",pch=19)
text(400,120,labels = "Hop4,pval=0.01",cex = 0.75)
points(400,108,col="green",pch=19)

random_genes.Hops=lapply(random_genes.Hops,function(x){x<-x[-1];x})
random_genes_protective.Hops=lapply(random_genes_protective.Hops,function(x){x<-x[-1];x})
random_genes_damaging.Hops=lapply(random_genes_damaging.Hops,function(x){x<-x[-1];x})











rosmap_DCN_AD_SPtoTanziGenes=shortest.paths(graph = rosmap_DCN.AD2,v = which(V(rosmap_DCN.AD2)$symbol%in%V(rosmap_DCN.AD2)$symbol),to = which(V(rosmap_DCN.AD2)$symbol%in%tanzi_selectGenes))
rosmap_DCN_AD_top10Superhubs_SPtoTanziGenes=shortest.paths(graph = rosmap_DCN.AD2,v = which(V(rosmap_DCN.AD2)$symbol%in%rosmap_DCN_AD.top10_Superhubs),to = which(V(rosmap_DCN.AD2)$symbol%in%tanzi_selectGenes))
rosmap_DCN_AD_top10Superhubs_SPtoTanziDamaging=shortest.paths(graph = rosmap_DCN.AD2,v = which(V(rosmap_DCN.AD2)$symbol%in%rosmap_DCN_AD.top10_Superhubs),to = which(V(rosmap_DCN.AD2)$symbol%in%tanzi_selectGenes.Damaging))
rosmap_DCN_AD_top10Superhubs_SPtoTanziDamaging=rosmap_DCN_AD_top10Superhubs_SPtoTanziDamaging[,-which(colSums(rosmap_DCN_AD_top10Superhubs_SPtoTanziDamaging)=='Inf')]
rosmap_DCN_AD_top10Superhubs_SPtoTanziProtective=shortest.paths(graph = rosmap_DCN.AD2,v = which(V(rosmap_DCN.AD2)$symbol%in%rosmap_DCN_AD.top10_Superhubs),to = which(V(rosmap_DCN.AD2)$symbol%in%tanzi_selectGenes.Protective))
rosmap_DCN_AD_top10Superhubs_SPtoTanziProtective=rosmap_DCN_AD_top10Superhubs_SPtoTanziProtective[,-which(colSums(rosmap_DCN_AD_top10Superhubs_SPtoTanziProtective)=='Inf')]
rosmap_DCN_AD_top10Superhubs_SPtoTanziDamaging.df=data.frame(t(rosmap_DCN_AD_top10Superhubs_SPtoTanziDamaging),stringsAsFactors = F)
rosmap_DCN_AD_top10Superhubs_SPtoTanziProtective.df=data.frame(t(rosmap_DCN_AD_top10Superhubs_SPtoTanziProtective),stringsAsFactors = F)
rosmap_DCN_AD_top10Superhubs_SPtoTanziDamaging.df$variantType='Damaging/Predisposed'
rosmap_DCN_AD_top10Superhubs_SPtoTanziProtective.df$variantType='Protective/Resistant'
rosmap_DCN_AD_top10Superhubs_SPtoTanziDamaging_Protective.df=rbind(rosmap_DCN_AD_top10Superhubs_SPtoTanziDamaging.df,rosmap_DCN_AD_top10Superhubs_SPtoTanziProtective.df)

rosmap_DCN_AD_SPtoTanziGenes2=apply(rosmap_DCN_AD_SPtoTanziGenes,1,function(x){x[which(x=='Inf')]<-0;x})
rosmap_DCN_AD_top10Superhubs_SPtoTanziGenes[,20]=0
colnames(rosmap_DCN_AD_SPtoTanziGenes2)=mapIds2(IDs = colnames(rosmap_DCN_AD_SPtoTanziGenes2),IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[1]][,2]
rownames(rosmap_DCN_AD_SPtoTanziGenes2)=mapIds2(IDs = rownames(rosmap_DCN_AD_SPtoTanziGenes2),IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[1]][,2]
############################################################################################################################



write.table(as_edgelist(graph = mayo_CER_DCN.AD2,names = F),"./MAYO/MAYO_ReSeq_CER_DCN_AD_edgelist.tsv",sep = "\t",col.names = F,row.names = F,quote = F)
write.table(data.frame(get.edgelist(graph = mayo_CER_DCN.AD2),stringsAsFactors = F),"./MAYO/MAYO_ReSeq_CER_DCN_AD_df.txt",sep = "\t",col.names = F,row.names = F,quote = F)
write.table(as_edgelist(graph = mayo_TCX_DCN.AD2,names = F),"./MAYO/MAYO_ReSeq_TCX_DCN_AD_edgelist.tsv",sep = "\t",col.names = F,row.names = F,quote = F)
write.table(data.frame(get.edgelist(graph = mayo_TCX_DCN.AD2),stringsAsFactors = F),"./MAYO/MAYO_ReSeq_TCX_DCN_AD_df.txt",sep = "\t",col.names = F,row.names = F,quote = F)
write.table(as_edgelist(graph = msmm_DCN$FP$AD,names = F),"./MSMM/MSMM_ReSeq_FP_DCN_AD_edgelist.tsv",sep = "\t",col.names = F,row.names = F,quote = F)
write.table(data.frame(get.edgelist(graph = mayo_CER_DCN.AD2),stringsAsFactors = F),"./MAYO/MAYO_ReSeq_CER_DCN_AD_df.txt",sep = "\t",col.names = F,row.names = F,quote = F)
write.table(as_edgelist(graph = msmm_DCN$IFG$AD,names = F),"./MSMM/MSMM_ReSeq_IFG_DCN_AD_edgelist.tsv",sep = "\t",col.names = F,row.names = F,quote = F)
write.table(data.frame(get.edgelist(graph = mayo_CER_DCN.AD2),stringsAsFactors = F),"./MAYO/MAYO_ReSeq_CER_DCN_AD_df.txt",sep = "\t",col.names = F,row.names = F,quote = F)
write.table(as_edgelist(graph = msmm_DCN$PHG$AD,names = F),"./MSMM/MSMM_ReSeq_PHG_DCN_AD_edgelist.tsv",sep = "\t",col.names = F,row.names = F,quote = F)
write.table(data.frame(get.edgelist(graph = mayo_CER_DCN.AD2),stringsAsFactors = F),"./MAYO/MAYO_ReSeq_CER_DCN_AD_df.txt",sep = "\t",col.names = F,row.names = F,quote = F)
write.table(as_edgelist(graph = msmm_DCN$STG$AD,names = F),"./MSMM/MSMM_ReSeq_STG_DCN_AD_edgelist.tsv",sep = "\t",col.names = F,row.names = F,quote = F)
write.table(data.frame(get.edgelist(graph = mayo_CER_DCN.AD2),stringsAsFactors = F),"./MAYO/MAYO_ReSeq_CER_DCN_AD_df.txt",sep = "\t",col.names = F,row.names = F,quote = F)
write.table(as_edgelist(graph = rosmap_DCN.AD2,names = F),"./ROSMAP/ROSMAP_ReSeq_DLPFC_DCN_AD_edgelist.tsv",sep = "\t",col.names = F,row.names = F,quote = F)
write.table(data.frame(get.edgelist(graph = rosmap_DCN.AD2),stringsAsFactors = F),"./MAYO/MAYO_ReSeq_CER_DCN_AD_df.txt",sep = "\t",col.names = F,row.names = F,quote = F)

write.table(data.frame(Indices=as.numeric(V(mayo_CER_DCN.AD2)[V(mayo_CER_DCN.AD2)$name]),Genes=V(mayo_CER_DCN.AD2)$name),"./MAYO/MAYO_ReSeq_CER_DCN_AD_gene_indices.txt",sep = "\t",col.names = F,row.names = F,quote = F)
write.table(data.frame(Indices=as.numeric(V(mayo_TCX_DCN.AD2)[V(mayo_TCX_DCN.AD2)$name]),Genes=V(mayo_TCX_DCN.AD2)$name),"./MAYO/MAYO_ReSeq_TCX_DCN_AD_gene_indices.txt",sep = "\t",col.names = F,row.names = F,quote = F)
write.table(data.frame(Indices=as.numeric(V(msmm_DCN$FP$AD)[V(msmm_DCN$FP$AD)$name]),Genes=V(msmm_DCN$FP$AD)$name),"./MSMM/MSMM_ReSeq_FP_DCN_AD_gene_indices.txt",sep = "\t",col.names = F,row.names = F,quote = F)
write.table(data.frame(Indices=as.numeric(V(msmm_DCN$IFG$AD)[V(msmm_DCN$IFG$AD)$name]),Genes=V(msmm_DCN$IFG$AD)$name),"./MSMM/MSMM_ReSeq_IFG_DCN_AD_gene_indices.txt",sep = "\t",col.names = F,row.names = F,quote = F)
write.table(data.frame(Indices=as.numeric(V(msmm_DCN$PHG$AD)[V(msmm_DCN$PHG$AD)$name]),Genes=V(msmm_DCN$PHG$AD)$name),"./MSMM/MSMM_ReSeq_PHG_DCN_AD_gene_indices.txt",sep = "\t",col.names = F,row.names = F,quote = F)
write.table(data.frame(Indices=as.numeric(V(msmm_DCN$STG$AD)[V(msmm_DCN$STG$AD)$name]),Genes=V(msmm_DCN$STG$AD)$name),"./MSMM/MSMM_ReSeq_STG_DCN_AD_gene_indices.txt",sep = "\t",col.names = F,row.names = F,quote = F)
write.table(data.frame(Indices=as.numeric(V(rosmap_DCN.AD2)[V(rosmap_DCN.AD2)$name]),Genes=V(rosmap_DCN.AD2)$name),"./ROSMAP/ROSMAP_ReSeq_DLPFC_DCN_AD_indices.txt",sep = "\t",col.names = F,row.names = F,quote = F)

#Determine functional significance of AMP-AD Consensus module overlap with DCNs
ConsensusModules=vector(mode = "list",length = 7)
ConsensusModules_df=vector(mode = "list",length = 7)
names(ConsensusModules)=names(ConsensusModules_df)=c("DLPFC","FP","IFG","PHG","STG","CBE","TCX")
ConsensusModules_df$DLPFC=read.csv("../ConsensusModules/consensusDLPFC.csv",header = T,as.is = T,stringsAsFactors = F)
ConsensusModules_df$CBE=read.csv("../ConsensusModules/consensusCBE.csv",header = T,as.is = T,stringsAsFactors = F)
ConsensusModules_df$TCX=read.csv("../ConsensusModules/consensusTCX.csv",header = T,as.is = T,stringsAsFactors = F)
ConsensusModules_df$FP=read.csv("../ConsensusModules/consensusFP.csv",header = T,as.is = T,stringsAsFactors = F)
ConsensusModules_df$IFG=read.csv("../ConsensusModules/consensusIFG.csv",header = T,as.is = T,stringsAsFactors = F)
ConsensusModules_df$PHG=read.csv("../ConsensusModules/consensusPHG.csv",header = T,as.is = T,stringsAsFactors = F)
ConsensusModules_df$STG=read.csv("../ConsensusModules/consensusSTG.csv",header = T,as.is = T,stringsAsFactors = F)



Consensus_UniqueModules=vector(mode = "list",length = 7)
names(Consensus_UniqueModules)=names(ConsensusModules)
Consensus_UniqueModules$DLPFC=unique(con_DLPFC$moduleNumber)
Consensus_UniqueModules$CBE=unique(con_CBE$moduleNumber)
Consensus_UniqueModules$TCX=unique(con_TCX$moduleNumber)
Consensus_UniqueModules$FP=unique(con_FP$moduleNumber)
Consensus_UniqueModules$IFG=unique(con_IFG$moduleNumber)
Consensus_UniqueModules$PHG=unique(con_PHG$moduleNumber)
Consensus_UniqueModules$STG=unique(con_STG$moduleNumber)

ConsensusModules$DLPFC=vector(mode = "list",length = length(modules_DPLFC))
ConsensusModules$CBE=vector(mode = "list",length = length(modules_CBE))
ConsensusModules$TCX=vector(mode = "list",length = length(modules_TCX))
ConsensusModules$FP=vector(mode = "list",length = length(modules_FP))
ConsensusModules$IFG=vector(mode = "list",length = length(modules_IFG))
ConsensusModules$PHG=vector(mode = "list",length = length(modules_PHG))
ConsensusModules$STG=vector(mode = "list",length = length(modules_STG))

for(i in 1:7){
  for(j in 1:length(Consensus_UniqueModules[[i]]))
  ConsensusModules[[i]][[j]]=ConsensusModules_df[[i]]$Gene.ID[which(ConsensusModules_df[[i]]$moduleNumber==Consensus_UniqueModules[[i]][j])]
  ConsensusModules[[i]]=lapply(lapply(lapply(ConsensusModules[[i]],function(x)mapIds2(IDs = x,IDFrom = "ENSEMBL",IDTo = "SYMBOL")),`[[`,1),`[[`,2)
  names(ConsensusModules[[i]])=Consensus_UniqueModules[[i]]
}

#Mt.Sinai DCNs with corr. p-values
ConsensusModules_Entrez=ConsensusModules_OlpPathways=ConsensusModules_OlpJaccard=vector(mode = "list",length = 7)
names(ConsensusModules_Entrez)=names(ConsensusModules_OlpJaccard)=names(ConsensusModules_OlpPathways)=names(ConsensusModules)
ConsensusModules_Entrez$DLPFC=ConsensusModules_Entrez$FP=ConsensusModules_Entrez$IFG=ConsensusModules_Entrez$PHG=ConsensusModules_Entrez$STG=ConsensusModules_Entrez$CBE=ConsensusModules_Entrez$TCX=list()
ConsensusModules_OlpPathways$DLPFC=ConsensusModules_OlpPathways$FP=ConsensusModules_OlpPathways$IFG=ConsensusModules_OlpPathways$PHG=ConsensusModules_OlpPathways$STG=ConsensusModules_OlpPathways$CBE=ConsensusModules_OlpPathways$TCX=list()
for(i in c(1:5,7)){
  ConsensusModules_Entrez[[i]]=lapply(ConsensusModules[[i]],function(x)mapIds2(IDs = x,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2])
  
  
}
#Module_entrez=lapply(ConsensusModules$DLPFC,function(x)mapIds2(IDs = x,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2])
ConsensusModules_OlpJaccard$DLPFC=lapply(ConsensusModules$DLPFC,jaccard,V(rosmap_DCN.AD2)$symbol)
ConsensusModules_OlpJaccard$FP=lapply(ConsensusModules$FP,jaccard,V(msmm_DCN$FP$AD)$name)
ConsensusModules_OlpJaccard$IFG=lapply(ConsensusModules$IFG,jaccard,V(msmm_DCN$IFG$AD)$name)
ConsensusModules_OlpJaccard$PHG=lapply(ConsensusModules$PHG,jaccard,V(msmm_DCN$PHG$AD)$name)
ConsensusModules_OlpJaccard$STG=lapply(ConsensusModules$STG,jaccard,V(msmm_DCN$STG$AD)$name)
ConsensusModules_OlpJaccard$TCX=lapply(ConsensusModules$TCX,jaccard,V(mayo_TCX_DCN.AD2)$name)

ConsensusModules_OlpPathways$DLPFC=lapply(ConsensusModules_Entrez$DLPFC,function(y)summary(enrichKEGG(gene = y,organism = "hsa",pvalueCutoff = 0.05,universe = V(rosmap_DCN.AD2)$entrez,pAdjustMethod = "BH"))[,2])
ConsensusModules_OlpPathways$FP=lapply(ConsensusModules_Entrez$FP,function(y)summary(enrichKEGG(gene = y,organism = "hsa",pvalueCutoff = 0.05,universe = V(msmm_DCN$FP$AD)$entrez,pAdjustMethod = "BH"))[,2])
ConsensusModules_OlpPathways$IFG=lapply(ConsensusModules_Entrez$IFG,function(y)summary(enrichKEGG(gene = y,organism = "hsa",pvalueCutoff = 0.05,universe = V(msmm_DCN$IFG$AD)$entrez,pAdjustMethod = "BH"))[,2])
ConsensusModules_OlpPathways$PHG=lapply(ConsensusModules_Entrez$PHG,function(y)summary(enrichKEGG(gene = y,organism = "hsa",pvalueCutoff = 0.05,universe = V(msmm_DCN$PHG$AD)$entrez,pAdjustMethod = "BH"))[,2])
ConsensusModules_OlpPathways$STG=lapply(ConsensusModules_Entrez$STG,function(y)summary(enrichKEGG(gene = y,organism = "hsa",pvalueCutoff = 0.05,universe = V(msmm_DCN$STG$AD)$entrez,pAdjustMethod = "BH"))[,2])
#ConsensusModules_OlpPathways$CBE=lapply(ConsensusModules_Entrez$CBE,function(y)summary(enrichKEGG(gene = y,organism = "hsa",pvalueCutoff = 0.05,universe = V()$entrez,pAdjustMethod = "BH"))[,2])
ConsensusModules_OlpPathways$TCX=lapply(ConsensusModules_Entrez$TCX,function(y)summary(enrichKEGG(gene = y,organism = "hsa",pvalueCutoff = 0.05,universe = V(mayo_TCX_DCN.AD2)$entrez,pAdjustMethod = "BH"))[,2])

ConsensusModules_OlpPathways_Concat=vector(mode = "list",length = 7)
names(ConsensusModules_OlpPathways_Concat)=names(ConsensusModules)
for(i in c(1:5,7)){
  ConsensusModules_OlpPathways_Concat[[i]]=lapply(ConsensusModules_OlpPathways[[i]],paste,collapse=";")
}

pval_vec=vector(mode = "list",length = 7)
names(pval_vec)=names(ConsensusModules)
pval_vec$DLPFC=pval_vec$FP=pval_vec$IFG=pval_vec$PHG=pval_vec$STG=pval_vec$CBE=pval_vec$TCX=list()
for(i in 1:length(ConsensusModules$DLPFC)){
  pval_vec$DLPFC[[i]]=phyper(lapply(ConsensusModules_Entrez$DLPFC,function(x)length(intersect(x,V(rosmap_DCN.AD2)$entrez)))[[i]],
                        lapply(ConsensusModules_Entrez$DLPFC,length)[[i]],
                        19000,
                        length(V(rosmap_DCN.AD2)$entrez),lower.tail = F)
  
}
for(i in 1:length(ConsensusModules$FP)){
  pval_vec$FP[[i]]=phyper(lapply(ConsensusModules_Entrez$FP,function(x)length(intersect(x,V(msmm_DCN$FP$AD)$entrez)))[[i]],
                     lapply(ConsensusModules_Entrez$FP,length)[[i]],
                     19000,
                     length(V(msmm_DCN$FP$AD)$entrez),lower.tail = F)
  
}
for(i in 1:length(ConsensusModules$IFG)){
  pval_vec$IFG[[i]]=phyper(lapply(ConsensusModules_Entrez$IFG,function(x)length(intersect(x,V(msmm_DCN$IFG$AD)$entrez)))[[i]],
                      lapply(ConsensusModules_Entrez$IFG,length)[[i]],
                      19000,
                      length(V(msmm_DCN$FP$AD)$entrez),lower.tail = F)
  
}
for(i in 1:length(ConsensusModules$PHG)){
  pval_vec$PHG[[i]]=phyper(lapply(ConsensusModules_Entrez$PHG,function(x)length(intersect(x,V(msmm_DCN$PHG$AD)$entrez)))[[i]],
                      lapply(ConsensusModules_Entrez$PHG,length)[[i]],
                      19000,
                      length(V(msmm_DCN$PHG$AD)$entrez),lower.tail = F)
  
}
for(i in 1:length(ConsensusModules$STG)){
  pval_vec$STG[[i]]=phyper(lapply(ConsensusModules_Entrez$STG,function(x)length(intersect(x,V(msmm_DCN$STG$AD)$entrez)))[[i]],
                      lapply(ConsensusModules_Entrez$STG,length)[[i]],
                      19000,
                      length(V(msmm_DCN$STG$AD)$entrez),lower.tail = F)
  
}
for(i in 1:length(ConsensusModules$TCX)){
  pval_vec$TCX[[i]]=phyper(lapply(ConsensusModules_Entrez$TCX,function(x)length(intersect(x,V(mayo_TCX_DCN.AD2)$entrez)))[[i]],
                      lapply(ConsensusModules_Entrez$TCX,length)[[i]],
                      19000,
                      length(V(mayo_TCX_DCN.AD2)$entrez),lower.tail = F)
  
}

ConsensusModules_Olp_Summary=vector(mode = "list",length=7)
names(ConsensusModules_Olp_Summary)=names(ConsensusModules)
for(i in c(1:5,7)){
  tmp=matrix(NA,nrow=length(Consensus_UniqueModules[[i]]),ncol = 5)
  for(j in 1:dim(tmp)[1]){
    tmp[j,]=cbind(ConsensusModuleNumber=Consensus_UniqueModules[[i]][j],
                  Consensus_Module_NrGenes=length(unname(unlist(ConsensusModules_Entrez[[i]][j]))),
                  if(length(unname(unlist(ConsensusModules_OlpPathways[[i]][j])))==0){
                    ConsensusModule_Olp_Pathways=0
                   }
                  else
                  ConsensusModule_Olp_Pathways=length(unname(unlist(ConsensusModules_OlpPathways[[i]][j]))),
                  ConsensusModule_Olp_PathwayNames=unname(unlist(ConsensusModules_OlpPathways_Concat[[i]][j])),
                  Jaccard=unname(unlist(ConsensusModules_OlpJaccard[[i]]))[j])
                  #Olp_pval=unlist(pval_vec[[i]][j]))
  }
  colnames(tmp)=c("Consensus Module Number","Genes in Consensus Module","Enriched KEGG pathways in DCN+Module overlap","KEGG pathway terms","Jaccard coefficient")
  tmp=data.frame(tmp,stringsAsFactors = F)
  tmp$Genes.in.Consensus.Module=as.numeric(tmp$Genes.in.Consensus.Module)
  tmp$Rank_Olp_Pathways=rank(-as.numeric(tmp$Enriched.KEGG.pathways.in.DCN.Module.overlap))
  tmp$Rank_Jaccard=rank(-as.numeric(tmp$Jaccard.coefficient))
  tmp$MeanRank=apply(tmp[,c(6:7)],1,mean)
  tmp$Olp_pval=unlist(pval_vec[[i]])
  ConsensusModules_Olp_Summary[[i]]=tmp
  write.table(ConsensusModules_Olp_Summary[[i]],paste(names(ConsensusModules_Olp_Summary)[i],"ConsensusModule_Overlap_Pathway.txt",sep = "_"),sep = "\t",col.names = T,row.names = F,quote = F)
}

#Determine functional significance of AMP-AD Consensus DEGs overlap with DCNs

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

pval_vec_deg=vector(mode = "list",length = 7)
names(pval_vec_deg)=names(ConsensusDEGs)
pval_vec_deg$DLPFC=pval_vec_deg$FP=pval_vec_deg$IFG=pval_vec_deg$PHG=pval_vec_deg$STG=pval_vec_deg$CER=pval_vec_deg$TCX=list()

pval_vec_deg$DLPFC=phyper(length(intersect(V(rosmap_DCN.AD2)$name,ConsensusDEGs$DLPFC$hgnc_symbol)),length(ConsensusDEGs$DLPFC$hgnc_symbol),19000,length(V(rosmap_DCN.AD2)$name),lower.tail = F)
pval_vec_deg$FP=phyper(length(intersect(V(msmm_DCN$FP$AD)$name,ConsensusDEGs$FP$hgnc_symbol)),length(ConsensusDEGs$FP$hgnc_symbol),19000,length(V(msmm_DCN$FP$AD)$name),lower.tail = F)
pval_vec_deg$IFG=phyper(length(intersect(V(msmm_DCN$IFG$AD)$name,ConsensusDEGs$IFG$hgnc_symbol)),length(ConsensusDEGs$IFG$hgnc_symbol),19000,length(V(msmm_DCN$IFG$AD)$name),lower.tail = F)
pval_vec_deg$PHG=phyper(length(intersect(V(msmm_DCN$PHG$AD)$name,ConsensusDEGs$PHG$hgnc_symbol)),length(ConsensusDEGs$PHG$hgnc_symbol),19000,length(V(msmm_DCN$PHG$AD)$name),lower.tail = F)
pval_vec_deg$STG=phyper(length(intersect(V(msmm_DCN$STG$AD)$name,ConsensusDEGs$STG$hgnc_symbol)),length(ConsensusDEGs$STG$hgnc_symbol),19000,length(V(msmm_DCN$STG$AD)$name),lower.tail = F)
pval_vec_deg$TCX=phyper(length(intersect(V(mayo_TCX_DCN.AD)$name,ConsensusDEGs$TCX$hgnc_symbol)),length(ConsensusDEGs$TCX$hgnc_symbol),19000,length(V(mayo_TCX_DCN.AD)$name),lower.tail = F)
pval_vec_deg$CER=NA

ConsensusDEGs_OlpPathways=vector(mode = "list",length = 7)
names(ConsensusDEGs_OlpPathways)=names(ConsensusDEGs)
ConsensusDEGs_OlpPathways=lapply(lapply(lapply(ConsensusDEGs,function(x)mapIds2(IDs = x$hgnc_symbol,IDFrom="SYMBOL",IDTo="ENTREZID")[[1]][,2]),function(y)data.frame(enrichKEGG(gene = y,organism = "hsa",pvalueCutoff = 0.05,pAdjustMethod = "BH"))),`[[`,2)

amp_ad.MetaDCN=vector(mode = "list",length = 7)
names(amp_ad.MetaDCN)=names(ConsensusDEGs)
amp_ad.MetaDCN$DLPFC=V(rosmap_DCN.AD2)$symbol
amp_ad.MetaDCN$FP=V(msmm_DCN$FP$AD)$name
amp_ad.MetaDCN$IFG=V(msmm_DCN$IFG$AD)$name
amp_ad.MetaDCN$PHG=V(msmm_DCN$PHG$AD)$name
amp_ad.MetaDCN$STG=V(msmm_DCN$STG$AD)$name
amp_ad.MetaDCN$CER=c("NO_DCN")
amp_ad.MetaDCN$TCX=V(mayo_TCX_DCN.AD)$name


tmp2=matrix(NA,nrow=7,ncol=8)
for(i in 1:7){
  tmp2[i,1]=names(ConsensusDEGs)[i]
  tmp2[i,2]=as.numeric(length(ConsensusDEGs[[i]]$hgnc_symbol))
  tmp2[i,3]=as.numeric(length(intersect(ConsensusDEGs[[i]]$hgnc_symbol,amp_ad.MetaDCN[[i]])))
  tmp2[i,4]=as.numeric(pval_vec_deg[[i]])
  if(length(ConsensusDEGs_OlpPathways[[i]])==0){
    tmp2[i,5]=0
  }
  else
  tmp2[i,5]=paste(ConsensusDEGs_OlpPathways[[i]],collapse = ";")
  tmp2[i,6]=as.numeric(jaccard(A = ConsensusDEGs[[i]]$hgnc_symbol,B = amp_ad.MetaDCN[[i]]))
  
}
tmp2=data.frame(tmp2,stringsAsFactors = F)
tmp2[,c(2:3)]=apply(tmp2[,c(2:3)],2,as.integer)
tmp2[,c(4,6)]=apply(tmp2[,c(4,6)],2,as.numeric)
colnames(tmp2)=c("BrainRegion","Number_Of_DEGs","Overlap_with_DEGs","p_Overlap","Overlap_EnrichedPathways","Jaccard_Overlap","Rank_Jaccard","Rank_EnrichedPathways")
tmp2$Rank_Jaccard=rank(-tmp2$Jaccard_Overlap,ties.method = "first")
tmp2$Rank_EnrichedPathways=rank(-unname(unlist(lapply(ConsensusDEGs_OlpPathways,length))),ties.method = "first")
write.table(tmp2,"../DEG_Analyses/ConsensusDEGs_DCN_OverlapSummary.txt",sep="\t",col.names = T,row.names = F,quote=F)

amp_ad_MetaDCN.GeneOverlap=matrix(NA,nrow=6,ncol=6)
rownames(amp_ad_MetaDCN.GeneOverlap)=colnames(amp_ad_MetaDCN.GeneOverlap)=names(amp_ad.MetaDCN[-6])
for(i in 1:6){
  amp_ad_MetaDCN.GeneOverlap[i,]=round(unlist(lapply(amp_ad.MetaDCN[-6],jaccard,amp_ad.MetaDCN[-6][[i]])),digits = 2)
}
### Determine TF overlap. TFs downloaded from TCoF v2
tcof_v2_db=read.table("../../TF_Databases/TCoF_V2_human_tf_interactions.csv",sep = ",",header = T,as.is = T)
phg_df_AD.DCEL075=msmm_phg_diffcorr[abs(msmm_phg_diffcorr$r.t)>=0.75,]
phg_df_Control.DCEL075=msmm_phg_diffcorr[abs(msmm_phg_diffcorr$r.c)>=0.75,]
ifg_df_AD.DCEL075=msmm_ifg_diffcorr[abs(msmm_ifg_diffcorr$r.t)>=0.75,]
ifg_df_Control.DCEL075=msmm_ifg_diffcorr[abs(msmm_ifg_diffcorr$r.c)>=0.75,]

phg_df_AD.DCEL075=msmm_phg_diffcorr[abs(msmm_phg_diffcorr$r.t)>=0.75,]
phg_df_Control.DCEL075=msmm_phg_diffcorr[abs(msmm_phg_diffcorr$r.c)>=0.75,]
ifg_df_AD.DCEL075=msmm_ifg_diffcorr[abs(msmm_ifg_diffcorr$r.t)>=0.75,]
ifg_df_Control.DCEL075=msmm_phg_diffcorr[abs(msmm_phg_diffcorr$r.c)>=0.75,]

phg_df_AD.DCEL075_graph=graph.data.frame(d = phg_df_AD.DCEL075[,c(2:3)],directed = F)
ifg_df_AD.DCEL075_graph=graph.data.frame(d = ifg_df_AD.DCEL075[,c(2:3)],directed = F)

mayo_df_TF_Control_AD_commonGenes=intersect(unique(tcof_v2_db$Symbol1),intersect(V(mayo_df_Control.DCE075_graph)$name,V(mayo_df_AD.DCE075_graph)$name))
mayo_df_Control_AD_DCE075.AD_InteractingGenes=lapply(lapply(lapply(mayo_df_TF_Control_AD_commonGenes,function(x)neighborhood(graph = mayo_df_AD.DCE075_graph,nodes = x,order = 1)),`[[`,1),function(y)paste(names(y),collapse = ","))
mayo_df_Control_AD_DCE075.Control_InteractingGenes=lapply(lapply(lapply(mayo_df_TF_Control_AD_commonGenes,function(x)neighborhood(graph = mayo_df_Control.DCE075_graph,nodes = x,order = 1)),`[[`,1),function(y)paste(names(y),collapse = ","))
names(mayo_df_Control_AD_DCE075.AD_InteractingGenes)=names(mayo_df_Control_AD_DCE075.Control_InteractingGenes)=paste("Diff_TF",mayo_df_TF_Control_AD_commonGenes,sep = "_")


mayo_df_Control_AD_DCE075.DiffDegree=matrix(NA,nrow=length(mayo_df_TF_Control_AD_commonGenes),ncol=4)
mayo_df_Control_AD_DCE075.DiffDegree[,1]=unname(degree(graph = mayo_df_Control.DCE075_graph,v = mayo_df_TF_Control_AD_commonGenes))
mayo_df_Control_AD_DCE075.DiffDegree[,2]=unname(degree(graph = mayo_df_AD.DCE075_graph,v = mayo_df_TF_Control_AD_commonGenes))
for(i in 1:dim(mayo_df_Control_AD_DCE075.DiffDegree)[1]){
  mayo_df_Control_AD_DCE075.DiffDegree[i,3]=mayo_df_Control_AD_DCE075.Control_InteractingGenes[[i]]
  mayo_df_Control_AD_DCE075.DiffDegree[i,4]=mayo_df_Control_AD_DCE075.AD_InteractingGenes[[i]]
}

mayo_df_Control_AD_DCE075.DiffDegree=data.frame(mayo_df_Control_AD_DCE075.DiffDegree,stringsAsFactors = F)
colnames(mayo_df_Control_AD_DCE075.DiffDegree)=c("Degree_Control","Degree_AD","Control_InteractingGenes","AD_InteractingGenes")

mayo_df_Control_AD_DCE075.DiffDegree$Degree_Control=as.numeric(mayo_df_Control_AD_DCE075.DiffDegree$Degree_Control)
mayo_df_Control_AD_DCE075.DiffDegree$Degree_AD=as.numeric(mayo_df_Control_AD_DCE075.DiffDegree$Degree_AD)
rownames(mayo_df_Control_AD_DCE075.DiffDegree)=mayo_df_TF_Control_AD_commonGenes
write.table(mayo_df_Control_AD_DCE075.DiffDegree,"MAYO/MAYO_TF_Control_AD_DiffDegree.txt",sep="\t",row.names = T,col.names = T,quote = F)


mayo_df_Control_AD_DCE075_AD_InteractingGenes.KEGG=lapply(lapply(lapply(lapply(mayo_df_Control_AD_DCE075.AD_InteractingGenes,strsplit,split=","),function(a)a[[1]]),function(b)gprofiler(query = b,organism = "hsapiens",max_p_value = 0.05,src_filter = "KEGG")),`[`,c("term.name"))
mayo_df_Control_AD_DCE075_Control_InteractingGenes.KEGG=lapply(lapply(lapply(lapply(mayo_df_Control_AD_DCE075.Control_InteractingGenes,strsplit,split=","),function(a)a[[1]]),function(b)gprofiler(query = b,organism = "hsapiens",max_p_value = 0.05,src_filter = "KEGG")),`[`,c("term.name","intersection"))

###########################################################################################
# Determine cell-type specific enrichment #
###########################################################################################
zhang_celltype_ADgenes=read.xls('/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/Zhang_19BrainRegions_Paper/Zhang_BrainCelltype_Markers.xlsx',skip=1,sheet=3,header=T,as.is=T)
zhang_celltype_ADgenes.list=zhang_celltype_PLQ_ADgenes.list=vector(mode = 'list',length = 5)
names(zhang_celltype_ADgenes.list)=names(zhang_celltype_PLQ_ADgenes.list)=sort(unique(zhang_celltype_ADgenes$Cell.type))
zhang_celltype_ADgenes.list$Astrocytes=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_ADgenes.list)[1],zhang_celltype_ADgenes$Cell.type)]
zhang_celltype_ADgenes.list$Endothelial=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_ADgenes.list)[2],zhang_celltype_ADgenes$Cell.type)]
zhang_celltype_ADgenes.list$Microglia=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_ADgenes.list)[3],zhang_celltype_ADgenes$Cell.type)]
zhang_celltype_ADgenes.list$Neurons=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_ADgenes.list)[4],zhang_celltype_ADgenes$Cell.type)]
zhang_celltype_ADgenes.list$Oligodendrocytes=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_ADgenes.list)[5],zhang_celltype_ADgenes$Cell.type)]
