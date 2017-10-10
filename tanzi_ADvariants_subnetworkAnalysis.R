library(Pi)
library(clusterProfiler)
library(data.table)
library(centiserve)

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
tanzi_variants=read.table("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Gene_lists_snp_analysis_Aug.txt",sep = "\t",header=T,as.is=T)
names(tanzi_selectGenes)=mapIds2(IDs = tanzi_selectGenes,IDFrom = "SYMBOL",IDTo = "ENSEMBL")[[1]][,2]

rosmap_DLPFC_DCN.AD=readRDS('./ROSMAP/rosmap_DCN_AD_pValCorr.RDS')
rosmap_DLPFC_DCN.AD2=delete_vertices(graph = rosmap_DLPFC_DCN.AD,v = grep(pattern = paste(mapIds2(IDs = V(rosmap_DLPFC_DCN.AD)$name,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[2]],collapse = "|"),x = V(rosmap_DLPFC_DCN.AD)$name))
V(rosmap_DLPFC_DCN.AD2)$symbol=mapIds2(IDs = V(rosmap_DLPFC_DCN.AD2)$name,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[1]][,2]
V(rosmap_DLPFC_DCN.AD2)$entrez=mapIds2(IDs = V(rosmap_DLPFC_DCN.AD2)$name,IDFrom = "ENSEMBL",IDTo = "ENTREZID")[[1]][,2]

mayo_TCX_DCN.AD=readRDS('./MAYO/MAYO_ReSeq_TCX_DCN_AD.RDS')
mayo_TCX_DCN.AD2=delete_vertices(graph = mayo_TCX_DCN.AD,v = grep(pattern = paste(mapIds2(IDs = V(mayo_TCX_DCN.AD)$name,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[2]],collapse = "|"),x = V(mayo_TCX_DCN.AD)$name))

V(mayo_TCX_DCN.AD2)$symbol=mapIds2(IDs = V(mayo_TCX_DCN.AD2)$name,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[1]][,2]
V(mayo_TCX_DCN.AD2)$entrez=mapIds2(IDs = V(mayo_TCX_DCN.AD2)$name,IDFrom = "ENSEMBL",IDTo = "ENTREZID")[[1]][,2]
mayo_TCX_DCN_AD2.SeedData=data.frame(V(mayo_TCX_DCN.AD2)$name,rep(0,length(V(mayo_TCX_DCN.AD2)$name)),stringsAsFactors = F)
mayo_TCX_DCN_AD2.SeedData[grep(pattern = paste(AD_genes,collapse = "|"),x = V(mayo_TCX_DCN.AD2)$symbol),2]=1
colnames(mayo_TCX_DCN_AD2.SeedData)=c()
mayo_TCX_DCN_AD2.pData=xPierGenes(data=mayo_TCX_DCN_AD2.SeedData,weighted=T,restart=0.65,network.customised=mayo_TCX_DCN.AD2,parallel = F)
mayo_TCX_DCN_AD2.pSubnet=xPierSubnet(pNode = mayo_TCX_DCN_AD2.pData,network.customised = mayo_TCX_DCN.AD2,subnet.significance = 0.05,verbose = T,priority.quantile = NA)
V(mayo_TCX_DCN_AD2.pSubnet)$name=mapIds2(IDs = V(mayo_TCX_DCN_AD2.pSubnet)$name,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[1]][,2]
V(mayo_TCX_DCN_AD2.pSubnet)$significance=as.numeric(V(mayo_TCX_DCN_AD2.pSubnet)$significance)
V(mayo_TCX_DCN_AD2.pSubnet)$score=as.numeric(V(mayo_TCX_DCN_AD2.pSubnet)$score)
V(mayo_TCX_DCN_AD2.pSubnet)$priority=as.numeric(V(mayo_TCX_DCN_AD2.pSubnet)$priority)
V(mayo_TCX_DCN_AD2.pSubnet)$colour[grep(pattern = paste(AD_genes,collapse = "|"),x = V(mayo_TCX_DCN_AD2.pSubnet)$name)]="#EE4000"
V(mayo_TCX_DCN_AD2.pSubnet)$colour[-grep(pattern = paste(AD_genes,collapse = "|"),x = V(mayo_TCX_DCN_AD2.pSubnet)$name)]="#8B8989"






























#################################################################################################
sp_rosmap_dcn=shortest.paths(graph=rosmap_DCN.AD2,v=V(rosmap_DCN.AD2),to=V(rosmap_DCN.AD2))
sp_rosmap_dcn2=upperTriangle(sp_rosmap_dcn)[which(upperTriangle(sp_rosmap_dcn)!='Inf')]
sp_rosmap_dcn_Control_LowP=shortest.paths(graph=rosmap_DCN.AD2,v=which(V(rosmap_DCN.AD2)$symbol%in%tanzi_variants$Control_LowP),to=which(V(rosmap_DCN.AD2)$symbol%in%tanzi_variants$Control_LowP))
sp_rosmap_dcn_Control_LowP2=upperTriangle(sp_rosmap_dcn_Control_LowP)[which(upperTriangle(sp_rosmap_dcn_Control_LowP)!='Inf')]
sp_rosmap_dcn_AD_LowP=shortest.paths(graph=rosmap_DCN.AD2,v=which(V(rosmap_DCN.AD2)$symbol%in%tanzi_variants$AD_LowP),to=which(V(rosmap_DCN.AD2)$symbol%in%tanzi_variants$AD_LowP))
sp_rosmap_dcn_AD_LowP2=upperTriangle(sp_rosmap_dcn_AD_LowP)[which(upperTriangle(sp_rosmap_dcn_AD_LowP)!='Inf')]
sp_rosmap_dcn_Control_HighP=shortest.paths(graph=rosmap_DCN.AD2,v=which(V(rosmap_DCN.AD2)$symbol%in%tanzi_variants$Control_HighP),to=which(V(rosmap_DCN.AD2)$symbol%in%tanzi_variants$Control_HighP))
sp_rosmap_dcn_Control_HighP2=upperTriangle(sp_rosmap_dcn_Control_HighP)[which(upperTriangle(sp_rosmap_dcn_Control_HighP)!='Inf')]
sp_rosmap_dcn_AD_HighP=shortest.paths(graph=rosmap_DCN.AD2,v=which(V(rosmap_DCN.AD2)$symbol%in%tanzi_variants$AD_HighP),to=which(V(rosmap_DCN.AD2)$symbol%in%tanzi_variants$AD_HighP))
sp_rosmap_dcn_AD_HighP2=upperTriangle(sp_rosmap_dcn_AD_HighP)[which(upperTriangle(sp_rosmap_dcn_AD_HighP)!='Inf')]

cls_rosmap_dcn_Control_LowP=closeness(graph=rosmap_DCN.AD2,vids=which(V(rosmap_DCN.AD2)$symbol%in%tanzi_variants$Control_LowP))
cls_rosmap_dcn_AD_LowP=closeness(graph=rosmap_DCN.AD2,vids=which(V(rosmap_DCN.AD2)$symbol%in%tanzi_variants$AD_LowP))
cls_rosmap_dcn_Control_HighP=closeness(graph=rosmap_DCN.AD2,vids=which(V(rosmap_DCN.AD2)$symbol%in%tanzi_variants$Control_HighP))
cls_rosmap_dcn_AD_HighP=closeness(graph=rosmap_DCN.AD2,vids=which(V(rosmap_DCN.AD2)$symbol%in%tanzi_variants$AD_HighP))

deg_rosmap_dcn_Control_LowP=degree(graph=rosmap_DCN.AD2,v=which(V(rosmap_DCN.AD2)$symbol%in%tanzi_variants$Control_LowP))
deg_rosmap_dcn_AD_LowP=degree(graph=rosmap_DCN.AD2,v=which(V(rosmap_DCN.AD2)$symbol%in%tanzi_variants$AD_LowP))
deg_rosmap_dcn_Control_HighP=degree(graph=rosmap_DCN.AD2,v=which(V(rosmap_DCN.AD2)$symbol%in%tanzi_variants$Control_HighP))
deg_rosmap_dcn_AD_HighP=degree(graph=rosmap_DCN.AD2,v=which(V(rosmap_DCN.AD2)$symbol%in%tanzi_variants$AD_HighP))
