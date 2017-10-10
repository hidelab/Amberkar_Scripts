library(igraph)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(Pi)
library(centiserve)

mapIds2<-function(IDs,IDFrom,IDTo){
require(org.Hs.eg.db)
idmap=mapIds(x = org.Hs.eg.db,keys = IDs,column = IDTo,keytype = IDFrom,multiVals = "first")
na_vec=names(idmap[is.na(idmap)==T])
idmap=idmap[is.na(idmap)==F]
idmap_df=data.frame("From"=names(idmap),"To"=unlist(unname(idmap)),stringsAsFactors = F)
return(list(map=idmap_df,noMap=na_vec))
}
cbindlist <- function(list) {
 	n <- length(list)
 	res <- NULL
 	for (i in seq(n)) res <- cbind(res, list[[i]])
 	return(res)
}
ipah_DCN05.HV=readRDS("IPAH_DCN05_HV.RDS")
ipah_DCN05.IPAH=readRDS("IPAH_DCN05_IPAH.RDS")
ipah_curated.variants=read.csv("PAH Targets Bertero 2014 and AL.csv",header = T,sep = ",",as.is = T,stringsAsFactors = F)
ipah_DCN05.HV=delete_vertices(graph = ipah_DCN05.HV,v = V(ipah_DCN05.HV)[mapIds2(IDs = V(ipah_DCN05.HV)$name,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[2]]])
ipah_DCN05.IPAH=delete_vertices(graph = ipah_DCN05.IPAH,v = V(ipah_DCN05.IPAH)[mapIds2(IDs = V(ipah_DCN05.IPAH)$name,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[2]]])

V(ipah_DCN05.HV)$symbol=mapIds2(IDs = V(ipah_DCN05.HV)$name,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[1]][,2]
V(ipah_DCN05.HV)$name=V(ipah_DCN05.HV)$symbol
V(ipah_DCN05.IPAH)$symbol=mapIds2(IDs = V(ipah_DCN05.IPAH)$name,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[1]][,2]
V(ipah_DCN05.IPAH)$name=V(ipah_DCN05.IPAH)$symbol

ipah_DCN_IPAH.seedData=data.frame(V(ipah_DCN05.HV)$name,rep(0,length(V(ipah_DCN05.HV)$name)),stringsAsFactors = F)
colnames(ipah_DCN_IPAH.seedData)=c()
ipah_DCN_IPAH.seedData[grep(pattern = paste(ipah_curated.variants$Gene,collapse = "|"),x = V(ipah_DCN05.IPAH)$symbol),2]=1
ipah_DCN_IPAH.pData=xPierGenes(data = ipah_DCN_IPAH.seedData,weighted=T,restart=0.95,network.customised=ipah_DCN05.IPAH)
ipah_DCN_IPAH.pSubnet=xPierSubnet(pNode = ipah_DCN_IPAH.pData,network.customised =ipah_DCN05.IPAH,subnet.significance = 0.01,verbose = T)
V(ipah_DCN_IPAH.pSubnet)$symbol=mapIds2(IDs = V(ipah_DCN_IPAH.pSubnet)$name,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[1]][,2]
V(ipah_DCN_IPAH.pSubnet)$color="#8B8989"
V(ipah_DCN_IPAH.pSubnet)$color[V(ipah_DCN_IPAH.pSubnet)$symbol%in%ipah_curated.variants$Gene]="#EE4000"
V(ipah_DCN_IPAH.pSubnet)$lobby=lobby(graph = ipah_DCN_IPAH.pSubnet,loops = F)
V(ipah_DCN_IPAH.pSubnet)$name=V(ipah_DCN_IPAH.pSubnet)$symbol
write.graph(ipah_DCN_IPAH.pSubnet,"IPAH_DCN_pSubnet.gml",format = "gml")


#Superhub analysis
hs_gene_table=toTable(org.Hs.egSYMBOL)
ipah.top10_Superhubs=as.list(V(ipah_DCN05.IPAH)$symbol[which(rank(-V(ipah_DCN05.IPAH)$lobby,ties.method = "first")%in%sort(rank(-V(ipah_DCN05.IPAH)$lobby,ties.method = "first"))[1:10])])
ipah.top10_Superhubs.Hops=matrix(NA,nrow=10,ncol = 6)
rownames(ipah.top10_Superhubs.Hops)=ipah.top10_Superhubs
colnames(ipah.top10_Superhubs.Hops)=paste("Hop",1:6,sep = "")
neighborhoodSignificance=function(superhubs,network,curated_genes,hs_gene_table,hop){
  neighborhood_overlaps=lapply(superhubs,function(x)neighborhood(graph = network,order=hop,which(V(network)$name%in%x)))
  rand_genes=sample(x = hs_gene_table$symbol,size = length(curated_genes),replace = F)
  neighborhood_overlaps.length=lapply(neighborhood_overlaps,function(y)length(intersect(names(y[[1]]),rand_genes)))
  neighborhood_size=lapply(ipah.top10_Superhubs,function(x)neighborhood.size(graph = ipah_DCN05.IPAH,order = 1,nodes = which(V(ipah_DCN05.IPAH)$name%in%x)))
  names(neighborhood_overlaps.length)=names(neighborhood_size)=superhubs
  return(list(neighborhood_overlaps.length,neighborhood_size))
}
for(j in 1:6){
  for(i in 1:10){
    ipah.top10_Superhubs.Hops[i,j]=length(intersect(names(neighborhood(graph = ipah_DCN05.IPAH,order = j,nodes = which(V(ipah_DCN05.IPAH)$name%in%ipah.top10_Superhubs[i]))[[1]]),ipah_curated.variants$Gene))
    }
}
ns.hop1=foreach(1:500,.packages="igraph") %dopar% {
  data.frame(hop1=neighborhoodSignificance(superhubs = ipah.top10_Superhubs,network = ipah_DCN05.IPAH,curated_genes = ipah_curated.variants$Gene,hs_gene_table = hs_gene_table,hop = 1))
}
ns.hop2=foreach(1:500,.packages="igraph") %dopar% {
  data.frame(hop2=neighborhoodSignificance(superhubs = ipah.top10_Superhubs,network = ipah_DCN05.IPAH,curated_genes = ipah_curated.variants$Gene,hs_gene_table = hs_gene_table,hop = 2))
}
ns.hop3=foreach(1:500,.packages="igraph") %dopar% {
  data.frame(hop3=neighborhoodSignificance(superhubs = ipah.top10_Superhubs,network = ipah_DCN05.IPAH,curated_genes = ipah_curated.variants$Gene,hs_gene_table = hs_gene_table,hop = 3))
}
ns.hop4=foreach(1:500,.packages="igraph") %dopar% {
  data.frame(hop4=neighborhoodSignificance(superhubs = ipah.top10_Superhubs,network = ipah_DCN05.IPAH,curated_genes = ipah_curated.variants$Gene,hs_gene_table = hs_gene_table,hop = 4))
}

lobby_hopSignificance=readRDS("lobby_hopSignificance_75bp.RDS")


