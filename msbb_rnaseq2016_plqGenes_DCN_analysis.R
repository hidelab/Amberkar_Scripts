library(igraph)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(gdata)

RewiredUniqueGenes=function(dcn,posThreshold,negThreshold){
  normal_poscoexp.genes=V(graph.data.frame(d=data.frame(ends(graph = dcn[[1]],es = which(E(dcn[[1]])$weight>posThreshold),names = T)),directed = F))$name
  disease_poscoexp.genes=V(graph.data.frame(d=data.frame(ends(graph = dcn[[2]],es = which(E(dcn[[2]])$weight>posThreshold),names = T)),directed = F))$name
  normal_negcoexp.genes=V(graph.data.frame(d=data.frame(ends(graph = dcn[[1]],es = which(E(dcn[[1]])$weight<negThreshold),names = T)),directed = F))$name
  disease_negcoexp.genes=V(graph.data.frame(d=data.frame(ends(graph = dcn[[2]],es = which(E(dcn[[2]])$weight<negThreshold),names = T)),directed = F))$name
  #Setup datastructure
  rewired_unique_genes=vector(mode = 'list',length = 2)
  names(rewired_unique_genes)=c(names(dcn)[1],names(dcn)[2])
  rewired_unique_genes[[1]]=rewired_unique_genes[[2]]=vector(mode = 'list',length = 4)
  names(rewired_unique_genes[[1]])=names(rewired_unique_genes[[2]])=c('PosCoexp_UniqGenes','NegCoexp_UniqGenes','PosCoexp_Degree','NegCoexp_Degree')
  #For Normal samples, always the first element of datastruct
  rewired_unique_genes[[1]]$PosCoexp_UniqGenes=V(dcn[[1]])$name[which(V(dcn[[1]])$name%in%setdiff(normal_poscoexp.genes,disease_poscoexp.genes))]
  rewired_unique_genes[[1]]$NegCoexp_UniqGenes=V(dcn[[1]])$name[which(V(dcn[[1]])$name%in%setdiff(normal_negcoexp.genes,disease_negcoexp.genes))]
  rewired_unique_genes[[1]]$PosCoexp_Degree=data.frame(Normal_PosCoexp_Unique_Genes=V(dcn[[1]])$name[which(V(dcn[[1]])$name%in%setdiff(normal_poscoexp.genes,disease_poscoexp.genes))],
                                                       Normal_PosCoexp_Unique_Degree=unname(degree(graph = dcn[[1]],v = which(V(dcn[[1]])$name%in%setdiff(normal_poscoexp.genes,disease_poscoexp.genes)))),
                                                       Disease_PosCoexp_Genes=0,
                                                       Diff=unname(degree(graph = dcn[[1]],v = which(V(dcn[[1]])$name%in%setdiff(normal_poscoexp.genes,disease_poscoexp.genes))))-0,stringsAsFactors = F)
  rewired_unique_genes[[1]]$PosCoexp_Degree=rewired_unique_genes[[1]]$PosCoexp_Degree[order(rewired_unique_genes[[1]]$PosCoexp_Degree$Diff,decreasing = T),]
  rewired_unique_genes[[1]]$NegCoexp_Degree=data.frame(Normal_NegCoexp_Unique_Genes=V(dcn[[1]])$name[which(V(dcn[[1]])$name%in%setdiff(normal_negcoexp.genes,disease_negcoexp.genes))],
                                                       Normal_NegCoexp_Unique_Degree=unname(degree(graph = dcn[[1]],v = which(V(dcn[[1]])$name%in%setdiff(normal_negcoexp.genes,disease_negcoexp.genes)))),
                                                       Disease_NegCoexp_Genes=0,
                                                       Diff=unname(degree(graph = dcn[[1]],v = which(V(dcn[[1]])$name%in%setdiff(normal_negcoexp.genes,disease_negcoexp.genes))))-0,stringsAsFactors = F)
  rewired_unique_genes[[1]]$NegCoexp_Degree=rewired_unique_genes[[1]]$NegCoexp_Degree[order(rewired_unique_genes[[1]]$NegCoexp_Degree$Diff,decreasing = T),]
  #For Disease samples, always the second element of datastruct
  rewired_unique_genes[[2]]$NegCoexp_UniqGenes=V(dcn[[2]])$name[which(V(dcn[[2]])$name%in%setdiff(disease_negcoexp.genes,normal_negcoexp.genes))]
  rewired_unique_genes[[2]]$PosCoexp_UniqGenes=V(dcn[[2]])$name[which(V(dcn[[2]])$name%in%setdiff(disease_poscoexp.genes,normal_poscoexp.genes))]
  rewired_unique_genes[[2]]$PosCoexp_Degree=data.frame(Disease_PosCoexp_Unique_Genes=V(dcn[[2]])$name[which(V(dcn[[2]])$name%in%setdiff(disease_poscoexp.genes,normal_poscoexp.genes))],
                                                       Disease_PosCoexp_Unique_Degree=unname(degree(graph = dcn[[2]],v = which(V(dcn[[2]])$name%in%setdiff(disease_poscoexp.genes,normal_poscoexp.genes)))),
                                                       Normal_PosCoexp_Genes=0,
                                                       Diff=unname(degree(graph = dcn[[2]],v = which(V(dcn[[2]])$name%in%setdiff(disease_poscoexp.genes,normal_poscoexp.genes))))-0,stringsAsFactors = F)
  rewired_unique_genes[[2]]$PosCoexp_Degree=rewired_unique_genes[[2]]$PosCoexp_Degree[order(rewired_unique_genes[[2]]$PosCoexp_Degree$Diff,decreasing = T),]
  rewired_unique_genes[[2]]$NegCoexp_Degree=data.frame(Disease_NegCoexp_Unique_Genes=V(dcn[[2]])$name[which(V(dcn[[2]])$name%in%setdiff(disease_negcoexp.genes,normal_negcoexp.genes))],
                                                       Disease_NegCoexp_Unique_Degree=unname(degree(graph = dcn[[2]],v = which(V(dcn[[2]])$name%in%setdiff(disease_negcoexp.genes,normal_negcoexp.genes)))),
                                                       Normal_NegCoexp_Genes=0,
                                                       Diff=unname(degree(graph = dcn[[2]],v = which(V(dcn[[2]])$name%in%setdiff(disease_negcoexp.genes,normal_negcoexp.genes))))-0,stringsAsFactors = F)
  rewired_unique_genes[[2]]$NegCoexp_Degree=rewired_unique_genes[[2]]$NegCoexp_Degree[order(rewired_unique_genes[[2]]$NegCoexp_Degree$Diff,decreasing = T),]
  return(rewired_unique_genes)
}

msbb_rnaseq2016.DCN=msbb_rnaseq2016_R2.DCN=msbb_rnaseq2016_PLQ.DCN=msbb_rnaseq2016_PLQGenes.cluster=msbb_dcn.clusters=msbb_diffcorr_fdr10=vector(mode = "list",length = 4)
names(msbb_rnaseq2016.DCN)=names(msbb_rnaseq2016_R2.DCN)=names(msbb_rnaseq2016_PLQ.DCN)=names(msbb_rnaseq2016_PLQGenes.cluster)=names(msbb_dcn.clusters)=names(msbb_diffcorr_fdr10)=c("FP","IFG","PHG","STG")

msbb_rnaseq2016.DCN$FP=msbb_rnaseq2016.DCN$IFG=msbb_rnaseq2016.DCN$PHG=msbb_rnaseq2016.DCN$STG=vector(mode="list",length=2)
msbb_rnaseq2016_PLQ.DCN$FP=msbb_rnaseq2016_PLQ.DCN$IFG=msbb_rnaseq2016_PLQ.DCN$PHG=msbb_rnaseq2016_PLQ.DCN$STG=vector(mode="list",length=2)
names(msbb_rnaseq2016.DCN$FP)=names(msbb_rnaseq2016.DCN$IFG)=names(msbb_rnaseq2016.DCN$PHG)=names(msbb_rnaseq2016.DCN$STG)=c("Low","High")
names(msbb_rnaseq2016_PLQ.DCN$FP)=names(msbb_rnaseq2016_PLQ.DCN$IFG)=names(msbb_rnaseq2016_PLQ.DCN$PHG)=names(msbb_rnaseq2016_PLQ.DCN$STG)=c("Low","High")

msbb_rnaseq2016.DCN$FP=readRDS("FP/DCN/FP_DCN_FDR01.RDS")
msbb_rnaseq2016.DCN$IFG=readRDS("IFG/DCN/IFG_DCN_FDR01.RDS")
msbb_rnaseq2016.DCN$PHG=readRDS("PHG/DCN/PHG_DCN_FDR01.RDS")
msbb_rnaseq2016.DCN$STG=readRDS("STG/DCN/STG_DCN_FDR01.RDS")
msbb_rnaseq2016_PLQ.DCN$FP=readRDS("FP/PLQ_DCN/FP_PLQGenes_FDR01.DCN.RDS")
msbb_rnaseq2016_PLQ.DCN$IFG=readRDS("IFG/PLQ_DCN/IFG_PLQGenes_FDR01.DCN.RDS")
msbb_rnaseq2016_PLQ.DCN$PHG=readRDS("PHG/PLQ_DCN/PHG_PLQGenes_FDR01.DCN.RDS")
msbb_rnaseq2016_PLQ.DCN$STG=readRDS("STG/PLQ_DCN/STG_PLQGenes_FDR01.DCN.RDS")

msbb_rnaseq2016_FP_uniq_rewired=RewiredUniqueGenes(dcn = msbb_rnaseq2016.DCN$FP,posThreshold = 1.5,negThreshold = 0.5)
msbb_rnaseq2016_IFG_uniq_rewired=RewiredUniqueGenes(dcn = msbb_rnaseq2016.DCN$IFG,posThreshold = 1.5,negThreshold = 0.5)
msbb_rnaseq2016_PHG_uniq_rewired=RewiredUniqueGenes(dcn = msbb_rnaseq2016.DCN$PHG,posThreshold = 1.5,negThreshold = 0.5)
msbb_rnaseq2016_STG_uniq_rewired=RewiredUniqueGenes(dcn = msbb_rnaseq2016.DCN$STG,posThreshold = 1.5,negThreshold = 0.5)

msbb_rnaseq2016_PLQ.PosCoexp=msbb_rnaseq2016_PLQ.NegCoexp=msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF=msbb_rnaseq2016_PLQGenes_NegCoexp.rewiredDF=msbb_rnaseq2016_PLQGenes=vector(mode = "list",length = 4)
names(msbb_rnaseq2016_PLQGenes)=names(msbb_rnaseq2016_PLQ.PosCoexp)=names(msbb_rnaseq2016_PLQ.NegCoexp)=names(msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF)=names(msbb_rnaseq2016_PLQGenes_NegCoexp.rewiredDF)=names(msbb_rnaseq2016_PLQGenes.cluster)
msbb_rnaseq2016_PLQ.PosCoexp$FP=msbb_rnaseq2016_PLQ.PosCoexp$IFG=msbb_rnaseq2016_PLQ.PosCoexp$PHG=msbb_rnaseq2016_PLQ.PosCoexp$STG=vector(mode = "list",length = 2)
msbb_rnaseq2016_PLQ.NegCoexp$FP=msbb_rnaseq2016_PLQ.NegCoexp$IFG=msbb_rnaseq2016_PLQ.NegCoexp$PHG=msbb_rnaseq2016_PLQ.NegCoexp$STG=vector(mode = "list",length = 2)
names(msbb_rnaseq2016_PLQ.PosCoexp$FP)=names(msbb_rnaseq2016_PLQ.NegCoexp$FP)=c("Low","High")
names(msbb_rnaseq2016_PLQ.PosCoexp$IFG)=names(msbb_rnaseq2016_PLQ.NegCoexp$IFG)=c("Low","High")
names(msbb_rnaseq2016_PLQ.PosCoexp$PHG)=names(msbb_rnaseq2016_PLQ.NegCoexp$PHG)=c("Low","High")
names(msbb_rnaseq2016_PLQ.PosCoexp$STG)=names(msbb_rnaseq2016_PLQ.NegCoexp$STG)=c("Low","High")

msbb_rnaseq2016_PLQ.PosCoexp$FP$Low=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016_PLQ.DCN$FP$Low,es = which(E(msbb_rnaseq2016_PLQ.DCN$FP$Low)$weight>1.5),names = T)),directed = F)
E(msbb_rnaseq2016_PLQ.PosCoexp$FP$Low)$weight=E(msbb_rnaseq2016_PLQ.DCN$FP$Low)$weight[which(E(msbb_rnaseq2016_PLQ.DCN$FP$Low)$weight>1.5)]
write.graph(graph = msbb_rnaseq2016_PLQ.PosCoexp$FP$Low,file = 'FP_PLQGenes_PosCoexp_Low.gml',format = 'gml')
msbb_rnaseq2016_PLQ.PosCoexp$IFG$Low=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016_PLQ.DCN$IFG$Low,es = which(E(msbb_rnaseq2016_PLQ.DCN$IFG$Low)$weight>1.5),names = T)),directed = F)
E(msbb_rnaseq2016_PLQ.PosCoexp$IFG$Low)$weight=E(msbb_rnaseq2016_PLQ.DCN$IFG$Low)$weight[which(E(msbb_rnaseq2016_PLQ.DCN$IFG$Low)$weight>1.5)]
write.graph(graph = msbb_rnaseq2016_PLQ.PosCoexp$IFG$Low,file = 'IFG_PLQGenes_PosCoexp_Low.gml',format = 'gml')
msbb_rnaseq2016_PLQ.PosCoexp$PHG$Low=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016_PLQ.DCN$PHG$Low,es = which(E(msbb_rnaseq2016_PLQ.DCN$PHG$Low)$weight>1.5),names = T)),directed = F)
E(msbb_rnaseq2016_PLQ.PosCoexp$PHG$Low)$weight=E(msbb_rnaseq2016_PLQ.DCN$PHG$Low)$weight[which(E(msbb_rnaseq2016_PLQ.DCN$PHG$Low)$weight>1.5)]
write.graph(graph = msbb_rnaseq2016_PLQ.PosCoexp$PHG$Low,file = 'PHG_PLQGenes_PosCoexp_Low.gml',format = 'gml')
msbb_rnaseq2016_PLQ.PosCoexp$STG$Low=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016_PLQ.DCN$STG$Low,es = which(E(msbb_rnaseq2016_PLQ.DCN$STG$Low)$weight>1.5),names = T)),directed = F)
E(msbb_rnaseq2016_PLQ.PosCoexp$STG$Low)$weight=E(msbb_rnaseq2016_PLQ.DCN$STG$Low)$weight[which(E(msbb_rnaseq2016_PLQ.DCN$STG$Low)$weight>1.5)]
write.graph(graph = msbb_rnaseq2016_PLQ.PosCoexp$STG$Low,file = 'STG_PLQGenes_PosCoexp_Low.gml',format = 'gml')

msbb_rnaseq2016_PLQ.PosCoexp$FP$High=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016_PLQ.DCN$FP$High,es = which(E(msbb_rnaseq2016_PLQ.DCN$FP$High)$weight>1.5),names = T)),directed = F)
E(msbb_rnaseq2016_PLQ.PosCoexp$FP$High)$weight=E(msbb_rnaseq2016_PLQ.DCN$FP$High)$weight[which(E(msbb_rnaseq2016_PLQ.DCN$FP$High)$weight>1.5)]
write.graph(graph = msbb_rnaseq2016_PLQ.PosCoexp$FP$High,file = 'FP_PLQGenes_PosCoexp_High.gml',format = 'gml')
msbb_rnaseq2016_PLQ.PosCoexp$IFG$High=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016_PLQ.DCN$IFG$High,es = which(E(msbb_rnaseq2016_PLQ.DCN$IFG$High)$weight>1.5),names = T)),directed = F)
E(msbb_rnaseq2016_PLQ.PosCoexp$IFG$High)$weight=E(msbb_rnaseq2016_PLQ.DCN$IFG$High)$weight[which(E(msbb_rnaseq2016_PLQ.DCN$IFG$High)$weight>1.5)]
write.graph(graph = msbb_rnaseq2016_PLQ.PosCoexp$IFG$High,file = 'IFG_PLQGenes_PosCoexp_High.gml',format = 'gml')
msbb_rnaseq2016_PLQ.PosCoexp$PHG$High=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016_PLQ.DCN$PHG$High,es = which(E(msbb_rnaseq2016_PLQ.DCN$PHG$High)$weight>1.5),names = T)),directed = F)
E(msbb_rnaseq2016_PLQ.PosCoexp$PHG$High)$weight=E(msbb_rnaseq2016_PLQ.DCN$PHG$High)$weight[which(E(msbb_rnaseq2016_PLQ.DCN$PHG$High)$weight>1.5)]
write.graph(graph = msbb_rnaseq2016_PLQ.PosCoexp$PHG$High,file = 'PHG_PLQGenes_PosCoexp_High.gml',format = 'gml')
msbb_rnaseq2016_PLQ.PosCoexp$STG$High=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016_PLQ.DCN$STG$High,es = which(E(msbb_rnaseq2016_PLQ.DCN$STG$High)$weight>1.5),names = T)),directed = F)
E(msbb_rnaseq2016_PLQ.PosCoexp$STG$High)$weight=E(msbb_rnaseq2016_PLQ.DCN$STG$High)$weight[which(E(msbb_rnaseq2016_PLQ.DCN$STG$High)$weight>1.5)]
write.graph(graph = msbb_rnaseq2016_PLQ.PosCoexp$STG$High,file = 'STG_PLQGenes_PosCoexp_High.gml',format = 'gml')

msbb_rnaseq2016_PLQ.NegCoexp$FP$High=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016_PLQ.DCN$FP$High,es = which(E(msbb_rnaseq2016_PLQ.DCN$FP$High)$weight<0.5),names = T)),directed = F)
E(msbb_rnaseq2016_PLQ.NegCoexp$FP$High)$weight=E(msbb_rnaseq2016_PLQ.DCN$FP$High)$weight[which(E(msbb_rnaseq2016_PLQ.DCN$FP$High)$weight<0.5)]
write.graph(graph = msbb_rnaseq2016_PLQ.NegCoexp$FP$High,file = 'FP_PLQGenes_NegCoexp_High.gml',format = 'gml')
msbb_rnaseq2016_PLQ.NegCoexp$IFG$High=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016_PLQ.DCN$IFG$High,es = which(E(msbb_rnaseq2016_PLQ.DCN$IFG$High)$weight<0.5),names = T)),directed = F)
E(msbb_rnaseq2016_PLQ.NegCoexp$IFG$High)$weight=E(msbb_rnaseq2016_PLQ.DCN$IFG$High)$weight[which(E(msbb_rnaseq2016_PLQ.DCN$IFG$High)$weight<0.5)]
write.graph(graph = msbb_rnaseq2016_PLQ.NegCoexp$IFG$High,file = 'IFG_PLQGenes_NegCoexp_High.gml',format = 'gml')
msbb_rnaseq2016_PLQ.NegCoexp$PHG$High=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016_PLQ.DCN$PHG$High,es = which(E(msbb_rnaseq2016_PLQ.DCN$PHG$High)$weight<0.5),names = T)),directed = F)
E(msbb_rnaseq2016_PLQ.NegCoexp$PHG$High)$weight=E(msbb_rnaseq2016_PLQ.DCN$PHG$High)$weight[which(E(msbb_rnaseq2016_PLQ.DCN$PHG$High)$weight<0.5)]
write.graph(graph = msbb_rnaseq2016_PLQ.NegCoexp$PHG$High,file = 'PHG_PLQGenes_NegCoexp_High.gml',format = 'gml')
msbb_rnaseq2016_PLQ.NegCoexp$STG$High=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016_PLQ.DCN$STG$High,es = which(E(msbb_rnaseq2016_PLQ.DCN$STG$High)$weight<0.5),names = T)),directed = F)
E(msbb_rnaseq2016_PLQ.NegCoexp$STG$High)$weight=E(msbb_rnaseq2016_PLQ.DCN$STG$High)$weight[which(E(msbb_rnaseq2016_PLQ.DCN$STG$High)$weight<0.5)]
write.graph(graph = msbb_rnaseq2016_PLQ.NegCoexp$STG$High,file = 'STG_PLQGenes_NegCoexp_High.gml',format = 'gml')

msbb_rnaseq2016_PLQ.NegCoexp$FP$Low=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016_PLQ.DCN$FP$Low,es = which(E(msbb_rnaseq2016_PLQ.DCN$FP$Low)$weight<0.5),names = T)),directed = F)
E(msbb_rnaseq2016_PLQ.NegCoexp$FP$Low)$weight=E(msbb_rnaseq2016_PLQ.DCN$FP$Low)$weight[which(E(msbb_rnaseq2016_PLQ.DCN$FP$Low)$weight<0.5)]
write.graph(graph = msbb_rnaseq2016_PLQ.NegCoexp$FP$Low,file = 'FP_PLQGenes_NegCoexp_Low.gml',format = 'gml')
msbb_rnaseq2016_PLQ.NegCoexp$IFG$Low=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016_PLQ.DCN$IFG$Low,es = which(E(msbb_rnaseq2016_PLQ.DCN$IFG$Low)$weight<0.5),names = T)),directed = F)
E(msbb_rnaseq2016_PLQ.NegCoexp$IFG$Low)$weight=E(msbb_rnaseq2016_PLQ.DCN$IFG$Low)$weight[which(E(msbb_rnaseq2016_PLQ.DCN$IFG$Low)$weight<0.5)]
write.graph(graph = msbb_rnaseq2016_PLQ.NegCoexp$IFG$Low,file = 'IFG_PLQGenes_NegCoexp_Low.gml',format = 'gml')
msbb_rnaseq2016_PLQ.NegCoexp$PHG$Low=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016_PLQ.DCN$PHG$Low,es = which(E(msbb_rnaseq2016_PLQ.DCN$PHG$Low)$weight<0.5),names = T)),directed = F)
E(msbb_rnaseq2016_PLQ.NegCoexp$PHG$Low)$weight=E(msbb_rnaseq2016_PLQ.DCN$PHG$Low)$weight[which(E(msbb_rnaseq2016_PLQ.DCN$PHG$Low)$weight<0.5)]
write.graph(graph = msbb_rnaseq2016_PLQ.NegCoexp$PHG$Low,file = 'PHG_PLQGenes_NegCoexp_Low.gml',format = 'gml')
msbb_rnaseq2016_PLQ.NegCoexp$STG$Low=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016_PLQ.DCN$STG$Low,es = which(E(msbb_rnaseq2016_PLQ.DCN$STG$Low)$weight<0.5),names = T)),directed = F)
E(msbb_rnaseq2016_PLQ.NegCoexp$STG$Low)$weight=E(msbb_rnaseq2016_PLQ.DCN$STG$Low)$weight[which(E(msbb_rnaseq2016_PLQ.DCN$STG$Low)$weight<0.5)]
write.graph(graph = msbb_rnaseq2016_PLQ.NegCoexp$STG$Low,file = 'STG_PLQGenes_NegCoexp_Low.gml',format = 'gml')

msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF$FP=data.frame(Common_Genes=intersect(V(msbb_rnaseq2016_PLQ.PosCoexp$FP$Low)$name,V(msbb_rnaseq2016_PLQ.PosCoexp$FP$High)$name),
                                                          Low=degree(graph = msbb_rnaseq2016_PLQ.PosCoexp$FP$Low,v = intersect(V(msbb_rnaseq2016_PLQ.PosCoexp$FP$Low)$name,V(msbb_rnaseq2016_PLQ.PosCoexp$FP$High)$name)),
                                                          High=degree(graph = msbb_rnaseq2016_PLQ.PosCoexp$FP$High,v = intersect(V(msbb_rnaseq2016_PLQ.PosCoexp$FP$Low)$name,V(msbb_rnaseq2016_PLQ.PosCoexp$FP$High)$name)),
                                                          Difference=degree(graph = msbb_rnaseq2016_PLQ.PosCoexp$FP$Low,v = intersect(V(msbb_rnaseq2016_PLQ.PosCoexp$FP$Low)$name,V(msbb_rnaseq2016_PLQ.PosCoexp$FP$High)$name))-degree(graph = msbb_rnaseq2016_PLQ.PosCoexp$FP$High,v = intersect(V(msbb_rnaseq2016_PLQ.PosCoexp$FP$Low)$name,V(msbb_rnaseq2016_PLQ.PosCoexp$FP$High)$name)),stringsAsFactors = F)
msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF$FP=msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF$FP[order(msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF$FP$Difference,decreasing = T),]
WriteXLS(x = msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF$FP,ExcelFileName = 'FP_PosCoexp_PLQGenes_rewiredDF.xls',SheetNames = 'FP_PosCoexp_PLQGenes',BoldHeaderRow = T)
msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF$IFG=data.frame(Common_Genes=intersect(V(msbb_rnaseq2016_PLQ.PosCoexp$IFG$Low)$name,V(msbb_rnaseq2016_PLQ.PosCoexp$IFG$High)$name),
                                                           Low=degree(graph = msbb_rnaseq2016_PLQ.PosCoexp$IFG$Low,v = intersect(V(msbb_rnaseq2016_PLQ.PosCoexp$IFG$Low)$name,V(msbb_rnaseq2016_PLQ.PosCoexp$IFG$High)$name)),
                                                           High=degree(graph = msbb_rnaseq2016_PLQ.PosCoexp$IFG$High,v = intersect(V(msbb_rnaseq2016_PLQ.PosCoexp$IFG$Low)$name,V(msbb_rnaseq2016_PLQ.PosCoexp$IFG$High)$name)),
                                                           Difference=degree(graph = msbb_rnaseq2016_PLQ.PosCoexp$IFG$Low,v = intersect(V(msbb_rnaseq2016_PLQ.PosCoexp$IFG$Low)$name,V(msbb_rnaseq2016_PLQ.PosCoexp$IFG$High)$name))-degree(graph = msbb_rnaseq2016_PLQ.PosCoexp$IFG$High,v = intersect(V(msbb_rnaseq2016_PLQ.PosCoexp$IFG$Low)$name,V(msbb_rnaseq2016_PLQ.PosCoexp$IFG$High)$name)),stringsAsFactors = F)
msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF$IFG=msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF$IFG[order(msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF$IFG$Difference,decreasing = T),]
WriteXLS(x = msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF$IFG,ExcelFileName = 'IFG_PosCoexp_PLQGenes_rewiredDF.xls',SheetNames = 'IFG_PosCoexp_PLQGenes',BoldHeaderRow = T)
msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF$PHG=data.frame(Common_Genes=intersect(V(msbb_rnaseq2016_PLQ.PosCoexp$PHG$Low)$name,V(msbb_rnaseq2016_PLQ.PosCoexp$PHG$High)$name),
                                                           Low=degree(graph = msbb_rnaseq2016_PLQ.PosCoexp$PHG$Low,v = intersect(V(msbb_rnaseq2016_PLQ.PosCoexp$PHG$Low)$name,V(msbb_rnaseq2016_PLQ.PosCoexp$PHG$High)$name)),
                                                           High=degree(graph = msbb_rnaseq2016_PLQ.PosCoexp$PHG$High,v = intersect(V(msbb_rnaseq2016_PLQ.PosCoexp$PHG$Low)$name,V(msbb_rnaseq2016_PLQ.PosCoexp$PHG$High)$name)),
                                                           Difference=degree(graph = msbb_rnaseq2016_PLQ.PosCoexp$PHG$Low,v = intersect(V(msbb_rnaseq2016_PLQ.PosCoexp$PHG$Low)$name,V(msbb_rnaseq2016_PLQ.PosCoexp$PHG$High)$name))-degree(graph = msbb_rnaseq2016_PLQ.PosCoexp$PHG$High,v = intersect(V(msbb_rnaseq2016_PLQ.PosCoexp$PHG$Low)$name,V(msbb_rnaseq2016_PLQ.PosCoexp$PHG$High)$name)),stringsAsFactors = F)
msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF$PHG=msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF$PHG[order(msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF$PHG$Difference,decreasing = T),]
WriteXLS(x = msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF$PHG,ExcelFileName = 'PHG_PosCoexp_PLQGenes_rewiredDF.xls',SheetNames = 'PHG_PosCoexp_PLQGenes',BoldHeaderRow = T)
msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF$STG=data.frame(Common_Genes=intersect(V(msbb_rnaseq2016_PLQ.PosCoexp$STG$Low)$name,V(msbb_rnaseq2016_PLQ.PosCoexp$STG$High)$name),
                                                           Low=degree(graph = msbb_rnaseq2016_PLQ.PosCoexp$STG$Low,v = intersect(V(msbb_rnaseq2016_PLQ.PosCoexp$STG$Low)$name,V(msbb_rnaseq2016_PLQ.PosCoexp$STG$High)$name)),
                                                           High=degree(graph = msbb_rnaseq2016_PLQ.PosCoexp$STG$High,v = intersect(V(msbb_rnaseq2016_PLQ.PosCoexp$STG$Low)$name,V(msbb_rnaseq2016_PLQ.PosCoexp$STG$High)$name)),
                                                           Difference=degree(graph = msbb_rnaseq2016_PLQ.PosCoexp$STG$Low,v = intersect(V(msbb_rnaseq2016_PLQ.PosCoexp$STG$Low)$name,V(msbb_rnaseq2016_PLQ.PosCoexp$STG$High)$name))-degree(graph = msbb_rnaseq2016_PLQ.PosCoexp$STG$High,v = intersect(V(msbb_rnaseq2016_PLQ.PosCoexp$STG$Low)$name,V(msbb_rnaseq2016_PLQ.PosCoexp$STG$High)$name)),stringsAsFactors = F)
msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF$STG=msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF$STG[order(msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF$STG$Difference,decreasing = T),]
WriteXLS(x = msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF$STG,ExcelFileName = 'STG_PosCoexp_PLQGenes_rewiredDF.xls',SheetNames = 'STG_PosCoexp_PLQGenes',BoldHeaderRow = T)

msbb_rnaseq2016_PLQGenes_NegCoexp.rewiredDF$FP=data.frame(Common_Genes=intersect(V(msbb_rnaseq2016_PLQ.NegCoexp$FP$Low)$name,V(msbb_rnaseq2016_PLQ.NegCoexp$FP$High)$name),
                                                          Low=degree(graph = msbb_rnaseq2016_PLQ.NegCoexp$FP$Low,v = intersect(V(msbb_rnaseq2016_PLQ.NegCoexp$FP$Low)$name,V(msbb_rnaseq2016_PLQ.NegCoexp$FP$High)$name)),
                                                          High=degree(graph = msbb_rnaseq2016_PLQ.NegCoexp$FP$High,v = intersect(V(msbb_rnaseq2016_PLQ.NegCoexp$FP$Low)$name,V(msbb_rnaseq2016_PLQ.NegCoexp$FP$High)$name)),
                                                          Difference=degree(graph = msbb_rnaseq2016_PLQ.NegCoexp$FP$Low,v = intersect(V(msbb_rnaseq2016_PLQ.NegCoexp$FP$Low)$name,V(msbb_rnaseq2016_PLQ.NegCoexp$FP$High)$name))-degree(graph = msbb_rnaseq2016_PLQ.NegCoexp$FP$High,v = intersect(V(msbb_rnaseq2016_PLQ.NegCoexp$FP$Low)$name,V(msbb_rnaseq2016_PLQ.NegCoexp$FP$High)$name)),stringsAsFactors = F)
msbb_rnaseq2016_PLQGenes_NegCoexp.rewiredDF$FP=msbb_rnaseq2016_PLQGenes_NegCoexp.rewiredDF$FP[order(msbb_rnaseq2016_PLQGenes_NegCoexp.rewiredDF$FP$Difference,decreasing = T),]
WriteXLS(x = msbb_rnaseq2016_PLQGenes_NegCoexp.rewiredDF$FP,ExcelFileName = 'FP_NegCoexp_PLQGenes_rewiredDF.xls',SheetNames = 'FP_NegCoexp_PLQGenes',BoldHeaderRow = T)
msbb_rnaseq2016_PLQGenes_NegCoexp.rewiredDF$IFG=data.frame(Common_Genes=intersect(V(msbb_rnaseq2016_PLQ.NegCoexp$IFG$Low)$name,V(msbb_rnaseq2016_PLQ.NegCoexp$IFG$High)$name),
                                                           Low=degree(graph = msbb_rnaseq2016_PLQ.NegCoexp$IFG$Low,v = intersect(V(msbb_rnaseq2016_PLQ.NegCoexp$IFG$Low)$name,V(msbb_rnaseq2016_PLQ.NegCoexp$IFG$High)$name)),
                                                           High=degree(graph = msbb_rnaseq2016_PLQ.NegCoexp$IFG$High,v = intersect(V(msbb_rnaseq2016_PLQ.NegCoexp$IFG$Low)$name,V(msbb_rnaseq2016_PLQ.NegCoexp$IFG$High)$name)),
                                                           Difference=degree(graph = msbb_rnaseq2016_PLQ.NegCoexp$IFG$Low,v = intersect(V(msbb_rnaseq2016_PLQ.NegCoexp$IFG$Low)$name,V(msbb_rnaseq2016_PLQ.NegCoexp$IFG$High)$name))-degree(graph = msbb_rnaseq2016_PLQ.NegCoexp$IFG$High,v = intersect(V(msbb_rnaseq2016_PLQ.NegCoexp$IFG$Low)$name,V(msbb_rnaseq2016_PLQ.NegCoexp$IFG$High)$name)),stringsAsFactors = F)
msbb_rnaseq2016_PLQGenes_NegCoexp.rewiredDF$IFG=msbb_rnaseq2016_PLQGenes_NegCoexp.rewiredDF$FP[order(msbb_rnaseq2016_PLQGenes_NegCoexp.rewiredDF$IFG$Difference,decreasing = T),]
WriteXLS(x = msbb_rnaseq2016_PLQGenes_NegCoexp.rewiredDF$IFG,ExcelFileName = 'IFG_NegCoexp_PLQGenes_rewiredDF.xls',SheetNames = 'IFG_NegCoexp_PLQGenes',BoldHeaderRow = T)
msbb_rnaseq2016_PLQGenes_NegCoexp.rewiredDF$PHG=data.frame(Common_Genes=intersect(V(msbb_rnaseq2016_PLQ.NegCoexp$PHG$Low)$name,V(msbb_rnaseq2016_PLQ.NegCoexp$PHG$High)$name),
                                                           Low=degree(graph = msbb_rnaseq2016_PLQ.NegCoexp$PHG$Low,v = intersect(V(msbb_rnaseq2016_PLQ.NegCoexp$PHG$Low)$name,V(msbb_rnaseq2016_PLQ.NegCoexp$PHG$High)$name)),
                                                           High=degree(graph = msbb_rnaseq2016_PLQ.NegCoexp$PHG$High,v = intersect(V(msbb_rnaseq2016_PLQ.NegCoexp$PHG$Low)$name,V(msbb_rnaseq2016_PLQ.NegCoexp$PHG$High)$name)),
                                                           Difference=degree(graph = msbb_rnaseq2016_PLQ.NegCoexp$PHG$Low,v = intersect(V(msbb_rnaseq2016_PLQ.NegCoexp$PHG$Low)$name,V(msbb_rnaseq2016_PLQ.NegCoexp$PHG$High)$name))-degree(graph = msbb_rnaseq2016_PLQ.NegCoexp$PHG$High,v = intersect(V(msbb_rnaseq2016_PLQ.NegCoexp$PHG$Low)$name,V(msbb_rnaseq2016_PLQ.NegCoexp$PHG$High)$name)),stringsAsFactors = F)
msbb_rnaseq2016_PLQGenes_NegCoexp.rewiredDF$PHG=msbb_rnaseq2016_PLQGenes_NegCoexp.rewiredDF$PHG[order(msbb_rnaseq2016_PLQGenes_NegCoexp.rewiredDF$PHG$Difference,decreasing = T),]
WriteXLS(x = msbb_rnaseq2016_PLQGenes_NegCoexp.rewiredDF$PHG,ExcelFileName = 'PHG_NegCoexp_PLQGenes_rewiredDF.xls',SheetNames = 'PHG_NegCoexp_PLQGenes',BoldHeaderRow = T)
msbb_rnaseq2016_PLQGenes_NegCoexp.rewiredDF$STG=data.frame(Common_Genes=intersect(V(msbb_rnaseq2016_PLQ.NegCoexp$STG$Low)$name,V(msbb_rnaseq2016_PLQ.NegCoexp$STG$High)$name),
                                                           Low=degree(graph = msbb_rnaseq2016_PLQ.NegCoexp$STG$Low,v = intersect(V(msbb_rnaseq2016_PLQ.NegCoexp$STG$Low)$name,V(msbb_rnaseq2016_PLQ.NegCoexp$STG$High)$name)),
                                                           High=degree(graph = msbb_rnaseq2016_PLQ.NegCoexp$STG$High,v = intersect(V(msbb_rnaseq2016_PLQ.NegCoexp$STG$Low)$name,V(msbb_rnaseq2016_PLQ.NegCoexp$STG$High)$name)),
                                                           Difference=degree(graph = msbb_rnaseq2016_PLQ.NegCoexp$STG$Low,v = intersect(V(msbb_rnaseq2016_PLQ.NegCoexp$STG$Low)$name,V(msbb_rnaseq2016_PLQ.NegCoexp$STG$High)$name))-degree(graph = msbb_rnaseq2016_PLQ.NegCoexp$STG$High,v = intersect(V(msbb_rnaseq2016_PLQ.NegCoexp$STG$Low)$name,V(msbb_rnaseq2016_PLQ.NegCoexp$STG$High)$name)),stringsAsFactors = F)
msbb_rnaseq2016_PLQGenes_NegCoexp.rewiredDF$STG=msbb_rnaseq2016_PLQGenes_NegCoexp.rewiredDF$STG[order(msbb_rnaseq2016_PLQGenes_NegCoexp.rewiredDF$STG$Difference,decreasing = T),]
WriteXLS(x = msbb_rnaseq2016_PLQGenes_NegCoexp.rewiredDF$STG,ExcelFileName = 'STG_NegCoexp_PLQGenes_rewiredDF.xls',SheetNames = 'STG_NegCoexp_PLQGenes',BoldHeaderRow = T)
# #For R2
# 
# 
# msbb_rnaseq2016_R2Genes_NegCoexp.rewiredDF$FP=data.frame(Common_Genes=intersect(V(msbb_rnaseq2016_R2.NegCoexp$FP$Low)$name,V(msbb_rnaseq2016_R2.NegCoexp$FP$High)$name),
#                                                          Low=degree(graph = msbb_rnaseq2016_R2.NegCoexp$FP$Low,v = intersect(V(msbb_rnaseq2016_R2.NegCoexp$FP$Low)$name,V(msbb_rnaseq2016_R2.NegCoexp$FP$High)$name)),
#                                                          High=degree(graph = msbb_rnaseq2016_R2.NegCoexp$FP$High,v = intersect(V(msbb_rnaseq2016_R2.NegCoexp$FP$Low)$name,V(msbb_rnaseq2016_R2.NegCoexp$FP$High)$name)),
#                                                          Difference=degree(graph = msbb_rnaseq2016_R2.NegCoexp$FP$Low,v = intersect(V(msbb_rnaseq2016_R2.NegCoexp$FP$Low)$name,V(msbb_rnaseq2016_R2.NegCoexp$FP$High)$name))-degree(graph = msbb_rnaseq2016_R2.NegCoexp$FP$High,v = intersect(V(msbb_rnaseq2016_R2.NegCoexp$FP$Low)$name,V(msbb_rnaseq2016_R2.NegCoexp$FP$High)$name)),stringsAsFactors = F)
# msbb_rnaseq2016_R2Genes_NegCoexp.rewiredDF$IFG=data.frame(Common_Genes=intersect(V(msbb_rnaseq2016_R2.NegCoexp$IFG$Low)$name,V(msbb_rnaseq2016_R2.NegCoexp$IFG$High)$name),
#                                                           Low=degree(graph = msbb_rnaseq2016_R2.NegCoexp$IFG$Low,v = intersect(V(msbb_rnaseq2016_R2.NegCoexp$IFG$Low)$name,V(msbb_rnaseq2016_R2.NegCoexp$IFG$High)$name)),
#                                                           High=degree(graph = msbb_rnaseq2016_R2.NegCoexp$IFG$High,v = intersect(V(msbb_rnaseq2016_R2.NegCoexp$IFG$Low)$name,V(msbb_rnaseq2016_R2.NegCoexp$IFG$High)$name)),
#                                                           Difference=degree(graph = msbb_rnaseq2016_R2.NegCoexp$IFG$Low,v = intersect(V(msbb_rnaseq2016_R2.NegCoexp$IFG$Low)$name,V(msbb_rnaseq2016_R2.NegCoexp$IFG$High)$name))-degree(graph = msbb_rnaseq2016_R2.NegCoexp$IFG$High,v = intersect(V(msbb_rnaseq2016_R2.NegCoexp$IFG$Low)$name,V(msbb_rnaseq2016_R2.NegCoexp$IFG$High)$name)),stringsAsFactors = F)
# msbb_rnaseq2016_R2Genes_NegCoexp.rewiredDF$PHG=data.frame(Common_Genes=intersect(V(msbb_rnaseq2016_R2.NegCoexp$PHG$Low)$name,V(msbb_rnaseq2016_R2.NegCoexp$PHG$High)$name),
#                                                           Low=degree(graph = msbb_rnaseq2016_R2.NegCoexp$PHG$Low,v = intersect(V(msbb_rnaseq2016_R2.NegCoexp$PHG$Low)$name,V(msbb_rnaseq2016_R2.NegCoexp$PHG$High)$name)),
#                                                           High=degree(graph = msbb_rnaseq2016_R2.NegCoexp$PHG$High,v = intersect(V(msbb_rnaseq2016_R2.NegCoexp$PHG$Low)$name,V(msbb_rnaseq2016_R2.NegCoexp$PHG$High)$name)),
#                                                           Difference=degree(graph = msbb_rnaseq2016_R2.NegCoexp$PHG$Low,v = intersect(V(msbb_rnaseq2016_R2.NegCoexp$PHG$Low)$name,V(msbb_rnaseq2016_R2.NegCoexp$PHG$High)$name))-degree(graph = msbb_rnaseq2016_R2.NegCoexp$PHG$High,v = intersect(V(msbb_rnaseq2016_R2.NegCoexp$PHG$Low)$name,V(msbb_rnaseq2016_R2.NegCoexp$PHG$High)$name)),stringsAsFactors = F)
# msbb_rnaseq2016_R2Genes_NegCoexp.rewiredDF$STG=data.frame(Common_Genes=intersect(V(msbb_rnaseq2016_R2.NegCoexp$STG$Low)$name,V(msbb_rnaseq2016_R2.NegCoexp$STG$High)$name),
#                                                           Low=degree(graph = msbb_rnaseq2016_R2.NegCoexp$STG$Low,v = intersect(V(msbb_rnaseq2016_R2.NegCoexp$STG$Low)$name,V(msbb_rnaseq2016_R2.NegCoexp$STG$High)$name)),
#                                                           High=degree(graph = msbb_rnaseq2016_R2.NegCoexp$STG$High,v = intersect(V(msbb_rnaseq2016_R2.NegCoexp$STG$Low)$name,V(msbb_rnaseq2016_R2.NegCoexp$STG$High)$name)),
#                                                           Difference=degree(graph = msbb_rnaseq2016_R2.NegCoexp$STG$Low,v = intersect(V(msbb_rnaseq2016_R2.NegCoexp$STG$Low)$name,V(msbb_rnaseq2016_R2.NegCoexp$STG$High)$name))-degree(graph = msbb_rnaseq2016_R2.NegCoexp$STG$High,v = intersect(V(msbb_rnaseq2016_R2.NegCoexp$STG$Low)$name,V(msbb_rnaseq2016_R2.NegCoexp$STG$High)$name)),stringsAsFactors = F)

#Read  unbiased DiffCorr to construct DCNs and perform rewiring
msbb_rnaseq2016.PosCoexp=msbb_rnaseq2016.NegCoexp=msbb_rnaseq2016Genes_PosCoexp.rewiredDF=msbb_rnaseq2016Genes_NegCoexp.rewiredDF=vector(mode = "list",length = 4)
names(msbb_rnaseq2016.PosCoexp)=names(msbb_rnaseq2016.NegCoexp)=names(msbb_rnaseq2016Genes_PosCoexp.rewiredDF)=names(msbb_rnaseq2016Genes_NegCoexp.rewiredDF)
msbb_rnaseq2016.PosCoexp$FP=msbb_rnaseq2016.PosCoexp$IFG=msbb_rnaseq2016.PosCoexp$PHG=msbb_rnaseq2016.PosCoexp$STG=vector(mode = "list",length = 2)
msbb_rnaseq2016.NegCoexp$FP=msbb_rnaseq2016.NegCoexp$IFG=msbb_rnaseq2016.NegCoexp$PHG=msbb_rnaseq2016.NegCoexp$STG=vector(mode = "list",length = 2)
names(msbb_rnaseq2016.PosCoexp$FP)=names(msbb_rnaseq2016.NegCoexp$FP)=c("Low","High")
names(msbb_rnaseq2016.PosCoexp$IFG)=names(msbb_rnaseq2016.NegCoexp$IFG)=c("Low","High")
names(msbb_rnaseq2016.PosCoexp$PHG)=names(msbb_rnaseq2016.NegCoexp$PHG)=c("Low","High")
names(msbb_rnaseq2016.PosCoexp$STG)=names(msbb_rnaseq2016.NegCoexp$STG)=c("Low","High")

msbb_rnaseq2016.PosCoexp$FP$Low=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016.DCN$FP$Low,es = which(E(msbb_rnaseq2016.DCN$FP$Low)$weight>1.5),names = T)),directed = F)
E(msbb_rnaseq2016.PosCoexp$FP$Low)$weight=E(msbb_rnaseq2016.DCN$FP$Low)$weight[which(E(msbb_rnaseq2016.DCN$FP$Low)$weight>1.5)]
write.graph(graph = msbb_rnaseq2016.PosCoexp$FP$Low,file = 'FPGenes_PosCoexp_Low.gml',format = 'gml')
msbb_rnaseq2016.PosCoexp$IFG$Low=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016.DCN$IFG$Low,es = which(E(msbb_rnaseq2016.DCN$IFG$Low)$weight>1.5),names = T)),directed = F)
E(msbb_rnaseq2016.PosCoexp$IFG$Low)$weight=E(msbb_rnaseq2016.DCN$IFG$Low)$weight[which(E(msbb_rnaseq2016.DCN$IFG$Low)$weight>1.5)]
write.graph(graph = msbb_rnaseq2016.PosCoexp$IFG$Low,file = 'IFGGenes_PosCoexp_Low.gml',format = 'gml')
msbb_rnaseq2016.PosCoexp$PHG$Low=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016.DCN$PHG$Low,es = which(E(msbb_rnaseq2016.DCN$PHG$Low)$weight>1.5),names = T)),directed = F)
E(msbb_rnaseq2016.PosCoexp$PHG$Low)$weight=E(msbb_rnaseq2016.DCN$PHG$Low)$weight[which(E(msbb_rnaseq2016.DCN$PHG$Low)$weight>1.5)]
write.graph(graph = msbb_rnaseq2016.PosCoexp$PHG$Low,file = 'PHGGenes_PosCoexp_Low.gml',format = 'gml')
msbb_rnaseq2016.PosCoexp$STG$Low=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016.DCN$STG$Low,es = which(E(msbb_rnaseq2016.DCN$STG$Low)$weight>1.5),names = T)),directed = F)
E(msbb_rnaseq2016.PosCoexp$STG$Low)$weight=E(msbb_rnaseq2016.DCN$STG$Low)$weight[which(E(msbb_rnaseq2016.DCN$STG$Low)$weight>1.5)]
write.graph(graph = msbb_rnaseq2016.PosCoexp$STG$Low,file = 'STGGenes_PosCoexp_Low.gml',format = 'gml')

msbb_rnaseq2016.PosCoexp$FP$High=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016.DCN$FP$High,es = which(E(msbb_rnaseq2016.DCN$FP$High)$weight>1.5),names = T)),directed = F)
E(msbb_rnaseq2016.PosCoexp$FP$High)$weight=E(msbb_rnaseq2016.DCN$FP$High)$weight[which(E(msbb_rnaseq2016.DCN$FP$High)$weight>1.5)]
write.graph(graph = msbb_rnaseq2016.PosCoexp$FP$High,file = 'FPGenes_PosCoexp_High.gml',format = 'gml')
msbb_rnaseq2016.PosCoexp$IFG$High=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016.DCN$IFG$High,es = which(E(msbb_rnaseq2016.DCN$IFG$High)$weight>1.5),names = T)),directed = F)
E(msbb_rnaseq2016.PosCoexp$IFG$High)$weight=E(msbb_rnaseq2016.DCN$IFG$High)$weight[which(E(msbb_rnaseq2016.DCN$IFG$High)$weight>1.5)]
write.graph(graph = msbb_rnaseq2016.PosCoexp$IFG$High,file = 'IFGGenes_PosCoexp_High.gml',format = 'gml')
msbb_rnaseq2016.PosCoexp$PHG$High=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016.DCN$PHG$High,es = which(E(msbb_rnaseq2016.DCN$PHG$High)$weight>1.5),names = T)),directed = F)
E(msbb_rnaseq2016.PosCoexp$PHG$High)$weight=E(msbb_rnaseq2016.DCN$PHG$High)$weight[which(E(msbb_rnaseq2016.DCN$PHG$High)$weight>1.5)]
write.graph(graph = msbb_rnaseq2016.PosCoexp$PHG$High,file = 'PHGGenes_PosCoexp_High.gml',format = 'gml')
msbb_rnaseq2016.PosCoexp$STG$High=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016.DCN$STG$High,es = which(E(msbb_rnaseq2016.DCN$STG$High)$weight>1.5),names = T)),directed = F)
E(msbb_rnaseq2016.PosCoexp$STG$High)$weight=E(msbb_rnaseq2016.DCN$STG$High)$weight[which(E(msbb_rnaseq2016.DCN$STG$High)$weight>1.5)]
write.graph(graph = msbb_rnaseq2016.PosCoexp$STG$High,file = 'STGGenes_PosCoexp_High.gml',format = 'gml')

msbb_rnaseq2016.NegCoexp$FP$High=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016.DCN$FP$High,es = which(E(msbb_rnaseq2016.DCN$FP$High)$weight<0.5),names = T)),directed = F)
E(msbb_rnaseq2016.NegCoexp$FP$High)$weight=E(msbb_rnaseq2016.DCN$FP$High)$weight[which(E(msbb_rnaseq2016.DCN$FP$High)$weight<0.5)]
write.graph(graph = msbb_rnaseq2016.NegCoexp$FP$High,file = 'FPGenes_NegCoexp_High.gml',format = 'gml')
msbb_rnaseq2016.NegCoexp$IFG$High=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016.DCN$IFG$High,es = which(E(msbb_rnaseq2016.DCN$IFG$High)$weight<0.5),names = T)),directed = F)
E(msbb_rnaseq2016.NegCoexp$IFG$High)$weight=E(msbb_rnaseq2016.DCN$IFG$High)$weight[which(E(msbb_rnaseq2016.DCN$IFG$High)$weight<0.5)]
write.graph(graph = msbb_rnaseq2016.NegCoexp$IFG$High,file = 'IFGGenes_NegCoexp_High.gml',format = 'gml')
msbb_rnaseq2016.NegCoexp$PHG$High=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016.DCN$PHG$High,es = which(E(msbb_rnaseq2016.DCN$PHG$High)$weight<0.5),names = T)),directed = F)
E(msbb_rnaseq2016.NegCoexp$PHG$High)$weight=E(msbb_rnaseq2016.DCN$PHG$High)$weight[which(E(msbb_rnaseq2016.DCN$PHG$High)$weight<0.5)]
write.graph(graph = msbb_rnaseq2016.NegCoexp$PHG$High,file = 'PHGGenes_NegCoexp_High.gml',format = 'gml')
msbb_rnaseq2016.NegCoexp$STG$High=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016.DCN$STG$High,es = which(E(msbb_rnaseq2016.DCN$STG$High)$weight<0.5),names = T)),directed = F)
E(msbb_rnaseq2016.NegCoexp$STG$High)$weight=E(msbb_rnaseq2016.DCN$STG$High)$weight[which(E(msbb_rnaseq2016.DCN$STG$High)$weight<0.5)]
write.graph(graph = msbb_rnaseq2016.NegCoexp$STG$High,file = 'STGGenes_NegCoexp_High.gml',format = 'gml')

msbb_rnaseq2016.NegCoexp$FP$Low=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016.DCN$FP$Low,es = which(E(msbb_rnaseq2016.DCN$FP$Low)$weight<0.5),names = T)),directed = F)
E(msbb_rnaseq2016.NegCoexp$FP$Low)$weight=E(msbb_rnaseq2016.DCN$FP$Low)$weight[which(E(msbb_rnaseq2016.DCN$FP$Low)$weight<0.5)]
write.graph(graph = msbb_rnaseq2016.NegCoexp$FP$Low,file = 'FPGenes_NegCoexp_Low.gml',format = 'gml')
msbb_rnaseq2016.NegCoexp$IFG$Low=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016.DCN$IFG$Low,es = which(E(msbb_rnaseq2016.DCN$IFG$Low)$weight<0.5),names = T)),directed = F)
E(msbb_rnaseq2016.NegCoexp$IFG$Low)$weight=E(msbb_rnaseq2016.DCN$IFG$Low)$weight[which(E(msbb_rnaseq2016.DCN$IFG$Low)$weight<0.5)]
write.graph(graph = msbb_rnaseq2016.NegCoexp$IFG$Low,file = 'IFGGenes_NegCoexp_Low.gml',format = 'gml')
msbb_rnaseq2016.NegCoexp$PHG$Low=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016.DCN$PHG$Low,es = which(E(msbb_rnaseq2016.DCN$PHG$Low)$weight<0.5),names = T)),directed = F)
E(msbb_rnaseq2016.NegCoexp$PHG$Low)$weight=E(msbb_rnaseq2016.DCN$PHG$Low)$weight[which(E(msbb_rnaseq2016.DCN$PHG$Low)$weight<0.5)]
write.graph(graph = msbb_rnaseq2016.NegCoexp$PHG$Low,file = 'PHGGenes_NegCoexp_Low.gml',format = 'gml')
msbb_rnaseq2016.NegCoexp$STG$Low=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016.DCN$STG$Low,es = which(E(msbb_rnaseq2016.DCN$STG$Low)$weight<0.5),names = T)),directed = F)
E(msbb_rnaseq2016.NegCoexp$STG$Low)$weight=E(msbb_rnaseq2016.DCN$STG$Low)$weight[which(E(msbb_rnaseq2016.DCN$STG$Low)$weight<0.5)]
write.graph(graph = msbb_rnaseq2016.NegCoexp$STG$Low,file = 'STGGenes_NegCoexp_Low.gml',format = 'gml')
msbb_rnaseq2016Genes_PosCoexp.rewiredDF$FP=data.frame(Common_Genes=intersect(V(msbb_rnaseq2016.PosCoexp$FP$Low)$name,V(msbb_rnaseq2016.PosCoexp$FP$High)$name),
                                                      Low=degree(graph = msbb_rnaseq2016.PosCoexp$FP$Low,v = intersect(V(msbb_rnaseq2016.PosCoexp$FP$Low)$name,V(msbb_rnaseq2016.PosCoexp$FP$High)$name)),
                                                      High=degree(graph = msbb_rnaseq2016.PosCoexp$FP$High,v = intersect(V(msbb_rnaseq2016.PosCoexp$FP$Low)$name,V(msbb_rnaseq2016.PosCoexp$FP$High)$name)),
                                                      Difference=degree(graph = msbb_rnaseq2016.PosCoexp$FP$Low,v = intersect(V(msbb_rnaseq2016.PosCoexp$FP$Low)$name,V(msbb_rnaseq2016.PosCoexp$FP$High)$name))-degree(graph = msbb_rnaseq2016.PosCoexp$FP$High,v = intersect(V(msbb_rnaseq2016.PosCoexp$FP$Low)$name,V(msbb_rnaseq2016.PosCoexp$FP$High)$name)),stringsAsFactors = F)
msbb_rnaseq2016Genes_PosCoexp.rewiredDF$FP=msbb_rnaseq2016Genes_PosCoexp.rewiredDF$FP[order(msbb_rnaseq2016Genes_PosCoexp.rewiredDF$FP$Difference,decreasing = T),]
WriteXLS(x = msbb_rnaseq2016Genes_PosCoexp.rewiredDF$FP,ExcelFileName = 'FP_PosCoexp_rewiredDF.xls',SheetNames = 'FP_PosCoexp',BoldHeaderRow = T)
msbb_rnaseq2016Genes_PosCoexp.rewiredDF$IFG=data.frame(Common_Genes=intersect(V(msbb_rnaseq2016.PosCoexp$IFG$Low)$name,V(msbb_rnaseq2016.PosCoexp$IFG$High)$name),
                                                       Low=degree(graph = msbb_rnaseq2016.PosCoexp$IFG$Low,v = intersect(V(msbb_rnaseq2016.PosCoexp$IFG$Low)$name,V(msbb_rnaseq2016.PosCoexp$IFG$High)$name)),
                                                       High=degree(graph = msbb_rnaseq2016.PosCoexp$IFG$High,v = intersect(V(msbb_rnaseq2016.PosCoexp$IFG$Low)$name,V(msbb_rnaseq2016.PosCoexp$IFG$High)$name)),
                                                       Difference=degree(graph = msbb_rnaseq2016.PosCoexp$IFG$Low,v = intersect(V(msbb_rnaseq2016.PosCoexp$IFG$Low)$name,V(msbb_rnaseq2016.PosCoexp$IFG$High)$name))-degree(graph = msbb_rnaseq2016.PosCoexp$IFG$High,v = intersect(V(msbb_rnaseq2016.PosCoexp$IFG$Low)$name,V(msbb_rnaseq2016.PosCoexp$IFG$High)$name)),stringsAsFactors = F)
msbb_rnaseq2016Genes_PosCoexp.rewiredDF$IFG=msbb_rnaseq2016Genes_PosCoexp.rewiredDF$IFG[order(msbb_rnaseq2016Genes_PosCoexp.rewiredDF$IFG$Difference,decreasing = T),]
WriteXLS(x = msbb_rnaseq2016Genes_PosCoexp.rewiredDF$IFG,ExcelFileName = 'IFG_PosCoexp_rewiredDF.xls',SheetNames = 'IFG_PosCoexp',BoldHeaderRow = T)
msbb_rnaseq2016Genes_PosCoexp.rewiredDF$PHG=data.frame(Common_Genes=intersect(V(msbb_rnaseq2016.PosCoexp$PHG$Low)$name,V(msbb_rnaseq2016.PosCoexp$PHG$High)$name),
                                                       Low=degree(graph = msbb_rnaseq2016.PosCoexp$PHG$Low,v = intersect(V(msbb_rnaseq2016.PosCoexp$PHG$Low)$name,V(msbb_rnaseq2016.PosCoexp$PHG$High)$name)),
                                                       High=degree(graph = msbb_rnaseq2016.PosCoexp$PHG$High,v = intersect(V(msbb_rnaseq2016.PosCoexp$PHG$Low)$name,V(msbb_rnaseq2016.PosCoexp$PHG$High)$name)),
                                                       Difference=degree(graph = msbb_rnaseq2016.PosCoexp$PHG$Low,v = intersect(V(msbb_rnaseq2016.PosCoexp$PHG$Low)$name,V(msbb_rnaseq2016.PosCoexp$PHG$High)$name))-degree(graph = msbb_rnaseq2016.PosCoexp$PHG$High,v = intersect(V(msbb_rnaseq2016.PosCoexp$PHG$Low)$name,V(msbb_rnaseq2016.PosCoexp$PHG$High)$name)),stringsAsFactors = F)
msbb_rnaseq2016Genes_PosCoexp.rewiredDF$PHG=msbb_rnaseq2016Genes_PosCoexp.rewiredDF$PHG[order(msbb_rnaseq2016Genes_PosCoexp.rewiredDF$PHG$Difference,decreasing = T),]
WriteXLS(x = msbb_rnaseq2016Genes_PosCoexp.rewiredDF$PHG,ExcelFileName = 'PHG_PosCoexp_rewiredDF.xls',SheetNames = 'PHG_PosCoexp',BoldHeaderRow = T)
msbb_rnaseq2016Genes_PosCoexp.rewiredDF$STG=data.frame(Common_Genes=intersect(V(msbb_rnaseq2016.PosCoexp$STG$Low)$name,V(msbb_rnaseq2016.PosCoexp$STG$High)$name),
                                                       Low=degree(graph = msbb_rnaseq2016.PosCoexp$STG$Low,v = intersect(V(msbb_rnaseq2016.PosCoexp$STG$Low)$name,V(msbb_rnaseq2016.PosCoexp$STG$High)$name)),
                                                       High=degree(graph = msbb_rnaseq2016.PosCoexp$STG$High,v = intersect(V(msbb_rnaseq2016.PosCoexp$STG$Low)$name,V(msbb_rnaseq2016.PosCoexp$STG$High)$name)),
                                                       Difference=degree(graph = msbb_rnaseq2016.PosCoexp$STG$Low,v = intersect(V(msbb_rnaseq2016.PosCoexp$STG$Low)$name,V(msbb_rnaseq2016.PosCoexp$STG$High)$name))-degree(graph = msbb_rnaseq2016.PosCoexp$STG$High,v = intersect(V(msbb_rnaseq2016.PosCoexp$STG$Low)$name,V(msbb_rnaseq2016.PosCoexp$STG$High)$name)),stringsAsFactors = F)
msbb_rnaseq2016Genes_PosCoexp.rewiredDF$STG=msbb_rnaseq2016Genes_PosCoexp.rewiredDF$STG[order(msbb_rnaseq2016Genes_PosCoexp.rewiredDF$STG$Difference,decreasing = T),]
WriteXLS(x = msbb_rnaseq2016Genes_PosCoexp.rewiredDF$STG,ExcelFileName = 'STG_PosCoexp_rewiredDF.xls',SheetNames = 'STG_PosCoexp',BoldHeaderRow = T)

msbb_rnaseq2016Genes_NegCoexp.rewiredDF$FP=data.frame(Common_Genes=intersect(V(msbb_rnaseq2016.NegCoexp$FP$Low)$name,V(msbb_rnaseq2016.NegCoexp$FP$High)$name),
                                                      Low=degree(graph = msbb_rnaseq2016.NegCoexp$FP$Low,v = intersect(V(msbb_rnaseq2016.NegCoexp$FP$Low)$name,V(msbb_rnaseq2016.NegCoexp$FP$High)$name)),
                                                      High=degree(graph = msbb_rnaseq2016.NegCoexp$FP$High,v = intersect(V(msbb_rnaseq2016.NegCoexp$FP$Low)$name,V(msbb_rnaseq2016.NegCoexp$FP$High)$name)),
                                                      Difference=degree(graph = msbb_rnaseq2016.NegCoexp$FP$Low,v = intersect(V(msbb_rnaseq2016.NegCoexp$FP$Low)$name,V(msbb_rnaseq2016.NegCoexp$FP$High)$name))-degree(graph = msbb_rnaseq2016.NegCoexp$FP$High,v = intersect(V(msbb_rnaseq2016.NegCoexp$FP$Low)$name,V(msbb_rnaseq2016.NegCoexp$FP$High)$name)),stringsAsFactors = F)
msbb_rnaseq2016Genes_NegCoexp.rewiredDF$FP=msbb_rnaseq2016Genes_NegCoexp.rewiredDF$FP[order(msbb_rnaseq2016Genes_NegCoexp.rewiredDF$FP$Difference,decreasing = T),]
WriteXLS(x = msbb_rnaseq2016_NegCoexp.rewiredDF$FP,ExcelFileName = 'FP_NegCoexp_rewiredDF.xls',SheetNames = 'FP_NegCoexp',BoldHeaderRow = T)
msbb_rnaseq2016Genes_NegCoexp.rewiredDF$IFG=data.frame(Common_Genes=intersect(V(msbb_rnaseq2016.NegCoexp$IFG$Low)$name,V(msbb_rnaseq2016.NegCoexp$IFG$High)$name),
                                                       Low=degree(graph = msbb_rnaseq2016.NegCoexp$IFG$Low,v = intersect(V(msbb_rnaseq2016.NegCoexp$IFG$Low)$name,V(msbb_rnaseq2016.NegCoexp$IFG$High)$name)),
                                                       High=degree(graph = msbb_rnaseq2016.NegCoexp$IFG$High,v = intersect(V(msbb_rnaseq2016.NegCoexp$IFG$Low)$name,V(msbb_rnaseq2016.NegCoexp$IFG$High)$name)),
                                                       Difference=degree(graph = msbb_rnaseq2016.NegCoexp$IFG$Low,v = intersect(V(msbb_rnaseq2016.NegCoexp$IFG$Low)$name,V(msbb_rnaseq2016.NegCoexp$IFG$High)$name))-degree(graph = msbb_rnaseq2016.NegCoexp$IFG$High,v = intersect(V(msbb_rnaseq2016.NegCoexp$IFG$Low)$name,V(msbb_rnaseq2016.NegCoexp$IFG$High)$name)),stringsAsFactors = F)
msbb_rnaseq2016Genes_NegCoexp.rewiredDF$IFG=msbb_rnaseq2016Genes_NegCoexp.rewiredDF$FP[order(msbb_rnaseq2016Genes_NegCoexp.rewiredDF$IFG$Difference,decreasing = T),]
WriteXLS(x = msbb_rnaseq2016Genes_NegCoexp.rewiredDF$IFG,ExcelFileName = 'IFG_NegCoexp_rewiredDF.xls',SheetNames = 'IFG_NegCoexp',BoldHeaderRow = T)
msbb_rnaseq2016Genes_NegCoexp.rewiredDF$PHG=data.frame(Common_Genes=intersect(V(msbb_rnaseq2016.NegCoexp$PHG$Low)$name,V(msbb_rnaseq2016.NegCoexp$PHG$High)$name),
                                                       Low=degree(graph = msbb_rnaseq2016.NegCoexp$PHG$Low,v = intersect(V(msbb_rnaseq2016.NegCoexp$PHG$Low)$name,V(msbb_rnaseq2016.NegCoexp$PHG$High)$name)),
                                                       High=degree(graph = msbb_rnaseq2016.NegCoexp$PHG$High,v = intersect(V(msbb_rnaseq2016.NegCoexp$PHG$Low)$name,V(msbb_rnaseq2016.NegCoexp$PHG$High)$name)),
                                                       Difference=degree(graph = msbb_rnaseq2016.NegCoexp$PHG$Low,v = intersect(V(msbb_rnaseq2016.NegCoexp$PHG$Low)$name,V(msbb_rnaseq2016.NegCoexp$PHG$High)$name))-degree(graph = msbb_rnaseq2016.NegCoexp$PHG$High,v = intersect(V(msbb_rnaseq2016.NegCoexp$PHG$Low)$name,V(msbb_rnaseq2016.NegCoexp$PHG$High)$name)),stringsAsFactors = F)
msbb_rnaseq2016Genes_NegCoexp.rewiredDF$PHG=msbb_rnaseq2016Genes_NegCoexp.rewiredDF$PHG[order(msbb_rnaseq2016Genes_NegCoexp.rewiredDF$PHG$Difference,decreasing = T),]
WriteXLS(x = msbb_rnaseq2016Genes_NegCoexp.rewiredDF$PHG,ExcelFileName = 'PHG_NegCoexp_rewiredDF.xls',SheetNames = 'PHG_NegCoexp',BoldHeaderRow = T)
msbb_rnaseq2016Genes_NegCoexp.rewiredDF$STG=data.frame(Common_Genes=intersect(V(msbb_rnaseq2016.NegCoexp$STG$Low)$name,V(msbb_rnaseq2016.NegCoexp$STG$High)$name),
                                                       Low=degree(graph = msbb_rnaseq2016.NegCoexp$STG$Low,v = intersect(V(msbb_rnaseq2016.NegCoexp$STG$Low)$name,V(msbb_rnaseq2016.NegCoexp$STG$High)$name)),
                                                       High=degree(graph = msbb_rnaseq2016.NegCoexp$STG$High,v = intersect(V(msbb_rnaseq2016.NegCoexp$STG$Low)$name,V(msbb_rnaseq2016.NegCoexp$STG$High)$name)),
                                                       Difference=degree(graph = msbb_rnaseq2016.NegCoexp$STG$Low,v = intersect(V(msbb_rnaseq2016.NegCoexp$STG$Low)$name,V(msbb_rnaseq2016.NegCoexp$STG$High)$name))-degree(graph = msbb_rnaseq2016.NegCoexp$STG$High,v = intersect(V(msbb_rnaseq2016.NegCoexp$STG$Low)$name,V(msbb_rnaseq2016.NegCoexp$STG$High)$name)),stringsAsFactors = F)
msbb_rnaseq2016Genes_NegCoexp.rewiredDF$STG=msbb_rnaseq2016Genes_NegCoexp.rewiredDF$STG[order(msbb_rnaseq2016Genes_NegCoexp.rewiredDF$STG$Difference,decreasing = T),]
WriteXLS(x = msbb_rnaseq2016Genes_NegCoexp.rewiredDF$STG,ExcelFileName = 'STG_NegCoexp_rewiredDF.xls',SheetNames = 'STG_NegCoexp',BoldHeaderRow = T)


#Read  in DEG results
msbb_rnaseq2016_PLQGenes.DEG=msbb_rnaseq2016_PLQGenes=msbb_rnaseq2016_AllAnalyses=vector(mode = "list",length = 4)
names(msbb_rnaseq2016_PLQGenes.DEG)=names(msbb_rnaseq2016_PLQGenes)=names(msbb_rnaseq2016_AllAnalyses)=names(msbb_rnaseq2016_PLQ.DCN)
msbb_rnaseq2016_PLQGenes.DEG$FP=read.table("DEG_Genes/MSBB_RNAseq2016_FP_DEG.txt",sep = "\t",header = T,as.is = T)
msbb_rnaseq2016_PLQGenes.DEG$IFG=read.table("DEG_Genes/MSBB_RNAseq2016_IFG_DEG.txt",sep = "\t",header = T,as.is = T)
msbb_rnaseq2016_PLQGenes.DEG$PHG=read.table("DEG_Genes/MSBB_RNAseq2016_PHG_DEG.txt",sep = "\t",header = T,as.is = T)
msbb_rnaseq2016_PLQGenes.DEG$STG=read.table("DEG_Genes/MSBB_RNAseq2016_STG_DEG.txt",sep = "\t",header = T,as.is = T)

msbb_rnaseq2016_PLQGenes$FP=read.table("PLQ_Assoc_Genes/FP_allPLQ_AssocGenes.txt",sep = "\t",header = T,as.is = T)
msbb_rnaseq2016_PLQGenes$FP$FDR=p.adjust(p=msbb_rnaseq2016_PLQGenes$FP$Rho.p,method = "fdr")
msbb_rnaseq2016_PLQGenes$IFG=read.table("PLQ_Assoc_Genes/IFG_allPLQ_AssocGenes.txt",sep = "\t",header = T,as.is = T)
msbb_rnaseq2016_PLQGenes$IFG$FDR=p.adjust(p=msbb_rnaseq2016_PLQGenes$IFG$Rho.p,method = "fdr")
msbb_rnaseq2016_PLQGenes$PHG=read.table("PLQ_Assoc_Genes/PHG_allPLQ_AssocGenes.txt",sep = "\t",header = T,as.is = T)
msbb_rnaseq2016_PLQGenes$PHG$FDR=p.adjust(p=msbb_rnaseq2016_PLQGenes$PHG$Rho.p,method = "fdr")
msbb_rnaseq2016_PLQGenes$STG=read.table("PLQ_Assoc_Genes/STG_allPLQ_AssocGenes.txt",sep = "\t",header = T,as.is = T)
msbb_rnaseq2016_PLQGenes$STG$FDR=p.adjust(p=msbb_rnaseq2016_PLQGenes$STG$Rho.p,method = "fdr")

#Perform comparative enrichment analysis using clusterProfiler
msbb_rnaseq2016_PLQGenes.cluster$FP=select(x=org.Hs.eg.db,keys = V(msbb_rnaseq2016_PLQ.DCN$FP$Low)$name,keytype = "SYMBOL",columns = "ENTREZID")[,2]
msbb_rnaseq2016_PLQGenes.cluster$IFG=select(x=org.Hs.eg.db,keys = V(msbb_rnaseq2016_PLQ.DCN$IFG$Low)$name,keytype = "SYMBOL",columns = "ENTREZID")[,2]
msbb_rnaseq2016_PLQGenes.cluster$PHG=select(x=org.Hs.eg.db,keys = V(msbb_rnaseq2016_PLQ.DCN$PHG$Low)$name,keytype = "SYMBOL",columns = "ENTREZID")[,2]
msbb_rnaseq2016_PLQGenes.cluster$STG=select(x=org.Hs.eg.db,keys = V(msbb_rnaseq2016_PLQ.DCN$STG$Low)$name,keytype = "SYMBOL",columns = "ENTREZID")[,2]

msbb_rnaseq2016_AllAnalyses$FP=list(PLQ_Genes=msbb_rnaseq2016_PLQGenes$FP$Genes[which(msbb_rnaseq2016_PLQGenes$FP$FDR<=0.05)],
                                    PLQ_DEG=msbb_rnaseq2016_PLQGenes.DEG$FP$row[which(abs(msbb_rnaseq2016_PLQGenes.DEG$FP$log2FoldChange)>1&msbb_rnaseq2016_PLQGenes.DEG$FP$padj<=0.05)],
                                    PLQ_DCN=V(msbb_rnaseq2016_PLQ.DCN$FP$High)$name,
                                    DCN=V(msbb_rnaseq2016.DCN$FP$High)$name,
                                    PLQ_PosCoexp_Loss=msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF$FP$Common_Genes[which(msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF$FP$Difference>0)],
                                    PLQ_PosCoexp_Gain=msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF$FP$Common_Genes[which(msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF$FP$Difference<0)],
                                    PLQ_NegCoexp_Loss=msbb_rnaseq2016_PLQGenes_NegCoexp.rewiredDF$FP$Common_Genes[which(msbb_rnaseq2016_PLQGenes_NegCoexp.rewiredDF$FP$Difference>0)],
                                    PLQ_NegCoexp_Gain=msbb_rnaseq2016_PLQGenes_NegCoexp.rewiredDF$FP$Common_Genes[which(msbb_rnaseq2016_PLQGenes_NegCoexp.rewiredDF$FP$Difference<0)],
                                    PosCoexp_Loss=msbb_rnaseq2016Genes_PosCoexp.rewiredDF$FP$Common_Genes[which(msbb_rnaseq2016Genes_PosCoexp.rewiredDF$FP$Difference>0)],
                                    PosCoexp_Gain=msbb_rnaseq2016Genes_PosCoexp.rewiredDF$FP$Common_Genes[which(msbb_rnaseq2016Genes_PosCoexp.rewiredDF$FP$Difference<0)],
                                    NegCoexp_Loss=msbb_rnaseq2016Genes_NegCoexp.rewiredDF$FP$Common_Genes[which(msbb_rnaseq2016Genes_NegCoexp.rewiredDF$FP$Difference>0)],
                                    NegCoexp_Gain=msbb_rnaseq2016Genes_NegCoexp.rewiredDF$FP$Common_Genes[which(msbb_rnaseq2016Genes_NegCoexp.rewiredDF$FP$Difference<0)],
                                    PosCoexp_Low_UniqGenes=msbb_rnaseq2016_FP_uniq_rewired$Low$PosCoexp_UniqGenes,
                                    PosCoexp_High_UniqGenes=msbb_rnaseq2016_FP_uniq_rewired$High$PosCoexp_UniqGenes,
                                    NegCoexp_Low_UniqGenes=msbb_rnaseq2016_FP_uniq_rewired$Low$NegCoexp_UniqGenes,
                                    NegCoexp_High_UniqGenes=msbb_rnaseq2016_FP_uniq_rewired$High$NegCoexp_UniqGenes)
msbb_rnaseq2016_AllAnalyses$IFG=list(PLQ_Genes=msbb_rnaseq2016_PLQGenes$IFG$Genes[which(msbb_rnaseq2016_PLQGenes$IFG$FDR<=0.01)],
                                     PLQ_DEG=msbb_rnaseq2016_PLQGenes.DEG$IFG$row[which(abs(msbb_rnaseq2016_PLQGenes.DEG$IFG$log2FoldChange)>1&msbb_rnaseq2016_PLQGenes.DEG$IFG$padj<=0.05)],
                                     PLQ_DCN=V(msbb_rnaseq2016_PLQ.DCN$IFG$High)$name,
                                     DCN=V(msbb_rnaseq2016.DCN$IFG$High)$name,
                                     PLQ_PosCoexp_Loss=msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF$IFG$Common_Genes[which(msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF$IFG$Difference>0)],
                                     PLQ_PosCoexp_Gain=msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF$IFG$Common_Genes[which(msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF$IFG$Difference<0)],
                                     PLQ_NegCoexp_Loss=msbb_rnaseq2016_PLQGenes_NegCoexp.rewiredDF$IFG$Common_Genes[which(msbb_rnaseq2016_PLQGenes_NegCoexp.rewiredDF$IFG$Difference>0)],
                                     PLQ_NegCoexp_Gain=msbb_rnaseq2016_PLQGenes_NegCoexp.rewiredDF$IFG$Common_Genes[which(msbb_rnaseq2016_PLQGenes_NegCoexp.rewiredDF$IFG$Difference<0)],
                                     PosCoexp_Loss=msbb_rnaseq2016Genes_PosCoexp.rewiredDF$IFG$Common_Genes[which(msbb_rnaseq2016Genes_PosCoexp.rewiredDF$IFG$Difference>0)],
                                     PosCoexp_Gain=msbb_rnaseq2016Genes_PosCoexp.rewiredDF$IFG$Common_Genes[which(msbb_rnaseq2016Genes_PosCoexp.rewiredDF$IFG$Difference<0)],
                                     NegCoexp_Loss=msbb_rnaseq2016Genes_NegCoexp.rewiredDF$IFG$Common_Genes[which(msbb_rnaseq2016Genes_NegCoexp.rewiredDF$IFG$Difference>0)],
                                     NegCoexp_Gain=msbb_rnaseq2016Genes_NegCoexp.rewiredDF$IFG$Common_Genes[which(msbb_rnaseq2016Genes_NegCoexp.rewiredDF$IFG$Difference<0)],
                                     PosCoexp_Low_UniqGenes=msbb_rnaseq2016_IFG_uniq_rewired$Low$PosCoexp_UniqGenes,
                                     PosCoexp_High_UniqGenes=msbb_rnaseq2016_IFG_uniq_rewired$High$PosCoexp_UniqGenes,
                                     NegCoexp_Low_UniqGenes=msbb_rnaseq2016_IFG_uniq_rewired$Low$NegCoexp_UniqGenes,
                                     NegCoexp_High_UniqGenes=msbb_rnaseq2016_IFG_uniq_rewired$High$NegCoexp_UniqGenes)
msbb_rnaseq2016_AllAnalyses$PHG=list(PLQ_Genes=msbb_rnaseq2016_PLQGenes$PHG$Genes[which(msbb_rnaseq2016_PLQGenes$PHG$FDR<=0.01)],
                                     PLQ_DEG=msbb_rnaseq2016_PLQGenes.DEG$PHG$row[which(abs(msbb_rnaseq2016_PLQGenes.DEG$PHG$log2FoldChange)>1&msbb_rnaseq2016_PLQGenes.DEG$PHG$padj<=0.05)],
                                     PLQ_DCN=V(msbb_rnaseq2016_PLQ.DCN$PHG$High)$name,
                                     DCN=V(msbb_rnaseq2016.DCN$PHG$High)$name,
                                     PLQ_PosCoexp_Loss=msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF$PHG$Common_Genes[which(msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF$PHG$Difference>0)],
                                     PLQ_PosCoexp_Gain=msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF$PHG$Common_Genes[which(msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF$PHG$Difference<0)],
                                     PLQ_NegCoexp_Loss=msbb_rnaseq2016_PLQGenes_NegCoexp.rewiredDF$PHG$Common_Genes[which(msbb_rnaseq2016_PLQGenes_NegCoexp.rewiredDF$PHG$Difference>0)],
                                     PLQ_NegCoexp_Gain=msbb_rnaseq2016_PLQGenes_NegCoexp.rewiredDF$PHG$Common_Genes[which(msbb_rnaseq2016_PLQGenes_NegCoexp.rewiredDF$PHG$Difference<0)],
                                     PosCoexp_Loss=msbb_rnaseq2016Genes_PosCoexp.rewiredDF$PHG$Common_Genes[which(msbb_rnaseq2016Genes_PosCoexp.rewiredDF$PHG$Difference>0)],
                                     PosCoexp_Gain=msbb_rnaseq2016Genes_PosCoexp.rewiredDF$PHG$Common_Genes[which(msbb_rnaseq2016Genes_PosCoexp.rewiredDF$PHG$Difference<0)],
                                     NegCoexp_Loss=msbb_rnaseq2016Genes_NegCoexp.rewiredDF$PHG$Common_Genes[which(msbb_rnaseq2016Genes_NegCoexp.rewiredDF$PHG$Difference>0)],
                                     NegCoexp_Gain=msbb_rnaseq2016Genes_NegCoexp.rewiredDF$PHG$Common_Genes[which(msbb_rnaseq2016Genes_NegCoexp.rewiredDF$PHG$Difference<0)],
                                     PosCoexp_Low_UniqGenes=msbb_rnaseq2016_PHG_uniq_rewired$Low$PosCoexp_UniqGenes,
                                     PosCoexp_High_UniqGenes=msbb_rnaseq2016_PHG_uniq_rewired$High$PosCoexp_UniqGenes,
                                     NegCoexp_Low_UniqGenes=msbb_rnaseq2016_PHG_uniq_rewired$Low$NegCoexp_UniqGenes,
                                     NegCoexp_High_UniqGenes=msbb_rnaseq2016_PHG_uniq_rewired$High$NegCoexp_UniqGenes)
msbb_rnaseq2016_AllAnalyses$STG=list(PLQ_Genes=msbb_rnaseq2016_PLQGenes$STG$Genes[which(msbb_rnaseq2016_PLQGenes$STG$FDR<=0.01)],
                                     PLQ_DEG=msbb_rnaseq2016_PLQGenes.DEG$STG$row[which(abs(msbb_rnaseq2016_PLQGenes.DEG$STG$log2FoldChange)>1&msbb_rnaseq2016_PLQGenes.DEG$STG$padj<=0.05)],
                                     PLQ_DCN=V(msbb_rnaseq2016_PLQ.DCN$STG$High)$name,
                                     DCN=V(msbb_rnaseq2016.DCN$STG$High)$name,
                                     PLQ_PosCoexp_Loss=msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF$STG$Common_Genes[which(msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF$STG$Difference>0)],
                                     PLQ_PosCoexp_Gain=msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF$STG$Common_Genes[which(msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF$STG$Difference<0)],
                                     PLQ_NegCoexp_Loss=msbb_rnaseq2016_PLQGenes_NegCoexp.rewiredDF$STG$Common_Genes[which(msbb_rnaseq2016_PLQGenes_NegCoexp.rewiredDF$STG$Difference>0)],
                                     PLQ_NegCoexp_Gain=msbb_rnaseq2016_PLQGenes_NegCoexp.rewiredDF$STG$Common_Genes[which(msbb_rnaseq2016_PLQGenes_NegCoexp.rewiredDF$STG$Difference<0)],
                                     PosCoexp_Loss=msbb_rnaseq2016Genes_PosCoexp.rewiredDF$STG$Common_Genes[which(msbb_rnaseq2016Genes_PosCoexp.rewiredDF$STG$Difference>0)],
                                     PosCoexp_Gain=msbb_rnaseq2016Genes_PosCoexp.rewiredDF$STG$Common_Genes[which(msbb_rnaseq2016Genes_PosCoexp.rewiredDF$STG$Difference<0)],
                                     NegCoexp_Loss=msbb_rnaseq2016Genes_NegCoexp.rewiredDF$STG$Common_Genes[which(msbb_rnaseq2016Genes_NegCoexp.rewiredDF$STG$Difference>0)],
                                     NegCoexp_Gain=msbb_rnaseq2016Genes_NegCoexp.rewiredDF$STG$Common_Genes[which(msbb_rnaseq2016Genes_NegCoexp.rewiredDF$STG$Difference<0)],
                                     PosCoexp_Low_UniqGenes=msbb_rnaseq2016_STG_uniq_rewired$Low$PosCoexp_UniqGenes,
                                     PosCoexp_High_UniqGenes=msbb_rnaseq2016_STG_uniq_rewired$High$PosCoexp_UniqGenes,
                                     NegCoexp_Low_UniqGenes=msbb_rnaseq2016_STG_uniq_rewired$Low$NegCoexp_UniqGenes,
                                     NegCoexp_High_UniqGenes=msbb_rnaseq2016_STG_uniq_rewired$High$NegCoexp_UniqGenes)

#PLQ rewiring only
PLQ_Diff_Rewiring=list(FP_PLQ_PosCoexp_Gain=msbb_rnaseq2016_AllAnalyses$FP$PLQ_PosCoexp_Gain,
                       FP_PLQ_PosCoexp_Loss=msbb_rnaseq2016_AllAnalyses$FP$PLQ_PosCoexp_Loss,
                       FP_PLQ_NegCoexp_Gain=msbb_rnaseq2016_AllAnalyses$FP$PLQ_NegCoexp_Gain,
                       FP_PLQ_NegCoexp_Loss=msbb_rnaseq2016_AllAnalyses$FP$PLQ_NegCoexp_Loss,
                       IFG_PLQ_PosCoexp_Gain=msbb_rnaseq2016_AllAnalyses$IFG$PLQ_PosCoexp_Gain,
                       IFG_PLQ_PosCoexp_Loss=msbb_rnaseq2016_AllAnalyses$IFG$PLQ_PosCoexp_Loss,
                       IFG_PLQ_NegCoexp_Gain=msbb_rnaseq2016_AllAnalyses$IFG$PLQ_NegCoexp_Gain,
                       IFG_PLQ_NegCoexp_Loss=msbb_rnaseq2016_AllAnalyses$IFG$PLQ_NegCoexp_Loss,
                       PHG_PLQ_PosCoexp_Gain=msbb_rnaseq2016_AllAnalyses$PHG$PLQ_PosCoexp_Gain,
                       PHG_PLQ_PosCoexp_Loss=msbb_rnaseq2016_AllAnalyses$PHG$PLQ_PosCoexp_Loss,
                       PHG_PLQ_NegCoexp_Gain=msbb_rnaseq2016_AllAnalyses$PHG$PLQ_NegCoexp_Gain,
                       PHG_PLQ_NegCoexp_Loss=msbb_rnaseq2016_AllAnalyses$PHG$PLQ_NegCoexp_Loss,
                       STG_PLQ_PosCoexp_Gain=msbb_rnaseq2016_AllAnalyses$STG$PLQ_PosCoexp_Gain,
                       STG_PLQ_PosCoexp_Loss=msbb_rnaseq2016_AllAnalyses$STG$PLQ_PosCoexp_Loss,
                       STG_PLQ_NegCoexp_Gain=msbb_rnaseq2016_AllAnalyses$STG$PLQ_NegCoexp_Gain,
                       STG_PLQ_NegCoexp_Loss=msbb_rnaseq2016_AllAnalyses$STG$PLQ_NegCoexp_Loss)

#Cell-type marker analysis
zhang_celltype_ADgenes=read.xls('../../../BrainExpression_Datasets/Zhang_19BrainRegions_Paper/13073_2016_355_MOESM1_ESM.xlsx',skip=1,,sheet=3,header=T,as.is=T)
zhang_celltype_ADgenes.list=zhang_celltype_PLQ_ADgenes.list=vector(mode = 'list',length = 5)
names(zhang_celltype_ADgenes.list)=names(zhang_celltype_PLQ_ADgenes.list)=sort(unique(zhang_celltype_ADgenes$Cell.type))
zhang_celltype_ADgenes.list$Astrocytes=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_ADgenes.list)[1],zhang_celltype_ADgenes$Cell.type)]
zhang_celltype_ADgenes.list$Endothelial=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_ADgenes.list)[2],zhang_celltype_ADgenes$Cell.type)]
zhang_celltype_ADgenes.list$Microglia=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_ADgenes.list)[3],zhang_celltype_ADgenes$Cell.type)]
zhang_celltype_ADgenes.list$Neurons=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_ADgenes.list)[4],zhang_celltype_ADgenes$Cell.type)]
zhang_celltype_ADgenes.list$Oligodendrocytes=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_ADgenes.list)[5],zhang_celltype_ADgenes$Cell.type)]

zhang_celltype_PLQ_ADgenes.list$Astrocytes=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_PLQ_ADgenes.list)[1],zhang_celltype_ADgenes$Cell.type)]
zhang_celltype_PLQ_ADgenes.list$Endothelial=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_PLQ_ADgenes.list)[2],zhang_celltype_ADgenes$Cell.type)]
zhang_celltype_PLQ_ADgenes.list$Microglia=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_PLQ_ADgenes.list)[3],zhang_celltype_ADgenes$Cell.type)]
zhang_celltype_PLQ_ADgenes.list$Neurons=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_PLQ_ADgenes.list)[4],zhang_celltype_ADgenes$Cell.type)]
zhang_celltype_PLQ_ADgenes.list$Oligodendrocytes=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_PLQ_ADgenes.list)[5],zhang_celltype_ADgenes$Cell.type)]

msbb_dcn.Zhange_celltype_markers=cbind(FP_celltype=lapply(zhang_celltype_ADgenes.list,intersect,msbb_dcn.clusters$FP),
                                       IFG_celltype=lapply(zhang_celltype_ADgenes.list,intersect,msbb_dcn.clusters$IFG),
                                       PHG_celltype=lapply(zhang_celltype_ADgenes.list,intersect,msbb_dcn.clusters$PHG),
                                       STG_celltype=lapply(zhang_celltype_ADgenes.list,intersect,msbb_dcn.clusters$STG))
msbb_dcn.Zhange_celltype_markers2=data.frame(msbb_dcn.Zhange_celltype_markers)

msbb_dcn.Zhange_PLQ_celltype_markers=cbind(FP_PLQ_celltype=lapply(zhang_celltype_PLQ_ADgenes.list,intersect,msbb_plq_dcn.clusters$FP),
                                           IFG_PLQ_celltype=lapply(zhang_celltype_PLQ_ADgenes.list,intersect,msbb_plq_dcn.clusters$IFG),
                                           PHG_PLQ_celltype=lapply(zhang_celltype_PLQ_ADgenes.list,intersect,msbb_plq_dcn.clusters$PHG),
                                           STG_PLQ_celltype=lapply(zhang_celltype_PLQ_ADgenes.list,intersect,msbb_plq_dcn.clusters$STG))
msbb_dcn.Zhange_PLQ_celltype_markers2=data.frame(msbb_dcn.Zhange_PLQ_celltype_markers)

msbb_dcn.Zhange_celltype_markers.combined=data.frame(cbind(FP_celltype=lapply(zhang_celltype_ADgenes.list,intersect,msbb_dcn.clusters$FP),
                                                           IFG_celltype=lapply(zhang_celltype_ADgenes.list,intersect,msbb_dcn.clusters$IFG),
                                                           PHG_celltype=lapply(zhang_celltype_ADgenes.list,intersect,msbb_dcn.clusters$PHG),
                                                           STG_celltype=lapply(zhang_celltype_ADgenes.list,intersect,msbb_dcn.clusters$STG),
                                                           FP_PLQ_celltype=lapply(zhang_celltype_PLQ_ADgenes.list,intersect,msbb_plq_dcn.clusters$FP),
                                                           IFG_PLQ_celltype=lapply(zhang_celltype_PLQ_ADgenes.list,intersect,msbb_plq_dcn.clusters$IFG),
                                                           PHG_PLQ_celltype=lapply(zhang_celltype_PLQ_ADgenes.list,intersect,msbb_plq_dcn.clusters$PHG),
                                                           STG_PLQ_celltype=lapply(zhang_celltype_PLQ_ADgenes.list,intersect,msbb_plq_dcn.clusters$STG)))

plot(compareCluster(geneClusters = lapply(lapply(lapply(msbb_dcn.Zhange_celltype_markers2,`[[`,3),select,x=org.Hs.eg.db,columns = 'ENTREZID',keytype = 'SYMBOL'),`[[`,2),fun = 'enrichKEGG',organism='hsa',pvalueCutoff=0.05))
plot(compareCluster(geneClusters = lapply(lapply(lapply(msbb_dcn.Zhange_celltype_markers2,`[[`,3),select,x=org.Hs.eg.db,columns = 'ENTREZID',keytype = 'SYMBOL'),`[[`,2),fun = 'enrichPathway',organism='human',pvalueCutoff=0.05))

#Rewiring threshold optimisation
pos_weight=seq(1.25,1.75,by=0.01)
neg_weight=seq(0.25,0.75,by=0.01)
poscoexp_pval_vector=negcoexp_pval_vector=poscoexp_degree=negcoexp_degree=vector(mode = "list",length = 4)
names(poscoexp_pval_vector)=names(negcoexp_pval_vector)=names(poscoexp_degree)=names(negcoexp_degree)=c("FP","IFG","PHG","STG")
poscoexp_pval_vector$FP=poscoexp_pval_vector$IFG=poscoexp_pval_vector$PHG=poscoexp_pval_vector$STG=vector(mode = "list",length = length(pos_weight))
negcoexp_pval_vector$FP=negcoexp_pval_vector$IFG=negcoexp_pval_vector$PHG=negcoexp_pval_vector$STG=vector(mode = "list",length = length(neg_weight))
names(poscoexp_pval_vector$FP)=names(poscoexp_pval_vector$IFG)=names(poscoexp_pval_vector$PHG)=names(poscoexp_pval_vector$STG)=pos_weight
names(negcoexp_pval_vector$FP)=names(negcoexp_pval_vector$IFG)=names(negcoexp_pval_vector$PHG)=names(negcoexp_pval_vector$STG)=neg_weight

poscoexp_degree$FP=poscoexp_degree$IFG=poscoexp_degree$PHG=poscoexp_degree$STG=vector(mode = "list",length = 2)
negcoexp_degree$FP=negcoexp_degree$IFG=negcoexp_degree$PHG=negcoexp_degree$STG=vector(mode = "list",length = 2)
names(poscoexp_degree$FP)=names(poscoexp_degree$IFG)=names(poscoexp_degree$PHG)=names(poscoexp_degree$STG)=c("Low","High")
names(negcoexp_degree$FP)=names(negcoexp_degree$IFG)=names(negcoexp_degree$PHG)=names(negcoexp_degree$STG)=c("Low","High")

for(r in 1:4){
  for (w in 1:length(pos_weight))
  {
    poscoexp_degree[[r]]$Low=degree(graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016.DCN[[r]]$Low,es = which(E(msbb_rnaseq2016.DCN[[r]]$Low)$weight<pos_weight[w]),names = T)),directed = F))
    poscoexp_degree[[r]]$High=degree(graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016.DCN[[r]]$High,es = which(E(msbb_rnaseq2016.DCN[[r]]$High)$weight>pos_weight[w]),names = T)),directed = F))
    negcoexp_degree[[r]]$Low=degree(graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016.DCN[[r]]$Low,es = which(E(msbb_rnaseq2016.DCN[[r]]$Low)$weight<neg_weight[w]),names = T)),directed = F))
    negcoexp_degree[[r]]$High=degree(graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016.DCN[[r]]$High,es = which(E(msbb_rnaseq2016.DCN[[r]]$High)$weight>neg_weight[w]),names = T)),directed = F))
    poscoexp_pval_vector[[r]][[w]]=-log10(wilcox.test(x=poscoexp_degree[[r]]$Low,y=poscoexp_degree[[r]]$High)$p.value)
    negcoexp_pval_vector[[r]][[w]]=-log10(wilcox.test(x=negcoexp_degree[[r]]$Low,y=negcoexp_degree[[r]]$High)$p.value)
  }
  
}
df1=cbind(FP_PosCoexp_pval=poscoexp_pval_vector$FP,
          IFG_PosCoexp_pval=poscoexp_pval_vector$IFG,
          PHG_PosCoexp_pval=poscoexp_pval_vector$PHG,
          STG_PosCoexp_pval=poscoexp_pval_vector$STG,
          FP_NegCoexp_pval=negcoexp_pval_vector$FP,
          IFG_NegCoexp_pval=negcoexp_pval_vector$IFG,
          PHG_NegCoexp_pval=negcoexp_pval_vector$PHG,
          STG_NegCoexp_pval=negcoexp_pval_vector$STG)
threshold_df=data.frame(df1)
rownames(threshold_df)=as.numeric(rownames(df1))-1

plot(compareCluster(geneClusters = lapply(lapply(msbb_rnaseq2016_AllAnalyses$IFG[c(5:7)],select,x = org.Hs.eg.db,keytype = 'SYMBOL',columns = 'ENTREZID'),`[[`,2),fun = 'enrichKEGG',organism='hsa',pvalueCutoff=0.1,pAdjustMethod='BH'),title='FP Rewiring,KEGG pathway enrichment')



graph.union(induced_subgraph(graph = msbb_rnaseq2016.PosCoexp$STG$Low,vids = msbb_rnaseq2016Genes_PosCoexp.rewiredDF$STG$Common_Genes[(msbb_rnaseq2016Genes_PosCoexp.rewiredDF$STG$Difference>0)]),induced_subgraph(graph = ttrust_tf_network,vids = which(V(ttrust_tf_network)$name%in%msbb_rnaseq2016Genes_PosCoexp.rewiredDF$STG$Common_Genes[(msbb_rnaseq2016Genes_PosCoexp.rewiredDF$STG$Difference>0)]%in%ttrust_tf_data$TF)))