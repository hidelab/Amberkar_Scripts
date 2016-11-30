library(org.Hs.eg.db)
library(igraph)
library(clusterProfiler)
library(cocor)

setwd("/Users/sandeepamberkar/Work/Data/ROSMAP/")
rosmap_cocor=readRDS("ROSMAPsigEdge.RDS")
rosmap_DCN$NCI=rosmap_DCN$AD=graph.data.frame(d = rosmap_cocor[,c(13,15)],directed = F)
E(rosmap_DCN$NCI)$weight=rosmap_cocor$r.c+1
E(rosmap_DCN$AD)$weight=rosmap_cocor$r.t+1

rosmap_DCN2=vector(mode = "list",length = 2)
rosmap_DCN2$NCI=rosmap_DCN2$AD=vector(mode = "list",length = 2)
names(rosmap_DCN2$NCI)=names(rosmap_DCN2$AD)=c("PosCoexp","NegCoexp")

rosmap_DCN2$NCI$PosCoexp=graph.data.frame(d = data.frame(ends(graph = rosmap_DCN$NCI,es = which(E(rosmap_DCN$NCI)$weight>1.25),names = T)),directed = F)
E(rosmap_DCN2$NCI$PosCoexp)$weight=E(rosmap_DCN$NCI)$weight[which(E(rosmap_DCN$NCI)$weight>1.25)]
rosmap_DCN2$AD$PosCoexp=graph.data.frame(d = data.frame(ends(graph = rosmap_DCN$AD,es = which(E(rosmap_DCN$AD)$weight>1.25),names = T)),directed = F)
E(rosmap_DCN2$AD$PosCoexp)$weight=E(rosmap_DCN$AD)$weight[which(E(rosmap_DCN$AD)$weight>1.25)]
rosmap_DCN2$NCI$NegCoexp=graph.data.frame(d = data.frame(ends(graph = rosmap_DCN$NCI,es = which(E(rosmap_DCN$NCI)$weight<0.75),names = T)),directed = F)
E(rosmap_DCN2$NCI$NegCoexp)$weight=E(rosmap_DCN$NCI)$weight[which(E(rosmap_DCN$NCI)$weight<0.75)]
rosmap_DCN2$AD$NegCoexp=graph.data.frame(d = data.frame(ends(graph = rosmap_DCN$AD,es = which(E(rosmap_DCN$AD)$weight<0.75),names = T)),directed = F)
E(rosmap_DCN2$AD$NegCoexp)$weight=E(rosmap_DCN$AD)$weight[which(E(rosmap_DCN$AD)$weight<0.75)]

rosmap_DCN.rewire=vector(mode = "list",length = 2)
names(rosmap_DCN.rewire)=c("PosCoexp","NegCoexp")
rosmap_DCN.rewire$PosCoexp=data.frame(Genes=intersect(V(rosmap_DCN2$AD$PosCoexp)$name,V(rosmap_DCN2$NCI$PosCoexp)$name),
                                      NCI=degree(graph = rosmap_DCN2$NCI$PosCoexp,v = intersect(V(rosmap_DCN2$AD$PosCoexp)$name,V(rosmap_DCN2$NCI$PosCoexp)$name)),
                                      AD=degree(graph = rosmap_DCN2$AD$PosCoexp,v = intersect(V(rosmap_DCN2$AD$PosCoexp)$name,V(rosmap_DCN2$NCI$PosCoexp)$name)),
                                      Difference=degree(graph = rosmap_DCN2$NCI$PosCoexp,v = intersect(V(rosmap_DCN2$AD$PosCoexp)$name,V(rosmap_DCN2$NCI$PosCoexp)$name))-degree(graph = rosmap_DCN2$AD$PosCoexp,v = intersect(V(rosmap_DCN2$AD$PosCoexp)$name,V(rosmap_DCN2$NCI$PosCoexp)$name)),stringsAsFactors = F)
rosmap_DCN.rewire$NegCoexp=data.frame(Genes=intersect(V(rosmap_DCN2$AD$NegCoexp)$name,V(rosmap_DCN2$NCI$NegCoexp)$name),
                                      NCI=degree(graph = rosmap_DCN2$NCI$NegCoexp,v = intersect(V(rosmap_DCN2$AD$NegCoexp)$name,V(rosmap_DCN2$NCI$NegCoexp)$name)),
                                      AD=degree(graph = rosmap_DCN2$AD$NegCoexp,v = intersect(V(rosmap_DCN2$AD$Neg)$name,V(rosmap_DCN2$NCI$NegCoexp)$name)),
                                      Difference=degree(graph = rosmap_DCN2$NCI$NegCoexp,v = intersect(V(rosmap_DCN2$AD$NegCoexp)$name,V(rosmap_DCN2$NCI$NegCoexp)$name))-degree(graph = rosmap_DCN2$AD$Neg,v = intersect(V(rosmap_DCN2$AD$Neg)$name,V(rosmap_DCN2$NCI$Neg)$name)),stringsAsFactors = F)
write.table(rosmap_DCN.rewire$PosCoexp,"rosmap_DCN_rewire_PosCoexp.txt",sep = "\t",col.names = T,row.names = F,quote = F)
write.table(rosmap_DCN.rewire$NegCoexp,"rosmap_DCN_rewire_NegCoexp.txt",sep = "\t",col.names = T,row.names = F,quote = F)

rosmap_DCN.rewire_cluster=vector(mode="list",length = 4)
names(rosmap_DCN.rewire_cluster)=c("PosCoexp_Gain","PosCoexp_Loss","NegCoexp_Gain","NegCoexp_Loss")
rosmap_DCN.rewire_cluster$PosCoexp_Gain=select(x = org.Hs.eg.db,keys = rosmap_DCN.rewire$PosCoexp$Genes[which(rosmap_DCN.rewire$PosCoexp$Difference<0)],columns = "ENTREZID",keytype = "SYMBOL")[,2]
rosmap_DCN.rewire_cluster$PosCoexp_Loss=select(x = org.Hs.eg.db,keys = rosmap_DCN.rewire$PosCoexp$Genes[which(rosmap_DCN.rewire$PosCoexp$Difference>0)],columns = "ENTREZID",keytype = "SYMBOL")[,2]
rosmap_DCN.rewire_cluster$NegCoexp_Gain=select(x = org.Hs.eg.db,keys = rosmap_DCN.rewire$NegCoexp$Genes[which(rosmap_DCN.rewire$NegCoexp$Difference<0)],columns = "ENTREZID",keytype = "SYMBOL")[,2]
rosmap_DCN.rewire_cluster$NegCoexp_Loss=select(x = org.Hs.eg.db,keys = rosmap_DCN.rewire$NegCoexp$Genes[which(rosmap_DCN.rewire$NegCoexp$Difference>0)],columns = "ENTREZID",keytype = "SYMBOL")[,2]

plot(compareCluster(geneClusters = rosmap_DCN.rewire_cluster,fun = "enrichKEGG",pvalueCutoff = 0.1,pAdjustMethod = "BH"))
 
write(rosmap_DCN.rewire$NegCoexp$Genes[which(rosmap_DCN.rewire$NegCoexp$Difference>0)],"rosmap_DCN_rewire_NegCoexp_Loss.txt",sep = "\n")
write(rosmap_DCN.rewire$NegCoexp$Genes[which(rosmap_DCN.rewire$NegCoexp$Difference<0)],"rosmap_DCN_rewire_NegCoexp_Gain.txt",sep = "\n")
write(rosmap_DCN.rewire$PosCoexp$Genes[which(rosmap_DCN.rewire$PosCoexp$Difference>0)],"rosmap_DCN_rewire_PosCoexp_Loss.txt",sep = "\n")
write(rosmap_DCN.rewire$PosCoexp$Genes[which(rosmap_DCN.rewire$PosCoexp$Difference<0)],"rosmap_DCN_rewire_PosCoexp_Gain.txt",sep = "\n")
write(select(x = org.Hs.eg.db,keys = c("2900","22941","2917","2785","2903","5567","9229","2890","57030"),keytype = "ENTREZID",columns = "SYMBOL")[,2],"rosmap_DCN_rewire_PosCoexp_Loss_Glutamatergic_synapse_genes.txt",sep = "\n")
write(select(x = org.Hs.eg.db,keys = c("2900","22941","2917","9229","775","23236","2893","2903"),keytype = "ENTREZID",columns = "SYMBOL")[,2],"rosmap_DCN_rewire_NegCoexp_Gain_Glutamatergic_synapse_genes.txt",sep = "\n")