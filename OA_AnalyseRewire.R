library(igraph)
library(org.Hs.eg.db)

OA_analysis=readRDS("~/Work/Collaborations/OA_EleZeggini/OA_D_vs_C_cpm_corr_change_allResults.RDS")
OA_DCN=vector(mode = "list",length = 2)
names(OA_DCN)=c("Healthy","Diseased")
OA_DCN$Healthy=graph.data.frame(d = OA_analysis[which((OA_analysis$FDR<=0.1)&(OA_analysis$p.cocor<=0.05)),c(2:3)],directed = F)
OA_DCN$Diseased=graph.data.frame(d = OA_analysis[which((OA_analysis$FDR<=0.1)&(OA_analysis$p.cocor<=0.05)),c(2:3)],directed = F)

E(OA_DCN$Healthy)$weight=OA_analysis[which((OA_analysis$FDR<=0.1)&(OA_analysis$p.cocor<=0.05)),]$r.c+1
E(OA_DCN$Diseased)$weight=OA_analysis[which((OA_analysis$FDR<=0.1)&(OA_analysis$p.cocor<=0.05)),]$r.t+1

g1=graph.data.frame(d = ends(graph = OA_DCN$Healthy,es = which(E(OA_DCN$Healthy)$weight>=1.5),names = T),directed = F)
g2=graph.data.frame(d = ends(graph = OA_DCN$Diseased,es = which(E(OA_DCN$Diseased)$weight>=1.5),names = T),directed = F)
g3=graph.data.frame(d = ends(graph = OA_DCN$Healthy,es = which(E(OA_DCN$Healthy)$weight<=0.5),names = T),directed = F)
g4=graph.data.frame(d = ends(graph = OA_DCN$Diseased,es = which(E(OA_DCN$Diseased)$weight<=0.5),names = T),directed = F)

OA_DCN.rewired=vector(mode = "list",length = 2)
names(OA_DCN.rewired)=c("PosCorr","NegCorr")
OA_DCN.rewired$PosCorr=data.frame(Genes=names(degree(graph = g1,v = intersect(V(g1)$name,V(g2)$name))),OA_Healthy=degree(graph = g1,v = intersect(V(g1)$name,V(g2)$name)),OA_Diseased=degree(graph = g2,v = intersect(V(g1)$name,V(g2)$name)),Diff=degree(graph = g1,v = intersect(V(g1)$name,V(g2)$name))-degree(graph = g2,v = intersect(V(g1)$name,V(g2)$name)),stringsAsFactors = F)
OA_DCN.rewired$NegCorr=data.frame(Genes=names(degree(graph = g3,v = intersect(V(g3)$name,V(g4)$name))),OA_Healthy=degree(graph = g3,v = intersect(V(g3)$name,V(g4)$name)),OA_Diseased=degree(graph = g4,v = intersect(V(g3)$name,V(g4)$name)),Diff=degree(graph = g3,v = intersect(V(g3)$name,V(g4)$name))-degree(graph = g4,v = intersect(V(g3)$name,V(g4)$name)),stringsAsFactors = F)

OA_DCN.rewired_cluster=vector(mode = "list",length = 4)
names(OA_DCN.rewired_cluster)=c("PosCoexp_Gain","PosCoexp_Loss","NegCoexp_Gain","NegCoexp_Loss")
OA_DCN.rewired_cluster$PosCoexp_Gain=select(x = org.Hs.eg.db,keys=OA_DCN.rewired$PosCorr$Genes[which(OA_DCN.rewired$PosCorr$Diff>0)],columns = "ENTREZID",keytype = "SYMBOL")[,2]
OA_DCN.rewired_cluster$PosCoexp_Loss=select(x = org.Hs.eg.db,keys=OA_DCN.rewired$PosCorr$Genes[which(OA_DCN.rewired$PosCorr$Diff<0)],columns = "ENTREZID",keytype = "SYMBOL")[,2]
OA_DCN.rewired_cluster$NegCoexp_Gain=select(x = org.Hs.eg.db,keys=OA_DCN.rewired$NegCorr$Genes[which(OA_DCN.rewired$PosCorr$Diff>0)],columns = "ENTREZID",keytype = "SYMBOL")[,2]
OA_DCN.rewired_cluster$NegCoexp_Loss=select(x = org.Hs.eg.db,keys=OA_DCN.rewired$NegCorr$Genes[which(OA_DCN.rewired$PosCorr$Diff<0)],columns = "ENTREZID",keytype = "SYMBOL")[,2]

OA_DCN.rewired_cluster2=vector(mode = "list",length = 4)
names(OA_DCN.rewired_cluster2)=c("PosCoexp_Gain","PosCoexp_Loss","NegCoexp_Gain","NegCoexp_Loss")
OA_DCN.rewired_cluster2$PosCoexp_Gain=OA_DCN.rewired$PosCorr$Genes[which(OA_DCN.rewired$PosCorr$Diff>0)]
OA_DCN.rewired_cluster2$PosCoexp_Loss=OA_DCN.rewired$PosCorr$Genes[which(OA_DCN.rewired$PosCorr$Diff<0)]
OA_DCN.rewired_cluster2$NegCoexp_Gain=OA_DCN.rewired$NegCorr$Genes[which(OA_DCN.rewired$PosCorr$Diff>0)]
OA_DCN.rewired_cluster2$NegCoexp_Loss=OA_DCN.rewired$NegCorr$Genes[which(OA_DCN.rewired$PosCorr$Diff<0)]

write.table(OA_DCN.rewired$PosCorr,"OA_DCN.rewired_PosCorr_Genes.txt",sep = "\t",col.names = T,row.names = F,quote=F)
write.table(OA_DCN.rewired$NegCorr,"OA_DCN.rewired_NegCorr_Genes.txt",sep = "\t",col.names = T,row.names = F,quote=F)
ge