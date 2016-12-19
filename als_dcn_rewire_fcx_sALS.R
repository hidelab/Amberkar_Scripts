library(data.table)
library(org.Hs.eg.db)
library(igraph)
library(parallel)
library(gtools)

setwd("/shared/hidelab2/user/md4zsa/Work/Data/ALS/sALS_DiffCorr/txt")
system("awk 'FNR==1 && NR!=1 { while (/^<header>/) getline; }    1 {print}' FCX*.txt >FCX_sALS_DiffCorr_Results.txt")
allResults$FDR=p.adjust(p = allResults$p.cocor,method = "fdr")
allResults$FDR.c=p.adjust(p = allResults$p.c,method = "fdr")
allResults$FDR.t=p.adjust(p = allResults$p.t,method = "fdr")

write.table(allResults,"FCX_ALS_allResults_DiffCorr.txt",sep="\t",col.names = T,row.names = F,quote=F)
fcx_als_allResults=read.table("FCX_ALS_allResults_DiffCorr.txt",sep="\t",header = T,as.is = T)

fcx_als_DCN=vector(mode = "list",length = 2)
names(fcx_als_DCN)=c("Control","sALS")
fcx_als_DCN$Control=fcx_als_DCN$sALS=graph.data.frame(d=fcx_als_allResults[which(fcx_als_allResults$p.c<=0.05&fcx_als_allResults$p.t<=0.05&fcx_als_allResults$FDR<=0.1),c(1:2)],directed = F)
E(fcx_als_DCN$Control)$weight=fcx_als_allResults$r.c+1
E(fcx_als_DCN$sALS)$weight=fcx_als_allResults$r.t+1

fcx_als_PosCoexp=fcx_als_NegCoexp=vector(mode = "list",length = 2)
names(fcx_als_PosCoexp)=names(fcx_als_NegCoexp)=names(fcx_als_DCN)
fcx_als_PosCoexp$Control=graph.data.frame(d=data.frame(ends(graph = fcx_als_DCN$Control,es = which(E(fcx_als_DCN$Control)$weight>1.5),names = T)),directed = F)
# E(fcx_als_PosCoexp$Control)$weight=E(fcx_als_DCN$Control)$weight[which(E(fcx_als_DCN$Control)$weight>1.5)]
fcx_als_PosCoexp$sALS=graph.data.frame(d=data.frame(ends(graph = fcx_als_DCN$sALS,es = which(E(fcx_als_DCN$sALS)$weight>1.5),names = T)),directed = F)
# E(fcx_als_PosCoexp$sALS)$weight=E(fcx_als_DCN$sALS)$weight[which(E(fcx_als_DCN$sALS)$weight>1.5)]
fcx_als_NegCoexp$Control=graph.data.frame(d=data.frame(ends(graph = fcx_als_DCN$Control,es = which(E(fcx_als_DCN$Control)$weight<0.5),names = T)),directed = F)
fcx_als_NegCoexp$sALS=graph.data.frame(d=data.frame(ends(graph = fcx_als_DCN$sALS,es = which(E(fcx_als_DCN$sALS)$weight<0.5),names = T)),directed = F)

fcx_als_PosCoexp_RewiredGenes=data.frame(Common_Genes=intersect(V(fcx_als_PosCoexp$Control)$name,V(fcx_als_PosCoexp$sALS)$name),
                                         Control=degree(graph = fcx_als_PosCoexp$Control,v = intersect(V(fcx_als_PosCoexp$Control)$name,V(fcx_als_PosCoexp$sALS)$name)),
                                         sALS=degree(graph = fcx_als_PosCoexp$sALS,v = intersect(V(fcx_als_PosCoexp$Control)$name,V(fcx_als_PosCoexp$sALS)$name)),
                                         Difference=degree(graph = fcx_als_PosCoexp$Control,v = intersect(V(fcx_als_PosCoexp$Control)$name,V(fcx_als_PosCoexp$sALS)$name))-degree(graph = fcx_als_PosCoexp$sALS,v = intersect(V(fcx_als_PosCoexp$Control)$name,V(fcx_als_PosCoexp$sALS)$name)),stringsAsFactors = F)
fcx_als_NegCoexp_RewiredGenes=data.frame(Common_Genes=intersect(V(fcx_als_NegCoexp$Control)$name,V(fcx_als_NegCoexp$sALS)$name),
                                         Control=degree(graph = fcx_als_NegCoexp$Control,v = intersect(V(fcx_als_NegCoexp$Control)$name,V(fcx_als_NegCoexp$sALS)$name)),
                                         sALS=degree(graph = fcx_als_NegCoexp$sALS,v = intersect(V(fcx_als_NegCoexp$Control)$name,V(fcx_als_NegCoexp$sALS)$name)),
                                         Difference=degree(graph = fcx_als_NegCoexp$Control,v = intersect(V(fcx_als_NegCoexp$Control)$name,V(fcx_als_NegCoexp$sALS)$name))-degree(graph = fcx_als_NegCoexp$sALS,v = intersect(V(fcx_als_NegCoexp$Control)$name,V(fcx_als_NegCoexp$sALS)$name)),stringsAsFactors = F)

write.table(fcx_als_PosCoexp_RewiredGenes,"FCX_sALS_PosCoexp_RewiredGenes.txt",sep="\t",col.names = T,row.names = F,quote = F)
write.table(fcx_als_NegCoexp_RewiredGenes,"FCX_sALS_NegCoexp_RewiredGenes.txt",sep="\t",col.names = T,row.names = F,quote = F)

fcx_als_DCN_RewiredGenes.cluster=vector(mode = "list",length = 4)
names(fcx_als_DCN_RewiredGenes.cluster)=c("PosCoexp_Loss","PosCoexp_Gain","NegCoexp_Loss","NegCoexp_Gain")
fcx_als_DCN_RewiredGenes.cluster$PosCoexp_Loss=select(x = org.Hs.eg.db,keys = fcx_als_PosCoexp_RewiredGenes$Common_Genes[which(fcx_als_PosCoexp_RewiredGenes$Difference>0)],columns = "ENTREZID",keytype = "SYMBOL")[,2]
fcx_als_DCN_RewiredGenes.cluster$PosCoexp_Gain=select(x = org.Hs.eg.db,keys = fcx_als_PosCoexp_RewiredGenes$Common_Genes[which(fcx_als_PosCoexp_RewiredGenes$Difference<0)],columns = "ENTREZID",keytype = "SYMBOL")[,2]
fcx_als_DCN_RewiredGenes.cluster$NegCoexp_Loss=select(x = org.Hs.eg.db,keys = fcx_als_NegCoexp_RewiredGenes$Common_Genes[which(fcx_als_NegCoexp_RewiredGenes$Difference>0)],columns = "ENTREZID",keytype = "SYMBOL")[,2]
fcx_als_DCN_RewiredGenes.cluster$NegCoexp_Gain=select(x = org.Hs.eg.db,keys = fcx_als_NegCoexp_RewiredGenes$Common_Genes[which(fcx_als_NegCoexp_RewiredGenes$Difference<0)],columns = "ENTREZID",keytype = "SYMBOL")[,2]

plot(compareCluster(geneClusters = fcx_als_DCN_RewiredGenes.cluster,fun = "enrichKEGG",pvalueCutoff=0.05,organism="hsa"))
plot(compareCluster(geneClusters = fcx_als_DCN_RewiredGenes.cluster,fun = "enrichPathway",pvalueCutoff=0.05,organism="human"))
 
