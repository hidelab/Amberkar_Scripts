library(igraph)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ReactomePA)


msbb_rnaseq2016_R2.DCN=vector(mode = "list",length = 4)
names(msbb_rnaseq2016_R2.DCN)=c("FP","IFG","PHG","STG")
msbb_rnaseq2016_R2.DCN$FP=msbb_rnaseq2016_R2.DCN$IFG=msbb_rnaseq2016_R2.DCN$PHG=msbb_rnaseq2016_R2.DCN$STG=vector(mode="list",length=2)
msbb_rnaseq2016_R2.DCN$FP=readRDS("DCN/FP_R2_DCN_FDR01.RDS")
msbb_rnaseq2016_R2.DCN$IFG=readRDS("DCN/IFG_R2_DCN_FDR01.RDS")
msbb_rnaseq2016_R2.DCN$PHG=readRDS("DCN/PHG_R2_DCN_FDR01.RDS")
msbb_rnaseq2016_R2.DCN$STG=readRDS("DCN/STG_R2_DCN_FDR01.RDS")

msbb_rnaseq2016_R2Genes.cluster=msbb_rnaseq2016_R2.PosCoexp=msbb_rnaseq2016_R2.NegCoexp=msbb_rnaseq2016_R2Genes_PosCoexp.rewiredDF=msbb_rnaseq2016_R2Genes_NegCoexp.rewiredDF=vector(mode = "list",length = 4)
names(msbb_rnaseq2016_R2.PosCoexp)=names(msbb_rnaseq2016_R2.NegCoexp)=names(msbb_rnaseq2016_R2Genes_PosCoexp.rewiredDF)=names(msbb_rnaseq2016_R2Genes_NegCoexp.rewiredDF)=names(msbb_rnaseq2016_R2Genes.cluster)
msbb_rnaseq2016_R2.PosCoexp$FP=msbb_rnaseq2016_R2.PosCoexp$IFG=msbb_rnaseq2016_R2.PosCoexp$PHG=msbb_rnaseq2016_R2.PosCoexp$STG=vector(mode = "list",length = 2)
msbb_rnaseq2016_R2.NegCoexp$FP=msbb_rnaseq2016_R2.NegCoexp$IFG=msbb_rnaseq2016_R2.NegCoexp$PHG=msbb_rnaseq2016_R2.NegCoexp$STG=vector(mode = "list",length = 2)
names(msbb_rnaseq2016_R2.PosCoexp$FP)=names(msbb_rnaseq2016_R2.NegCoexp$FP)=c("Low","High")
names(msbb_rnaseq2016_R2.PosCoexp$IFG)=names(msbb_rnaseq2016_R2.NegCoexp$IFG)=c("Low","High")
names(msbb_rnaseq2016_R2.PosCoexp$PHG)=names(msbb_rnaseq2016_R2.NegCoexp$PHG)=c("Low","High")
names(msbb_rnaseq2016_R2.PosCoexp$STG)=names(msbb_rnaseq2016_R2.NegCoexp$STG)=c("Low","High")

msbb_rnaseq2016_R2.PosCoexp$FP$Low=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016_R2.DCN$FP$Low,es = which(E(msbb_rnaseq2016_R2.DCN$FP$Low)$weight>1.5),names = T)),directed = F)
msbb_rnaseq2016_R2.PosCoexp$IFG$Low=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016_R2.DCN$IFG$Low,es = which(E(msbb_rnaseq2016_R2.DCN$IFG$Low)$weight>1.5),names = T)),directed = F)
msbb_rnaseq2016_R2.PosCoexp$PHG$Low=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016_R2.DCN$PHG$Low,es = which(E(msbb_rnaseq2016_R2.DCN$PHG$Low)$weight>1.5),names = T)),directed = F)
msbb_rnaseq2016_R2.PosCoexp$STG$Low=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016_R2.DCN$STG$Low,es = which(E(msbb_rnaseq2016_R2.DCN$STG$Low)$weight>1.5),names = T)),directed = F)

msbb_rnaseq2016_R2.PosCoexp$FP$High=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016_R2.DCN$FP$High,es = which(E(msbb_rnaseq2016_R2.DCN$FP$High)$weight>1.5),names = T)),directed = F)
msbb_rnaseq2016_R2.PosCoexp$IFG$High=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016_R2.DCN$IFG$High,es = which(E(msbb_rnaseq2016_R2.DCN$IFG$High)$weight>1.5),names = T)),directed = F)
msbb_rnaseq2016_R2.PosCoexp$PHG$High=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016_R2.DCN$PHG$High,es = which(E(msbb_rnaseq2016_R2.DCN$PHG$High)$weight>1.5),names = T)),directed = F)
msbb_rnaseq2016_R2.PosCoexp$STG$High=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016_R2.DCN$STG$High,es = which(E(msbb_rnaseq2016_R2.DCN$STG$High)$weight>1.5),names = T)),directed = F)

msbb_rnaseq2016_R2.NegCoexp$FP$Low=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016_R2.DCN$FP$Low,es = which(E(msbb_rnaseq2016_R2.DCN$FP$Low)$weight<0.5),names = T)),directed = F)
msbb_rnaseq2016_R2.NegCoexp$IFG$Low=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016_R2.DCN$IFG$Low,es = which(E(msbb_rnaseq2016_R2.DCN$IFG$Low)$weight<0.5),names = T)),directed = F)
msbb_rnaseq2016_R2.NegCoexp$PHG$Low=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016_R2.DCN$PHG$Low,es = which(E(msbb_rnaseq2016_R2.DCN$PHG$Low)$weight<0.5),names = T)),directed = F)
msbb_rnaseq2016_R2.NegCoexp$STG$Low=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016_R2.DCN$STG$Low,es = which(E(msbb_rnaseq2016_R2.DCN$STG$Low)$weight<0.5),names = T)),directed = F)

msbb_rnaseq2016_R2.NegCoexp$FP$High=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016_R2.DCN$FP$High,es = which(E(msbb_rnaseq2016_R2.DCN$FP$High)$weight<0.5),names = T)),directed = F)
msbb_rnaseq2016_R2.NegCoexp$IFG$High=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016_R2.DCN$IFG$High,es = which(E(msbb_rnaseq2016_R2.DCN$IFG$High)$weight<0.5),names = T)),directed = F)
msbb_rnaseq2016_R2.NegCoexp$PHG$High=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016_R2.DCN$PHG$High,es = which(E(msbb_rnaseq2016_R2.DCN$PHG$High)$weight<0.5),names = T)),directed = F)
msbb_rnaseq2016_R2.NegCoexp$STG$High=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016_R2.DCN$STG$High,es = which(E(msbb_rnaseq2016_R2.DCN$STG$High)$weight<0.5),names = T)),directed = F)

msbb_rnaseq2016_PLQGenes.DEG=msbb_rnaseq2016_PLQGenes=msbb_rnaseq2016_AllAnalyses=vector(mode = "list",length = 4)
names(msbb_rnaseq2016_PLQGenes.DEG)=names(msbb_rnaseq2016_PLQGenes)=names(msbb_rnaseq2016_AllAnalyses)
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

msbb_rnaseq2016_R2Genes_PosCoexp.rewiredDF$FP=data.frame(Common_Genes=intersect(V(msbb_rnaseq2016_R2.PosCoexp$FP$Low)$name,V(msbb_rnaseq2016_R2.PosCoexp$FP$High)$name),
                                                         Low=degree(graph = msbb_rnaseq2016_R2.PosCoexp$FP$Low,v = intersect(V(msbb_rnaseq2016_R2.PosCoexp$FP$Low)$name,V(msbb_rnaseq2016_R2.PosCoexp$FP$High)$name)),
                                                         High=degree(graph = msbb_rnaseq2016_R2.PosCoexp$FP$High,v = intersect(V(msbb_rnaseq2016_R2.PosCoexp$FP$Low)$name,V(msbb_rnaseq2016_R2.PosCoexp$FP$High)$name)),
                                                         Difference=degree(graph = msbb_rnaseq2016_R2.PosCoexp$FP$Low,v = intersect(V(msbb_rnaseq2016_R2.PosCoexp$FP$Low)$name,V(msbb_rnaseq2016_R2.PosCoexp$FP$High)$name))-degree(graph = msbb_rnaseq2016_R2.PosCoexp$FP$High,v = intersect(V(msbb_rnaseq2016_R2.PosCoexp$FP$Low)$name,V(msbb_rnaseq2016_R2.PosCoexp$FP$High)$name)),stringsAsFactors = F)
msbb_rnaseq2016_R2Genes_PosCoexp.rewiredDF$IFG=data.frame(Common_Genes=intersect(V(msbb_rnaseq2016_R2.PosCoexp$IFG$Low)$name,V(msbb_rnaseq2016_R2.PosCoexp$IFG$High)$name),
                                                          Low=degree(graph = msbb_rnaseq2016_R2.PosCoexp$IFG$Low,v = intersect(V(msbb_rnaseq2016_R2.PosCoexp$IFG$Low)$name,V(msbb_rnaseq2016_R2.PosCoexp$IFG$High)$name)),
                                                          High=degree(graph = msbb_rnaseq2016_R2.PosCoexp$IFG$High,v = intersect(V(msbb_rnaseq2016_R2.PosCoexp$IFG$Low)$name,V(msbb_rnaseq2016_R2.PosCoexp$IFG$High)$name)),
                                                          Difference=degree(graph = msbb_rnaseq2016_R2.PosCoexp$IFG$Low,v = intersect(V(msbb_rnaseq2016_R2.PosCoexp$IFG$Low)$name,V(msbb_rnaseq2016_R2.PosCoexp$IFG$High)$name))-degree(graph = msbb_rnaseq2016_R2.PosCoexp$IFG$High,v = intersect(V(msbb_rnaseq2016_R2.PosCoexp$IFG$Low)$name,V(msbb_rnaseq2016_R2.PosCoexp$IFG$High)$name)),stringsAsFactors = F)
msbb_rnaseq2016_R2Genes_PosCoexp.rewiredDF$PHG=data.frame(Common_Genes=intersect(V(msbb_rnaseq2016_R2.PosCoexp$PHG$Low)$name,V(msbb_rnaseq2016_R2.PosCoexp$PHG$High)$name),
                                                          Low=degree(graph = msbb_rnaseq2016_R2.PosCoexp$PHG$Low,v = intersect(V(msbb_rnaseq2016_R2.PosCoexp$PHG$Low)$name,V(msbb_rnaseq2016_R2.PosCoexp$PHG$High)$name)),
                                                          High=degree(graph = msbb_rnaseq2016_R2.PosCoexp$PHG$High,v = intersect(V(msbb_rnaseq2016_R2.PosCoexp$PHG$Low)$name,V(msbb_rnaseq2016_R2.PosCoexp$PHG$High)$name)),
                                                          Difference=degree(graph = msbb_rnaseq2016_R2.PosCoexp$PHG$Low,v = intersect(V(msbb_rnaseq2016_R2.PosCoexp$PHG$Low)$name,V(msbb_rnaseq2016_R2.PosCoexp$PHG$High)$name))-degree(graph = msbb_rnaseq2016_R2.PosCoexp$PHG$High,v = intersect(V(msbb_rnaseq2016_R2.PosCoexp$PHG$Low)$name,V(msbb_rnaseq2016_R2.PosCoexp$PHG$High)$name)),stringsAsFactors = F)
msbb_rnaseq2016_R2Genes_PosCoexp.rewiredDF$STG=data.frame(Common_Genes=intersect(V(msbb_rnaseq2016_R2.PosCoexp$STG$Low)$name,V(msbb_rnaseq2016_R2.PosCoexp$STG$High)$name),
                                                          Low=degree(graph = msbb_rnaseq2016_R2.PosCoexp$STG$Low,v = intersect(V(msbb_rnaseq2016_R2.PosCoexp$STG$Low)$name,V(msbb_rnaseq2016_R2.PosCoexp$STG$High)$name)),
                                                          High=degree(graph = msbb_rnaseq2016_R2.PosCoexp$STG$High,v = intersect(V(msbb_rnaseq2016_R2.PosCoexp$STG$Low)$name,V(msbb_rnaseq2016_R2.PosCoexp$STG$High)$name)),
                                                          Difference=degree(graph = msbb_rnaseq2016_R2.PosCoexp$STG$Low,v = intersect(V(msbb_rnaseq2016_R2.PosCoexp$STG$Low)$name,V(msbb_rnaseq2016_R2.PosCoexp$STG$High)$name))-degree(graph = msbb_rnaseq2016_R2.PosCoexp$STG$High,v = intersect(V(msbb_rnaseq2016_R2.PosCoexp$STG$Low)$name,V(msbb_rnaseq2016_R2.PosCoexp$STG$High)$name)),stringsAsFactors = F)


msbb_rnaseq2016_AllAnalyses$FP=list(PLQ_Genes=msbb_rnaseq2016_PLQGenes$FP$Genes[which(msbb_rnaseq2016_PLQGenes$FP$FDR<=0.05)],
                                    PLQ_DEG=msbb_rnaseq2016_PLQGenes.DEG$FP$row[which(abs(msbb_rnaseq2016_PLQGenes.DEG$FP$log2FoldChange)>1&msbb_rnaseq2016_PLQGenes.DEG$FP$padj<=0.05)],
                                    R2_DCN=V(msbb_rnaseq2016_R2.DCN$FP$High)$name,
                                    PosCoexp_Loss=msbb_rnaseq2016_R2Genes_PosCoexp.rewiredDF$FP$Common_Genes[which(msbb_rnaseq2016_R2Genes_PosCoexp.rewiredDF$FP$Difference>0)],
                                    PosCoexp_Gain=msbb_rnaseq2016_R2Genes_PosCoexp.rewiredDF$FP$Common_Genes[which(msbb_rnaseq2016_R2Genes_PosCoexp.rewiredDF$FP$Difference<0)],
                                    NegCoexp_Loss=msbb_rnaseq2016_R2Genes_NegCoexp.rewiredDF$FP$Common_Genes[which(msbb_rnaseq2016_R2Genes_NegCoexp.rewiredDF$FP$Difference>0)],
                                    NegCoexp_Gain=msbb_rnaseq2016_R2Genes_NegCoexp.rewiredDF$FP$Common_Genes[which(msbb_rnaseq2016_R2Genes_NegCoexp.rewiredDF$FP$Difference<0)])
msbb_rnaseq2016_AllAnalyses$IFG=list(PLQ_Genes=msbb_rnaseq2016_PLQGenes$IFG$Genes[which(msbb_rnaseq2016_PLQGenes$IFG$FDR<=0.01)],
                                     PLQ_DEG=msbb_rnaseq2016_PLQGenes.DEG$IFG$row[which(abs(msbb_rnaseq2016_PLQGenes.DEG$IFG$log2FoldChange)>1&msbb_rnaseq2016_PLQGenes.DEG$IFG$padj<=0.05)],
                                     R2_DCN=V(msbb_rnaseq2016_R2.DCN$IFG$High)$name,
                                     PosCoexp_Loss=msbb_rnaseq2016_R2Genes_PosCoexp.rewiredDF$IFG$Common_Genes[which(msbb_rnaseq2016_R2Genes_PosCoexp.rewiredDF$IFG$Difference>0)],
                                     PosCoexp_Gain=msbb_rnaseq2016_R2Genes_PosCoexp.rewiredDF$IFG$Common_Genes[which(msbb_rnaseq2016_R2Genes_PosCoexp.rewiredDF$IFG$Difference<0)],
                                     NegCoexp_Loss=msbb_rnaseq2016_R2Genes_NegCoexp.rewiredDF$IFG$Common_Genes[which(msbb_rnaseq2016_R2Genes_NegCoexp.rewiredDF$IFG$Difference>0)],
                                     NegCoexp_Gain=msbb_rnaseq2016_R2Genes_NegCoexp.rewiredDF$IFG$Common_Genes[which(msbb_rnaseq2016_R2Genes_NegCoexp.rewiredDF$IFG$Difference<0)])
msbb_rnaseq2016_AllAnalyses$PHG=list(PLQ_Genes=msbb_rnaseq2016_PLQGenes$PHG$Genes[which(msbb_rnaseq2016_PLQGenes$PHG$FDR<=0.01)],
                                     PLQ_DEG=msbb_rnaseq2016_PLQGenes.DEG$PHG$row[which(abs(msbb_rnaseq2016_PLQGenes.DEG$PHG$log2FoldChange)>1&msbb_rnaseq2016_PLQGenes.DEG$PHG$padj<=0.05)],
                                     R2_DCN=V(msbb_rnaseq2016_R2.DCN$PHG$High)$name,
                                     PosCoexp_Loss=msbb_rnaseq2016_R2Genes_PosCoexp.rewiredDF$PHG$Common_Genes[which(msbb_rnaseq2016_R2Genes_PosCoexp.rewiredDF$PHG$Difference>0)],
                                     PosCoexp_Gain=msbb_rnaseq2016_R2Genes_PosCoexp.rewiredDF$PHG$Common_Genes[which(msbb_rnaseq2016_R2Genes_PosCoexp.rewiredDF$PHG$Difference<0)],
                                     NegCoexp_Loss=msbb_rnaseq2016_R2Genes_NegCoexp.rewiredDF$PHG$Common_Genes[which(msbb_rnaseq2016_R2Genes_NegCoexp.rewiredDF$PHG$Difference>0)],
                                     NegCoexp_Gain=msbb_rnaseq2016_R2Genes_NegCoexp.rewiredDF$PHG$Common_Genes[which(msbb_rnaseq2016_R2Genes_NegCoexp.rewiredDF$PHG$Difference<0)])
msbb_rnaseq2016_AllAnalyses$STG=list(PLQ_Genes=msbb_rnaseq2016_PLQGenes$STG$Genes[which(msbb_rnaseq2016_PLQGenes$STG$FDR<=0.01)],
                                     PLQ_DEG=msbb_rnaseq2016_PLQGenes.DEG$STG$row[which(abs(msbb_rnaseq2016_PLQGenes.DEG$STG$log2FoldChange)>1&msbb_rnaseq2016_PLQGenes.DEG$STG$padj<=0.05)],
                                     R2_DCN=V(msbb_rnaseq2016_R2.DCN$STG$High)$name,
                                     PosCoexp_Loss=msbb_rnaseq2016_R2Genes_PosCoexp.rewiredDF$STG$Common_Genes[which(msbb_rnaseq2016_R2Genes_PosCoexp.rewiredDF$STG$Difference>0)],
                                     PosCoexp_Gain=msbb_rnaseq2016_R2Genes_PosCoexp.rewiredDF$STG$Common_Genes[which(msbb_rnaseq2016_R2Genes_PosCoexp.rewiredDF$STG$Difference<0)],
                                     NegCoexp_Loss=msbb_rnaseq2016_R2Genes_NegCoexp.rewiredDF$STG$Common_Genes[which(msbb_rnaseq2016_R2Genes_NegCoexp.rewiredDF$STG$Difference>0)],
                                     NegCoexp_Gain=msbb_rnaseq2016_R2Genes_NegCoexp.rewiredDF$STG$Common_Genes[which(msbb_rnaseq2016_R2Genes_NegCoexp.rewiredDF$STG$Difference<0)])


#Test for determining edge weight threshold
pos_weight=seq(1.25,1.75,by=0.01)
neg_weight=seq(0.25,0.75,by=0.01)
poscoexp_pval_vector=negcoexp_pval_vector=poscoexp_degree=negcoexp_degree=vector(mode = "list",length = 4)
names(poscoexp_pval_vector)=names(negcoexp_pval_vector)=names(pval_vector)=names(poscoexp_degree)=names(negcoexp_degree)=c("FP","IFG","PHG","STG")
poscoexp_pval_vector$FP=poscoexp_pval_vector$IFG=poscoexp_pval_vector$PHG=poscoexp_pval_vector$STG=vector(mode = "list",length = length(pos_weight))
negcoexp_pval_vector$FP=negcoexp_pval_vector$IFG=negcoexp_pval_vector$PHG=negcoexp_pval_vector$STG=vector(mode = "list",length = length(neg_weight))
names(poscoexp_pval_vector$FP)=names(poscoexp_pval_vector$IFG)=names(poscoexp_pval_vector$PHG)=names(poscoexp_pval_vector$STG)=pos_weight
names(negcoexp_pval_vector$FP)=names(negcoexp_pval_vector$IFG)=names(negcoexp_pval_vector$PHG)=names(negcoexp_pval_vector$STG)=neg_weight

poscoexp_degree$FP=poscoexp_degree$IFG=poscoexp_degree$PHG=poscoexp_degree$STG=vector(mode = "list",length = 2)
negcoexp_degree$FP=negcoexp_degree$IFG=negcoexp_degree$PHG=negcoexp_degree$STG=vector(mode = "list",length = 2)
names(poscoexp_degree$FP)=names(poscoexp_degree$IFG)=names(poscoexp_degree$PHG)=names(poscoexp_degree$STG)=c("Low","High")
names(negcoexp_degree$FP)=names(negcoexp_degree$IFG)=names(negcoexp_degree$PHG)=names(negcoexp_degree$STG)=c("Low","High")

for(r in 1:4)
{
  for (w in 1:length(pos_weight))
  {
    poscoexp_degree[[r]]$Low=degree(graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016_R2.DCN[[r]]$Low,es = which(E(msbb_rnaseq2016_R2.DCN[[r]]$Low)$weight<pos_weight[w]),names = T)),directed = F))
    poscoexp_degree[[r]]$High=degree(graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016_R2.DCN[[r]]$High,es = which(E(msbb_rnaseq2016_R2.DCN[[r]]$High)$weight>pos_weight[w]),names = T)),directed = F))
    negcoexp_degree[[r]]$Low=degree(graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016_R2.DCN[[r]]$Low,es = which(E(msbb_rnaseq2016_R2.DCN[[r]]$Low)$weight<neg_weight[w]),names = T)),directed = F))
    negcoexp_degree[[r]]$High=degree(graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016_R2.DCN[[r]]$High,es = which(E(msbb_rnaseq2016_R2.DCN[[r]]$High)$weight>neg_weight[w]),names = T)),directed = F))
    poscoexp_pval_vector[[r]][[w]]=-log10(wilcox.test(x=poscoexp_degree[[r]]$Low,y=poscoexp_degree[[r]]$High)$p.value)
    negcoexp_pval_vector[[r]][[w]]=-log10(wilcox.test(x=negcoexp_degree[[r]]$Low,y=negcoexp_degree[[r]]$High)$p.value)
  }
  
}


