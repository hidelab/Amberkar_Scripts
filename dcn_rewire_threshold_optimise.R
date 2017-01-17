library(igraph)
library(org.Hs.eg.db)
library(clusterProfiler)


msbb_rnaseq2016.DCN=msbb_rnaseq2016_R2.DCN=msbb_rnaseq2016_PLQ.DCN=msbb_rnaseq2016_PLQGenes.cluster=vector(mode = "list",length = 4)
names(msbb_rnaseq2016.DCN)=names(msbb_rnaseq2016_R2.DCN)=names(msbb_rnaseq2016_PLQ.DCN)=names(msbb_rnaseq2016_PLQGenes.cluster)=c("FP","IFG","PHG","STG")

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

pos_weight=seq(1.25,1.75,by=0.05)
neg_weight=seq(0.25,0.75,by=0.05)
poscoexp_pval_vector=negcoexp_pval_vector=vector(mode = 'list',length = 4)
names(poscoexp_pval_vector)=names(negcoexp_pval_vector)=c("FP","IFG","PHG","STG")
poscoexp_pval_vector$FP=poscoexp_pval_vector$IFG=poscoexp_pval_vector$PHG=poscoexp_pval_vector$STG=vector(mode = 'list',length=length(pos_weight))
names(poscoexp_pval_vector$FP)=names(poscoexp_pval_vector$IFG)=names(poscoexp_pval_vector$PHG)=names(poscoexp_pval_vector$STG)=pos_weight
negcoexp_pval_vector$FP=negcoexp_pval_vector$IFG=negcoexp_pval_vector$PHG=negcoexp_pval_vector$STG=vector(mode = 'list',length=length(pos_weight))
names(negcoexp_pval_vector$FP)=names(negcoexp_pval_vector$IFG)=names(negcoexp_pval_vector$PHG)=names(negcoexp_pval_vector$STG)=neg_weight
for(i in 1:length(pos_weight)){
  poscoexp_pval_vector$FP[[i]]=vector(mode = 'list',length=6)
  names(poscoexp_pval_vector$FP[[i]])=c('Degree','Closeness','PageRank','GO.BP_SemSim','GO.CC_SemSim','GO.MF_SemSim')
  poscoexp_pval_vector$IFG[[i]]=vector(mode = 'list',length=6)
  names(poscoexp_pval_vector$IFG[[i]])=c('Degree','Closeness','PageRank','GO.BP_SemSim','GO.CC_SemSim','GO.MF_SemSim')
  poscoexp_pval_vector$PHG[[i]]=vector(mode = 'list',length=6)
  names(poscoexp_pval_vector$PHG[[i]])=c('Degree','Closeness','PageRank','GO.BP_SemSim','GO.CC_SemSim','GO.MF_SemSim')
  poscoexp_pval_vector$STG[[i]]=vector(mode = 'list',length=6)
  names(poscoexp_pval_vector$STG[[i]])=c('Degree','Closeness','PageRank','GO.BP_SemSim','GO.CC_SemSim','GO.MF_SemSim')
  negcoexp_pval_vector$FP[[i]]=vector(mode = 'list',length=6)
  names(negcoexp_pval_vector$FP[[i]])=c('Degree','Closeness','PageRank','GO.BP_SemSim','GO.CC_SemSim','GO.MF_SemSim')
  negcoexp_pval_vector$IFG[[i]]=vector(mode = 'list',length=6)
  names(negcoexp_pval_vector$IFG[[i]])=c('Degree','Closeness','PageRank','GO.BP_SemSim','GO.CC_SemSim','GO.MF_SemSim')
  negcoexp_pval_vector$PHG[[i]]=vector(mode = 'list',length=6)
  names(negcoexp_pval_vector$PHG[[i]])=c('Degree','Closeness','PageRank','GO.BP_SemSim','GO.CC_SemSim','GO.MF_SemSim')
  negcoexp_pval_vector$STG[[i]]=vector(mode = 'list',length=6)
  names(negcoexp_pval_vector$STG[[i]])=c('Degree','Closeness','PageRank','GO.BP_SemSim','GO.CC_SemSim','GO.MF_SemSim')
}

poscoexp_degree=negcoexp_degree=poscoexp_pgr=negcoexp_pgr=poscoexp_cls=negcoexp_cls=vector(mode = "list",length = 4)
names(poscoexp_pval_vector)=names(negcoexp_pval_vector)=names(poscoexp_degree)=names(negcoexp_degree)=names(poscoexp_pgr)=names(negcoexp_pgr)=names(poscoexp_cls)=names(negcoexp_cls)=c("FP","IFG","PHG","STG")
poscoexp_pval_vector$FP=poscoexp_pval_vector$IFG=poscoexp_pval_vector$PHG=poscoexp_pval_vector$STG=vector(mode = "list",length = length(pos_weight))
negcoexp_pval_vector$FP=negcoexp_pval_vector$IFG=negcoexp_pval_vector$PHG=negcoexp_pval_vector$STG=vector(mode = "list",length = length(neg_weight))
names(poscoexp_pval_vector$FP)=names(poscoexp_pval_vector$IFG)=names(poscoexp_pval_vector$PHG)=names(poscoexp_pval_vector$STG)=pos_weight
names(negcoexp_pval_vector$FP)=names(negcoexp_pval_vector$IFG)=names(negcoexp_pval_vector$PHG)=names(negcoexp_pval_vector$STG)=neg_weight

poscoexp_degree$FP=poscoexp_degree$IFG=poscoexp_degree$PHG=poscoexp_degree$STG=vector(mode = "list",length = 2)
negcoexp_degree$FP=negcoexp_degree$IFG=negcoexp_degree$PHG=negcoexp_degree$STG=vector(mode = "list",length = 2)
poscoexp_pgr$FP=poscoexp_pgr$IFG=poscoexp_pgr$PHG=poscoexp_pgr$STG=vector(mode = "list",length = 2)
negcoexp_pgr$FP=negcoexp_pgr$IFG=negcoexp_pgr$PHG=negcoexp_pgr$STG=vector(mode = "list",length = 2)
poscoexp_cls$FP=poscoexp_cls$IFG=poscoexp_cls$PHG=poscoexp_cls$STG=vector(mode = "list",length = 2)
negcoexp_cls$FP=negcoexp_cls$IFG=negcoexp_cls$PHG=negcoexp_cls$STG=vector(mode = "list",length = 2)
names(poscoexp_degree$FP)=names(poscoexp_degree$IFG)=names(poscoexp_degree$PHG)=names(poscoexp_degree$STG)=c("Low","High")
names(negcoexp_degree$FP)=names(negcoexp_degree$IFG)=names(negcoexp_degree$PHG)=names(negcoexp_degree$STG)=c("Low","High")
names(poscoexp_pgr$FP)=names(poscoexp_pgr$IFG)=names(poscoexp_pgr$PHG)=names(poscoexp_pgr$STG)=c("Low","High")
names(negcoexp_pgr$FP)=names(negcoexp_pgr$IFG)=names(negcoexp_pgr$PHG)=names(negcoexp_pgr$STG)=c("Low","High")
names(poscoexp_cls$FP)=names(poscoexp_cls$IFG)=names(poscoexp_cls$PHG)=names(poscoexp_cls$STG)=c("Low","High")
names(negcoexp_cls$FP)=names(negcoexp_cls$IFG)=names(negcoexp_cls$PHG)=names(negcoexp_cls$STG)=c("Low","High")


for(r in 1:4){
  for (w in 1:length(pos_weight)){
    for(c in 1:3){
      poscoexp_degree[[r]]$Low=degree(graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016.DCN[[r]]$Low,es = which(E(msbb_rnaseq2016.DCN[[r]]$Low)$weight<pos_weight[w]),names = T)),directed = F))
      poscoexp_degree[[r]]$High=degree(graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016.DCN[[r]]$High,es = which(E(msbb_rnaseq2016.DCN[[r]]$High)$weight>pos_weight[w]),names = T)),directed = F))
      negcoexp_degree[[r]]$Low=degree(graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016.DCN[[r]]$Low,es = which(E(msbb_rnaseq2016.DCN[[r]]$Low)$weight<neg_weight[w]),names = T)),directed = F))
      negcoexp_degree[[r]]$High=degree(graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016.DCN[[r]]$High,es = which(E(msbb_rnaseq2016.DCN[[r]]$High)$weight>neg_weight[w]),names = T)),directed = F))
      poscoexp_cls[[r]]$Low=closeness(graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016.DCN[[r]]$Low,es = which(E(msbb_rnaseq2016.DCN[[r]]$Low)$weight<pos_weight[w]),names = T)),directed = F))
      poscoexp_cls[[r]]$High=closeness(graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016.DCN[[r]]$High,es = which(E(msbb_rnaseq2016.DCN[[r]]$High)$weight>pos_weight[w]),names = T)),directed = F))
      negcoexp_cls[[r]]$Low=closeness(graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016.DCN[[r]]$Low,es = which(E(msbb_rnaseq2016.DCN[[r]]$Low)$weight<neg_weight[w]),names = T)),directed = F))
      negcoexp_cls[[r]]$High=closeness(graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016.DCN[[r]]$High,es = which(E(msbb_rnaseq2016.DCN[[r]]$High)$weight>neg_weight[w]),names = T)),directed = F))
      poscoexp_pgr[[r]]$Low=page_rank(graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016.DCN[[r]]$Low,es = which(E(msbb_rnaseq2016.DCN[[r]]$Low)$weight<pos_weight[w]),names = T)),directed = F))$vector
      poscoexp_pgr[[r]]$High=page_rank(graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016.DCN[[r]]$High,es = which(E(msbb_rnaseq2016.DCN[[r]]$High)$weight>pos_weight[w]),names = T)),directed = F))$vector
      negcoexp_pgr[[r]]$Low=page_rank(graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016.DCN[[r]]$Low,es = which(E(msbb_rnaseq2016.DCN[[r]]$Low)$weight<neg_weight[w]),names = T)),directed = F))$vector
      negcoexp_pgr[[r]]$High=page_rank(graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016.DCN[[r]]$High,es = which(E(msbb_rnaseq2016.DCN[[r]]$High)$weight>neg_weight[w]),names = T)),directed = F))$vector
      
      poscoexp_pval_vector[[r]][[w]][[1]]=-log10(wilcox.test(x=poscoexp_degree[[r]]$Low,y=poscoexp_degree[[r]]$High)$p.value)
      negcoexp_pval_vector[[r]][[w]][[1]]=-log10(wilcox.test(x=negcoexp_degree[[r]]$Low,y=negcoexp_degree[[r]]$High)$p.value)
      poscoexp_pval_vector[[r]][[w]][[2]]=-log10(wilcox.test(x=poscoexp_pgr[[r]]$Low,y=poscoexp_pgr[[r]]$High)$p.value)
      negcoexp_pval_vector[[r]][[w]][[2]]=-log10(wilcox.test(x=negcoexp_pgr[[r]]$Low,y=negcoexp_pgr[[r]]$High)$p.value)
      poscoexp_pval_vector[[r]][[w]][[3]]=-log10(wilcox.test(x=poscoexp_cls[[r]]$Low,y=poscoexp_cls[[r]]$High)$p.value)
      negcoexp_pval_vector[[r]][[w]][[3]]=-log10(wilcox.test(x=negcoexp_cls[[r]]$Low,y=negcoexp_cls[[r]]$High)$p.value)
      
    }
  }
}
saveRDS(poscoexp_pval_vector,'MSBB_PosCoexp_PvalVector.RDS')
saveRDS(negcoexp_pval_vector,'MSBB_NegCoexp_PvalVector.RDS')
# msbb_dcn.clusters=msbb_plq_dcn.clusters=list()
# msbb_dcn.clusters[[1]]=V(msbb_rnaseq2016.DCN$FP$High)$name
# msbb_dcn.clusters[[2]]=V(msbb_rnaseq2016.DCN$IFG$High)$name
# msbb_dcn.clusters[[3]]=V(msbb_rnaseq2016.DCN$PHG$High)$name
# msbb_dcn.clusters[[4]]=V(msbb_rnaseq2016.DCN$STG$High)$name
# 
# msbb_plq_dcn.clusters[[1]]=V(msbb_rnaseq2016_PLQ.DCN$FP$High)$name
# msbb_plq_dcn.clusters[[2]]=V(msbb_rnaseq2016_PLQ.DCN$IFG$High)$name
# msbb_plq_dcn.clusters[[3]]=V(msbb_rnaseq2016_PLQ.DCN$PHG$High)$name
# msbb_plq_dcn.clusters[[4]]=V(msbb_rnaseq2016_PLQ.DCN$STG$High)$name
# 
# nc=8
# blocksize=50
# for (r in 2:4){
#   for (w in 1:length(neg_weight)){
#     for (ont in 4:6){
#       
#       network_genes.normal=V(graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016.DCN[[r]]$Low,es = which(E(msbb_rnaseq2016.DCN[[r]]$Low)$weight<pos_weight[w]),names = T)),directed = F))$name
#       network_genes.AD=V(graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016.DCN[[r]]$High,es = which(E(msbb_rnaseq2016.DCN[[r]]$High)$weight>pos_weight[w]),names = T)),directed = F))$name
#       #common_genes=sort(intersect(network_genes.AD,network_genes.normal))
#       i<-0
#       start<-i*blocksize+1
#       end<-min((i+1)*blocksize, common_genes)
#       go_semsim.normal=go_semsim.AD=c()
#       while(start < common_genes){
#         cat(paste("In brain region,",names(msbb_dcn.clusters)[r],"computing semantic similarity for",names(hsGO)[ont],"processing block ...",i,"\n",sep=" "))
#         input1=network_genes.normal[start:end]
#         input2=network_genes.AD[start:end]
#         res.normal = mgeneSim(genes = input1,semData = hsGO[[ont-3]],measure = 'Wang',drop = 'IEA',combine = 'BMA',verbose = T)
#         res2=rowSums(res.normal)
#         go_semsim.normal=append(go_semsim.normal,res2,after = length(go_semsim.normal))
#         res.AD = mgeneSim(genes = input2,semData = hsGO[[ont-3]],measure = 'Wang',drop = 'IEA',combine = 'BMA',verbose = T)
#         res3=rowSums(res.AD)
#         go_semsim.AD=append(go_semsim.AD,res3,after = length(go_semsim.AD))
#         i<-i+1
#         start<-i*blocksize+1
#         end<-min((i+1)*blocksize, common_genes)
#       }
#       poscoexp_pval_vector[[r]][[w]][[ont-3]]=go_semsim
#     }
#   }
# }
#     
#       
#          
# 
# df1=cbind(FP_PosCoexp_pval=poscoexp_pval_vector$FP,
#           IFG_PosCoexp_pval=poscoexp_pval_vector$IFG,
#           PHG_PosCoexp_pval=poscoexp_pval_vector$PHG,
#           STG_PosCoexp_pval=poscoexp_pval_vector$STG,
#           FP_NegCoexp_pval=negcoexp_pval_vector$FP,
#           IFG_NegCoexp_pval=negcoexp_pval_vector$IFG,
#           PHG_NegCoexp_pval=negcoexp_pval_vector$PHG,
#           STG_NegCoexp_pval=negcoexp_pval_vector$STG)
# threshold_df=data.frame(df1)
# rownames(threshold_df)=as.numeric(rownames(df1))-1