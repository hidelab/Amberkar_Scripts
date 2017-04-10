library(igraph)
library(org.Hs.eg.db)
library(clusterProfiler)

setwd('/shared/hidelab2/user/md4zsa/Work/Data/MSMM_RNAseq/MSMM_RNAseq_FinalRelease2')
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
nc=8
blocksize=50
for (r in 2:4){
  for (w in 1:length(neg_weight)){
    for (ont in 4:6){

      poscoexp_network_genes.normal=V(graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016.DCN[[r]]$Low,es = which(E(msbb_rnaseq2016.DCN[[r]]$Low)$weight<pos_weight[w]),names = T)),directed = F))$name
      poscoexp_network_genes.AD=V(graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016.DCN[[r]]$High,es = which(E(msbb_rnaseq2016.DCN[[r]]$High)$weight>pos_weight[w]),names = T)),directed = F))$name
      negcoexp_network_genes.normal=V(graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016.DCN[[r]]$Low,es = which(E(msbb_rnaseq2016.DCN[[r]]$Low)$weight<neg_weight[w]),names = T)),directed = F))$name
      negcoexp_network_genes.AD=V(graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016.DCN[[r]]$High,es = which(E(msbb_rnaseq2016.DCN[[r]]$High)$weight>neg_weight[w]),names = T)),directed = F))$name
      #common_genes=sort(intersect(network_genes.AD,network_genes.normal))
      i<-0
      start<-i*blocksize+1
      end<-min((i+1)*blocksize, length(poscoexp_network_genes.normal))
      poscoexp_go_semsim.normal=poscoexp_go_semsim.AD=c()
      while(start < length(poscoexp_network_genes.normal)){
        cat(paste("In brain region,",names(msbb_rnaseq2016.DCN)[r],"computing semantic similarity for normal genes and ",names(hsGO)[ont-3],"processing block ...",i,"\n",sep=" "))
        input1=poscoexp_network_genes.normal[start:end]
        poscoexp_res.normal = mgeneSim(genes = input1,semData = hsGO[[ont-3]],measure = 'Wang',drop = 'IEA',combine = 'BMA',verbose = T)
        res2=rowSums(poscoexp_res.normal)
        poscoexp_go_semsim.normal=append(poscoexp_network_genes.normal,res2,after = length(poscoexp_network_genes.normal))
        i<-i+1
        start<-i*blocksize+1
        end<-min((i+1)*blocksize, length(poscoexp_network_genes.normal))
      }
      i<-0
      start<-i*blocksize+1
      end<-min((i+1)*blocksize, length(poscoexp_network_genes.AD))
      while(start < length(poscoexp_network_genes.AD)){
        cat(paste("In brain region,",names(msbb_rnaseq2016.DCN)[r],"computing semantic similarity for AD genes and ",names(hsGO)[ont-3],"processing block ...",i,"\n",sep=" "))
        input2=poscoexp_network_genes.AD[start:end]
        poscoexp_res.AD = mgeneSim(genes = input2,semData = hsGO[[ont-3]],measure = 'Wang',drop = 'IEA',combine = 'BMA',verbose = T)
        res3=rowSums(poscoexp_res.AD)
        poscoexp_go_semsim.AD=append(poscoexp_go_semsim.AD,res3,after = length(poscoexp_go_semsim.AD))
        i<-i+1
        start<-i*blocksize+1
        end<-min((i+1)*blocksize, length(poscoexp_network_genes.AD))
      }
      #For NegCoexp
      i<-0
      start<-i*blocksize+1
      end<-min((i+1)*blocksize, length(negcoexp_network_genes.normal))
      negcoexp_go_semsim.normal=negcoexp_go_semsim.AD=c()
      while(start < length(negcoexp_network_genes.normal)){
        cat(paste("In brain region,",names(msbb_rnaseq2016.DCN)[r],"computing semantic similarity for normal genes and ",names(hsGO)[ont-3],"processing block ...",i,"\n",sep=" "))
        input3=negcoexp_network_genes.normal[start:end]
        negcoexp_res.normal = mgeneSim(genes = input3,semData = hsGO[[ont-3]],measure = 'Wang',drop = 'IEA',combine = 'BMA',verbose = T)
        res2=rowSums(negcoexp_res.normal)
        negcoexp_go_semsim.normal=append(negcoexp_network_genes.normal,res2,after = length(negcoexp_network_genes.normal))
        i<-i+1
        start<-i*blocksize+1
        end<-min((i+1)*blocksize, length(negcoexp_network_genes.normal))
      }
      i<-0
      start<-i*blocksize+1
      end<-min((i+1)*blocksize, length(negcoexp_network_genes.AD))
      while(start < length(negcoexp_network_genes.AD)){
        cat(paste("In brain region,",names(msbb_rnaseq2016.DCN)[r],"computing semantic similarity for AD genes and ",names(hsGO)[ont-3],"processing block ...",i,"\n",sep=" "))
        input4=negcoexp_network_genes.AD[start:end]
        negcoexp_res.AD = mgeneSim(genes = input4,semData = hsGO[[ont-3]],measure = 'Wang',drop = 'IEA',combine = 'BMA',verbose = T)
        res3=rowSums(negcoexp_res.AD)
        negcoexp_go_semsim.AD=append(negcoexp_go_semsim.AD,res3,after = length(negcoexp_go_semsim.AD))
        i<-i+1
        start<-i*blocksize+1
        end<-min((i+1)*blocksize, length(negcoexp_network_genes.AD))
      }
      poscoexp_pval_vector[[r]][[w]][[ont]]=-log10(wilcox.test(x=poscoexp_go_semsim.normal,y=poscoexp_go_semsim.AD)$p.val)
      negcoexp_pval_vector[[r]][[w]][[ont]]=-log10(wilcox.test(x=negcoexp_go_semsim.normal,y=negcoexp_go_semsim.AD)$p.val)
    }
  }
}




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