library(igraph)
library(org.Hs.eg.db)
library(GOSemSim)
libr

OA_analysis=readRDS("~/Work/Collaborations/OA_EleZeggini/OA_D_vs_C_cpm_corr_change_allResults.RDS")
OA_DCN=vector(mode = "list",length = 2)
names(OA_DCN)=c("Healthy","Diseased")
OA_DCN$Healthy=graph.data.frame(d = OA_analysis[which(OA_analysis$FDR<=0.1),c(2:3)],directed = F)
OA_DCN$Diseased=graph.data.frame(d = OA_analysis[which(OA_analysis$FDR<=0.1),c(2:3)],directed = F)
oa_healthy_poscoexp=subset(OA_analysis,subset = (FDR<0.1&r.c>0.25))
oa_disease_poscoexp=subset(OA_analysis,subset = (FDR<0.1&r.t>0.25))
oa_healthy_negcoexp=subset(OA_analysis,subset = (FDR<0.1&r.c< -0.25))
oa_disease_negcoexp=subset(OA_analysis,subset = (FDR<0.1&r.t< -0.25))
oa_healthy_poscoexp.genes=union(oa_healthy_poscoexp$Gene.A,oa_healthy_poscoexp$Gene.B)
oa_disease_poscoexp.genes=union(oa_disease_poscoexp$Gene.A,oa_disease_poscoexp$Gene.B)
oa_healthy_negcoexp.genes=union(oa_healthy_negcoexp$Gene.A,oa_healthy_negcoexp$Gene.B)
oa_disease_negcoexp.genes=union(oa_disease_negcoexp$Gene.A,oa_disease_negcoexp$Gene.B)

E(OA_DCN$Healthy)$weight=OA_analysis[which((OA_analysis$FDR<=0.1)),]$r.c+1
E(OA_DCN$Diseased)$weight=OA_analysis[which((OA_analysis$FDR<=0.1)),]$r.t+1

g1=graph.data.frame(d = ends(graph = OA_DCN$Healthy,es = which(E(OA_DCN$Healthy)$weight>=1.5),names = T),directed = F)
g2=graph.data.frame(d = ends(graph = OA_DCN$Diseased,es = which(E(OA_DCN$Diseased)$weight>=1.5),names = T),directed = F)
g3=graph.data.frame(d = ends(graph = OA_DCN$Healthy,es = which(E(OA_DCN$Healthy)$weight<=0.5),names = T),directed = F)
g4=graph.data.frame(d = ends(graph = OA_DCN$Diseased,es = which(E(OA_DCN$Diseased)$weight<=0.5),names = T),directed = F)

OA_DCN.rewired=vector(mode = "list",length = 6)
names(OA_DCN.rewired)=c("PosCorr","NegCorr","Healthy_Uniq_PosCorr","Disease_Uniq_PosCorr","Healthy_Uniq_NegCorr","Disease_Uniq_NegCorr")
OA_DCN.rewired$PosCorr=data.frame(Genes=names(degree(graph = g1,v = intersect(V(g1)$name,V(g2)$name))),
                                  OA_Healthy=degree(graph = g1,v = intersect(V(g1)$name,V(g2)$name)),
                                  OA_Diseased=degree(graph = g2,v = intersect(V(g1)$name,V(g2)$name)),
                                  Diff=degree(graph = g1,v = intersect(V(g1)$name,V(g2)$name))-degree(graph = g2,v = intersect(V(g1)$name,V(g2)$name)),stringsAsFactors = F)

OA_DCN.rewired$NegCorr=data.frame(Genes=names(degree(graph = g3,v = intersect(V(g3)$name,V(g4)$name))),
                                  OA_Healthy=degree(graph = g3,v = intersect(V(g3)$name,V(g4)$name)),
                                  OA_Diseased=degree(graph = g4,v = intersect(V(g3)$name,V(g4)$name)),
                                  Diff=degree(graph = g3,v = intersect(V(g3)$name,V(g4)$name))-degree(graph = g4,v = intersect(V(g3)$name,V(g4)$name)),stringsAsFactors = F)
OA_DCN.rewired$Healthy_Uniq_PosCorr=data.frame(OA_Healthy_UniqGenes=names(degree(graph = OA_DCN$Healthy,v = which(V(OA_DCN$Healthy)$name%in%setdiff(oa_healthy_poscoexp.genes,oa_disease_poscoexp.genes)))),
                                               Degree_OA_Healthy_UniqGenes=unname(degree(graph = OA_DCN$Healthy,v = which(V(OA_DCN$Healthy)$name%in%setdiff(oa_healthy_poscoexp.genes,oa_disease_poscoexp.genes)))),
                                               Degree_OA_Disease_UniqGenes=0,
                                               Diff=unname(degree(graph = OA_DCN$Healthy,v = which(V(OA_DCN$Healthy)$name%in%setdiff(oa_healthy_poscoexp.genes,oa_disease_poscoexp.genes))))-0,stringsAsFactors = F)
OA_DCN.rewired$Disease_Uniq_PosCorr=data.frame(OA_Disease_UniqGenes=names(degree(graph = OA_DCN$Disease,v = which(V(OA_DCN$Disease)$name%in%setdiff(oa_disease_poscoexp.genes,oa_healthy_poscoexp.genes)))),
                                               Degree_OA_Disease_UniqGenes=unname(degree(graph = OA_DCN$Disease,v = which(V(OA_DCN$Disease)$name%in%setdiff(oa_disease_poscoexp.genes,oa_healthy_poscoexp.genes)))),
                                               Degree_OA_Healthy_UniqGenes=0,
                                               Diff=unname(degree(graph = OA_DCN$Disease,v = which(V(OA_DCN$Disease)$name%in%setdiff(oa_disease_poscoexp.genes,oa_healthy_poscoexp.genes))))-0,stringsAsFactors = F)

OA_DCN.rewired$Healthy_Uniq_NegCorr=data.frame(OA_Healthy_UniqGenes=names(degree(graph = OA_DCN$Healthy,v = which(V(OA_DCN$Healthy)$name%in%setdiff(oa_healthy_negcoexp.genes,oa_disease_negcoexp.genes)))),
                                               Degree_OA_Healthy_UniqGenes=unname(degree(graph = OA_DCN$Healthy,v = which(V(OA_DCN$Healthy)$name%in%setdiff(oa_healthy_negcoexp.genes,oa_disease_negcoexp.genes)))),
                                               Degree_OA_Disease_UniqGenes=0,
                                               Diff=unname(degree(graph = OA_DCN$Healthy,v = which(V(OA_DCN$Healthy)$name%in%setdiff(oa_healthy_negcoexp.genes,oa_disease_negcoexp.genes))))-0,stringsAsFactors = F)

OA_DCN.rewired$Disease_Uniq_NegCorr=data.frame(OA_Disease_UniqGenes=names(degree(graph = OA_DCN$Disease,v = which(V(OA_DCN$Disease)$name%in%setdiff(oa_disease_negcoexp.genes,oa_healthy_negcoexp.genes)))),
                                               Degree_OA_Disease_UniqGenes=unname(degree(graph = OA_DCN$Disease,v = which(V(OA_DCN$Disease)$name%in%setdiff(oa_disease_negcoexp.genes,oa_healthy_negcoexp.genes)))),
                                               Degree_OA_Healthy_UniqGenes=0,
                                               Diff=unname(degree(graph = OA_DCN$Disease,v = which(V(OA_DCN$Disease)$name%in%setdiff(oa_disease_negcoexp.genes,oa_healthy_negcoexp.genes))))-0,stringsAsFactors = F)

OA_DCN.rewired_cluster=vector(mode = "list",length = 9)
names(OA_DCN.rewired_cluster)=c("PosCoexp_Gain","PosCoexp_Loss","NegCoexp_Gain","NegCoexp_Loss","PosCoexp_UniqHealthy","PosCoexp_UniqDisease","NegCoexp_UniqHealthy","NegCoexp_UniqDisease","DEG")
OA_DCN.rewired_cluster$PosCoexp_Gain=select(x = org.Hs.eg.db,keys=OA_DCN.rewired$PosCorr$Genes[which(OA_DCN.rewired$PosCorr$Diff>0)],columns = "ENTREZID",keytype = "SYMBOL")[,2]
OA_DCN.rewired_cluster$PosCoexp_Loss=select(x = org.Hs.eg.db,keys=OA_DCN.rewired$PosCorr$Genes[which(OA_DCN.rewired$PosCorr$Diff<0)],columns = "ENTREZID",keytype = "SYMBOL")[,2]
OA_DCN.rewired_cluster$NegCoexp_Gain=select(x = org.Hs.eg.db,keys=OA_DCN.rewired$NegCorr$Genes[which(OA_DCN.rewired$PosCorr$Diff>0)],columns = "ENTREZID",keytype = "SYMBOL")[,2]
OA_DCN.rewired_cluster$NegCoexp_Loss=select(x = org.Hs.eg.db,keys=OA_DCN.rewired$NegCorr$Genes[which(OA_DCN.rewired$PosCorr$Diff<0)],columns = "ENTREZID",keytype = "SYMBOL")[,2]
OA_DCN.rewired_cluster$PosCoexp_UniqHealthy=select(x = org.Hs.eg.db,keys=OA_DCN.rewired$Healthy_Uniq_PosCorr$OA_Healthy_UniqGenes,columns = "ENTREZID",keytype = "SYMBOL")[,2]
OA_DCN.rewired_cluster$PosCoexp_UniqDisease=select(x = org.Hs.eg.db,keys=OA_DCN.rewired$Disease_Uniq_PosCorr$OA_Disease_UniqGenes,columns = "ENTREZID",keytype = "SYMBOL")[,2]
OA_DCN.rewired_cluster$NegCoexp_UniqHealthy=select(x = org.Hs.eg.db,keys=OA_DCN.rewired$Healthy_Uniq_NegCorr$OA_Healthy_UniqGenes,columns = "ENTREZID",keytype = "SYMBOL")[,2]
OA_DCN.rewired_cluster$NegCoexp_UniqDisease=select(x = org.Hs.eg.db,keys=OA_DCN.rewired$Disease_Uniq_NegCorr$OA_Disease_UniqGenes,columns = "ENTREZID",keytype = "SYMBOL")[,2]
OA_DCN.rewired_cluster$DEG=select(x = org.Hs.eg.db,keys = oa_deg_df$genes[which(abs(oa_deg_df$logFC)>log(2)&oa_deg_df$FDR<=0.05)],columns = "ENTREZID",keytype = "SYMBOL")[,2]

OA_DCN.rewired_cluster2=vector(mode = "list",length = 8)
names(OA_DCN.rewired_cluster2)=c("PosCoexp_Gain","PosCoexp_Loss","NegCoexp_Gain","NegCoexp_Loss","PosCoexp_UniqHealthy","PosCoexp_UniqDisease","NegCoexp_UniqHealthy","NegCoexp_UniqDisease")
OA_DCN.rewired_cluster2$PosCoexp_Gain=OA_DCN.rewired$PosCorr$Genes[which(OA_DCN.rewired$PosCorr$Diff>0)]
OA_DCN.rewired_cluster2$PosCoexp_Loss=OA_DCN.rewired$PosCorr$Genes[which(OA_DCN.rewired$PosCorr$Diff<0)]
OA_DCN.rewired_cluster2$NegCoexp_Gain=OA_DCN.rewired$NegCorr$Genes[which(OA_DCN.rewired$PosCorr$Diff>0)]
OA_DCN.rewired_cluster2$NegCoexp_Loss=OA_DCN.rewired$NegCorr$Genes[which(OA_DCN.rewired$PosCorr$Diff<0)]
OA_DCN.rewired_cluster2$PosCoexp_UniqHealthy=OA_DCN.rewired$Healthy_Uniq_PosCorr$OA_Healthy_UniqGenes
OA_DCN.rewired_cluster2$PosCoexp_UniqDisease=OA_DCN.rewired$Disease_Uniq_PosCorr$OA_Disease_UniqGenes
OA_DCN.rewired_cluster2$NegCoexp_UniqHealthy=OA_DCN.rewired$Healthy_Uniq_NegCorr$OA_Healthy_UniqGenes
OA_DCN.rewired_cluster2$NegCoexp_UniqDisease=OA_DCN.rewired$Disease_Uniq_NegCorr$OA_Disease_UniqGenes


write.table(OA_DCN.rewired$PosCorr,"OA_DCN.rewired_PosCorr_Genes.txt",sep = "\t",col.names = T,row.names = F,quote=F)
write.table(OA_DCN.rewired$NegCorr,"OA_DCN.rewired_NegCorr_Genes.txt",sep = "\t",col.names = T,row.names = F,quote=F)
write.table(OA_DCN.rewired$Healthy_Uniq_PosCorr,'OA_DCN.rewired_Healthy_Uniq_PosCorr_Genes.txt',sep = "\t",col.names = T,row.names = F,quote=F)
write.table(OA_DCN.rewired$Disease_Uniq_PosCorr,'OA_DCN.rewired_Disease_Uniq_PosCorr_Genes.txt',sep = "\t",col.names = T,row.names = F,quote=F)
write.table(OA_DCN.rewired$Healthy_Uniq_NegCorr,'OA_DCN.rewired_Healthy_Uniq_NegCorr_Genes.txt',sep = "\t",col.names = T,row.names = F,quote=F)
write.table(OA_DCN.rewired$Disease_Uniq_NegCorr,'OA_DCN.rewired_Disease_Uniq_NegCorr_Genes.txt',sep = "\t",col.names = T,row.names = F,quote=F)
write(OA_DCN.rewired$Healthy_Uniq_PosCorr$OA_Healthy_UniqGenes,'OA_DCN.rewired_Healthy_Uniq_PosCorr_Genelist.txt',sep = '\n')
write(OA_DCN.rewired$Disease_Uniq_PosCorr$OA_Disease_UniqGenes,'OA_DCN.rewired_Disease_Uniq_PosCorr_Genelist.txt',sep = '\n')
write(OA_DCN.rewired$Healthy_Uniq_NegCorr$OA_Healthy_UniqGenes,'OA_DCN.rewired_Healthy_Uniq_NegCorr_Genelist.txt',sep = '\n')
write(OA_DCN.rewired$Disease_Uniq_NegCorr$OA_Disease_UniqGenes,'OA_DCN.rewired_Disease_Uniq_NegCorr_Genelist.txt',sep = '\n')
write.graph(graph = OA_DCN$Diseased,file = 'OA_DCN_Diseased.gml',format = 'gml')
write.graph(graph = OA_DCN$Healthy,file = 'OA_DCN_Healthy.gml',format = 'gml')

OA_AllAnalyses=list()
OA_AllAnalyses[[1]]=select(x = org.Hs.eg.db,keys = oa_deg_df$genes[which(abs(oa_deg_df$logFC)>log(2)&oa_deg_df$FDR<=0.05)])[,2]
OA_AllAnalyses[[1]]=


pos_weight=seq(1.25,1.75,by=0.05)
neg_weight=seq(0.25,0.75,by=0.05)
poscoexp_pval_vector=vector(mode = 'list',length=length(pos_weight))
negcoexp_pval_vector=vector(mode = 'list',length=length(pos_weight))
names(poscoexp_pval_vector)=pos_weight
names(negcoexp_pval_vector)=neg_weight
for(i in 1:length(pos_weight)){
  poscoexp_pval_vector[[i]]=vector(mode = 'list',length=8)
  names(poscoexp_pval_vector[[i]])=c('Degree','Closeness','PageRank','GO.BP_SemSim','GO.CC_SemSim','GO.MF_SemSim','Unique_HealthyPos_Genes','Unique_DiseasePos_Genes')
  negcoexp_pval_vector[[i]]=vector(mode = 'list',length=8)
  names(negcoexp_pval_vector[[i]])=c('Degree','Closeness','PageRank','GO.BP_SemSim','GO.CC_SemSim','GO.MF_SemSim','Unique_HealthyNeg_Genes','Unique_DiseaseNeg_Genes')
  
}

poscoexp_degree=negcoexp_degree=poscoexp_pgr=negcoexp_pgr=poscoexp_cls=negcoexp_cls=poscoexp_GO=negcoexp_GO=vector(mode = 'list',length = 2)
names(poscoexp_degree)=names(negcoexp_degree)=names(poscoexp_pgr)=names(negcoexp_pgr)=names(poscoexp_cls)=names(negcoexp_cls)=names(poscoexp_GO)=names(negcoexp_GO)=c('OA_Healthy','OA_Disease')
poscoexp_GO$OA_Healthy=poscoexp_GO$OA_Disease=negcoexp_GO$OA_Healthy=negcoexp_GO$OA_Disease=vector(mode = 'list',length=3)
names(poscoexp_GO$OA_Healthy)=names(poscoexp_GO$OA_Disease)=names(negcoexp_GO$OA_Healthy)=names(negcoexp_GO$OA_Disease)=c('BP','CC','MF')

hsGO=vector(mode = "list",length = 3)
names(hsGO)=c('BP','CC','MF')
hsGO$BP=godata(OrgDb='org.Hs.eg.db',keytype='SYMBOL',ont='BP')
hsGO$CC=godata(OrgDb='org.Hs.eg.db',keytype='SYMBOL',ont='CC')
hsGO$MF=godata(OrgDb='org.Hs.eg.db',keytype='SYMBOL',ont='MF')

for (w in 1:length(pos_weight)){
    
      healthypos<-subset(OA_analysis, subset=(FDR < 0.1 & r.c > pos_weight[w]-1))
      diseasepos<-subset(OA_analysis, subset=(FDR < 0.1 & r.t > pos_weight[w]-1))
      healthyneg<-subset(OA_analysis, subset=(FDR < 0.1 & r.c > neg_weight[w]-1))
      diseaseneg<-subset(OA_analysis, subset=(FDR < 0.1 & r.t > neg_weight[w]-1))
      healthypos_genes<-union(healthypos$Gene.A, healthypos$Gene.B)
      diseasepos_genes<-union(diseasepos$Gene.A, diseasepos$Gene.B)
      healthyneg_genes<-union(healthyneg$Gene.A, healthyneg$Gene.B)
      diseaseneg_genes<-union(diseaseneg$Gene.A, diseaseneg$Gene.B)
      
      poscoexp_degree$OA_Healthy=degree(graph.data.frame(d=data.frame(ends(graph = OA_DCN$Healthy,es = which(E(OA_DCN$Healthy)$weight<pos_weight[w]),names = T)),directed = F))
      poscoexp_degree$OA_Disease=degree(graph.data.frame(d=data.frame(ends(graph = OA_DCN$Diseased,es = which(E(OA_DCN$Diseased)$weight>pos_weight[w]),names = T)),directed = F))
      
      negcoexp_degree$OA_Healthy=degree(graph.data.frame(d=data.frame(ends(graph = OA_DCN$Healthy,es = which(E(OA_DCN$Healthy)$weight<neg_weight[w]),names = T)),directed = F))
      negcoexp_degree$OA_Disease=degree(graph.data.frame(d=data.frame(ends(graph = OA_DCN$Diseased,es = which(E(OA_DCN$Diseased)$weight>neg_weight[w]),names = T)),directed = F))
      
      poscoexp_cls$OA_Healthy=closeness(graph.data.frame(d=data.frame(ends(graph = OA_DCN$Healthy,es = which(E(OA_DCN$Healthy)$weight<pos_weight[w]),names = T)),directed = F))
      poscoexp_cls$OA_Disease=closeness(graph.data.frame(d=data.frame(ends(graph = OA_DCN$Diseased,es = which(E(OA_DCN$Diseased)$weight>pos_weight[w]),names = T)),directed = F))
      negcoexp_cls$OA_Healthy=closeness(graph.data.frame(d=data.frame(ends(graph = OA_DCN$Healthy,es = which(E(OA_DCN$Healthy)$weight<neg_weight[w]),names = T)),directed = F))
      negcoexp_cls$OA_Disease=closeness(graph.data.frame(d=data.frame(ends(graph = OA_DCN$Diseased,es = which(E(OA_DCN$Diseased)$weight>neg_weight[w]),names = T)),directed = F))
      
      poscoexp_pgr$OA_Healthy=page_rank(graph.data.frame(d=data.frame(ends(graph = OA_DCN$Healthy,es = which(E(OA_DCN$Healthy)$weight<pos_weight[w]),names = T)),directed = F))$vector
      poscoexp_pgr$OA_Disease=page_rank(graph.data.frame(d=data.frame(ends(graph = OA_DCN$Diseased,es = which(E(OA_DCN$Diseased)$weight>pos_weight[w]),names = T)),directed = F))$vector
      negcoexp_pgr$OA_Healthy=page_rank(graph.data.frame(d=data.frame(ends(graph = OA_DCN$Healthy,es = which(E(OA_DCN$Healthy)$weight<neg_weight[w]),names = T)),directed = F))$vector
      negcoexp_pgr$OA_Disease=page_rank(graph.data.frame(d=data.frame(ends(graph = OA_DCN$Diseased,es = which(E(OA_DCN$Diseased)$weight>neg_weight[w]),names = T)),directed = F))$vector
      
      # poscoexp_GO$OA_Healthy$BP=mgeneSim(genes = V(graph.data.frame(d=data.frame(ends(graph = OA_DCN$Healthy,es = which(E(OA_DCN$Healthy)$weight<pos_weight[w]),names = T)),directed = F))$name,semData = hsGO$BP,measure = 'Wang',drop = 'IEA',combine = 'BMA',verbose = T)
      # poscoexp_GO$OA_Disease$BP=mgeneSim(genes = V(graph.data.frame(d=data.frame(ends(graph = OA_DCN$Diseased,es = which(E(OA_DCN$Diseased)$weight>pos_weight[w]),names = T)),directed = F))$name,semData = hsGO$BP,measure = 'Wang',drop = 'IEA',combine = 'BMA',verbose = T)
      # negcoexp_GO$OA_Healthy$BP=mgeneSim(genes = V(graph.data.frame(d=data.frame(ends(graph = OA_DCN$Healthy,es = which(E(OA_DCN$Healthy)$weight<neg_weight[w]),names = T)),directed = F))$name,semData = hsGO$BP,measure = 'Wang',drop = 'IEA',combine = 'BMA',verbose = T)
      # negcoexp_GO$OA_Disease$BP=mgeneSim(genes = V(graph.data.frame(d=data.frame(ends(graph = OA_DCN$Diseased,es = which(E(OA_DCN$Diseased)$weight>neg_weight[w]),names = T)),directed = F))$name,semData = hsGO$BP,measure = 'Wang',drop = 'IEA',combine = 'BMA',verbose = T)
      # 
      # poscoexp_GO$OA_Healthy$CC=mgeneSim(genes = V(graph.data.frame(d=data.frame(ends(graph = OA_DCN$Healthy,es = which(E(OA_DCN$Healthy)$weight<pos_weight[w]),names = T)),directed = F))$name,semData = hsGO$CC,measure = 'Wang',drop = 'IEA',combine = 'BMA',verbose = T)
      # poscoexp_GO$OA_Disease$CC=mgeneSim(genes = V(graph.data.frame(d=data.frame(ends(graph = OA_DCN$Diseased,es = which(E(OA_DCN$Diseased)$weight>pos_weight[w]),names = T)),directed = F))$name,semData = hsGO$CC,measure = 'Wang',drop = 'IEA',combine = 'BMA',verbose = T)
      # negcoexp_GO$OA_Healthy$CC=mgeneSim(genes = V(graph.data.frame(d=data.frame(ends(graph = OA_DCN$Healthy,es = which(E(OA_DCN$Healthy)$weight<neg_weight[w]),names = T)),directed = F))$name,semData = hsGO$CC,measure = 'Wang',drop = 'IEA',combine = 'BMA',verbose = T)
      # negcoexp_GO$OA_Disease$CC=mgeneSim(genes = V(graph.data.frame(d=data.frame(ends(graph = OA_DCN$Diseased,es = which(E(OA_DCN$Diseased)$weight>neg_weight[w]),names = T)),directed = F))$name,semData = hsGO$CC,measure = 'Wang',drop = 'IEA',combine = 'BMA',verbose = T)
      # 
      # poscoexp_GO$OA_Healthy$MF=mgeneSim(genes = V(graph.data.frame(d=data.frame(ends(graph = OA_DCN$Healthy,es = which(E(OA_DCN$Healthy)$weight<pos_weight[w]),names = T)),directed = F))$name,semData = hsGO$MF,measure = 'Wang',drop = 'IEA',combine = 'BMA',verbose = T)
      # poscoexp_GO$OA_Disease$MF=mgeneSim(genes = V(graph.data.frame(d=data.frame(ends(graph = OA_DCN$Diseased,es = which(E(OA_DCN$Diseased)$weight>pos_weight[w]),names = T)),directed = F))$name,semData = hsGO$MF,measure = 'Wang',drop = 'IEA',combine = 'BMA',verbose = T)
      # negcoexp_GO$OA_Healthy$MF=mgeneSim(genes = V(graph.data.frame(d=data.frame(ends(graph = OA_DCN$Healthy,es = which(E(OA_DCN$Healthy)$weight<neg_weight[w]),names = T)),directed = F))$name,semData = hsGO$MF,measure = 'Wang',drop = 'IEA',combine = 'BMA',verbose = T)
      # negcoexp_GO$OA_Disease$MF=mgeneSim(genes = V(graph.data.frame(d=data.frame(ends(graph = OA_DCN$Diseased,es = which(E(OA_DCN$Diseased)$weight>neg_weight[w]),names = T)),directed = F))$name,semData = hsGO$MF,measure = 'Wang',drop = 'IEA',combine = 'BMA',verbose = T)
      # 
      poscoexp_pval_vector[[w]][[1]]=-log10(wilcox.test(x=poscoexp_degree$OA_Healthy,y=poscoexp_degree$OA_Disease)$p.value)
      negcoexp_pval_vector[[w]][[1]]=-log10(wilcox.test(x=negcoexp_degree$OA_Healthy,y=negcoexp_degree$OA_Disease)$p.value)
      poscoexp_pval_vector[[w]][[2]]=-log10(wilcox.test(x=poscoexp_pgr$OA_Healthy,y=poscoexp_pgr$OA_Disease)$p.value)
      negcoexp_pval_vector[[w]][[2]]=-log10(wilcox.test(x=negcoexp_pgr$OA_Healthy,y=negcoexp_pgr$OA_Disease)$p.value)
      poscoexp_pval_vector[[w]][[3]]=-log10(wilcox.test(x=poscoexp_cls$OA_Healthy,y=poscoexp_cls$OA_Disease)$p.value)
      negcoexp_pval_vector[[w]][[3]]=-log10(wilcox.test(x=negcoexp_cls$OA_Healthy,y=negcoexp_cls$OA_Disease)$p.value)
      # poscoexp_pval_vector[[w]][[4]]=-log10(wilcox.test(x=poscoexp_GO$OA_Healthy$BP,y=poscoexp_GO$OA_Disease$BP)$p.value)
      # negcoexp_pval_vector[[w]][[4]]=-log10(wilcox.test(x=negcoexp_GO$OA_Healthy$BP,y=negcoexp_GO$OA_Disease$BP)$p.value)
      # poscoexp_pval_vector[[w]][[5]]=-log10(wilcox.test(x=poscoexp_GO$OA_Healthy$CC,y=poscoexp_GO$OA_Disease$CC)$p.value)
      # negcoexp_pval_vector[[w]][[5]]=-log10(wilcox.test(x=negcoexp_GO$OA_Healthy$CC,y=negcoexp_GO$OA_Disease$MF)$p.value)
      # poscoexp_pval_vector[[w]][[6]]=-log10(wilcox.test(x=poscoexp_GO$OA_Healthy$MF,y=poscoexp_GO$OA_Disease$MF)$p.value)
      # negcoexp_pval_vector[[w]][[6]]=-log10(wilcox.test(x=negcoexp_GO$OA_Healthy$MF,y=negcoexp_GO$OA_Disease$MF)$p.value)
      poscoexp_pval_vector[[w]][[7]]=length(degree(graph = OA_DCN$Healthy,v = V(OA_DCN$Healthy)$name%in%setdiff(healthypos_genes, diseasepos_genes)))
      poscoexp_pval_vector[[w]][[8]]=length(degree(graph = OA_DCN$Diseased,v = V(OA_DCN$Diseased)$name%in%setdiff(diseasepos_genes,healthypos_genes)))
      negcoexp_pval_vector[[w]][[7]]=length(degree(graph = OA_DCN$Healthy,v = V(OA_DCN$Healthy)$name%in%setdiff(healthyneg_genes, diseaseneg_genes)))
      negcoexp_pval_vector[[w]][[8]]=length(degree(graph = OA_DCN$Diseased,v = V(OA_DCN$Diseased)$name%in%setdiff(diseaseneg_genes, healthyneg_genes)))
}

