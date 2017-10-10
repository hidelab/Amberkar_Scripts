library(igraph)
library(gdata)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
source("http://peterhaschke.com/Code/multiplot.R")
source("https://raw.githubusercontent.com/SamBuckberry/RUN-WGCNA/master/WGCNA_functions.R")


setwd('/Users/sandeepamberkar/Work/Data/MSMM-MSBB-HBTRC-Data/MSMM_RNAseq/MSMM_RNAseq_FinalRelease2/')
#Coexpression matrices for PLQ associated genes
msmm_rnaseq_PLQgenes_coexp_adjMatrix=readRDS('msmm_rnaseq_PLQgenes_coexp_adjMatrix.RDS')
msmm_rnaseq_PLQgenes.CoexpNet=lapply(msmm_rnaseq_PLQgenes_coexp_adjMatrix,graph_from_adjacency_matrix,mode = 'undirected',diag = F,weighted = T)
msmm_rnaseq_PLQgenes.CoexpNet_edgeWeight=lapply(msmm_rnaseq_PLQgenes.CoexpNet,function(x)E(x)$weight<-E(x)$weight+1)
for (i in 1:4){
  E(msmm_rnaseq_PLQgenes.CoexpNet[[i]])$weight=msmm_rnaseq_PLQgenes.CoexpNet_edgeWeight[[i]]
}
#Read edited variants
tanzi_damaging_snps=as.data.frame(read.xls('../../../../Collaborations/Tanzi_WGS/Damaging_exonic_SNPs.xlsx',sheet = 1),stringsAsFactors = F)
tanzi_protective_snps=as.data.frame(read.xls('../../../../Collaborations/Tanzi_WGS/Protective_exonic_SNPs.xlsx',sheet = 1),stringsAsFactors = F)
#Read curated PPI table
hs_ppi_table=fread('/Users/sandeepamberkar/Work/Data/PPI-Data/iref14_Human_UP_noDup_table.txt',sep = '\t',header = T,data.table = T,showProgress = T)
#Remove isoforms and collapse to gene level
hs_ppi.Isoform_A=grep(pattern = '-',x = hs_ppi_table$V1,value = F)
names(hs_ppi.Isoform_A)=grep(pattern = '-',x = hs_ppi_table$V1,value = T)
hs_ppi.Isoform_B=grep(pattern = '-',x = hs_ppi_table$V2,value = F)
names(hs_ppi.Isoform_B)=grep(pattern = '-',x = hs_ppi_table$V2,value = T)
hs_ppi_table$V1[hs_ppi.Isoform_A]=gsub(pattern = '-\\d+',replacement='',x = grep(pattern = '-',x = names(hs_ppi.Isoform_A),value = T))
hs_ppi_table$V2[hs_ppi.Isoform_B]=gsub(pattern = '-\\d+',replacement='',x = grep(pattern = '-',x = names(hs_ppi.Isoform_B),value = T))
hs_ppi_table$Gene.A=unname(mapIds(x = org.Hs.eg.db,keys = hs_ppi_table$V1,column = 'SYMBOL',keytype = 'UNIPROT'))
hs_ppi_table$Gene.B=unname(mapIds(x = org.Hs.eg.db,keys = hs_ppi_table$V2,column = 'SYMBOL',keytype = 'UNIPROT'))
hs_ppi=graph.data.frame(d = hs_ppi_table[which((is.na(hs_ppi_table$Gene.A)==F)&(is.na(hs_ppi_table$Gene.B)==F)),c(4:5)],directed = F)
#Build datastructure to store subnetworks and centralities
hs_ppi.subnet=hs_ppi_subnet.Centralities=vector(mode = 'list',length = 2)
names(hs_ppi.subnet)=names(hs_ppi_subnet.Centralities)=c('TanziProtective','TanziDamaging')
hs_ppi.subnet$TanziProtective=induced.subgraph(graph = hs_ppi,vids = unique(tanzi_protective_snps$Gene))
hs_ppi.subnet$TanziDamaging=induced.subgraph(graph = hs_ppi,vids = unique(tanzi_damaging_snps$Gene))
hs_ppi_subnet.Centralities$TanziProtective=hs_ppi_subnet.Centralities$TanziDamaging=vector(mode = 'list',length = 3)
names(hs_ppi_subnet.Centralities$TanziProtective)=names(hs_ppi_subnet.Centralities$TanziDamaging)=c('Degree','Betweenness','Short.Paths')
hs_ppi_subnet.Centralities$TanziProtective$Degree=degree(graph = hs_ppi.subnet$TanziProtective)
hs_ppi_subnet.Centralities$TanziProtective$Betweenness=betweenness(graph = hs_ppi.subnet$TanziProtective)
#Remove nodes which have Inf pathlength
tmp_sp1=upperTriangle(shortest.paths(graph = hs_ppi.subnet$TanziProtective,v = unique(tanzi_protective_snps$Gene)))
hs_ppi_subnet.Centralities$TanziProtective$Short.Paths=tmp_sp1[which(tmp_sp1!='Inf')]
hs_ppi_subnet.Centralities$TanziDamaging$Degree=degree(graph = hs_ppi.subnet$TanziDamaging)
hs_ppi_subnet.Centralities$TanziDamaging$Betweenness=betweenness(graph = hs_ppi.subnet$TanziDamaging)
tmp_sp2=upperTriangle(shortest.paths(graph = hs_ppi.subnet$TanziDamaging,v = unique(tanzi_damaging_snps$Gene)))
hs_ppi_subnet.Centralities$TanziDamaging$Short.Paths=tmp_sp2[which(tmp_sp2!='Inf')]

#Read cell type specific markers
zhang_celltype_ADgenes=read.xls('../../../BrainExpression_Datasets/Zhang_19BrainRegions_Paper/Zhang_BrainCelltype_Markers.xlsx',skip=1,sheet=3,header=T,as.is=T)
zhang_celltype_ADgenes.list=zhang_celltype_PLQ_ADgenes.list=vector(mode = 'list',length = 5)
names(zhang_celltype_ADgenes.list)=names(zhang_celltype_PLQ_ADgenes.list)=sort(unique(zhang_celltype_ADgenes$Cell.type))
zhang_celltype_ADgenes.list$Astrocytes=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_ADgenes.list)[1],zhang_celltype_ADgenes$Cell.type)]
zhang_celltype_ADgenes.list$Endothelial=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_ADgenes.list)[2],zhang_celltype_ADgenes$Cell.type)]
zhang_celltype_ADgenes.list$Microglia=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_ADgenes.list)[3],zhang_celltype_ADgenes$Cell.type)]
zhang_celltype_ADgenes.list$Neurons=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_ADgenes.list)[4],zhang_celltype_ADgenes$Cell.type)]
zhang_celltype_ADgenes.list$Oligodendrocytes=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_ADgenes.list)[5],zhang_celltype_ADgenes$Cell.type)]

#Build datastructure to store PLQ-genes DCN
msbb_rnaseq2016_PLQ.DCN=msbb_rnaseq2016_PLQ_TanziProtective.DCN=msbb_rnaseq2016_PLQ_TanziDamaging.DCN=vector(mode = "list",length = 4)
hs_ppi_TanziProtective.subnet.Degree=hs_ppi_TanziProtective.subnet.Betweenness=hs_ppi_TanziProtective.subnet.ShortestPaths=vector(mode = "list",length = 4)
hs_ppi_TanziDamaging.subnet.Degree=hs_ppi_TanziDamaging.subnet.Betweenness=hs_ppi_TanziDamaging.subnet.ShortestPaths=vector(mode = "list",length = 4)
msbb_rnaseq2016_PLQ_DCN.Degree=msbb_rnaseq2016_PLQ_DCN.Betweenness=msbb_rnaseq2016_PLQ_DCN.ShortestPaths=vector(mode = "list",length = 4)
msbb_rnaseq2016_PLQ_TanziProtective_DCN.Degree=msbb_rnaseq2016_PLQ_TanziProtective_DCN.Betweenness=msbb_rnaseq2016_PLQ_TanziProtective_DCN.ShortestPaths=vector(mode = "list",length = 4)
msbb_rnaseq2016_PLQ_TanziDamaging_DCN.Degree=msbb_rnaseq2016_PLQ_TanziDamaging_DCN.Betweenness=msbb_rnaseq2016_PLQ_TanziDamaging_DCN.ShortestPaths=vector(mode = "list",length = 4)
names(msbb_rnaseq2016_PLQ.DCN)=names(msbb_rnaseq2016_PLQ_TanziProtective.DCN)=names(msbb_rnaseq2016_PLQ_TanziDamaging.DCN)=c("FP","IFG","PHG","STG")
names(msbb_rnaseq2016_PLQ_DCN.Degree)=names(msbb_rnaseq2016_PLQ_DCN.Betweenness)=names(msbb_rnaseq2016_PLQ_DCN.ShortestPaths)=c("FP","IFG","PHG","STG")
names(msbb_rnaseq2016_PLQ_TanziProtective_DCN.Degree)=names(msbb_rnaseq2016_PLQ_TanziProtective_DCN.Betweenness)=names(msbb_rnaseq2016_PLQ_TanziProtective_DCN.ShortestPaths)=c("FP","IFG","PHG","STG")
names(msbb_rnaseq2016_PLQ_TanziDamaging_DCN.Degree)=names(msbb_rnaseq2016_PLQ_TanziDamaging_DCN.Betweenness)=names(msbb_rnaseq2016_PLQ_TanziDamaging_DCN.ShortestPaths)=c("FP","IFG","PHG","STG")

msbb_rnaseq2016_PLQ.DCN$FP=msbb_rnaseq2016_PLQ.DCN$IFG=msbb_rnaseq2016_PLQ.DCN$PHG=msbb_rnaseq2016_PLQ.DCN$STG=vector(mode="list",length=2)
msbb_rnaseq2016_PLQ_DCN.Degree$FP=msbb_rnaseq2016_PLQ_DCN.Degree$IFG=msbb_rnaseq2016_PLQ_DCN.Degree$PHG=msbb_rnaseq2016_PLQ_DCN.Degree$STG=vector(mode="list",length=2)
msbb_rnaseq2016_PLQ_DCN.Betweenness$FP=msbb_rnaseq2016_PLQ_DCN.Betweenness$IFG=msbb_rnaseq2016_PLQ_DCN.Betweenness$PHG=msbb_rnaseq2016_PLQ_DCN.Betweenness$STG=vector(mode="list",length=2)
msbb_rnaseq2016_PLQ_DCN.ShortestPaths$FP=msbb_rnaseq2016_PLQ_DCN.ShortestPaths$IFG=msbb_rnaseq2016_PLQ_DCN.ShortestPaths$PHG=msbb_rnaseq2016_PLQ_DCN.ShortestPaths$STG=vector(mode="list",length=2)
msbb_rnaseq2016_PLQ_TanziProtective.DCN$FP=msbb_rnaseq2016_PLQ_TanziProtective.DCN$IFG=msbb_rnaseq2016_PLQ_TanziProtective.DCN$PHG=msbb_rnaseq2016_PLQ_TanziProtective.DCN$STG=vector(mode="list",length=2)
msbb_rnaseq2016_PLQ_TanziProtective_DCN.Degree$FP=msbb_rnaseq2016_PLQ_TanziProtective_DCN.Degree$IFG=msbb_rnaseq2016_PLQ_TanziProtective_DCN.Degree$PHG=msbb_rnaseq2016_PLQ_TanziProtective_DCN.Degree$STG=vector(mode="list",length=2)
msbb_rnaseq2016_PLQ_TanziProtective_DCN.Betweenness$FP=msbb_rnaseq2016_PLQ_TanziProtective_DCN.Betweenness$IFG=msbb_rnaseq2016_PLQ_TanziProtective_DCN.Betweenness$PHG=msbb_rnaseq2016_PLQ_TanziProtective_DCN.Betweenness$STG=vector(mode="list",length=2)
msbb_rnaseq2016_PLQ_TanziProtective_DCN.ShortestPaths$FP=msbb_rnaseq2016_PLQ_TanziProtective_DCN.ShortestPaths$IFG=msbb_rnaseq2016_PLQ_TanziProtective_DCN.ShortestPaths$PHG=msbb_rnaseq2016_PLQ_TanziProtective_DCN.ShortestPaths$STG=vector(mode="list",length=2)
msbb_rnaseq2016_PLQ_TanziDamaging.DCN$FP=msbb_rnaseq2016_PLQ_TanziDamaging.DCN$IFG=msbb_rnaseq2016_PLQ_TanziDamaging.DCN$PHG=msbb_rnaseq2016_PLQ_TanziDamaging.DCN$STG=vector(mode="list",length=2)
msbb_rnaseq2016_PLQ_TanziDamaging_DCN.Degree$FP=msbb_rnaseq2016_PLQ_TanziDamaging_DCN.Degree$IFG=msbb_rnaseq2016_PLQ_TanziDamaging_DCN.Degree$PHG=msbb_rnaseq2016_PLQ_TanziDamaging_DCN.Degree$STG=vector(mode="list",length=2)
msbb_rnaseq2016_PLQ_TanziDamaging_DCN.Betweenness$FP=msbb_rnaseq2016_PLQ_TanziDamaging_DCN.Betweenness$IFG=msbb_rnaseq2016_PLQ_TanziDamaging_DCN.Betweenness$PHG=msbb_rnaseq2016_PLQ_TanziDamaging_DCN.Betweenness$STG=vector(mode="list",length=2)
msbb_rnaseq2016_PLQ_TanziDamaging_DCN.ShortestPaths$FP=msbb_rnaseq2016_PLQ_TanziDamaging_DCN.ShortestPaths$IFG=msbb_rnaseq2016_PLQ_TanziDamaging_DCN.ShortestPaths$PHG=msbb_rnaseq2016_PLQ_TanziDamaging_DCN.ShortestPaths$STG=vector(mode="list",length=2)

names(msbb_rnaseq2016_PLQ.DCN$FP)=names(msbb_rnaseq2016_PLQ.DCN$IFG)=names(msbb_rnaseq2016_PLQ.DCN$PHG)=names(msbb_rnaseq2016_PLQ.DCN$STG)=c("Low","High")
names(msbb_rnaseq2016_PLQ_DCN.Degree$FP)=names(msbb_rnaseq2016_PLQ_DCN.Degree$IFG)=names(msbb_rnaseq2016_PLQ_DCN.Degree$PHG)=names(msbb_rnaseq2016_PLQ_DCN.Degree$STG)=c("Low","High")
names(msbb_rnaseq2016_PLQ_DCN.Betweenness$FP)=names(msbb_rnaseq2016_PLQ_DCN.Betweenness$IFG)=names(msbb_rnaseq2016_PLQ_DCN.Betweenness$PHG)=names(msbb_rnaseq2016_PLQ_DCN.Betweenness$STG)=c("Low","High")
names(msbb_rnaseq2016_PLQ_DCN.ShortestPaths$FP)=names(msbb_rnaseq2016_PLQ_DCN.ShortestPaths$IFG)=names(msbb_rnaseq2016_PLQ_DCN.ShortestPaths$PHG)=names(msbb_rnaseq2016_PLQ_DCN.ShortestPaths$STG)=c("Low","High")
names(msbb_rnaseq2016_PLQ_TanziProtective.DCN$FP)=names(msbb_rnaseq2016_PLQ_TanziProtective.DCN$IFG)=names(msbb_rnaseq2016_PLQ_TanziProtective.DCN$PHG)=names(msbb_rnaseq2016_PLQ_TanziProtective.DCN$STG)=c("Low","High")
names(msbb_rnaseq2016_PLQ_TanziDamaging.DCN$FP)=names(msbb_rnaseq2016_PLQ_TanziDamaging.DCN$IFG)=names(msbb_rnaseq2016_PLQ_TanziDamaging.DCN$PHG)=names(msbb_rnaseq2016_PLQ_TanziDamaging.DCN$STG)=c("Low","High")
names(msbb_rnaseq2016_PLQ_TanziDamaging_DCN.Degree$FP)=names(msbb_rnaseq2016_PLQ_TanziDamaging_DCN.Degree$IFG)=names(msbb_rnaseq2016_PLQ_TanziDamaging_DCN.Degree$PHG)=names(msbb_rnaseq2016_PLQ_TanziDamaging_DCN.Degree$STG)=c("Low","High")
names(msbb_rnaseq2016_PLQ_TanziDamaging_DCN.Betweenness$FP)=names(msbb_rnaseq2016_PLQ_TanziDamaging_DCN.Betweenness$IFG)=names(msbb_rnaseq2016_PLQ_TanziDamaging_DCN.Betweenness$PHG)=names(msbb_rnaseq2016_PLQ_TanziDamaging_DCN.Betweenness$STG)=c("Low","High")
names(msbb_rnaseq2016_PLQ_TanziDamaging_DCN.ShortestPaths$FP)=names(msbb_rnaseq2016_PLQ_TanziDamaging_DCN.ShortestPaths$IFG)=names(msbb_rnaseq2016_PLQ_TanziDamaging_DCN.ShortestPaths$PHG)=names(msbb_rnaseq2016_PLQ_TanziDamaging_DCN.ShortestPaths$STG)=c("Low","High")
names(msbb_rnaseq2016_PLQ_TanziProtective_DCN.Degree$FP)=names(msbb_rnaseq2016_PLQ_TanziProtective_DCN.Degree$IFG)=names(msbb_rnaseq2016_PLQ_TanziProtective_DCN.Degree$PHG)=names(msbb_rnaseq2016_PLQ_TanziProtective_DCN.Degree$STG)=c("Low","High")
names(msbb_rnaseq2016_PLQ_TanziProtective_DCN.Betweenness$FP)=names(msbb_rnaseq2016_PLQ_TanziProtective_DCN.Betweenness$IFG)=names(msbb_rnaseq2016_PLQ_TanziProtective_DCN.Betweenness$PHG)=names(msbb_rnaseq2016_PLQ_TanziProtective_DCN.Betweenness$STG)=c("Low","High")
names(msbb_rnaseq2016_PLQ_TanziProtective_DCN.ShortestPaths$FP)=names(msbb_rnaseq2016_PLQ_TanziProtective_DCN.ShortestPaths$IFG)=names(msbb_rnaseq2016_PLQ_TanziProtective_DCN.ShortestPaths$PHG)=names(msbb_rnaseq2016_PLQ_TanziProtective_DCN.ShortestPaths$STG)=c("Low","High")

#Read DCNs
msbb_rnaseq2016_PLQ.DCN$FP=readRDS("FP/PLQ_DCN/FP_PLQGenes_FDR01.DCN.RDS")
msbb_rnaseq2016_PLQ.DCN$IFG=readRDS("IFG/PLQ_DCN/IFG_PLQGenes_FDR01.DCN.RDS")
msbb_rnaseq2016_PLQ.DCN$PHG=readRDS("PHG/PLQ_DCN/PHG_PLQGenes_FDR01.DCN.RDS")
msbb_rnaseq2016_PLQ.DCN$PHG=msbb_rnaseq2016_PLQ.DCN$PHG[-c(1:2)]
msbb_rnaseq2016_PLQ.DCN$STG=readRDS("STG/PLQ_DCN/STG_PLQGenes_FDR01.DCN.RDS")

#Compute network centralities for all subnetworks
for (t in 1:4){
  msbb_rnaseq2016_PLQ_TanziProtective.DCN[[t]]=lapply(msbb_rnaseq2016_PLQ.DCN[[t]],induced_subgraph,vids = unique(tanzi_protective_snps$Gene))
  msbb_rnaseq2016_PLQ_TanziDamaging.DCN[[t]]=lapply(msbb_rnaseq2016_PLQ.DCN[[t]],induced_subgraph,vids = unique(tanzi_damaging_snps$Gene))
  msbb_rnaseq2016_PLQ_DCN.Degree[[t]]=lapply(msbb_rnaseq2016_PLQ.DCN[[t]],degree)
  msbb_rnaseq2016_PLQ_DCN.Betweenness[[t]]=lapply(msbb_rnaseq2016_PLQ.DCN[[t]],function(x)betweenness(graph = x,weights = E(x)$weight))
  msbb_rnaseq2016_PLQ_DCN.ShortestPaths[[t]]=lapply(lapply(msbb_rnaseq2016_PLQ.DCN$STG,function(x)upperTriangle(shortest.paths(graph = x,weights = E(x)$weight))),function(y)y[which(y!='Inf')])
  
  msbb_rnaseq2016_PLQ_TanziProtective_DCN.Degree[[t]]=lapply(msbb_rnaseq2016_PLQ_TanziProtective.DCN[[t]],degree)
  msbb_rnaseq2016_PLQ_TanziProtective_DCN.Betweenness[[t]]=lapply(msbb_rnaseq2016_PLQ_TanziProtective.DCN[[t]],function(x)betweenness(graph = x,weights = E(x)$weight))
  msbb_rnaseq2016_PLQ_TanziProtective_DCN.ShortestPaths[[t]]=lapply(lapply(msbb_rnaseq2016_PLQ_TanziProtective.DCN[[t]],function(x)upperTriangle(shortest.paths(graph = x,weights = E(x)$weight))),function(y)y[which(y!='Inf')])
  
  msbb_rnaseq2016_PLQ_TanziDamaging_DCN.Degree[[t]]=lapply(msbb_rnaseq2016_PLQ_TanziDamaging.DCN[[t]],degree)
  msbb_rnaseq2016_PLQ_TanziDamaging_DCN.Betweenness[[t]]=lapply(msbb_rnaseq2016_PLQ_TanziDamaging.DCN[[t]],function(x)betweenness(graph = x,weights = E(x)$weight))
  msbb_rnaseq2016_PLQ_TanziDamaging_DCN.ShortestPaths[[t]]=lapply(lapply(msbb_rnaseq2016_PLQ_TanziDamaging.DCN[[t]],function(x)upperTriangle(shortest.paths(graph = x,weights = E(x)$weight))),function(y)y[which(y!='Inf')])
  #Compute shortes path for Low PLQ DCNs
  sp_df_low_plq_dcn=data.frame(sp=log10(msbb_rnaseq2016_PLQ_DCN.ShortestPaths[[t]]$Low),stringsAsFactors = F)
  sp_df_low_protectiveSNPs=data.frame(sp=log10(msbb_rnaseq2016_PLQ_TanziProtective_DCN.ShortestPaths[[t]]$Low),stringsAsFactors = F)
  sp_df_low_damagingSNPs=data.frame(sp=log10(msbb_rnaseq2016_PLQ_TanziDamaging_DCN.ShortestPaths[[t]]$Low),stringsAsFactors = F)
  sp_df_hs_ppi_protectiveSNPs=data.frame(sp=log10(hs_ppi_subnet.Centralities$TanziProtective$Short.Paths),stringsAsFactors = F)
  sp_df_hs_ppi_damagingSNPs=data.frame(sp=log10(hs_ppi_subnet.Centralities$TanziDamaging$Short.Paths),stringsAsFactors = F)
  sp_df_low_plq_dcn$type=paste(names(msbb_rnaseq2016_PLQ_DCN.ShortestPaths)[t],'LowPLQ-genes.DCN',sep = '_')
  sp_df_low_protectiveSNPs$type=paste(names(msbb_rnaseq2016_PLQ_DCN.ShortestPaths)[t],'LowPLQ','w/Protective_SNPs',sep = '_')
  sp_df_low_damagingSNPs$type=paste(names(msbb_rnaseq2016_PLQ_DCN.ShortestPaths)[t],'LowPLQ','w/Damaging_SNPs',sep = '_')
  sp_df_hs_ppi_protectiveSNPs$type='hs_ppi_w/Protective_SNPs'
  sp_df_hs_ppi_damagingSNPs$type='hs_ppi_w/Damaging_SNPs'
  sp_df_low=rbind(sp_df_low_plq_dcn,sp_df_low_damagingSNPs,sp_df_low_protectiveSNPs,sp_df_hs_ppi_protectiveSNPs,sp_df_hs_ppi_damagingSNPs)
  g1=ggplot(sp_df_low,aes(sp,fill=type))+geom_density(alpha=0.5)+ggtitle(label = 'ShortestPath distribution',subtitle = paste(names(msbb_rnaseq2016_PLQ.DCN)[t],'LowPLQ_DCN w/ PLQ-assoc.genes,Tanzi variants',sep = '_'))
  #Compute shortest path separately for High PLQ DCNs
  sp_df_high_plq_dcn=data.frame(sp=log10(msbb_rnaseq2016_PLQ_DCN.ShortestPaths[[t]]$High),stringsAsFactors = F)
  sp_df_high_protectiveSNPs=data.frame(sp=log10(msbb_rnaseq2016_PLQ_TanziProtective_DCN.ShortestPaths[[t]]$High),stringsAsFactors = F)
  sp_df_high_damagingSNPs=data.frame(sp=log10(msbb_rnaseq2016_PLQ_TanziDamaging_DCN.ShortestPaths[[t]]$High),stringsAsFactors = F)
  sp_df_hs_ppi_protectiveSNPs=data.frame(sp=log10(hs_ppi_subnet.Centralities$TanziProtective$Short.Paths),stringsAsFactors = F)
  sp_df_hs_ppi_damagingSNPs=data.frame(sp=log10(hs_ppi_subnet.Centralities$TanziDamaging$Short.Paths),stringsAsFactors = F)
  sp_df_high_plq_dcn$type=paste(names(msbb_rnaseq2016_PLQ_DCN.ShortestPaths)[t],'HighPLQ-genes.DCN',sep = '_')
  sp_df_high_protectiveSNPs$type=paste(names(msbb_rnaseq2016_PLQ_DCN.ShortestPaths)[t],'HighPLQ','w/Protective_SNPs',sep = '_')
  sp_df_high_damagingSNPs$type=paste(names(msbb_rnaseq2016_PLQ_DCN.ShortestPaths)[t],'HighPLQ','w/Damaging_SNPs',sep = '_')
  sp_df_hs_ppi_protectiveSNPs$type='hs_ppi_w/Protective_SNPs'
  sp_df_hs_ppi_damagingSNPs$type='hs_ppi_w/Damaging_SNPs'
  sp_df_high=rbind(sp_df_high_plq_dcn,sp_df_high_protectiveSNPs,sp_df_high_damagingSNPs,sp_df_hs_ppi_protectiveSNPs,sp_df_hs_ppi_damagingSNPs)
  g2=ggplot(sp_df_high,aes(sp,fill=type))+geom_density(alpha=0.5)+ggtitle(label = 'ShortestPath distribution',subtitle = paste(names(msbb_rnaseq2016_PLQ.DCN)[t],'HighPLQ_DCN w/ PLQ-assoc.genes,Tanzi variants',sep = '_'))
  #sp_df=rbind(sp_df_high_plq_dcn,sp_df_high_damagingSNPs,sp_df_high_protectiveSNPs,sp_df_low_plq_dcn,sp_df_low_protectiveSNPs,sp_df_low_damagingSNPs)
  
  #Compute betweenness for Low PLQ DCNs
  betweenness_df_low_stg=data.frame(betweenness=log10(msbb_rnaseq2016_PLQ_DCN.Betweenness[[t]]$Low),stringsAsFactors = F)
  betweenness_df_low_stg_protectiveSNPs=data.frame(betweenness=log10(msbb_rnaseq2016_PLQ_TanziProtective_DCN.Betweenness[[t]]$Low),stringsAsFactors = F)
  betweenness_df_low_stg_damagingSNPs=data.frame(betweenness=log10(msbb_rnaseq2016_PLQ_TanziDamaging_DCN.Betweenness[[t]]$Low),stringsAsFactors = F)
  betweenness_hs_ppi_protectiveSNPs=data.frame(betweenness=log10(hs_ppi_subnet.Centralities$TanziProtective$Betweenness),stringsAsFactors = F)
  betweenness_hs_ppi_damagingSNPs=data.frame(betweenness=log10(hs_ppi_subnet.Centralities$TanziDamaging$Betweenness),stringsAsFactors = F)
  betweenness_df_low_stg$type=paste(names(msbb_rnaseq2016_PLQ_DCN.Betweenness)[t],'LowPLQ-genes.DCN',sep = '_')
  betweenness_df_low_stg_protectiveSNPs$type=paste(names(msbb_rnaseq2016_PLQ_DCN.ShortestPaths)[t],'LowPLQ','w/Protective_SNPs',sep = '_')
  betweenness_df_low_stg_damagingSNPs$type=paste(names(msbb_rnaseq2016_PLQ_DCN.ShortestPaths)[t],'LowPLQ','w/Damaging_SNPs',sep = '_')
  betweenness_hs_ppi_protectiveSNPs$type='hs_ppi_w/Protective_SNPs'
  betweenness_hs_ppi_damagingSNPs$type='hs_ppi_w/Damaging_SNPs'
  betweenness_df_low=rbind(betweenness_hs_ppi_protectiveSNPs,betweenness_hs_ppi_damagingSNPs,betweenness_df_low_stg,betweenness_df_low_stg_damagingSNPs,betweenness_df_low_stg_protectiveSNPs)#
  g3=ggplot(betweenness_df_low,aes(betweenness,fill=type))+geom_density(alpha=0.5)+ggtitle(label = 'Betweenness distribution',subtitle = paste(names(msbb_rnaseq2016_PLQ.DCN)[t],'DCN w/ Low-PLQ-assoc.genes,Tanzi variants',sep = '_'))+ylab('log10 density')
  
  #Compute betweenness separately for High PLQ DCNs
  betweenness_df_high_stg=data.frame(betweenness=log10(msbb_rnaseq2016_PLQ_DCN.Betweenness[[t]]$High),stringsAsFactors = F)
  betweenness_df_high_stg_protectiveSNPs=data.frame(betweenness=log10(msbb_rnaseq2016_PLQ_TanziProtective_DCN.Betweenness[[t]]$High),stringsAsFactors = F)
  betweenness_df_high_stg_damagingSNPs=data.frame(betweenness=log10(msbb_rnaseq2016_PLQ_TanziDamaging_DCN.Betweenness[[t]]$High),stringsAsFactors = F)
  betweenness_hs_ppi_protectiveSNPs=data.frame(betweenness=log10(hs_ppi_subnet.Centralities$TanziProtective$Betweenness),stringsAsFactors = F)
  betweenness_hs_ppi_damagingSNPs=data.frame(betweenness=log10(hs_ppi_subnet.Centralities$TanziDamaging$Betweenness),stringsAsFactors = F)
  betweenness_df_high_stg$type=paste(names(msbb_rnaseq2016_PLQ_DCN.Betweenness)[t],'HighPLQ-genes.DCN',sep = '_')
  betweenness_df_high_stg_protectiveSNPs$type=paste(names(msbb_rnaseq2016_PLQ_DCN.ShortestPaths)[t],'HighPLQ','w/Protective_SNPs',sep = '_')
  betweenness_df_high_stg_damagingSNPs$type=paste(names(msbb_rnaseq2016_PLQ_DCN.ShortestPaths)[t],'HighPLQ','w/Damaging_SNPs',sep = '_')
  betweenness_hs_ppi_protectiveSNPs$type='hs_ppi_w/Protective_SNPs'
  betweenness_hs_ppi_damagingSNPs$type='hs_ppi_w/Damaging_SNPs'
  betweenness_df_high=rbind(betweenness_hs_ppi_protectiveSNPs,betweenness_hs_ppi_damagingSNPs,betweenness_df_high_stg,betweenness_df_high_stg_damagingSNPs,betweenness_df_high_stg_protectiveSNPs)#
  g4=ggplot(betweenness_df_high,aes(betweenness,fill=type))+geom_density(alpha=0.5)+ggtitle(label = 'Betweenness distribution',subtitle = paste(names(msbb_rnaseq2016_PLQ.DCN)[t],'DCN w/ High-PLQ-assoc.genes,Tanzi variants',sep = '_'))+ylab('log10 density')
  
  degree_df=data.frame(degree=log10(msbb_rnaseq2016_PLQ_TanziProtective_DCN.Degree[[t]]$Low),stringsAsFactors = F)
  degree_df_protectiveSNPs=data.frame(degree=log10(msbb_rnaseq2016_PLQ_TanziProtective_DCN.Degree[[t]]$Low),stringsAsFactors = F)
  degree_df_damagingSNPs=data.frame(degree=log10(msbb_rnaseq2016_PLQ_TanziDamaging_DCN.Degree[[t]]$Low),stringsAsFactors = F)
  degree_df_hs_ppi_protectiveSNPs=data.frame(degree=log10(hs_ppi_subnet.Centralities$TanziProtective$Degree))
  degree_df_hs_ppi_damagingSNPs=data.frame(degree=log10(hs_ppi_subnet.Centralities$TanziDamaging$Degree))
  degree_df$type=paste(names(msbb_rnaseq2016_PLQ_DCN.Degree)[t],'PLQ-genes.DCN',sep = '_')
  degree_df_protectiveSNPs$type=paste(names(msbb_rnaseq2016_PLQ_DCN.Degree)[t],'PLQ','w/Protective_SNPs',sep = '_')
  degree_df_damagingSNPs$type=paste(names(msbb_rnaseq2016_PLQ_DCN.Degree)[t],'PLQ','w/Damaging_SNPs',sep = '_')
  degree_df_hs_ppi_protectiveSNPs$type='hs_ppi_w/Protective_SNPs'
  degree_df_hs_ppi_damagingSNPs$type='hs_ppi_w/Damaging_SNPs'
  degree_df=rbind(degree_df,degree_df_damagingSNPs,degree_df_protectiveSNPs,degree_df_hs_ppi_protectiveSNPs,degree_df_hs_ppi_damagingSNPs)
  g5=ggplot(degree_df,aes(degree,fill=type))+geom_density(alpha=0.5)+ggtitle(label = 'Degree distribution',subtitle = paste(names(msbb_rnaseq2016_PLQ.DCN)[t],'DCN w/ PLQ-assoc.genes,Tanzi variants',sep = '_'))+ylab('log10 density')
  
  jpeg(filename = paste(names(msbb_rnaseq2016_PLQ.DCN)[t],'NetCentralityPlots.jpeg',sep = "_"),width = 4500,height = 1500,units = 'px',quality = 150)
  multiplot(g1,g2,g3,g4,g5,cols = 5)
  dev.off()
}

sp_stg_df_low_plq_dcn=data.frame(sp=log10(msbb_rnaseq2016_PLQ_DCN.ShortestPaths$STG$Low),stringsAsFactors = F)
sp_stg_df_high_plq_dcn=data.frame(sp=log10(msbb_rnaseq2016_PLQ_DCN.ShortestPaths$STG$High),stringsAsFactors = F)
sp_stg_low_damagingSNPs=data.frame(sp=log10(msbb_rnaseq2016_PLQ_TanziDamaging_DCN.ShortestPaths$STG$Low),stringsAsFactors = F)
sp_stg_high_damagingSNPs=data.frame(sp=log10(msbb_rnaseq2016_PLQ_TanziDamaging_DCN.ShortestPaths$STG$High),stringsAsFactors = F)
sp_stg_low_protectiveSNPs=data.frame(sp=log10(msbb_rnaseq2016_PLQ_TanziProtective_DCN.ShortestPaths$STG$Low),stringsAsFactors = F)
sp_stg_high_protectiveSNPs=data.frame(sp=log10(msbb_rnaseq2016_PLQ_TanziProtective_DCN.ShortestPaths$STG$High),stringsAsFactors = F)
sp_stg_df_low_plq_dcn$type='Low-PLQ.DCN'
sp_stg_df_high_plq_dcn$type='High-PLQ.DCN'
sp_stg_low_damagingSNPs$type='Low-PLQ.DCN w/DamagingSNPs'
sp_stg_high_damagingSNPs$type='High-PLQ.DCN w/DamagingSNPs'
sp_stg_low_protectiveSNPs$type='Low-PLQ.DCN w/ProtectiveSNPs'
sp_stg_high_protectiveSNPs$type='High-PLQ.DCN w/ProtectiveSNPs'


deg_stg_df_low_plq_dcn=data.frame(sp=log10(msbb_rnaseq2016_PLQ_DCN.Degree$STG$Low),stringsAsFactors = F)
deg_stg_df_high_plq_dcn=data.frame(sp=log10(msbb_rnaseq2016_PLQ_DCN.Degree$STG$High),stringsAsFactors = F)
deg_stg_low_damagingSNPs=data.frame(sp=log10(msbb_rnaseq2016_PLQ_TanziDamaging_DCN.Degree$STG$Low),stringsAsFactors = F)
deg_stg_high_damagingSNPs=data.frame(sp=log10(msbb_rnaseq2016_PLQ_TanziDamaging_DCN.Degree$STG$High),stringsAsFactors = F)
deg_stg_low_protectiveSNPs=data.frame(sp=log10(msbb_rnaseq2016_PLQ_TanziProtective_DCN.Degree$STG$Low),stringsAsFactors = F)
deg_stg_high_protectiveSNPs=data.frame(sp=log10(msbb_rnaseq2016_PLQ_TanziProtective_DCN.Degree$STG$High),stringsAsFactors = F)
deg_stg_df_low_plq_dcn$type='Low-PLQ.DCN'
deg_stg_df_high_plq_dcn$type='High-PLQ.DCN'
deg_stg_low_damagingSNPs$type='Low-PLQ.DCN w/DamagingSNPs'
deg_stg_high_damagingSNPs$type='High-PLQ.DCN w/DamagingSNPs'
deg_stg_low_protectiveSNPs$type='Low-PLQ.DCN w/ProtectiveSNPs'
deg_stg_high_protectiveSNPs$type='High-PLQ.DCN w/ProtectiveSNPs'

deg_stg_df=rbind(deg_stg_df_low_plq_dcn,deg_stg_df_high_plq_dcn,deg_stg_low_damagingSNPs,deg_stg_high_damagingSNPs,deg_stg_low_protectiveSNPs,deg_stg_high_protectiveSNPs)
sp_stg_df=rbind(sp_stg_df_low_plq_dcn,sp_stg_df_high_plq_dcn,sp_stg_low_damagingSNPs,sp_stg_high_damagingSNPs,sp_stg_low_protectiveSNPs,sp_stg_high_protectiveSNPs)


plot(compareCluster(geneClusters = list(hs_ppi_subnet_TanziProtective=select(x = org.Hs.eg.db,keys = V(hs_ppi.subnet$TanziProtective)$name,columns = 'ENTREZID',keytype = 'SYMBOL')[,2],
                                        hs_ppi_subnet_TanziDamaging=select(x = org.Hs.eg.db,keys = V(hs_ppi.subnet$TanziDamaging)$name,columns = 'ENTREZID',keytype = 'SYMBOL')[,2],
                                        STG_PLQ_TanziDamaging_DCN=select(x = org.Hs.eg.db,keys = V(msbb_rnaseq2016_PLQ_TanziDamaging.DCN$PHG$Low)$name,columns = 'ENTREZID',keytype = 'SYMBOL')[,2],
                                        STG_PLQ_TanziProtective_DCN=select(x = org.Hs.eg.db,keys = V(msbb_rnaseq2016_PLQ_TanziProtective.DCN$PHG$Low)$name,columns = 'ENTREZID',keytype = 'SYMBOL')[,2]),fun = 'enrichPathway',pvalueCutoff=0.1,pAdjustMethod='BH'))
#Test diff edge weight cutoff thresholds
low_cutoffs=seq(from=-0.25+1,to=-0.5+1,by = -0.01)
high_cutoffs=seq(from=0.25+1,to=0.5+1,by = 0.01)
rewired_edges=rewired_degree=rewired_betweenness=vector(mode = 'list',length = 4)
names(rewired_edges)=names(rewired_degree)=names(rewired_betweenness)=c('FP','IFG','PHG','STG')
rewired_edges$FP=rewired_edges$IFG=rewired_edges$PHG=rewired_edges$STG=vector(mode = 'list')
rewired_degree$FP=rewired_degree$IFG=rewired_degree$PHG=rewired_degree$STG=vector(mode = 'list')
rewired_betweenness$FP=rewired_betweenness$IFG=rewired_betweenness$PHG=rewired_betweenness$STG=vector(mode = 'list')

for (t in 1:4){
  for (c in 1:length(high_cutoffs)){
    msbb_rnaseq2016_PLQ.PosCoexp[[t]]$Low=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016_PLQ.DCN[[t]]$Low,es = which(E(msbb_rnaseq2016_PLQ.DCN[[t]]$Low)$weight>high_cutoffs[c]),names = T)),directed = F)
    E(msbb_rnaseq2016_PLQ.PosCoexp[[t]]$Low)$weight=E(msbb_rnaseq2016_PLQ.DCN[[t]]$Low)$weight[which(E(msbb_rnaseq2016_PLQ.DCN[[t]]$Low)$weight>high_cutoffs[c])]
    
    msbb_rnaseq2016_PLQ.PosCoexp[[t]]$High=graph.data.frame(d=data.frame(ends(graph = msbb_rnaseq2016_PLQ.DCN[[t]]$High,es = which(E(msbb_rnaseq2016_PLQ.DCN[[t]]$High)$weight>high_cutoffs[c]),names = T)),directed = F)
    E(msbb_rnaseq2016_PLQ.PosCoexp[[t]]$High)$weight=E(msbb_rnaseq2016_PLQ.DCN[[t]]$High)$weight[which(E(msbb_rnaseq2016_PLQ.DCN[[t]]$High)$weight>high_cutoffs[c])]
    msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF[[t]]=data.frame(Common_Genes=intersect(V(msbb_rnaseq2016_PLQ.PosCoexp[[t]]$Low)$name,V(msbb_rnaseq2016_PLQ.PosCoexp[[t]]$High)$name),
    Low=degree(graph = msbb_rnaseq2016_PLQ.PosCoexp[[t]]$Low,v = intersect(V(msbb_rnaseq2016_PLQ.PosCoexp[[t]]$Low)$name,V(msbb_rnaseq2016_PLQ.PosCoexp[[t]]$High)$name)),
    High=degree(graph = msbb_rnaseq2016_PLQ.PosCoexp[[t]]$High,v = intersect(V(msbb_rnaseq2016_PLQ.PosCoexp[[t]]$Low)$name,V(msbb_rnaseq2016_PLQ.PosCoexp[[t]]$High)$name)),
    Difference=degree(graph = msbb_rnaseq2016_PLQ.PosCoexp[[t]]$Low,v = intersect(V(msbb_rnaseq2016_PLQ.PosCoexp[[t]]$Low)$name,V(msbb_rnaseq2016_PLQ.PosCoexp[[t]]$High)$name))-degree(graph = msbb_rnaseq2016_PLQ.PosCoexp[[t]]$High,v = intersect(V(msbb_rnaseq2016_PLQ.PosCoexp[[t]]$Low)$name,V(msbb_rnaseq2016_PLQ.PosCoexp[[t]]$High)$name)),
    Difference_BWN=betweenness(graph = msbb_rnaseq2016_PLQ.PosCoexp[[t]]$Low,v = intersect(V(msbb_rnaseq2016_PLQ.PosCoexp[[t]]$Low)$name,V(msbb_rnaseq2016_PLQ.PosCoexp[[t]]$High)$name))-betweenness(graph = msbb_rnaseq2016_PLQ.PosCoexp[[t]]$High,v = intersect(V(msbb_rnaseq2016_PLQ.PosCoexp[[t]]$Low)$name,V(msbb_rnaseq2016_PLQ.PosCoexp[[t]]$High)$name)),stringsAsFactors = F)
    msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF[[t]]=msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF[[t]][order(msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF[[t]]$Difference,decreasing = T),]
    rewired_edges[[t]][[c]]=length(E(msbb_rnaseq2016_PLQ.PosCoexp[[t]]$Low))-length(E(msbb_rnaseq2016_PLQ.PosCoexp[[t]]$High))
    rewired_degree[[t]][[c]]=msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF[[t]]$Difference[msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF[[t]]$Difference>0]
    rewired_betweenness[[t]][[c]]=msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF[[t]]$Difference_BWN[msbb_rnaseq2016_PLQGenes_PosCoexp.rewiredDF[[t]]$Difference_BWN>0]
  }  
}
names(rewired_edges$FP)=names(rewired_edges$IFG)=names(rewired_edges$PHG)=names(rewired_edges$STG)=high_cutoffs
names(rewired_degree$FP)=names(rewired_degree$IFG)=names(rewired_degree$PHG)=names(rewired_degree$STG)=high_cutoffs
names(rewired_betweenness$FP)=names(rewired_betweenness$IFG)=names(rewired_betweenness$PHG)=names(rewired_betweenness$STG)=high_cutoffs

rw_bwn_fp=data.frame(Cutoff_threshold=high_cutoffs-1,rw_betweenness=unname(unlist(lapply(rewired_betweenness$FP,sum))),stringsAsFactors = F)
rw_deg_fp=data.frame(Cutoff_threshold=high_cutoffs-1,rw_degree=unname(unlist(lapply(rewired_degree$FP,sum))),stringsAsFactors = F)
rw_edge_fp=data.frame(Cutoff_threshold=high_cutoffs-1,rw_edges=unlist(rewired_edges$FP),stringsAsFactors = F)
g1=ggplot(rw_edge_fp,aes(x=Cutoff_threshold,y=rw_edges))+geom_point()+geom_smooth(method = 'loess')+ggtitle(label = 'Rewired Edges - FP',subtitle = '#rewired edges across edge cutoff thresholds')+theme(plot.title = element_text(face = 'bold',size = 15),axis.title = element_text(face = 'bold',size = 15))
g2=ggplot(rw_deg_fp,aes(x=Cutoff_threshold,y=rw_degree))+geom_point()+geom_smooth(method = 'loess')+ggtitle(label = 'Rewired Degrees - FP',subtitle = '#rewired genes across edge cutoff thresholds')+theme(plot.title = element_text(face = 'bold',size = 15),axis.title.x = element_text(face = 'bold',size = 15),axis.title.y = element_text(face = 'bold',size = 12))
multiplot(g1,g2,cols=2)

rw_bwn_ifg=data.frame(Cutoff_threshold=high_cutoffs-1,rw_betweenness=unname(unlist(lapply(rewired_betweenness$IFG,sum))),stringsAsFactors = F)
rw_deg_ifg=data.frame(Cutoff_threshold=high_cutoffs-1,rw_degree=unname(unlist(lapply(rewired_degree$IFG,sum))),stringsAsFactors = F)
rw_edge_ifg=data.frame(Cutoff_threshold=high_cutoffs-1,rw_edges=unlist(rewired_edges$IFG),stringsAsFactors = F)
g3=ggplot(rw_edge_ifg,aes(x=Cutoff_threshold,y=rw_edges))+geom_point()+geom_smooth(method = 'loess')+ggtitle(label = 'Rewired Edges - IFG',subtitle = '#rewired edges across edge cutoff thresholds')+theme(plot.title = element_text(face = 'bold',size = 15),axis.title = element_text(face = 'bold',size = 15))
g4=ggplot(rw_deg_ifg,aes(x=Cutoff_threshold,y=rw_degree))+geom_point()+geom_smooth(method = 'loess')+ggtitle(label = 'Rewired Degrees - IFG',subtitle = '#rewired genes across edge cutoff thresholds')+theme(plot.title = element_text(face = 'bold',size = 15),axis.title.x = element_text(face = 'bold',size = 15),axis.title.y = element_text(face = 'bold',size = 12))
multiplot(g3,g4,cols=2)

rw_bwn_phg=data.frame(Cutoff_threshold=high_cutoffs-1,rw_betweenness=unname(unlist(lapply(rewired_betweenness$PHG,sum))),stringsAsFactors = F)
rw_deg_phg=data.frame(Cutoff_threshold=high_cutoffs-1,rw_degree=unname(unlist(lapply(rewired_degree$PHG,sum))),stringsAsFactors = F)
rw_edge_phg=data.frame(Cutoff_threshold=high_cutoffs-1,rw_edges=unlist(rewired_edges$PHG),stringsAsFactors = F)
g5=ggplot(rw_edge_phg,aes(x=Cutoff_threshold,y=rw_edges))+geom_point()+geom_smooth(method = 'loess')+ggtitle(label = 'Rewired Edges - PHG',subtitle = '#rewired edges across edge cutoff thresholds')+theme(plot.title = element_text(face = 'bold',size = 15),axis.title = element_text(face = 'bold',size = 15))
g6=ggplot(rw_deg_phg,aes(x=Cutoff_threshold,y=rw_degree))+geom_point()+geom_smooth(method = 'loess')+ggtitle(label = 'Rewired Degrees - PHG',subtitle = '#rewired genes across edge cutoff thresholds')+theme(plot.title = element_text(face = 'bold',size = 15),axis.title.x = element_text(face = 'bold',size = 15),axis.title.y = element_text(face = 'bold',size = 12))
multiplot(g5,g6,cols=2)

rw_bwn_stg=data.frame(Cutoff_threshold=high_cutoffs-1,rw_betweenness=unname(unlist(lapply(rewired_betweenness$STG,sum))),stringsAsFactors = F)
rw_deg_stg=data.frame(Cutoff_threshold=high_cutoffs-1,rw_degree=unname(unlist(lapply(rewired_degree$STG,sum))),stringsAsFactors = F)
rw_edge_stg=data.frame(Cutoff_threshold=high_cutoffs-1,rw_edges=unlist(rewired_edges$STG),stringsAsFactors = F)
g7=ggplot(rw_edge_stg,aes(x=Cutoff_threshold,y=rw_edges))+geom_point()+geom_smooth(method = 'loess')+ggtitle(label = 'Rewired Edges - STG',subtitle = '#rewired edges across edge cutoff thresholds')+theme(plot.title = element_text(face = 'bold',size = 15),axis.title = element_text(face = 'bold',size = 15))
g8=ggplot(rw_deg_stg,aes(x=Cutoff_threshold,y=rw_degree))+geom_point()+geom_smooth(method = 'loess')+ggtitle(label = 'Rewired Degrees - STG',subtitle = '#rewired genes across edge cutoff thresholds')+theme(plot.title = element_text(face = 'bold',size = 15),axis.title.x = element_text(face = 'bold',size = 15),axis.title.y = element_text(face = 'bold',size = 12))
multiplot(g7,g8,cols=2)

rw_fp=data.frame(Cutoff_threshold=high_cutoffs,rw_edges=unlist(unname(rewired_edges$FP)),rw_degree=unname(unlist(lapply(rewired_degree$IFG,sum))),rw_betweenness=unname(unlist(lapply(rewired_betweenness$FP,sum))))
rw_fp_df=melt(rw_fp,id='Cutoff_threshold',Centralities=c("rw_edges","rw_degree"))
rw_ifg=data.frame(Cutoff_threshold=high_cutoffs,rw_edges=unlist(unname(rewired_edges$IFG)),rw_degree=unlist(lapply(rewired_degree$IFG,sum)),rw_betweenness=unname(unlist(lapply(rewired_betweenness$IFG,sum))))
rw_ifg_df=melt(rw_ifg,id='Cutoff_threshold',Centralities=c("rw_edges","rw_degree","rw_betweenness"))
ggplot(rw_ifg_df,aes(Cutoff_threshold,value,colour=variable))+geom_point()+geom_smooth(method = 'loess')

rw_phg=data.frame(Cutoff_threshold=high_cutoffs,rw_edges=unlist(unname(rewired_edges$PHG)),rw_degree=unlist(lapply(rewired_degree$PHG,sum)))
rw_stg=data.frame(Cutoff_threshold=high_cutoffs,rw_edges=unlist(unname(rewired_edges$STG)),rw_degree=unlist(lapply(rewired_degree$STG,sum)))
rw_fp$region='FP'
rw_ifg$region='IFG'
rw_phg$region='PHG'
rw_stg$region='STG'
rw=rbind(rw_ifg,rw_phg,rw_stg)
