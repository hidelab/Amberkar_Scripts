library(org.Hs.eg.db)
library(igraph)
library(cocor)
library(gdata)

setwd("/Users/sandeepamberkar/Work/Data/MSMM-MSBB-HBTRC-Data/MSMM_RNAseq/MSMM_RNAseq_FinalRelease/")
load("msmm_rnaseq_contScore_cocor_p005.RData")

msmm_rnaseq_Plaque.Corr005.Genes_cluster=vector(mode = "list",length = 3)
names(msmm_rnaseq_Plaque.Corr005.Genes_cluster)=names(msmm_rnaseq)
msmm_rnaseq_Plaque.Corr005.Genes_cluster$BM_10=select(x = org.Hs.eg.db,keys = msmm_rnaseq_Plaque.Corr005.Genes$BM_10,columns = "ENTREZID",keytype = "SYMBOL")[,2]
msmm_rnaseq_Plaque.Corr005.Genes_cluster$BM_10=msmm_rnaseq_Plaque.Corr005.Genes_cluster$BM_10[which(is.na(msmm_rnaseq_Plaque.Corr005.Genes_cluster$BM_10)==F)]
msmm_rnaseq_Plaque.Corr005.Genes_cluster$BM_22=select(x = org.Hs.eg.db,keys = msmm_rnaseq_Plaque.Corr005.Genes$BM_22,columns = "ENTREZID",keytype = "SYMBOL")[,2]
msmm_rnaseq_Plaque.Corr005.Genes_cluster$BM_22=msmm_rnaseq_Plaque.Corr005.Genes_cluster$BM_22[which(is.na(msmm_rnaseq_Plaque.Corr005.Genes_cluster$BM_22)==F)]
msmm_rnaseq_Plaque.Corr005.Genes_cluster$BM_36=select(x = org.Hs.eg.db,keys = msmm_rnaseq_Plaque.Corr005.Genes$BM_36,columns = "ENTREZID",keytype = "SYMBOL")[,2]
msmm_rnaseq_Plaque.Corr005.Genes_cluster$BM_36=msmm_rnaseq_Plaque.Corr005.Genes_cluster$BM_36[which(is.na(msmm_rnaseq_Plaque.Corr005.Genes_cluster$BM_36)==F)]


ffpe_neurons=read.xls("../../../BrainExpression_Datasets/FFPE-Proteomics-Wiesniewski/srep15456-s6.xls",sheet = 2,header=T,as.is=T)
#ffpe_neurons=ffpe_neurons$Accession[which((ffpe_neurons$Neuronal..more.stringent.database.=="Yes")&(ffpe_neurons$Alzheimer.s.associated.protein=="Yes"))]
ffpe_neurons2=gsub(pattern = "\\-[0-9]",replacement = "",x = ffpe_neurons$Accession[which((ffpe_neurons$Neuronal..more.stringent.database.=="Yes")&(ffpe_neurons$Alzheimer.s.associated.protein=="Yes"))])
#Resolve dual IDs manually
ffpe_neurons2[c(6,17,56,72,80,121,191)]=c("P68366","P09104","P38606","P30048","P18669","P48735","P27338")
ffpe_neurons_genes=select(x = org.Hs.eg.db,keys = ffpe_neurons2,columns = "SYMBOL",keytype = "UNIPROT")[,2]

hippocampal_proteome_ADdown=read.xls("../../../BrainExpression_Datasets/Hippocampus_Proteome_AD_Hondius_etal/mmc2.xlsx",sheet = 2,header=T,as.is=T)
hippocampal_proteome_ADup=read.xls("../../../BrainExpression_Datasets/Hippocampus_Proteome_AD_Hondius_etal/mmc2.xlsx",sheet = 1,header=T,as.is=T)

mouse_human_ortholog=read.table("../../../BrainExpression_Datasets/Mouse_Human_Orthologs_EnsemblGRC38.txt",header = T,sep="\t",as.is=T)
mm_BrainProteome_sharma.4x=read.xls("../../../BrainExpression_Datasets/Mouse_Brain_Proteome_Sharma_etal/nn.4160-S8.xlsx",sheet = 1,skip=1,header=T,as.is=T)

  
hippocampal_proteome_AD_PHG_DCN=c(c(intersect(hippocampal_proteome_ADdown$Genes,V(msmm_rnaseq.DCN$BM_36$Std)$name),c("TUBB4A","NPTN","NAP1L2","NAP1L3")),c(intersect(hippocampal_proteome_ADup$Genes,V(msmm_rnaseq.DCN$BM_36$Std)$name)))
lcm_neurons_DCN_PHG_AD=intersect(select(x = org.Hs.eg.db,keys = ffpe_neurons2,columns = "SYMBOL",keytype = "UNIPROT")[,2],V(msmm_rnaseq.DCN$BM_36$Std)$name)
mm_BrainProteome_humanOrtholog=vector(mode = "list",length = 4)
mm_BrainProteome_sharma.CellType=vector(mode = "list",length=4)
names(mm_BrainProteome_sharma.CellType)=names(mm_BrainProteome_humanOrtholog)=c("Microglia","Astrocytes","Oligodendrocytes","Neurons")
for (c in 1:4){
  mm_BrainProteome_sharma.CellType[[c]]=c(unlist(strsplit(x = grep(pattern = ";",mm_BrainProteome_sharma$Gene.names.1[which(mm_BrainProteome_sharma[,c+11]=="+")],value = T),split = ";")),grep(pattern = ";",mm_BrainProteome_sharma$Gene.names.1[which(mm_BrainProteome_sharma[,c+11]=="+")],invert = T,value = T))
  mm_BrainProteome_humanOrtholog[[c]]=mouse_human_ortholog$Associated.Gene.Name.1[which(mm_BrainProteome_sharma.CellType[[c]]%in%mouse_human_ortholog$Associated.Gene.Name)]
}



human_genes_orgDB=toTable(org.Hs.egSYMBOL)
overlap_dist=vector(mode = "list",length = 3)
overlap_dist2=vector(mode = "list",length = 1)
names(overlap_dist)=c("lcm_neurons","hipp_proteome","mm_proteome")
overlap_dist[[1]]=overlap_dist[[2]]=overlap_dist[[3]]=list()

for (o in 1:10000){
  hippocampal_proteome_AD.random=sample(x=human_genes_orgDB$symbol,size = length(union(hippocampal_proteome_ADdown$Genes,hippocampal_proteome_ADup$Genes)),replace = F)
  lcm_neurons_AD.random=sample(x=human_genes_orgDB$symbol,size = length(ffpe_neurons_genes),replace = F)
  mm_proteome_AD.random=sample(x = human_genes_orgDB$symbol,size = length(mm_BrainProteome_4x_humanOrtholog),replace = F)
  #hippocampal_proteome_hondius_ADup.random=sample(x=human_genes_orgDB$symbol,size = length(hippocampal_proteome_hondius_ADup$Gene..leading.protein.),replace = F)
  overlap_dist$hipp_proteome[[o]]=length(intersect(hippocampal_proteome_AD.random,msmm_rnaseq.DCN_PHG_rewired$Common.Genes))
  overlap_dist$lcm_neurons[[o]]=length(intersect(lcm_neurons_AD.random,msmm_rnaseq.DCN_PHG_rewired$Common.Genes))
  overlap_dist$mm_proteome[[o]]=length(intersect(mm_proteome_AD.random,msmm_rnaseq.DCN_PHG_rewired$Common.Genes))
  #overlap_dist[[o]]$down=length(intersect(hippocampal_proteome_hondius_ADdown.random,V(msmm_rnaseq.DCN$BM_36$Std)$name))
}
for (o in 1:10000){
  db_intersect=sample(x=human_genes_orgDB$symbol,size = length(intersect(union(hippocampal_proteome_ADdown$Genes,hippocampal_proteome_ADup$Genes),
                                                                         mm_BrainProteome_4x_humanOrtholog)),replace = F)
  overlap_dist2[[o]]=length(intersect(db_intersect,msmm_rnaseq.DCN_PHG_rewired$Common.Genes))
}

obj1_phg_lcmNeurons=newGeneOverlap(listA = ffpe_neurons_genes,listB = msmm_rnaseq.DCN_PHG_rewired$Common.Genes,spec = "hg19.gene")
obj2_phg_hippProteome=newGeneOverlap(listA = union(hippocampal_proteome_ADdown$Genes,hippocampal_proteome_ADup$Genes),listB = msmm_rnaseq.DCN_PHG_rewired$Common.Genes,spec = "hg19.gene")
obj3_phg_mmProteome=newGeneOverlap(listA = mm_BrainProteome_4x_humanOrtholog,listB = msmm_rnaseq.DCN_PHG_rewired$Common.Genes,spec = "hg19.gene")
obj4_phg_mm_hipp_Proteome=newGeneOverlap(listA = intersect(union(hippocampal_proteome_ADdown$Genes,hippocampal_proteome_ADup$Genes),mm_BrainProteome_4x_humanOrtholog),listB = msmm_rnaseq.DCN_PHG_rewired$Common.Genes,spec = "hg19.gene")

par(mfrow=c(3,1))
hist(unlist(overlap_dist$lcm_neurons),xlab = "#Overlaps",main = "#Overlap distribution for Drummund etal LCM Neurons, rewired PHG",col="blue4",xlim = c(0,60),breaks = 50)
abline(v=length(intersect(V(msmm_rnaseq.DCN$BM_36$Std)$name,ffpe_neurons_genes)),col="red4",lwd=2)
text(length(intersect(V(msmm_rnaseq.DCN$BM_36$Std)$name,ffpe_neurons_genes)),0, paste("pval=",testGeneOverlap(obj1_phg_lcmNeurons)@pval,sep=""), col = "black",adj = c(-.1, -.1))
hist(unlist(overlap_dist$hipp_proteome),xlab = "#Overlaps",main = "#Overlap distribution for Hondius etal Hippocampal proteome, rewired PHG",col="olivedrab3",xlim = c(0,60),breaks = 50)
abline(v=length(intersect(V(msmm_rnaseq.DCN$BM_36$Std)$name,union(hippocampal_proteome_ADdown$Genes,hippocampal_proteome_ADup$Genes))),col="red4",lwd=2)
text(length(intersect(V(msmm_rnaseq.DCN$BM_36$Std)$name,union(hippocampal_proteome_ADdown$Genes,hippocampal_proteome_ADup$Genes))),0, paste("pval=",testGeneOverlap(obj2_phg_hippProteome)@pval,sep=""), col = "black",adj = c(-.1, -.1))
hist(unlist(overlap_dist$mm_proteome),xlab = "#Overlaps",main = "#Overlap distribution for Sharma etal cellular proteome, rewired PHG",col="firebrick1",xlim = c(0,250),breaks = 50)
abline(v=length(intersect(V(msmm_rnaseq.DCN$BM_36$Std)$name,mm_BrainProteome_4x_humanOrtholog)),col="red4",lwd=2)
text(length(intersect(V(msmm_rnaseq.DCN$BM_36$Std)$name,mm_BrainProteome_4x_humanOrtholog)),0, paste("pval=",testGeneOverlap(obj3_phg_mmProteome)@pval,sep=""), col = "black",adj = c(-.1, -.1))


msmm_rnaseq.DCN_PHG_Filt1.5_Low=read.graph("BM_36_Low_p005_DiffCoexp_Subgraph_Filt1.5DCorr.graphml",format = "graphml")
msmm_rnaseq.DCN_PHG_Filt1.5_High=read.graph("BM_36_High_p005_DiffCoexp_Subgraph_Filt1.5DCorr.graphml",format = "graphml")

ebc=cluster_edge_betweenness(graph = msmm_rnaseq.DCN_PHG_Filt1.5_Low,directed = F,modularity = T,weights = E(msmm_rnaseq.DCN$BM_36$Std)$weight)
comms=communities(ebc)
indices=names(unlist(lapply(communities(x = ebc),function(x)which(length(x)>=10))))


