library(org.Hs.eg.db)
library(igraph)
library(cocor)
library(gdata)

setwd("/Users/sandeepamberkar/Work/Data/MSMM-MSBB-HBTRC-Data/MSMM_RNAseq/MSMM_RNAseq_FinalRelease/")
load("msmm_rnaseq_contScore_cocor_p005.RData")

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
mm_BrainProteome_4x_humanOrtholog=mouse_human_ortholog$Associated.Gene.Name.1[which(mouse_human_ortholog$Associated.Gene.Name%in%mm_BrainProteome_sharma.4x$Gene.names[which(mm_BrainProteome_sharma.4x$X.four.fold.expressed.in.atleast.one.cell.type=="+")])]

human_genes_orgDB=toTable(org.Hs.egSYMBOL)
overlap_dist=vector(mode = "list",length = 3)
names(overlap_dist)=c("lcm_neurons","hipp_proteome","mm_proteome")
overlap_dist[[1]]=overlap_dist[[2]]=overlap_dist[[3]]=list()

for (o in 1:10000){
  hippocampal_proteome_AD.random=sample(x=human_genes_orgDB$symbol,size = length(union(hippocampal_proteome_ADdown$Genes,hippocampal_proteome_ADup$Genes)),replace = F)
  lcm_neurons_AD.random=sample(x=human_genes_orgDB$symbol,size = length(ffpe_neurons_genes),replace = F)
  mm_proteome_AD.random=sample(x = human_genes_orgDB$symbol,size = length(mm_BrainProteome_4x_humanOrtholog),replace = F)
  #hippocampal_proteome_hondius_ADup.random=sample(x=human_genes_orgDB$symbol,size = length(hippocampal_proteome_hondius_ADup$Gene..leading.protein.),replace = F)
  overlap_dist$hipp_proteome[[o]]=length(intersect(hippocampal_proteome_AD.random,V(msmm_rnaseq.DCN$BM_36$Std)$name))
  overlap_dist$lcm_neurons[[o]]=length(intersect(lcm_neurons_AD.random,V(msmm_rnaseq.DCN$BM_36$Std)$name))
  overlap_dist$mm_proteome[[o]]=length(intersect(mm_proteome_AD.random,V(msmm_rnaseq.DCN$BM_36$Std)$name))
  #overlap_dist[[o]]$down=length(intersect(hippocampal_proteome_hondius_ADdown.random,V(msmm_rnaseq.DCN$BM_36$Std)$name))
}

par(mfrow=c(3,1))
hist(unlist(overlap_dist$lcm_neurons),xlab = "#Overlaps",main = "#Overlap distribution for Drummund etal LCM Neurons",col="blue4",xlim = c(0,50),breaks = 50)
abline(v=length(intersect(V(msmm_rnaseq.DCN$BM_36$Std)$name,ffpe_neurons_genes)),col="red4",lwd=2)
text(length(intersect(V(msmm_rnaseq.DCN$BM_36$Std)$name,ffpe_neurons_genes))+0.5,0, "pval=7.4e-12", col = "black",adj = c(-.1, -.1))
hist(unlist(overlap_dist$hipp_proteome),xlab = "#Overlaps",main = "#Overlap distribution for Hondius etal Hippocampal proteome",col="blue4",xlim = c(0,50),breaks = 50)
abline(v=length(intersect(V(msmm_rnaseq.DCN$BM_36$Std)$name,union(hippocampal_proteome_ADdown$Genes,hippocampal_proteome_ADup$Genes))),col="red4",lwd=2)
text(length(intersect(V(msmm_rnaseq.DCN$BM_36$Std)$name,union(hippocampal_proteome_ADdown$Genes,hippocampal_proteome_ADup$Genes)))+0.5,0, "pval=3.9e-15", col = "black",adj = c(-.1, -.1))
hist(unlist(overlap_dist$mm_proteome),xlab = "#Overlaps",main = "#Overlap distribution for Sharma etal cellular proteome",col="blue4",xlim = c(0,250),breaks = 50)
abline(v=length(intersect(V(msmm_rnaseq.DCN$BM_36$Std)$name,mm_BrainProteome_4x_humanOrtholog)),col="red4",lwd=2)
text(length(intersect(V(msmm_rnaseq.DCN$BM_36$Std)$name,mm_BrainProteome_4x_humanOrtholog))+0.5,0, "pval=2.5e-16", col = "black",adj = c(-.1, -.1))

ebc=cluster_edge_betweenness(graph = msmm_rnaseq.DCN$BM_36$Std,directed = F,modularity = T,weights = E(msmm_rnaseq.DCN$BM_36$Std)$weight)
comms=communities(ebc)
indices=names(unlist(lapply(communities(x = ebc),function(x)which(length(x)>=10))))


