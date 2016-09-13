library(org.Hs.eg.db)
library(igraph)
library(cocor)
library(gdata)
library(dataframes2xls)
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

#Rewiring analysis
g1=graph.data.frame(d = ends(graph = msmm_rnaseq.DCN$BM_36$High,es = which(E(msmm_rnaseq.DCN$BM_36$High)$weight>=1),names = T),directed = F)
g2=graph.data.frame(d = ends(graph = msmm_rnaseq.DCN$BM_36$Low,es = which(E(msmm_rnaseq.DCN$BM_36$Low)$weight>=1),names = T),directed = F)
g3=graph.data.frame(d = ends(graph = msmm_rnaseq.DCN$BM_22$High,es = which(E(msmm_rnaseq.DCN$BM_22$High)$weight>=1),names = T),directed = F)
g4=graph.data.frame(d = ends(graph = msmm_rnaseq.DCN$BM_22$Low,es = which(E(msmm_rnaseq.DCN$BM_22$Low)$weight>=1),names = T),directed = F)
g11=induced_subgraph(graph = g1,vids = intersect(V(g1)$name,V(g2)$name))
g22=induced_subgraph(graph = g2,vids = intersect(V(g1)$name,V(g2)$name))
g33=induced_subgraph(graph = g3,vids = intersect(V(g3)$name,V(g4)$name))
g44=induced_subgraph(graph = g4,vids = intersect(V(g3)$name,V(g4)$name))
msmm_rnaseq.DCN_rewired=vector(mode = "list",length = 2)
names(msmm_rnaseq.DCN_rewired)=c("BM_22","BM_36")
msmm_rnaseq.DCN_rewired$BM_22=data.frame(Common.Genes=names(degree(graph = g33,v = sort(names(degree(g33))))),BM22_LowPLQ=unname(degree(graph = g33,v = sort(names(degree(g33))))),BM22_HighPLQ=unname(degree(graph = g44,v = sort(names(degree(g44))))),Difference=abs(unname(degree(graph = g33,v = sort(names(degree(g33)))))-unname(degree(graph = g44,v = sort(names(degree(g44)))))),stringsAsFactors = F)
msmm_rnaseq.DCN_rewired$BM_36=data.frame(Common.Genes=names(degree(graph = g11,v = sort(names(degree(g11))))),BM36_LowPLQ=unname(degree(graph = g11,v = sort(names(degree(g11))))),BM36_HighPLQ=unname(degree(graph = g22,v = sort(names(degree(g22))))),Difference=abs(unname(degree(graph = g11,v = sort(names(degree(g11)))))-unname(degree(graph = g22,v = sort(names(degree(g22)))))),stringsAsFactors = F)
write.table(msmm_rnaseq.DCN_rewired$BM_22,"BM22_DCN_Filt1.0_RewiredGenes.txt",sep="\t",col.names = T,row.names = T,quote = F)
write.table(msmm_rnaseq.DCN_rewired$BM_36,"BM36_DCN_Filt1.0_RewiredGenes.txt",sep="\t",col.names = T,row.names = T,quote = F)

#Validation and Hypothesis generation
ffpe_neurons=read.xls("../../../BrainExpression_Datasets/FFPE-Proteomics-Wiesniewski/srep15456-s6.xls",sheet = 2,header=T,as.is=T)
#ffpe_neurons=ffpe_neurons$Accession[which((ffpe_neurons$Neuronal..more.stringent.database.=="Yes")&(ffpe_neurons$Alzheimer.s.associated.protein=="Yes"))]
ffpe_neurons2=gsub(pattern = "\\-[0-9]",replacement = "",x = ffpe_neurons$Accession[which((ffpe_neurons$Neuronal..more.stringent.database.=="Yes")&(ffpe_neurons$Alzheimer.s.associated.protein=="Yes"))])
#Resolve dual IDs manually
ffpe_neurons2[c(6,17,56,72,80,121,191)]=c("P68366","P09104","P38606","P30048","P18669","P48735","P27338")
ffpe_neurons_genes=select(x = org.Hs.eg.db,keys = ffpe_neurons2,columns = "SYMBOL",keytype = "UNIPROT")[,2]

hippocampal_proteome_ADdown=read.xls("../../../BrainExpression_Datasets/Hippocampus_Proteome_AD_Hondius_etal/mmc2.xlsx",sheet = 2,header=T,as.is=T)
hippocampal_proteome_ADup=read.xls("../../../BrainExpression_Datasets/Hippocampus_Proteome_AD_Hondius_etal/mmc2.xlsx",sheet = 1,header=T,as.is=T)
hippocampal_proteome_AD=union(c(grep(pattern = ";",hippocampal_proteome_ADup$Genes,value = T,invert = T),unlist(lapply(strsplit(grep(pattern = ";",hippocampal_proteome_ADup$Genes,value = T),split = ";"),`[[`,1))),
                          c(grep(pattern = ";",hippocampal_proteome_ADdown$Genes,value = T,invert = T),unlist(lapply(strsplit(grep(pattern = ";",hippocampal_proteome_ADdown$Genes,value = T),split = ";"),`[[`,1))))
mouse_human_ortholog=read.table("../../../BrainExpression_Datasets/Mouse_Human_Orthologs_EnsemblGRC38.txt",header = T,sep="\t",as.is=T)
mm_BrainProteome_sharma=read.xls("../../../BrainExpression_Datasets/Mouse_Brain_Proteome_Sharma_etal/nn.4160-S8.xlsx",sheet = 1,skip=1,header=T,as.is=T)


mm_BrainProteome_humanOrtholog=vector(mode = "list",length = 4)
mm_BrainProteome_sharma.CellType=vector(mode = "list",length=4)
names(mm_BrainProteome_sharma.CellType)=names(mm_BrainProteome_humanOrtholog)=c("Microglia","Astrocytes","Oligodendrocytes","Neurons")
for (c in 1:4){
  mm_BrainProteome_sharma.CellType[[c]]=c(unlist(strsplit(x = grep(pattern = ";",mm_BrainProteome_sharma$Gene.names.1[which((mm_BrainProteome_sharma[,c+11]=="+")&(mm_BrainProteome_sharma[,5]=="+"))],value = T),split = ";")),grep(pattern = ";",mm_BrainProteome_sharma$Gene.names.1[which((mm_BrainProteome_sharma[,c+11]=="+")&(mm_BrainProteome_sharma[,5]=="+"))],invert = T,value = T))
  mm_BrainProteome_humanOrtholog[[c]]=mouse_human_ortholog$Human.associated.gene.name[which((mouse_human_ortholog$Associated.Gene.Name%in%mm_BrainProteome_sharma.CellType[[c]]))]
}

library(biomaRt)
ensembl_85=useMart(biomart = "ENSEMBL_MART_ENSEMBL",host = "www.ensembl.org",dataset = "hsapiens_gene_ensembl")
tanzi_select_variants=read.xls("../../../../Collaborations/Tanzi_WGS/wking.top1000_samplebased_HET_METAL.xlsx",header=T,sheet=1,as.is=T)
tanzi_select_variants.chr_regions=paste(tanzi_select_variants$MarkerName,"-",as.integer(unlist(lapply(strsplit(x = tanzi_select_variants$MarkerName,split = ":"),`[[`,2)))+1,sep = "")
write(paste("chr",tanzi_select_variants$MarkerName,"-",as.integer(unlist(lapply(strsplit(x = tanzi_select_variants$MarkerName,split = ":"),`[[`,2)))+1,sep = ""),file = "../../../../Collaborations/Tanzi_WGS/Tanzi_HET_METAL_select_variants_chr_regions.txt",sep="\n")
tanzi_select_variants_geneNames=unique(select(x = org.Hs.eg.db,keys = tanzi_select_variants_Refseq$name,keytype = "REFSEQ",columns = "SYMBOL")[,2])

human_genes_orgDB=toTable(org.Hs.egSYMBOL)
overlap_dist=vector(mode = "list",length = 4)
overlap_dist2=vector(mode = "list",length = 1)
names(overlap_dist)=c("lcm_neurons","hipp_proteome","tanzi_select","mm_proteome")
overlap_dist[[1]]=overlap_dist[[2]]=overlap_dist[[3]]=overlap_dist[[4]]=list()
overlap_dist[[4]]=vector(mode = "list",length = 4)
names(overlap_dist[[4]])=names(mm_BrainProteome_humanOrtholog)

for (o in 1:10000){
  hippocampal_proteome_AD.random=sample(x=human_genes_orgDB$symbol,size = length(hippocampal_proteome_AD),replace = F)
  lcm_neurons_AD.random=sample(x=human_genes_orgDB$symbol,size = length(ffpe_neurons_genes),replace = F)
  # mm_proteome_AD.random=sample(x = human_genes_orgDB$symbol,size = length(mm_BrainProteome_humanOrtholog),replace = F)
  tanzi_select_variants.random=sample(x=human_genes_orgDB$symbol,size = length(tanzi_select_variants_geneNames),replace = F)
  
  overlap_dist$hipp_proteome[[o]]=length(intersect(hippocampal_proteome_AD.random,V(msmm_rnaseq.DCN$BM_36$Std)$name))
  overlap_dist$lcm_neurons[[o]]=length(intersect(lcm_neurons_AD.random,V(msmm_rnaseq.DCN$BM_36$Std)$name))
  # overlap_dist$mm_proteome[[o]]=length(intersect(mm_proteome_AD.random,V(msmm_rnaseq.DCN$BM_36$Std)$name))
  overlap_dist$tanzi_select[[o]]=length(intersect(tanzi_select_variants.random,V(msmm_rnaseq.DCN$BM_36$Std)$name))
}
for (c in 1:length(names(mm_BrainProteome_humanOrtholog))){
  for (o in 1:10000){
    mm_proteome.random=sample(x=human_genes_orgDB$symbol,size = length(mm_BrainProteome_humanOrtholog[[c]]),replace = F)
    overlap_dist$mm_proteome[[c]][[o]]=length(intersect(mm_proteome.random,tanzi_select_variants_geneNames))
  }
}


obj1_phg_lcmNeurons=newGeneOverlap(listA = ffpe_neurons_genes,listB = V(msmm_rnaseq.DCN$p001$BM_36$Std)$name,spec = "hg19.gene")
obj2_phg_hippProteome=newGeneOverlap(listA = hippocampal_proteome_AD,listB = V(msmm_rnaseq.DCN$p001$BM_36$Std)$name,spec = "hg19.gene")

obj4_phg_mm_hipp_Proteome=newGeneOverlap(listA = intersect(hippocampal_proteome_AD,mm_BrainProteome_humanOrtholog),listB = V(msmm_rnaseq.DCN$BM_36$Std)$name,spec = "hg19.gene")
obj5_phg_arr_rnaseq=newGeneOverlap(listA = select(x = org.Hs.eg.db,keys = V(msbb_array19.PHG_DCN$Std)$name,keytype = "ENTREZID",columns = "SYMBOL")[,2],listB = V(msmm_rnaseq.DCN$p001$BM_36$Std)$name,spec = "hg19.gene")
obj6_phg_tanzi=newGeneOverlap(listA = unique(select(x = org.Hs.eg.db,keys = tanzi_select_variants_Refseq$name,keytype = "REFSEQ",columns = "SYMBOL")[,2]),listB = V(msmm_rnaseq.DCN$p001$BM_36$Std)$name,spec = "hg19.gene")

obj31_phg_mmProteome=newGeneOverlap(listA = mm_BrainProteome_humanOrtholog$Microglia,listB = V(msmm_rnaseq.DCN$p001$BM_36$Std)$name,spec = "hg19.gene")
obj32_phg_mmProteome=newGeneOverlap(listA = mm_BrainProteome_humanOrtholog$Astrocytes,listB = V(msmm_rnaseq.DCN$p001$BM_36$Std)$name,spec = "hg19.gene")
obj33_phg_mmProteome=newGeneOverlap(listA = mm_BrainProteome_humanOrtholog$Oligodendrocytes,listB = V(msmm_rnaseq.DCN$p001$BM_36$Std)$name,spec = "hg19.gene")
obj34_phg_mmProteome=newGeneOverlap(listA = mm_BrainProteome_humanOrtholog$Neurons,listB = V(msmm_rnaseq.DCN$p001$BM_36$Std)$name,spec = "hg19.gene")

obj41_phg_mmProteome=newGeneOverlap(listA = mm_BrainProteome_humanOrtholog$Microglia,listB = tanzi_select_variants_geneNames,spec = "hg19.gene")
obj42_phg_mmProteome=newGeneOverlap(listA = mm_BrainProteome_humanOrtholog$Astrocytes,listB = tanzi_select_variants_geneNames,spec = "hg19.gene")
obj43_phg_mmProteome=newGeneOverlap(listA = mm_BrainProteome_humanOrtholog$Oligodendrocytes,listB = tanzi_select_variants_geneNames,spec = "hg19.gene")
obj44_phg_mmProteome=newGeneOverlap(listA = mm_BrainProteome_humanOrtholog$Neurons,listB = tanzi_select_variants_geneNames,spec = "hg19.gene")

par(mfrow=c(3,1))
hist(unlist(overlap_dist$lcm_neurons),xlab = "#Overlaps",main = "#Overlap distribution for Drummund etal LCM Neurons, DCN_PHG",col="blue4",xlim = c(0,55),breaks = 50)
abline(v=length(intersect(V(msmm_rnaseq.DCN$p001$BM_36$Std)$name,ffpe_neurons_genes)),col="red4",lwd=2)
text(length(intersect(V(msmm_rnaseq.DCN$p001$BM_36$Std)$name,ffpe_neurons_genes)),0, paste("pval=",testGeneOverlap(obj1_phg_lcmNeurons)@pval,sep=""), col = "black",adj = c(-.1, -.1))
hist(unlist(overlap_dist$hipp_proteome),xlab = "#Overlaps",main = "#Overlap distribution for Hondius etal Hippocampal proteome, DCN_PHG",col="olivedrab3",xlim = c(0,55),breaks = 50)
abline(v=length(intersect(V(msmm_rnaseq.DCN$p001$BM_36$Std)$name,hippocampal_proteome_AD)),col="red4",lwd=2)
text(length(intersect(V(msmm_rnaseq.DCN$p001$BM_36$Std)$name,hippocampal_proteome_AD)),0, paste("pval=",testGeneOverlap(obj2_phg_hippProteome)@pval,sep=""), col = "black",adj = c(-.1, -.1))
# hist(unlist(overlap_dist$mm_proteome),xlab = "#Overlaps",main = "#Overlap distribution for Sharma etal cellular proteome, rewired PHG",col="firebrick1",xlim = c(0,250),breaks = 50)
# abline(v=length(intersect(V(msmm_rnaseq.DCN$BM_36$Std)$name,mm_BrainProteome_humanOrtholog)),col="red4",lwd=2)
# text(length(intersect(V(msmm_rnaseq.DCN$BM_36$Std)$name,mm_BrainProteome_humanOrtholog)),0, paste("pval=",testGeneOverlap(obj3_phg_mmProteome)@pval,sep=""), col = "black",adj = c(-.1, -.1))
hist(unlist(overlap_dist$tanzi_select),xlab = "#Overlaps",main = "#Overlap distribution for Tanzi select 1k",col="orange1",xlim = c(0,15),breaks = 50)
abline(v=length(intersect(V(msmm_rnaseq.DCN$p001$BM_36$Std)$name,tanzi_select_variants_geneNames)),col="red4",lwd=2)
text(length(intersect(V(msmm_rnaseq.DCN$p001$BM_36$Std)$name,tanzi_select_variants_geneNames)),0, paste("pval=",testGeneOverlap(obj6_phg_tanzi)@pval,sep=""), col = "black",adj = c(-.1, -.1))

par(mfrow=c(2,2))
hist(unlist(as.numeric(overlap_dist$mm_proteome$Microglia)),xlab = "#Overlaps",main = "#Overlap distribution Sharma et al, Microglia w/ Tanzi select 1k",col="blue4",xlim = c(0,15),ylim=c(0,2000),breaks = 50)
abline(v=length(intersect(tanzi_select_variants_geneNames,mm_BrainProteome_humanOrtholog$Microglia)),col="red4",lwd=2)
text(length(intersect(tanzi_select_variants_geneNames,mm_BrainProteome_humanOrtholog$Microglia)),0, paste("pval=",round(testGeneOverlap(obj41_phg_mmProteome)@pval,digits = 3),sep=""), col = "black",adj = c(-.1, -.1))
hist(unlist(as.numeric(overlap_dist$mm_proteome$Astrocytes)),xlab = "#Overlaps",main = "#Overlap distribution Sharma et al, Astrocytes w/ Tanzi select 1k",col="olivedrab3",xlim = c(0,15),ylim=c(0,2000),breaks = 50)
abline(v=length(intersect(tanzi_select_variants_geneNames,mm_BrainProteome_humanOrtholog$Astrocytes)),col="red4",lwd=2)
text(length(intersect(tanzi_select_variants_geneNames,mm_BrainProteome_humanOrtholog$Astrocytes)),0, paste("pval=",round(testGeneOverlap(obj42_phg_mmProteome)@pval,digits = 3),sep=""), col = "black",adj = c(-.1, -.1))
hist(unlist(as.numeric(overlap_dist$mm_proteome$Oligodendrocytes)),xlab = "#Overlaps",main = "#Overlap distribution Sharma et al, Oligodendrocytes w/ Tanzi select 1k",col="darkorchid4",xlim = c(0,15),ylim=c(0,2000),breaks = 50)
abline(v=length(intersect(tanzi_select_variants_geneNames,mm_BrainProteome_humanOrtholog$Oligodendrocytes)),col="red4",lwd=2)
text(length(intersect(tanzi_select_variants_geneNames,mm_BrainProteome_humanOrtholog$Oligodendrocytes)),0, paste("pval=",round(testGeneOverlap(obj43_phg_mmProteome)@pval,digits = 3),sep=""), col = "black",adj = c(-.1, -.1))
hist(unlist(as.numeric(overlap_dist$mm_proteome$Neurons)),xlab = "#Overlaps",main = "#Overlap distribution Sharma et al, Neurons w/ Tanzi select 1k",col="orange2",xlim = c(0,15),ylim=c(0,2000),breaks = 50)
abline(v=length(intersect(tanzi_select_variants_geneNames,mm_BrainProteome_humanOrtholog$Neurons)),col="red4",lwd=2)
text(length(intersect(tanzi_select_variants_geneNames,mm_BrainProteome_humanOrtholog$Neurons)),0, paste("pval=",round(testGeneOverlap(obj44_phg_mmProteome)@pval,digits = 3),sep=""), col = "black",adj = c(-.1, -.1))

par(mfrow=c(2,2))
hist(unlist(as.numeric(overlap_dist$mm_proteome$Microglia)),xlab = "#Overlaps",main = "#Overlap distribution Sharma et al, Microglia w/ DCN_PHG",col="blue4",xlim = c(0,100),ylim=c(0,2000),breaks = 50)
abline(v=length(intersect(V(msmm_rnaseq.DCN$p001$BM_36$Std)$name,mm_BrainProteome_humanOrtholog$Microglia)),col="red4",lwd=2)
text(length(intersect(V(msmm_rnaseq.DCN$p001$BM_36$Std)$name,mm_BrainProteome_humanOrtholog$Microglia)),0, paste("pval=",testGeneOverlap(obj31_phg_mmProteome)@pval,sep=""), col = "black",adj = c(-.1, -.1))
hist(unlist(as.numeric(overlap_dist$mm_proteome$Astrocytes)),xlab = "#Overlaps",main = "#Overlap distribution Sharma et al, Astrocytes w/ DCN_PHG",col="olivedrab3",xlim = c(0,100),ylim=c(0,2000),breaks = 50)
abline(v=length(intersect(V(msmm_rnaseq.DCN$p001$BM_36$Std)$name,mm_BrainProteome_humanOrtholog$Astrocytes)),col="red4",lwd=2)
text(length(intersect(V(msmm_rnaseq.DCN$p001$BM_36$Std)$name,mm_BrainProteome_humanOrtholog$Astrocytes)),0, paste("pval=",testGeneOverlap(obj32_phg_mmProteome)@pval,sep=""), col = "black",adj = c(-.1, -.1))
hist(unlist(as.numeric(overlap_dist$mm_proteome$Oligodendrocytes)),xlab = "#Overlaps",main = "#Overlap distribution Sharma et al, Oligodendrocytes w/ DCN_PHG",col="darkorchid4",xlim = c(0,100),ylim=c(0,2000),breaks = 50)
abline(v=length(intersect(V(msmm_rnaseq.DCN$p001$BM_36$Std)$name,mm_BrainProteome_humanOrtholog$Oligodendrocytes)),col="red4",lwd=2)
text(length(intersect(V(msmm_rnaseq.DCN$p001$BM_36$Std)$name,mm_BrainProteome_humanOrtholog$Oligodendrocytes)),0, paste("pval=",testGeneOverlap(obj33_phg_mmProteome)@pval,sep=""), col = "black",adj = c(-.1, -.1))
hist(unlist(as.numeric(overlap_dist$mm_proteome$Neurons)),xlab = "#Overlaps",main = "#Overlap distribution Sharma et al, Neurons w/ DCN_PHG",col="orange2",xlim = c(0,130),ylim=c(0,2000),breaks = 50)
abline(v=length(intersect(V(msmm_rnaseq.DCN$p001$BM_36$Std)$name,mm_BrainProteome_humanOrtholog$Neurons)),col="red4",lwd=2)
text(length(intersect(V(msmm_rnaseq.DCN$p001$BM_36$Std)$name,mm_BrainProteome_humanOrtholog$Neurons)),0, paste("pval=",testGeneOverlap(obj34_phg_mmProteome)@pval,sep=""), col = "black",adj = c(-.1, -.1))


ebc=cluster_edge_betweenness(graph = msmm_rnaseq.DCN_PHG_Filt1.5_Low,directed = F,modularity = T,weights = E(msmm_rnaseq.DCN$BM_36$Std)$weight)
comms=communities(ebc)
indices=names(unlist(lapply(communities(x = ebc),function(x)which(length(x)>=10))))



###################################################################################
msbb_array19.PHG=read.table("../../MSBB_Array19/Normalised_Data/msbb_array19.corr_TrackGenes_PLQ_Genes005_PHG.txt",header = T,sep = "\t",as.is=T)
msbb_array19.PHG_covariates=read.table("../../MSBB_Array19/Normalised_Data/msbb_array19_covariates_PHG.txt",sep="\t",header = T,as.is = T)
msbb_array19.PHG_LowPLQ_samples=paste("X",msbb_array19.PHG_covariates$BrainBank[which(msbb_array19.PHG_covariates$PLQ_Mn<=1)],sep = "")
msbb_array19.PHG_HighPLQ_samples=paste("X",msbb_array19.PHG_covariates$BrainBank[which(msbb_array19.PHG_covariates$PLQ_Mn>=15)],sep = "")
