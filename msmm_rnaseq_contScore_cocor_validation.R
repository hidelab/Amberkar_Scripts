library(org.Hs.eg.db)
library(igraph)
library(cocor)
library(gdata)
library(dataframes2xls)
library(GeneOverlap)

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

msmm_rnaseq.cocor_filtered=vector(mode = "list",length = 2)
names(msmm_rnaseq.cocor_filtered)=c("p005","p001")
msmm_rnaseq.cocor_filtered$p005=msmm_rnaseq.cocor_filtered$p001=vector(mode = "list",length = 3)
names(msmm_rnaseq.cocor_filtered$p005)=names(msmm_rnaseq.cocor_filtered$p001)=names(msmm_rnaseq.cocor)
for (t in 1:length(names(msmm_rnaseq))){
  msmm_rnaseq.cocor_filtered$p005[[t]]=msmm_rnaseq.cocor[[t]][((msmm_rnaseq.cocor[[t]]$p.cocor<=0.05)&(msmm_rnaseq.cocor[[t]]$FDR<=0.1)),]
  msmm_rnaseq.cocor_filtered$p001[[t]]=msmm_rnaseq.cocor[[t]][((msmm_rnaseq.cocor[[t]]$p.cocor<=0.01)&(msmm_rnaseq.cocor[[t]]$FDR<=0.1)),]
  write.table(msmm_rnaseq.cocor_filtered$p005[[t]],file = paste(names(msmm_rnaseq.cocor_filtered$p005)[t],"p005_cocor_filtered_interactions.txt",sep="_"),sep = "\t",col.names = T,row.names = T,quote=F)
  write.table(msmm_rnaseq.cocor_filtered$p001[[t]],file = paste(names(msmm_rnaseq.cocor_filtered$p001)[t],"p001_cocor_filtered_interactions.txt",sep="_"),sep = "\t",col.names = T,row.names = T,quote=F)
}


msmm_rnaseq.DCN=vector(mode = "list",length = 2)
names(msmm_rnaseq.DCN)=c("p005","p001")
msmm_rnaseq.DCN[[1]]=msmm_rnaseq.DCN[[2]]=vector(mode = "list",length = 3)
names(msmm_rnaseq.DCN$p005)=names(msmm_rnaseq.DCN$p001)=names(msmm_rnaseq.cocor)
msmm_rnaseq.DCN$p005$BM_10=msmm_rnaseq.DCN$p005$BM_22=msmm_rnaseq.DCN$p005$BM_36=vector(mode = "list",length = 3)
msmm_rnaseq.DCN$p001$BM_10=msmm_rnaseq.DCN$p001$BM_22=msmm_rnaseq.DCN$p001$BM_36=vector(mode = "list",length = 3)
names(msmm_rnaseq.DCN$p001$BM_10)=names(msmm_rnaseq.DCN$p001$BM_22)=names(msmm_rnaseq.DCN$p001$BM_36)=names(msmm_rnaseq.DCN$p005$BM_10)=names(msmm_rnaseq.DCN$p005$BM_22)=names(msmm_rnaseq.DCN$p005$BM_36)=c("Low","High","Std")

for (t in 1:length(msmm_rnaseq)){
  msmm_rnaseq.DCN$p005[[t]][[1]]=msmm_rnaseq.DCN$p005[[t]][[2]]=msmm_rnaseq.DCN$p005[[t]][[3]]=graph.data.frame(d = msmm_rnaseq.cocor_filtered$p005[[t]][,c(1:2)],directed = F)
  msmm_rnaseq.DCN$p001[[t]][[1]]=msmm_rnaseq.DCN$p001[[t]][[2]]=msmm_rnaseq.DCN$p001[[t]][[3]]=graph.data.frame(d = msmm_rnaseq.cocor_filtered$p001[[t]][,c(1:2)],directed = F)
  E(msmm_rnaseq.DCN$p005[[t]][[1]])$weight=msmm_rnaseq.cocor_filtered$p005[[t]]$r.c+1
  E(msmm_rnaseq.DCN$p005[[t]][[2]])$weight=msmm_rnaseq.cocor_filtered$p005[[t]]$r.t+1
  E(msmm_rnaseq.DCN$p005[[t]][[3]])$weight=msmm_rnaseq.cocor_filtered$p005[[t]]$abs.corr.change
  E(msmm_rnaseq.DCN$p001[[t]][[1]])$weight=msmm_rnaseq.cocor_filtered$p001[[t]]$r.c+1
  E(msmm_rnaseq.DCN$p001[[t]][[2]])$weight=msmm_rnaseq.cocor_filtered$p001[[t]]$r.t+1
  E(msmm_rnaseq.DCN$p001[[t]][[3]])$weight=msmm_rnaseq.cocor_filtered$p001[[t]]$abs.corr.change
  write.graph(msmm_rnaseq.DCN$p001[[t]][[1]],paste(names(msmm_rnaseq.DCN$p001)[t],names(msmm_rnaseq.DCN$p001[[t]])[1],"p001","DiffCoexp","Subgraph.gml",sep = "_"),format = "gml")
  write.graph(msmm_rnaseq.DCN$p001[[t]][[2]],paste(names(msmm_rnaseq.DCN$p001)[t],names(msmm_rnaseq.DCN$p001[[t]])[2],"p001","DiffCoexp","Subgraph.gml",sep = "_"),format = "gml")
  write.graph(msmm_rnaseq.DCN$p001[[t]][[3]],paste(names(msmm_rnaseq.DCN$p001)[t],names(msmm_rnaseq.DCN$p001[[t]])[3],"p001","DiffCoexp","Subgraph.gml",sep = "_"),format = "gml")
  
}
  
  
#Rewiring analysis
g1=graph.data.frame(d = ends(graph = msmm_rnaseq.DCN$p001$BM_36$High,es = which(E(msmm_rnaseq.DCN$p001$BM_36$High)$weight>=1),names = T),directed = F)
g2=graph.data.frame(d = ends(graph = msmm_rnaseq.DCN$p001$BM_36$Low,es = which(E(msmm_rnaseq.DCN$p001$BM_36$Low)$weight>=1),names = T),directed = F)
g3=graph.data.frame(d = ends(graph = msmm_rnaseq.DCN$p001$BM_22$High,es = which(E(msmm_rnaseq.DCN$p001$BM_22$High)$weight>=1),names = T),directed = F)
g4=graph.data.frame(d = ends(graph = msmm_rnaseq.DCN$p001$BM_22$Low,es = which(E(msmm_rnaseq.DCN$p001$BM_22$Low)$weight>=1),names = T),directed = F)
g5=graph.data.frame(d = ends(graph = msmm_rnaseq.DCN$p001$BM_36$High,es = which(E(msmm_rnaseq.DCN$p001$BM_36$High)$weight<=0.5),names = T),directed = F)
g6=graph.data.frame(d = ends(graph = msmm_rnaseq.DCN$p001$BM_36$Low,es = which(E(msmm_rnaseq.DCN$p001$BM_36$Low)$weight<=0.5),names = T),directed = F)
g7=graph.data.frame(d = ends(graph = msmm_rnaseq.DCN$p001$BM_22$High,es = which(E(msmm_rnaseq.DCN$p001$BM_22$High)$weight<=0.5),names = T),directed = F)
g8=graph.data.frame(d = ends(graph = msmm_rnaseq.DCN$p001$BM_22$Low,es = which(E(msmm_rnaseq.DCN$p001$BM_22$Low)$weight<=0.5),names = T),directed = F)

g11=induced_subgraph(graph = g1,vids = intersect(V(g1)$name,V(g2)$name))
g22=induced_subgraph(graph = g2,vids = intersect(V(g1)$name,V(g2)$name))
g33=induced_subgraph(graph = g3,vids = intersect(V(g3)$name,V(g4)$name))
g44=induced_subgraph(graph = g4,vids = intersect(V(g3)$name,V(g4)$name))
g55=induced_subgraph(graph = g5,vids = intersect(V(g5)$name,V(g6)$name))
g66=induced_subgraph(graph = g6,vids = intersect(V(g5)$name,V(g6)$name))
g77=induced_subgraph(graph = g7,vids = intersect(V(g7)$name,V(g8)$name))
g88=induced_subgraph(graph = g8,vids = intersect(V(g7)$name,V(g8)$name))

msmm_rnaseq.DCN_rewired=vector(mode = "list",length = 2)
names(msmm_rnaseq.DCN_rewired)=c("BM_22","BM_36")
msmm_rnaseq.DCN_rewired$BM_22=vector(mode = "list",length = 2)
msmm_rnaseq.DCN_rewired$BM_36=vector(mode = "list",length = 2)
names(msmm_rnaseq.DCN_rewired$BM_22)=names(msmm_rnaseq.DCN_rewired$BM_36)=c("PLQ_PosCorr","PLQ_NegCorr")
msmm_rnaseq.DCN_rewired$BM_22$PLQ_PosCorr=data.frame(Common_Genes=names(degree(graph = g3,v = intersect(V(g3)$name,V(g4)$name))),Low_PLQ=degree(graph = g4,v = intersect(V(g3)$name,V(g4)$name)),High_PLQ=degree(graph = g3,v = intersect(V(g3)$name,V(g4)$name)),Difference=degree(graph = g4,v = intersect(V(g3)$name,V(g4)$name))-degree(graph = g3,v = intersect(V(g3)$name,V(g4)$name)),stringsAsFactors = F)
msmm_rnaseq.DCN_rewired$BM_22$PLQ_NegCorr=data.frame(Common_Genes=names(degree(graph = g7,v = intersect(V(g7)$name,V(g8)$name))),Low_PLQ=degree(graph = g8,v = intersect(V(g7)$name,V(g8)$name)),High_PLQ=degree(graph = g7,v = intersect(V(g7)$name,V(g8)$name)),Difference=degree(graph = g8,v = intersect(V(g7)$name,V(g8)$name))-degree(graph = g7,v = intersect(V(g7)$name,V(g8)$name)),stringsAsFactors = F)
msmm_rnaseq.DCN_rewired$BM_36$PLQ_PosCorr=data.frame(Common_Genes=names(degree(graph = g1,v = intersect(V(g1)$name,V(g2)$name))),Low_PLQ=degree(graph = g2,v = intersect(V(g1)$name,V(g2)$name)),High_PLQ=degree(graph = g1,v = intersect(V(g1)$name,V(g2)$name)),Difference=degree(graph = g2,v = intersect(V(g1)$name,V(g2)$name))-degree(graph = g1,v = intersect(V(g1)$name,V(g2)$name)),stringsAsFactors = F)
msmm_rnaseq.DCN_rewired$BM_36$PLQ_NegCorr=data.frame(Common_Genes=names(degree(graph = g5,v = intersect(V(g5)$name,V(g6)$name))),Low_PLQ=degree(graph = g6,v = intersect(V(g5)$name,V(g6)$name)),High_PLQ=degree(graph = g5,v = intersect(V(g5)$name,V(g6)$name)),Difference=degree(graph = g5,v = intersect(V(g5)$name,V(g6)$name))-degree(graph = g6,v = intersect(V(g5)$name,V(g6)$name)),stringsAsFactors = F)

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
mouse_human_ortholog2=read.table("../../../BrainExpression_Datasets/Mouse_Human_Orthologs_EnsemblGRC37.txt",header = T,sep="\t",as.is=T)
#mm_BrainProteome_sharma=read.xls("../../../BrainExpression_Datasets/Mouse_Brain_Proteome_Sharma_etal/nn.4160-S8.xlsx",sheet = 1,skip=1,header=T,as.is=T)
mm_BrainProteome_IC_files=list.files(path = "../../../BrainExpression_Datasets/MouseBrainProteome_MannLab/IsolatedCells",pattern = "IC",full.names = T)
mm_BrainProteome_CC_files=list.files(path = "../../../BrainExpression_Datasets/MouseBrainProteome_MannLab/CulturedCells",pattern = "CC",full.names = T)
mm_BrainProteome_concRNAseq=read.xls("../../../BrainExpression_Datasets/Mouse_Brain_Proteome_Sharma_etal/nn.4160-S5.xlsx",sheet = 1,header=1,as.is=T)
colnames(mm_BrainProteome_concRNAseq)=c("Oligodendrocytes.Div1.R1","Oligodendrocytes.Div1.R2","Oligodendrocytes.Div1.R3","Oligodendrocytes.div2.5.R1","Oligodendrocytes.div2.5.R2","Oligodendrocytes.div2.5.R3","Oligodendrocytes.div4.R1","Oligodendrocytes.div4.R2","Oligodendrocytes.div4.R3","adult.microglia.R1","adult.microglia.R2","adult.microglia.R3","Astrocytes.R1","Astrocytes.R2","Astrocytes.R3","cortical.neurons.div05.R1","cortical.neurons.div05.R2","cortical.neurons.div05.R3","cortical.neurons.div10.R1","cortical.neurons.div10.R2","cortical.neurons.div10.R3","GeneName")
cond=c(rep("Oligodendrocytes",9),rep("Microglia",3),rep("Astrocytes",3),rep("Neurons",6))
type=c(rep("paired-end",length(cond)))
mm_BrainProteome_concRNAseq.expSetup=data.frame(condition=cond,type=type)
nan_counts=unique(c(which(is.na(mm_BrainProteome_concRNAseq$Oligodendrocytes.Div1.R1)==T),
which(is.na(mm_BrainProteome_concRNAseq$Oligodendrocytes.Div1.R2)==T),
which(is.na(mm_BrainProteome_concRNAseq$Oligodendrocytes.Div1.R3)==T),
which(is.na(mm_BrainProteome_concRNAseq$Oligodendrocytes.div2.5.R1)==T),
which(is.na(mm_BrainProteome_concRNAseq$Oligodendrocytes.div2.5.R2)==T),
which(is.na(mm_BrainProteome_concRNAseq$Oligodendrocytes.div2.5.R3)==T),
which(is.na(mm_BrainProteome_concRNAseq$Oligodendrocytes.div4.R1)==T),
which(is.na(mm_BrainProteome_concRNAseq$Oligodendrocytes.div4.R2)==T),
which(is.na(mm_BrainProteome_concRNAseq$Oligodendrocytes.div4.R3)==T),
which(is.na(mm_BrainProteome_concRNAseq$adult.microglia.R1)==T),
which(is.na(mm_BrainProteome_concRNAseq$adult.microglia.R2)==T),
which(is.na(mm_BrainProteome_concRNAseq$adult.microglia.R3)==T),
which(is.na(mm_BrainProteome_concRNAseq$Astrocytes.R1)==T),
which(is.na(mm_BrainProteome_concRNAseq$Astrocytes.R2)==T),
which(is.na(mm_BrainProteome_concRNAseq$Astrocytes.R3)==T),
which(is.na(mm_BrainProteome_concRNAseq$cortical.neurons.div05.R1)==T),
which(is.na(mm_BrainProteome_concRNAseq$cortical.neurons.div05.R2)==T),
which(is.na(mm_BrainProteome_concRNAseq$cortical.neurons.div05.R3)==T),
which(is.na(mm_BrainProteome_concRNAseq$cortical.neurons.div10.R1)==T),
which(is.na(mm_BrainProteome_concRNAseq$cortical.neurons.div10.R2)==T),
which(is.na(mm_BrainProteome_concRNAseq$cortical.neurons.div10.R3)==T)))
mm_BrainProteome_concRNAseq2=sapply((mm_BrainProteome_concRNAseq[-c(nan_counts,1393,5652),-22]),as.integer)
rownames(mm_BrainProteome_concRNAseq2)=mm_BrainProteome_concRNAseq$GeneName[-c(nan_counts,1393,5652)]
mm_BrainProteome_concRNAseq.dds=DESeqDataSetFromMatrix(countData = mm_BrainProteome_concRNAseq2,colData = mm_BrainProteome_concRNAseq.expSetup,design = ~condition)
mm_BrainProteome_concRNAseq.dds=DESeq(mm_BrainProteome_concRNAseq.dds)
mm_BrainProteome_concRNAseq.counts=data.frame(log2(counts(mm_BrainProteome_concRNAseq.dds,normalized=T)))

mm_BrainProteome_humanOrtholog=vector(mode = "list",length = 4)
mm_BrainProteome_sharma.CellType=vector(mode = "list",length=2)
names(mm_BrainProteome_sharma.CellType)=c("Isolated_Cells","Cultered_Cells")
mm_BrainProteome_sharma.CellType$Isolated_Cells=mm_BrainProteome_sharma.CellType$Cultered_Cells=vector(mode = "list",length=4)
names(mm_BrainProteome_sharma.CellType$Isolated_Cells)=names(mm_BrainProteome_sharma.CellType$Cultered_Cells)=names(mm_BrainProteome_humanOrtholog)=c("Microglia","Astrocytes","Oligodendrocytes","Neurons")
for (c in 1:4){
  ic_df=read.xls(mm_BrainProteome_IC_files[c],sheet=1,skip=1,header=T,as.is=T)
  cc_df=read.xls(mm_BrainProteome_CC_files[c],sheet=1,header=T,as.is=T)
  mm_BrainProteome_sharma.CellType$Isolated_Cells[[c]]=c(grep(pattern = ";",ic_df$Gene.names,value = T,invert = T),unlist(strsplit(x = grep(pattern = ";",ic_df$Gene.names,value = T),split = ";")))
  mm_BrainProteome_sharma.CellType$Cultered_Cells[[c]]=c(grep(pattern = ";",cc_df$Gene.names,value = T,invert = T),unlist(strsplit(x = grep(pattern = ";",cc_df$Gene.names,value = T),split = ";")))
  mm_BrainProteome_humanOrtholog[[c]]=mouse_human_ortholog$Human.associated.gene.name[which(mouse_human_ortholog$Associated.Gene.Name%in%union(mm_BrainProteome_sharma.CellType$Isolated_Cells[[c]],mm_BrainProteome_sharma.CellType$Cultered_Cells[[c]]))]
}

load("~/Work/Data/MSMM-MSBB-HBTRC-Data/MSBB_Array19/Normalised_Data/msbb_corr_NTr_TrackGenes005.RData")
PHG_PP_NTr_genes=select(x = org.Hs.eg.db,keys = msbb_array19.corr_TrackGenes_NTr_Genes005$PHG,keytype = "ENTREZID",columns = "SYMBOL")[,2]
library(biomaRt)
ensembl_85=useMart(biomart = "ENSEMBL_MART_ENSEMBL",host = "www.ensembl.org",dataset = "hsapiens_gene_ensembl")
tanzi_select_variants=read.xls("../../../../Collaborations/Tanzi_WGS/wking.top1000_samplebased_HET_METAL.xlsx",header=T,sheet=1,as.is=T)
tanzi_select_variants.chr_regions=paste(tanzi_select_variants$MarkerName,"-",as.integer(unlist(lapply(strsplit(x = tanzi_select_variants$MarkerName,split = ":"),`[[`,2)))+1,sep = "")
write(paste("chr",tanzi_select_variants$MarkerName,"-",as.integer(unlist(lapply(strsplit(x = tanzi_select_variants$MarkerName,split = ":"),`[[`,2)))+1,sep = ""),file = "../../../../Collaborations/Tanzi_WGS/Tanzi_HET_METAL_select_variants_chr_regions.txt",sep="\n")
tanzi_select_variants_UCSC=read.xls("../../../../Collaborations/Tanzi_WGS/Tanzi_HET_METAL_select_variants_chr_regions_UCSC_Output_hg19.xlsx",sheet = 1,header=T,as.is=T)
tanzi_select_variants_geneNames=unique(tanzi_select_variants_UCSC$name2)

human_genes_orgDB=toTable(org.Hs.egSYMBOL)
overlap_dist=vector(mode = "list",length = 4)
overlap_dist_tanzi_mmProteome_celltype=vector(mode = "list",length = 4)
overlap_dist2=overlap_dist_PP_Ntr=vector(mode = "list",length = 1)
names(overlap_dist)=c("lcm_neurons","hipp_proteome","tanzi_select","mm_proteome")
overlap_dist[[1]]=overlap_dist[[2]]=overlap_dist[[3]]=overlap_dist[[4]]=list()
overlap_dist[[4]]=vector(mode = "list",length = 4)
names(overlap_dist[[4]])=names(overlap_dist_tanzi_mmProteome_celltype)=names(mm_BrainProteome_humanOrtholog)

for (o in 1:10000){
  hippocampal_proteome_AD.random=sample(x=human_genes_orgDB$symbol,size = length(hippocampal_proteome_AD),replace = F)
  lcm_neurons_AD.random=sample(x=human_genes_orgDB$symbol,size = length(ffpe_neurons_genes),replace = F)
  # mm_proteome_AD.random=sample(x = human_genes_orgDB$symbol,size = length(mm_BrainProteome_humanOrtholog),replace = F)
  tanzi_select_variants.random=sample(x=human_genes_orgDB$symbol,size = length(unique(tanzi_select_variants_UCSC$name2)),replace = F)
  
  overlap_dist$hipp_proteome[[o]]=length(intersect(hippocampal_proteome_AD.random,V(msmm_rnaseq.DCN$p001$BM_36$Std)$name))
  overlap_dist$lcm_neurons[[o]]=length(intersect(lcm_neurons_AD.random,V(msmm_rnaseq.DCN$p001$BM_36$Std)$name))
  # overlap_dist$mm_proteome[[o]]=length(intersect(mm_proteome_AD.random,V(msmm_rnaseq.DCN$BM_36$Std)$name))
  overlap_dist$tanzi_select[[o]]=length(intersect(tanzi_select_variants.random,V(msmm_rnaseq.DCN$p001$BM_36$Std)$name))
  PP_NTr_random=sample(x = human_genes_orgDB$symbol,size=length(PHG_PP_NTr_genes))
  overlap_dist_PP_Ntr[[o]]=length(intersect(PP_NTr_random,V(msmm_rnaseq.DCN$p001$BM_36$Std)$name))
}
for (c in 1:length(names(mm_BrainProteome_humanOrtholog))){
  for (o in 1:10000){
    mm_proteome.random=sample(x=human_genes_orgDB$symbol,size = length(mm_BrainProteome_humanOrtholog[[c]]),replace = F)
    overlap_dist$mm_proteome[[c]][[o]]=length(intersect(mm_proteome.random,V(msmm_rnaseq.DCN$p001$BM_36$Std)$name))
  }
}
for (c in 1:length(names(mm_BrainProteome_humanOrtholog))){
  for (o in 1:10000){
    mm_proteome.random=sample(x=human_genes_orgDB$symbol,size = length(mm_BrainProteome_humanOrtholog[[c]]),replace = F)
    overlap_dist_tanzi_mmProteome_celltype[[c]][[o]]=length(intersect(mm_proteome.random,tanzi_select_variants_geneNames))
  }
}

obj1_phg_lcmNeurons=newGeneOverlap(listA = ffpe_neurons_genes,listB = V(msmm_rnaseq.DCN$p001$BM_36$Std)$name,spec = "hg19.gene")
obj2_phg_hippProteome=newGeneOverlap(listA = hippocampal_proteome_AD,listB = V(msmm_rnaseq.DCN$p001$BM_36$Std)$name,spec = "hg19.gene")

obj4_phg_mm_hipp_Proteome=newGeneOverlap(listA = intersect(hippocampal_proteome_AD,mm_BrainProteome_humanOrtholog),listB = V(msmm_rnaseq.DCN$BM_36$Std)$name,spec = "hg19.gene")
obj5_phg_arr_rnaseq=newGeneOverlap(listA = select(x = org.Hs.eg.db,keys = V(msbb_array19.PHG_DCN$Std)$name,keytype = "ENTREZID",columns = "SYMBOL")[,2],listB = V(msmm_rnaseq.DCN$p001$BM_36$Std)$name,spec = "hg19.gene")
obj6_phg_tanzi=newGeneOverlap(listA = tanzi_select_variants_geneNames,listB = V(msmm_rnaseq.DCN$p001$BM_36$Std)$name,spec = "hg19.gene")

obj31_phg_mmProteome=newGeneOverlap(listA = mm_BrainProteome_humanOrtholog$Microglia,listB = V(msmm_rnaseq.DCN$p001$BM_36$Std)$name,spec = "hg19.gene")
obj32_phg_mmProteome=newGeneOverlap(listA = mm_BrainProteome_humanOrtholog$Astrocytes,listB = V(msmm_rnaseq.DCN$p001$BM_36$Std)$name,spec = "hg19.gene")
obj33_phg_mmProteome=newGeneOverlap(listA = mm_BrainProteome_humanOrtholog$Oligodendrocytes,listB = V(msmm_rnaseq.DCN$p001$BM_36$Std)$name,spec = "hg19.gene")
obj34_phg_mmProteome=newGeneOverlap(listA = mm_BrainProteome_humanOrtholog$Neurons,listB = V(msmm_rnaseq.DCN$p001$BM_36$Std)$name,spec = "hg19.gene")

obj41_phg_mmProteome=newGeneOverlap(listA = mm_BrainProteome_humanOrtholog$Microglia,listB = tanzi_select_variants_geneNames,spec = "hg19.gene")
obj42_phg_mmProteome=newGeneOverlap(listA = mm_BrainProteome_humanOrtholog$Astrocytes,listB = tanzi_select_variants_geneNames,spec = "hg19.gene")
obj43_phg_mmProteome=newGeneOverlap(listA = mm_BrainProteome_humanOrtholog$Oligodendrocytes,listB = tanzi_select_variants_geneNames,spec = "hg19.gene")
obj44_phg_mmProteome=newGeneOverlap(listA = mm_BrainProteome_humanOrtholog$Neurons,listB = tanzi_select_variants_geneNames,spec = "hg19.gene")

obj7_NTr_PP=newGeneOverlap(listA = PHG_PP_NTr_genes,listB = V(msmm_rnaseq.DCN$p001$BM_36$Std)$name,spec = "hg19.gene")

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
hist(unlist(as.numeric(overlap_dist_tanzi_mmProteome_celltype$Microglia)),xlab = "#Overlaps",main = "#Overlap distribution Sharma et al, Microglia w/ Tanzi select 1k",col="blue4",xlim = c(0,50),ylim=c(0,2000),breaks = 50)
abline(v=length(intersect(tanzi_select_variants_geneNames,overlap_dist_tanzi_mmProteome_celltype$Microglia)),col="red4",lwd=2)
text(length(intersect(tanzi_select_variants_geneNames,mm_BrainProteome_humanOrtholog$Microglia)),0, paste("pval=",testGeneOverlap(obj41_phg_mmProteome)@pval,sep=""), col = "black",adj = c(-.1, -.1))
hist(unlist(as.numeric(overlap_dist_tanzi_mmProteome_celltype$Astrocytes)),xlab = "#Overlaps",main = "#Overlap distribution Sharma et al, Astrocytes w/ Tanzi select 1k",col="olivedrab3",xlim = c(0,50),ylim=c(0,2000),breaks = 50)
abline(v=length(intersect(tanzi_select_variants_geneNames,mm_BrainProteome_humanOrtholog$Astrocytes)),col="red4",lwd=2)
text(length(intersect(tanzi_select_variants_geneNames,mm_BrainProteome_humanOrtholog$Astrocytes)),0, paste("pval=",testGeneOverlap(obj42_phg_mmProteome)@pval,sep=""), col = "black",adj = c(-.1, -.1))
hist(unlist(as.numeric(overlap_dist_tanzi_mmProteome_celltype$Oligodendrocytes)),xlab = "#Overlaps",main = "#Overlap distribution Sharma et al, Oligodendrocytes w/ Tanzi select 1k",col="darkorchid4",xlim = c(0,50),ylim=c(0,2000),breaks = 50)
abline(v=length(intersect(tanzi_select_variants_geneNames,mm_BrainProteome_humanOrtholog$Oligodendrocytes)),col="red4",lwd=2)
text(length(intersect(tanzi_select_variants_geneNames,mm_BrainProteome_humanOrtholog$Oligodendrocytes)),0, paste("pval=",testGeneOverlap(obj43_phg_mmProteome)@pval,sep=""), col = "black",adj = c(-.1, -.1))
hist(unlist(as.numeric(overlap_dist_tanzi_mmProteome_celltype$Neurons)),xlab = "#Overlaps",main = "#Overlap distribution Sharma et al, Neurons w/ Tanzi select 1k",col="orange2",xlim = c(0,50),ylim=c(0,2000),breaks = 50)
abline(v=length(intersect(tanzi_select_variants_geneNames,mm_BrainProteome_humanOrtholog$Neurons)),col="red4",lwd=2)
text(length(intersect(tanzi_select_variants_geneNames,mm_BrainProteome_humanOrtholog$Neurons)),0, paste("pval=",testGeneOverlap(obj44_phg_mmProteome)@pval,sep=""), col = "black",adj = c(-.1, -.1))


par(mfrow=c(2,2))
hist(unlist(as.numeric(overlap_dist$mm_proteome$Microglia)),xlab = "#Overlaps",main = "#Overlap distribution Sharma et al, Microglia w/ DCN_PHG",col="blue4",xlim = c(0,70),ylim=c(0,2000),breaks = 50)
abline(v=length(intersect(V(msmm_rnaseq.DCN$p001$BM_36$Std)$name,mm_BrainProteome_humanOrtholog$Microglia)),col="red4",lwd=2)
text(length(intersect(V(msmm_rnaseq.DCN$p001$BM_36$Std)$name,mm_BrainProteome_humanOrtholog$Microglia)),0, paste("pval=",testGeneOverlap(obj31_phg_mmProteome)@pval,sep=""), col = "black",adj = c(-.1, -.1))
hist(unlist(as.numeric(overlap_dist$mm_proteome$Astrocytes)),xlab = "#Overlaps",main = "#Overlap distribution Sharma et al, Astrocytes w/ DCN_PHG",col="olivedrab3",xlim = c(0,80),ylim=c(0,2000),breaks = 50)
abline(v=length(intersect(V(msmm_rnaseq.DCN$p001$BM_36$Std)$name,mm_BrainProteome_humanOrtholog$Astrocytes)),col="red4",lwd=2)
text(length(intersect(V(msmm_rnaseq.DCN$p001$BM_36$Std)$name,mm_BrainProteome_humanOrtholog$Astrocytes)),0, paste("pval=",testGeneOverlap(obj32_phg_mmProteome)@pval,sep=""), col = "black",adj = c(-.1, -.1))
hist(unlist(as.numeric(overlap_dist$mm_proteome$Oligodendrocytes)),xlab = "#Overlaps",main = "#Overlap distribution Sharma et al, Oligodendrocytes w/ DCN_PHG",col="darkorchid4",xlim = c(0,140),ylim=c(0,2000),breaks = 50)
abline(v=length(intersect(V(msmm_rnaseq.DCN$p001$BM_36$Std)$name,mm_BrainProteome_humanOrtholog$Oligodendrocytes)),col="red4",lwd=2)
text(length(intersect(V(msmm_rnaseq.DCN$p001$BM_36$Std)$name,mm_BrainProteome_humanOrtholog$Oligodendrocytes)),0, paste("pval=",testGeneOverlap(obj33_phg_mmProteome)@pval,sep=""), col = "black",adj = c(-.1, -.1))
hist(unlist(as.numeric(overlap_dist$mm_proteome$Neurons)),xlab = "#Overlaps",main = "#Overlap distribution Sharma et al, Neurons w/ DCN_PHG",col="orange2",xlim = c(0,50),ylim=c(0,2000),breaks = 50)
abline(v=length(intersect(V(msmm_rnaseq.DCN$p001$BM_36$Std)$name,mm_BrainProteome_humanOrtholog$Neurons)),col="red4",lwd=2)
text(length(intersect(V(msmm_rnaseq.DCN$p001$BM_36$Std)$name,mm_BrainProteome_humanOrtholog$Neurons)),0, paste("pval=",testGeneOverlap(obj34_phg_mmProteome)@pval,sep=""), col = "black",adj = c(-.1, -.1))

par(mfrow=c(1,1))
hist(unlist(as.numeric(overlap_dist_PP_Ntr)),xlab = "#Overlaps",main = "#Overlap distribution with PHG NTr PP",col="orange2",xlim = c(0,100),ylim=c(0,2000),breaks = 50)
abline(v=length(intersect(PHG_PP_NTr_genes,V(msmm_rnaseq.DCN$p001$BM_36$Std)$name)),col="red4",lwd=2)
text(length(intersect(PHG_PP_NTr_genes,V(msmm_rnaseq.DCN$p001$BM_36$Std)$name)),0, paste("pval=",round(testGeneOverlap(obj5_NTr_PP)@pval,digits = 3),sep=""), col = "black",adj = c(-.1, -.1))

msmm_rnaseq.DCN_PHG_CellType_Subnet=lapply(lapply(mm_BrainProteome_humanOrtholog,intersect,V(g11)$name),induced_subgraph,graph=g11)
#msmm_rnaseq.DCN_PHG_CellType_Subnet=lapply(lapply(mm_BrainProteome_humanOrtholog,intersect,V(g22)$name),induced_subgraph,graph=g22)
write.graph(msmm_rnaseq.DCN_PHG_CellType_Subnet$Microglia,"msmm_rnaseq_PHG_DCN_Microglia_Subnet.gml",format = "gml")
write.graph(msmm_rnaseq.DCN_PHG_CellType_Subnet$Astrocytes,"msmm_rnaseq_PHG_DCN_Astrocytes_Subnet.gml",format = "gml")
write.graph(msmm_rnaseq.DCN_PHG_CellType_Subnet$Oligodendrocytes,"msmm_rnaseq_PHG_DCN_Oligodendrocytes_Subnet.gml",format = "gml")
write.graph(msmm_rnaseq.DCN_PHG_CellType_Subnet$Neurons,"msmm_rnaseq_PHG_DCN_Neurons_Subnet.gml",format = "gml")

olp_dist_generic=vector(mode = "list",length=3)
names(olp_dist_generic)=c("FP_STG","STG_PHG","FP_PHG")
for (i in 1:10000){
  fp_random=sample(x=human_genes_orgDB$symbol,size = length(msmm_rnaseq_Plaque.Corr005.Genes$BM_10),replace = F)
  stg_random=sample(x=human_genes_orgDB$symbol,size = length(msmm_rnaseq_Plaque.Corr005.Genes$BM_22),replace = F)
  phg_random=sample(x=human_genes_orgDB$symbol,size = length(msmm_rnaseq_Plaque.Corr005.Genes$BM_36),replace = F)
  olp_dist_generic$FP_STG[[i]]=length(intersect(fp_random,stg_random))
  olp_dist_generic$STG_PHG[[i]]=length(intersect(stg_random,phg_random))
  olp_dist_generic$FP_PHG[[i]]=length(intersect(fp_random,phg_random))
}
obj81=newGeneOverlap(listA = msmm_rnaseq_Plaque.Corr005.Genes$BM_10,listB = msmm_rnaseq_Plaque.Corr005.Genes$BM_22,spec = "hg19.gene")
obj82=newGeneOverlap(listA = msmm_rnaseq_Plaque.Corr005.Genes$BM_22,listB = msmm_rnaseq_Plaque.Corr005.Genes$BM_36,spec = "hg19.gene")
obj83=newGeneOverlap(listA = msmm_rnaseq_Plaque.Corr005.Genes$BM_10,listB = msmm_rnaseq_Plaque.Corr005.Genes$BM_36,spec = "hg19.gene")
L_cluster=vector(mode = "list",length=3)
L_cluster[[1]]=na.omit(select(x=org.Hs.eg.db,keys=msmm_rnaseq_Plaque.Corr005.Genes$BM_36,keytype = "SYMBOL",columns = "ENTREZID")[,2])[1:1261]
L_cluster[[2]]=na.omit(select(x=org.Hs.eg.db,keys=msmm_rnaseq.DCN_rewired$BM_36$Common.Genes,keytype = "SYMBOL",columns = "ENTREZID")[,2])[1:361]
L_cluster[[3]]=na.omit(select(x=org.Hs.eg.db,keys=V(msmm_rnaseq.DCN$p001$BM_36$Std)$name,keytype = "SYMBOL",columns = "ENTREZID")[,2])[1:717]
names(L_cluster)=c("PHG_PLQ_Genes","PHG_RewiredGenes","PHG_DCN")
###################################################################################
msbb_array19.PHG=read.table("../../MSBB_Array19/Normalised_Data/msbb_array19.corr_TrackGenes_PLQ_Genes005_PHG.txt",header = T,sep = "\t",as.is=T)
msbb_array19.PHG_covariates=read.table("../../MSBB_Array19/Normalised_Data/msbb_array19_covariates_PHG.txt",sep="\t",header = T,as.is = T)
msbb_array19.PHG_LowPLQ_samples=paste("X",msbb_array19.PHG_covariates$BrainBank[which(msbb_array19.PHG_covariates$PLQ_Mn<=1)],sep = "")
msbb_array19.PHG_HighPLQ_samples=paste("X",msbb_array19.PHG_covariates$BrainBank[which(msbb_array19.PHG_covariates$PLQ_Mn>=15)],sep = "")
