library(igraph)
library(org.Hs.eg.db)
library(dplyr)
library(data.table)
library(foreach)
library(enrichR)
library(magrittr)

mapIds2<-function(IDs,IDFrom,IDTo){
  require(org.Hs.eg.db)
  idmap=mapIds(x = org.Hs.eg.db,keys = IDs,column = IDTo,keytype = IDFrom,multiVals = "first")
  na_vec=names(idmap[is.na(idmap)==T])
  idmap=idmap[is.na(idmap)==F]
  idmap_df=data.frame("From"=names(idmap),"To"=unlist(unname(idmap)),stringsAsFactors = F)
  return(list(map=idmap_df,noMap=na_vec))
}
jaccard=function(A,B){
  jc=set_cardinality(intersect(A,B))/set_cardinality(union(A,B))
  return(jc)
}
dbs <- listEnrichrDbs()
kegg_dbs=grep(pattern = "KEGG",x = dbs$libraryName,value = T)[3]
biocarta_dbs=grep(pattern = "BioCarta|Panther",x = dbs$libraryName,value = T)[4]
panther_dbs=grep(pattern = "BioCarta|Panther",x = dbs$libraryName,value = T)[5]
go_bp2017=grep(pattern = "GO_*",x = dbs$libraryName,value = T)[9]
go_cc2017=grep(pattern = "GO_*",x = dbs$libraryName,value = T)[7]
go_mf2017=grep(pattern = "GO_*",x = dbs$libraryName,value = T)[8]

# Read variant data
known_AD_genes=scan("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Gene_Lists/Known_AD_genes_V2.txt",what = "char",sep = "\n")
known_AD_genes.cons=scan("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Gene_Lists/Known_AD_genes_ConservativeList.txt",what = "char",sep = "\n")
tanzi_ranked_genes.great=read.table("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Gene_Lists/Ranked_genes_corrected_allele_top_4000_SNPs_genes_from_great.txt",header = T,sep = "\t",stringsAsFactors = F)
tanzi_ranked_regions.great=read.table("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Gene_Lists/Region_genes_top_1000.txt",header = T,sep = "\t",stringsAsFactors = F)
tanzi_ranked_genes.great$SNP_Rank=gsub(pattern = "SNP_rank_",replacement = "",x = tanzi_ranked_genes.great$SNP_Rank)
tanzi_ranked_genes.AD=unique(tanzi_ranked_genes.great%>%mutate(SNP_Rank=as.integer(gsub(pattern = "SNP_rank_",replacement = "",x = SNP_Rank)))%>%filter(Corrected_association=="Case"&SNP_Rank<=1000)%>%pull(Gene))
tanzi_ranked_genes.Control=unique(tanzi_ranked_genes.great%>%mutate(SNP_Rank=as.integer(gsub(pattern = "SNP_rank_",replacement = "",x = SNP_Rank)))%>%filter(Corrected_association=="Control"&SNP_Rank<=1000)%>%pull(Gene))
tanzi_ranked_genes.union=union(tanzi_ranked_genes.AD,tanzi_ranked_genes.Control)
tanzi_ranked_region_genes.union=scan("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Gene_Lists/Tanzi_Rank_Region_SNPgenes.txt",sep = "\n",what = "char")
tableXY_2018=read.table("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Gene_Lists/Genes_TableXY.txt",sep="\t",header = F,as.is = T)
# Read PPI data
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


#Read clusters from ClusterONE outputs
cl1_clusters=read.table("/Users/sandeepamberkar/Work/Data/PPI-Data/Cl1_Clusters/cl1_d06_s5_olp10_kAD_SNPgene_seeded_clusters.txt",sep = ",",header = T,row.names = 1,stringsAsFactors = F)
cl1_clusters.member_genes=strsplit(cl1_clusters$Members,split = " ")
cl1_clusters.size=lapply(strsplit(cl1_clusters$Members,split = " "),length)
cl1_overlap.variants=lapply(strsplit(cl1_clusters$Members,split = " "),function(x)intersect(x,c(known_AD_genes.cons,tanzi_ranked_genes.union)))
exp_variants=lapply(cl1_clusters.size,function(x)round(length(c(known_AD_genes.cons,tanzi_ranked_genes.union))/length(V(hs_ppi)$name)*x,digits = 1))
#exp_region_variants=lapply(cl1_clusters.size,function(x)round(length(tanzi_ranked_regions.great$Gene_GREAT)/20805*x,digits = 1))

names(cl1_overlap.variants)=names(exp_variants)=names(cl1_clusters.member_genes)=names(cl1_clusters.size)=rownames(cl1_clusters)
cl1_overlap.variants=cl1_overlap.variants[lapply(cl1_overlap.variants,length)>1]
cl1_clusters.size=cl1_clusters.size[names(cl1_overlap.variants)]
exp_variants=exp_variants[names(cl1_overlap.variants)]
#exp_region_variants=exp_region_variants[names(cl1_overlap.region_variants)]

# Module enrichment test with SNP genes
fisher_results.SNPs=vector(mode = "list",length = length(cl1_overlap.variants))
names(fisher_results.SNPs)=names(cl1_overlap.variants)
for(m in 1:length(cl1_overlap.variants)){
  fisher_results.SNPs[[m]]=fisher.test(matrix(c(length(cl1_overlap.variants[[m]]),
                                           length(c(known_AD_genes.cons,tanzi_ranked_genes.union))-length(cl1_overlap.variants[[m]]),
                                           cl1_clusters$Size[as.numeric(names(cl1_overlap.variants)[m])]-length(cl1_overlap.variants[[m]]),
                                           length(V(hs_ppi)$name)-length(c(known_AD_genes.cons,tanzi_ranked_genes.union))-length(cl1_overlap.variants[[m]])-cl1_clusters$Size[as.numeric(names(cl1_overlap.variants)[m])]-length(cl1_overlap.variants[[m]])),nrow = 2))
  
}
res.variants=foreach(i=1:length(cl1_overlap.variants))%do%{
  
    
  data.frame(Variants_in_Module=unlist(lapply(cl1_overlap.variants[i],paste,collapse=",")),
             Module_Size=unlist(cl1_clusters.size[i]),
             Expected_Variants=unname(unlist(exp_variants[i])),
             Overlap_Size=unlist(lapply(cl1_overlap.variants[i],length)),
             Fisher_statistic=unname(unlist(lapply(fisher_results.SNPs[i],function(s)s$estimate))),
             pval=unname(unlist(lapply(fisher_results.SNPs[i],function(p)p$p.value))),
             stringsAsFactors = F)
}

cl1_variant_enrichment.SNP=data.frame(rbindlist(res.variants),stringsAsFactors = F)
cl1_variant_enrichment.SNP$adj.p=p.adjust(p = cl1_variant_enrichment.SNP$pval,method = "BH")
rownames(cl1_variant_enrichment.SNP)=names(fisher_results.SNPs)
cl1_variant_enrichment.SNP=cl1_variant_enrichment.SNP%>%rownames_to_column("ModuleNumber")%>%dplyr::filter(adj.p<=0.1)
cl1_variant_enrichment.SNP=cl1_variant_enrichment.SNP%>%mutate(KnownAD_in_module=unlist(lapply(cl1_variant_enrichment.SNP$Variants_in_Module,function(a)paste(intersect(strsplit(x = a,split=",")[[1]],known_AD_genes.cons),collapse = ","))))
cl1_variant_enrichment.SNP=cl1_variant_enrichment.SNP%>%mutate(SNPgenes_in_module=unlist(lapply(cl1_variant_enrichment.SNP$Variants_in_Module,function(a)paste(intersect(strsplit(x = a,split=",")[[1]],tanzi_ranked_genes.union),collapse = ","))))
cl1_variant_enrichment_SNP.KEGG=lapply(cl1_clusters.member_genes[cl1_variant_enrichment.SNP$ModuleNumber],function(x)enrichr(genes = x,databases = kegg_dbs)[[1]]%>%filter(Adjusted.P.value<0.1)%>%pull(Term))
cl1_variant_enrichment_SNP.KEGG[which(unlist(lapply(cl1_variant_enrichment_SNP.KEGG,function(y)length(y)>0))==F)]="No enriched pathways"
cl1_variant_enrichment.SNP$KEGG=unlist(lapply(lapply(cl1_variant_enrichment_SNP.KEGG,intersect,tanzi_ranked_genes.KEGG$Term),function(x)paste(gsub(pattern = "_Homo sapiens_hsa\\d+",replacement = "",x = x),collapse = ";")))
cl1_variant_enrichment.SNP$KEGG[which(cl1_variant_enrichment.SNP$ModuleNumber%in%names(lapply(cl1_variant_enrichment_SNP.KEGG,function(x)intersect(x,tanzi_ranked_genes_KEGG.knownAD_mapped))[lapply(cl1_variant_enrichment_SNP.KEGG,function(x)length(intersect(x,tanzi_ranked_genes_KEGG.knownAD_mapped)))>0]))]=unlist(lapply(cl1_variant_enrichment_SNP.KEGG,function(x)paste(intersect(x,tanzi_ranked_genes_KEGG.knownAD_mapped),collapse = ";"))[lapply(cl1_variant_enrichment_SNP.KEGG,function(x)length(intersect(x,tanzi_ranked_genes_KEGG.knownAD_mapped)))>0])
cl1_variant_enrichment.SNP$KEGG_Module=unlist(lapply(lapply(cl1_clusters.member_genes[cl1_variant_enrichment.SNP$ModuleNumber],function(x)enrichr(genes = x,databases = kegg_dbs)[[1]]%>%filter(Adjusted.P.value<0.1)%>%pull(Term)),function(x)paste(gsub(pattern = "_Homo sapiens_hsa\\d+",replacement = "",x = x),collapse = ";")))
cl1_variant_enrichment.SNP_Mod=read.xls("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Cl1_SeededRank_SNPgene_EnrichedModules.xlsx",sheet = 1,method = "tab")
fwrite(cl1_variant_enrichment.SNP,"/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Cl1_SeededRank_SNPgene_EnrichedModules.txt",sep = "\t",col.names = T,row.names = T,quote = F)
# Check cluster enrichment with regional variants



sig_cluster_SNP.member_genes=cl1_clusters.member_genes[as.numeric(rownames(cl1_variant_enrichment.SNP[cl1_variant_enrichment.SNP$adj.p<=0.1,]))]
#sig_cluster_region_variants.member_genes=cl1_clusters.member_genes[as.numeric(rownames(cl1_variant_enrichment.region_variants[cl1_variant_enrichment.region_variants$pval<=0.05,]))]

sig_cluster_SNP.comp_matrix=matrix(NA,nrow = length(sig_cluster_SNP.member_genes),ncol=length(sig_cluster_SNP.member_genes))
rownames(sig_cluster_SNP.comp_matrix)=colnames(sig_cluster_SNP.comp_matrix)=names(sig_cluster_SNP.member_genes)
for(i in 1:dim(sig_cluster_SNP.comp_matrix)[1]){
  sig_cluster_SNP.comp_matrix[i,]=unlist(lapply(sig_cluster_SNP.member_genes,function(x)jaccard(x,sig_cluster_SNP.member_genes[[i]])))
}




# IGAP enriched pathways
igap_AD_pathways=read.xls("/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/IGAP_AD_Study/1-s2.0-S1552526014024923-mmc2.xls",sheet = 1)


# Read Peptide Atlas data
peptide_atlas=read.table("/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/Peptide-Atlas-HumanBrain-2015.tsv",sep = "\t",header = T,as.is = T,fill = T)
peptide_atlas.proteins=peptide_atlas%>%filter(organism_full_name!="")%>%pull(represented_by_biosequence_id)
peptide_atlas.proteins[grep(pattern = "-",peptide_atlas.proteins)]=unlist(lapply(strsplit(x = grep(pattern = "-",peptide_atlas.proteins,value = T),split = "-"),`[[`,1))
peptide_atlas.proteins_HGNC=mapIds2(IDs = peptide_atlas.proteins,IDFrom = "UNIPROT",IDTo = "SYMBOL")[[1]][,2]



# Build modular network
sig_cluster_SNP_intermod.subnet=graph_from_adjacency_matrix(adjmatrix = sig_cluster_SNP.comp_matrix,mode = "undirected",weighted = T,diag=F)
#sig_cluster_region_intermod.subnet=graph_from_adjacency_matrix(adjmatrix = sig_cluster_region_variants.comp_matrix,mode = "undirected",weighted = T,diag=F)
V(sig_cluster_SNP_intermod.subnet)$colour="#85929E"
V(sig_cluster_SNP_intermod.subnet)$colour[V(sig_cluster_SNP_intermod.subnet)$name%in%intersect(names(sig_cluster_region_variants.KEGG_gp),names(sig_cluster_SNP.KEGG_gp))]="#A569BD"
V(sig_cluster_SNP_intermod.subnet)$size=unlist(unname(lapply(lapply(lapply(sig_cluster_SNP.KEGG_gp,length),function(x)(x*500/sum(unlist(lapply(sig_cluster_SNP.KEGG_gp,length))))),function(y)round(y,digits=1))))
write.graph(graph = sig_cluster_SNP_intermod.subnet,file="/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Cl1_SigCluster_IntermodularNetwork.gml",format = "gml")
sig_cluster_SNP_intermod.cliques=cliques(graph = sig_cluster_SNP_intermod.subnet,min = 3)



lapply(lapply(cl1_clusters.member_genes[sig_cluster_SNP_intermod.cliques[[1]]$name],function(x)enrichr(genes = x,databases = kegg_dbs)[[1]]%>%filter(Adjusted.P.value<0.05)%>%pull(Term)),intersect,tanzi_ranked_genes.KEGG%>%filter(Adjusted.P.value<=0.1)%>%pull(Term))
lapply(lapply(cl1_clusters.member_genes[c("261","12","59","58")],function(x)enrichr(genes = x,databases = kegg_dbs)[[1]]%>%filter(Adjusted.P.value<0.05)%>%pull(Term)),intersect,tanzi_ranked_genes.KEGG%>%filter(Adjusted.P.value<=0.1)%>%pull(Term))

peptide_atlas=read.table("/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/Peptide-Atlas-HumanBrain-2015.tsv",sep = "\t",header = T,as.is = T,fill = T)
peptide_atlas.proteins=peptide_atlas%>%filter(organism_full_name!="")%>%pull(represented_by_biosequence_id)
peptide_atlas.proteins[grep(pattern = "-",peptide_atlas.proteins)]=unlist(lapply(strsplit(x = grep(pattern = "-",peptide_atlas.proteins,value = T),split = "-"),`[[`,1))
peptide_atlas.proteins_HGNC=mapIds2(IDs = peptide_atlas.proteins,IDFrom = "UNIPROT",IDTo = "SYMBOL")[[1]][,2]

brain_celltype=read.table("/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/cell_type_human_500.tsv",header = T,sep = "\t",as.is = T)
zhang_celltype_ADgenes=read.xls('/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/Zhang_19BrainRegions_Paper/Zhang_BrainCelltype_Markers.xlsx',skip=1,sheet=3,header=T,as.is=T)
zhang_tyrobp_modules=read.xls("/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/Zhang_TYROBP/mmc1.xlsx",sheet = 4,as.is=T)
zhang_celltype_ADgenes.list=vector(mode = 'list',length = 5)
names(zhang_celltype_ADgenes.list)=sort(unique(zhang_celltype_ADgenes$Cell.type))
zhang_celltype_ADgenes.list$Astrocytes=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_ADgenes.list)[1],zhang_celltype_ADgenes$Cell.type)]
zhang_celltype_ADgenes.list$Endothelial=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_ADgenes.list)[2],zhang_celltype_ADgenes$Cell.type)]
zhang_celltype_ADgenes.list$Microglia=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_ADgenes.list)[3],zhang_celltype_ADgenes$Cell.type)]
zhang_celltype_ADgenes.list$Neurons=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_ADgenes.list)[4],zhang_celltype_ADgenes$Cell.type)]
zhang_celltype_ADgenes.list$Oligodendrocytes=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_ADgenes.list)[5],zhang_celltype_ADgenes$Cell.type)]

tanzi_celltype_overlap.list=lapply(zhang_celltype_ADgenes.list,intersect,union(tanzi_ranked_genes.AD,tanzi_ranked_genes.Control))
tanzi_celltype_overlap_length.list=lapply(tanzi_celltype_overlap.list,length)
zhang_celltype_total.list=lapply(zhang_celltype_ADgenes.list,length)
tanzi_celltype_exp.list=lapply(zhang_celltype_total.list,function(a)round(a/20805*length(tanzi_ranked_genes.union),digits = 1))

fisher_results.tanzi_celltype=vector(mode = "list",length = length(tanzi_celltype_overlap.list))
names(fisher_results.tanzi_celltype)=names(zhang_celltype_ADgenes.list)
for(m in 1:length(tanzi_celltype_overlap.list)){
  fisher_results.tanzi_celltype[[m]]=fisher.test(matrix(c(tanzi_celltype_overlap_length.list[[m]],
                                                         length(tanzi_ranked_genes.union)-tanzi_celltype_overlap_length.list[[m]],
                                                         zhang_celltype_total.list[[m]]-tanzi_celltype_overlap_length.list[[m]],
                                                         20805-tanzi_celltype_overlap_length.list[[m]]-(length(tanzi_ranked_genes.union)-tanzi_celltype_overlap_length.list[[m]])-(zhang_celltype_total.list[[m]]-tanzi_celltype_overlap_length.list[[m]])),nrow = 2),alternative = "g")
                                                         
  
}
res.tanzi_celltype=foreach(i=1:length(tanzi_celltype_overlap.list))%do%{
  
  data.frame(CellType=names(tanzi_celltype_overlap.list)[i],
             SNPgenes_in_CelltypeMarkers=unlist(lapply(tanzi_celltype_overlap.list[i],paste,collapse=",")),
             Number_of_CellTypeMarkers=unlist(lapply(zhang_celltype_ADgenes.list[i],length)),
             Expected_CellTypeMarkers=unname(unlist(tanzi_celltype_exp.list[i])),
             Overlap_Size=unlist(tanzi_celltype_overlap_length.list[i]),
             Fisher_statistic=unname(unlist(lapply(fisher_results.tanzi_celltype[i],function(s)s$estimate))),
             pval=unname(unlist(lapply(fisher_results.tanzi_celltype[i],function(p)p$p.value))),
             stringsAsFactors = F)
}
res.tanzi_celltype.df=data.frame(rbindlist(res.tanzi_celltype),stringsAsFactors = F)
res.tanzi_celltype.df=res.tanzi_celltype.df%>%mutate(Representation=if_else(Overlap_Size>Expected_CellTypeMarkers,true = "Over",false = "Under"))%>%filter(Representation=="Over")
res.tanzi_celltype.df=res.tanzi_celltype.df%>%mutate(Adj.p=round(p.adjust(pval,method = "fdr"),digits = 5))
res.tanzi_celltype.KEGG=lapply(tanzi_celltype_overlap.list[res.tanzi_celltype.df$CellType],function(x)enrichr(genes = x,databases = kegg_dbs)[[1]]%>%filter(Adjusted.P.value<0.1)%>%pull(Term))
res.tanzi_celltype.df$KEGG=unlist(lapply(lapply(res.tanzi_celltype.KEGG,intersect,tanzi_ranked_genes.KEGG$Term),function(x)paste(gsub(pattern = "_Homo sapiens_hsa\\d+",replacement = "",x = x),collapse = ";")))
fwrite(res.tanzi_celltype.df,"/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Tanzi_RankedSNP_genes_CelltypeEnrichment.txt",sep="\t",col.names = T,row.names = F)

#1stN and direct interactions subnetwork
tanzi_fav_genes_1stN.subgraph=read.sif(sif.file = "/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Tanzi_FavGenes_1stN_Subgraph.sif",format = "igraph",directed = F)
#tanzi_fav_genes_1stN.subgraph=Reduce(graph.union,lapply(neighborhood(graph = hs_ppi,order = 1,nodes = V(hs_ppi)$name[which(V(hs_ppi)$name%in%tanzi_fav_genes)]),function(x)induced.subgraph(graph = hs_ppi,vids = V(hs_ppi)[which(V(hs_ppi)$name%in%x$name)])))
tanzi_fav_genes_knownAD_cons_1stN.subgraph=Reduce(graph.union,lapply(neighborhood(graph = hs_ppi,order = 1,nodes = V(hs_ppi)$name[which(V(hs_ppi)$name%in%c(tanzi_fav_genes,known_AD_genes.cons))]),function(x)induced.subgraph(graph = hs_ppi,vids = V(hs_ppi)[which(V(hs_ppi)$name%in%x$name)])))


tanzi_fav_genes.Seed_df=data.frame(V(hs_ppi)$name,rep(0,length(V(hs_ppi)$name)),stringsAsFactors=F)
tanzi_fav_genes.Seed_df[grep(pattern=paste(c(tanzi_fav_genes,known_AD_genes.cons),collapse="|"),x=tanzi_fav_genes.Seed_df[,1]),2]=1

colnames(tanzi_fav_genes.Seed_df)=c()
tanzi_fav_genes.pData=xPierGenes(data=tanzi_fav_genes.Seed_df,network.customised=hs_ppi,restart=0.9,verbose=T)
tanzi_fav_genes.pSubnet=xPierSubnet(pNode=tanzi_fav_genes.pData,network.customised=hs_ppi,subnet.significance=0.05,verbose=T)

tanzi_fav_genes.pSubnet=simplify(tanzi_fav_genes.pSubnet,remove.multiple = T,remove.loops = T)


g=igraph.to.graphNEL(tanzi_fav_genes.pSubnet)
g=initEdgeAttribute(graph = g,attribute.name = "weight",attribute.type = "numeric",default.value = 0.0)
g=initNodeAttribute(graph = g,attribute.name = "significance",attribute.type = "char",default.value = V(tanzi_fav_genes.pSubnet)$name)
g=initNodeAttribute(graph = g,attribute.name = "significance",attribute.type = "numeric",default.value = as.numeric(V(tanzi_fav_genes.pSubnet)$significance))
g=initNodeAttribute(graph = g,attribute.name = "score",attribute.type = "numeric",default.value = as.numeric(V(tanzi_fav_genes.pSubnet)$score))
g=initNodeAttribute(graph = g,attribute.name = "priority",attribute.type = "numeric",default.value = as.numeric(V(tanzi_fav_genes.pSubnet)$priority))

g=initNodeAttribute(graph = g,attribute.name = "label",attribute.type = "char",default.value = V(tanzi_fav_genes.pSubnet)$name)
cw=CytoscapeWindow(title = "tanzi_fav_genes_1stN_subgraph",graph = g,overwriteWindow = T)

tanzi_fav_genes_1stN.KEGG=enrichr(genes = V(tanzi_fav_genes.pSubnet)$name,databases = kegg_dbs)[[1]]%>%filter(Adjusted.P.value<=0.1)
tanzi_fav_genes_1stN.GO_BP=enrichr(genes = V(tanzi_fav_genes.pSubnet)$name,databases = go_bp2017)[[1]]%>%filter(Adjusted.P.value<=0.1)
tanzi_fav_genes_1stN.GO_CC=enrichr(genes = V(tanzi_fav_genes.pSubnet)$name,databases = go_cc2017)[[1]]%>%filter(Adjusted.P.value<=0.1)
tanzi_fav_genes_1stN.GO_MF=enrichr(genes = V(tanzi_fav_genes.pSubnet)$name,databases = go_mf2017)[[1]]%>%filter(Adjusted.P.value<=0.1)

msigdb_c7=read.gmt(gmtfile = "/Users/sandeepamberkar/Work/Software/MSigDB_Collections/c7.all.v6.1.entrez.gmt")
msigdb_c2=read.gmt(gmtfile = "/Users/sandeepamberkar/Work/Software/MSigDB_Collections/c2.all.v6.1.entrez.gmt")
msigdb_h=read.gmt(gmtfile = "/Users/sandeepamberkar/Work/Software/MSigDB_Collections/h.all.v6.1.entrez.gmt")
msigdb_c5.BP=read.gmt(gmtfile = "/Users/sandeepamberkar/Work/Software/MSigDB_Collections/c5.bp.v6.1.entrez.gmt")
msigdb_c5.CC=read.gmt(gmtfile = "/Users/sandeepamberkar/Work/Software/MSigDB_Collections/c5.cc.v6.1.entrez.gmt")
msigdb_c5.MF=read.gmt(gmtfile = "/Users/sandeepamberkar/Work/Software/MSigDB_Collections/c5.mf.v6.1.entrez.gmt")

tanzi_fav_genes_1stN.Msigdb_C2=enricher(gene = mapIds2(IDs = V(tanzi_fav_genes.pSubnet)$name,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2],TERM2GENE = msigdb_c2,pvalueCutoff = 0.1,pAdjustMethod = "BH")
tanzi_fav_genes_1stN.Msigdb_H=enricher(gene = mapIds2(IDs = V(tanzi_fav_genes.pSubnet)$name,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2],TERM2GENE = msigdb_h,pvalueCutoff = 0.1,pAdjustMethod = "BH")

#IPA output
tanzi_knownAD_XY_IPA.df=read.table("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/IPA-Output/amended.FINAL.TABLEXY.interactions/final.tablexy.regional.snp.consAD.direct.interactions.txt",header = T,sep = "\t",skip = 1,as.is = T)
tanzi_knownAD_XY_IPA.subgraph=graph.data.frame(d = tanzi_knownAD_XY_IPA.df[,c(1,3)],directed = F)
E(tanzi_knownAD_XY_IPA.subgraph)$name=tanzi_knownAD_XY_IPA.df$Relationship.Type
write.graph(tanzi_knownAD_XY_IPA.subgraph,file = "/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/tanzi_knownAD_XY_IPA_subgraph.gml",format = "gml")
V(tanzi_knownAD_XY_IPA.subgraph)$name
tanzi_knownAD_XY_IPA.msigdb_C2=enricher(gene = mapIds2(IDs = V(tanzi_knownAD_XY_IPA.subgraph)$name,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2],pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c2)@result
tanzi_knownAD_XY_IPA.msigdb_C2$GeneRatio=paste(" ",tanzi_knownAD_XY_IPA.msigdb_C2$GeneRatio,sep = "")
tanzi_knownAD_XY_IPA.msigdb_C2$BgRatio=paste(" ",tanzi_knownAD_XY_IPA.msigdb_C2$BgRatio,sep = "")
tanzi_knownAD_XY_IPA.msigdb_C2$genes_HGNC=unlist(lapply(tanzi_knownAD_XY_IPA.msigdb_C2$geneID,function(y)paste(mapIds2(IDs = strsplit(x = y,split = "/")[[1]],IDFrom = "ENTREZID",IDTo = "SYMBOL")[[1]][,2],collapse = ",")))
tanzi_knownAD_XY_IPA.msigdb_C2$knownAD=unlist(lapply(lapply(tanzi_knownAD_XY_IPA.msigdb_C2$genes_HGNC,function(x)strsplit(x = x,split=",")[[1]]),function(y)paste(intersect(y,known_AD_genes.cons),collapse = ",")))
tanzi_knownAD_XY_IPA.msigdb_C2$SAGs=unlist(lapply(lapply(tanzi_knownAD_XY_IPA.msigdb_C2$genes_HGNC,function(x)strsplit(x = x,split=",")[[1]]),function(y)paste(intersect(y,tanzi_ranked_genes.union),collapse = ",")))
tanzi_knownAD_XY_IPA.msigdb_C2$XY=unlist(lapply(lapply(tanzi_knownAD_XY_IPA.msigdb_C2$genes_HGNC,function(x)strsplit(x = x,split=",")[[1]]),function(y)paste(intersect(y,tableXY_GREAT_genes),collapse = ",")))
tanzi_knownAD_XY_IPA.msigdb_C2$Regional=unlist(lapply(lapply(tanzi_knownAD_XY_IPA.msigdb_C2$genes_HGNC,function(x)strsplit(x = x,split=",")[[1]]),function(y)paste(intersect(y,tanzi_region_genes),collapse = ",")))
fwrite(tanzi_knownAD_XY_IPA.msigdb_C2,"/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/IPA-Output/amended.FINAL.TABLEXY.interactions/tanzi_knownAD_XY_IPA.msigdb_C2.txt",sep = "\t",col.names = T,row.names = F,quote = F)

tanzi_knownAD_XY_IPA.msigdb_C7=enricher(gene = mapIds2(IDs = V(tanzi_knownAD_XY_IPA.subgraph)$name,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2],pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c7)@result
tanzi_knownAD_XY_IPA.msigdb_C7$GeneRatio=paste(" ",tanzi_knownAD_XY_IPA.msigdb_C7$GeneRatio,sep = "")
tanzi_knownAD_XY_IPA.msigdb_C7$BgRatio=paste(" ",tanzi_knownAD_XY_IPA.msigdb_C7$BgRatio,sep = "")
tanzi_knownAD_XY_IPA.msigdb_C7$genes_HGNC=unlist(lapply(tanzi_knownAD_XY_IPA.msigdb_C7$geneID,function(y)paste(mapIds2(IDs = strsplit(x = y,split = "/")[[1]],IDFrom = "ENTREZID",IDTo = "SYMBOL")[[1]][,2],collapse = ",")))
tanzi_knownAD_XY_IPA.msigdb_C7$knownAD=unlist(lapply(lapply(tanzi_knownAD_XY_IPA.msigdb_C7$genes_HGNC,function(x)strsplit(x = x,split=",")[[1]]),function(y)paste(intersect(y,known_AD_genes.cons),collapse = ",")))
tanzi_knownAD_XY_IPA.msigdb_C7$SAGs=unlist(lapply(lapply(tanzi_knownAD_XY_IPA.msigdb_C7$genes_HGNC,function(x)strsplit(x = x,split=",")[[1]]),function(y)paste(intersect(y,tanzi_ranked_genes.union),collapse = ",")))
tanzi_knownAD_XY_IPA.msigdb_C7$XY=unlist(lapply(lapply(tanzi_knownAD_XY_IPA.msigdb_C7$genes_HGNC,function(x)strsplit(x = x,split=",")[[1]]),function(y)paste(intersect(y,tableXY_GREAT_genes),collapse = ",")))
tanzi_knownAD_XY_IPA.msigdb_C7$Regional=unlist(lapply(lapply(tanzi_knownAD_XY_IPA.msigdb_C7$genes_HGNC,function(x)strsplit(x = x,split=",")[[1]]),function(y)paste(intersect(y,tanzi_region_genes),collapse = ",")))
fwrite(tanzi_knownAD_XY_IPA.msigdb_C7,"/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/IPA-Output/amended.FINAL.TABLEXY.interactions/tanzi_knownAD_XY_IPA.msigdb_C7.txt",sep = "\t",col.names = T,row.names = F,quote = F)

tanzi_knownAD_XY_IPA.msigdb_BP=enricher(gene = mapIds2(IDs = V(tanzi_knownAD_XY_IPA.subgraph)$name,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2],pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c5.BP)@result
tanzi_knownAD_XY_IPA.msigdb_BP$GeneRatio=paste(" ",tanzi_knownAD_XY_IPA.msigdb_BP$GeneRatio,sep = "")
tanzi_knownAD_XY_IPA.msigdb_BP$BgRatio=paste(" ",tanzi_knownAD_XY_IPA.msigdb_BP$BgRatio,sep = "")
tanzi_knownAD_XY_IPA.msigdb_BP$genes_HGNC=unlist(lapply(tanzi_knownAD_XY_IPA.msigdb_BP$geneID,function(y)paste(mapIds2(IDs = strsplit(x = y,split = "/")[[1]],IDFrom = "ENTREZID",IDTo = "SYMBOL")[[1]][,2],collapse = ",")))
tanzi_knownAD_XY_IPA.msigdb_BP$knownAD=unlist(lapply(lapply(tanzi_knownAD_XY_IPA.msigdb_BP$genes_HGNC,function(x)strsplit(x = x,split=",")[[1]]),function(y)paste(intersect(y,known_AD_genes.cons),collapse = ",")))
tanzi_knownAD_XY_IPA.msigdb_BP$SAGs=unlist(lapply(lapply(tanzi_knownAD_XY_IPA.msigdb_BP$genes_HGNC,function(x)strsplit(x = x,split=",")[[1]]),function(y)paste(intersect(y,tanzi_ranked_genes.union),collapse = ",")))
tanzi_knownAD_XY_IPA.msigdb_BP$XY=unlist(lapply(lapply(tanzi_knownAD_XY_IPA.msigdb_BP$genes_HGNC,function(x)strsplit(x = x,split=",")[[1]]),function(y)paste(intersect(y,tableXY_GREAT_genes),collapse = ",")))
tanzi_knownAD_XY_IPA.msigdb_BP$Regional=unlist(lapply(lapply(tanzi_knownAD_XY_IPA.msigdb_BP$genes_HGNC,function(x)strsplit(x = x,split=",")[[1]]),function(y)paste(intersect(y,tanzi_region_genes),collapse = ",")))
fwrite(tanzi_knownAD_XY_IPA.msigdb_BP,"/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/IPA-Output/amended.FINAL.TABLEXY.interactions/tanzi_knownAD_XY_IPA.msigdb_BP.txt",sep = "\t",col.names = T,row.names = F,quote = F)

tanzi_knownAD_XY_IPA.msigdb_CC=enricher(gene = mapIds2(IDs = V(tanzi_knownAD_XY_IPA.subgraph)$name,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2],pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c5.CC)@result
tanzi_knownAD_XY_IPA.msigdb_CC$GeneRatio=paste(" ",tanzi_knownAD_XY_IPA.msigdb_CC$GeneRatio,sep = "")
tanzi_knownAD_XY_IPA.msigdb_CC$BgRatio=paste(" ",tanzi_knownAD_XY_IPA.msigdb_CC$BgRatio,sep = "")
tanzi_knownAD_XY_IPA.msigdb_CC$genes_HGNC=unlist(lapply(tanzi_knownAD_XY_IPA.msigdb_CC$geneID,function(y)paste(mapIds2(IDs = strsplit(x = y,split = "/")[[1]],IDFrom = "ENTREZID",IDTo = "SYMBOL")[[1]][,2],collapse = ",")))
tanzi_knownAD_XY_IPA.msigdb_CC$knownAD=unlist(lapply(lapply(tanzi_knownAD_XY_IPA.msigdb_CC$genes_HGNC,function(x)strsplit(x = x,split=",")[[1]]),function(y)paste(intersect(y,known_AD_genes.cons),collapse = ",")))
tanzi_knownAD_XY_IPA.msigdb_CC$SAGs=unlist(lapply(lapply(tanzi_knownAD_XY_IPA.msigdb_CC$genes_HGNC,function(x)strsplit(x = x,split=",")[[1]]),function(y)paste(intersect(y,tanzi_ranked_genes.union),collapse = ",")))
tanzi_knownAD_XY_IPA.msigdb_CC$XY=unlist(lapply(lapply(tanzi_knownAD_XY_IPA.msigdb_CC$genes_HGNC,function(x)strsplit(x = x,split=",")[[1]]),function(y)paste(intersect(y,tableXY_GREAT_genes),collapse = ",")))
tanzi_knownAD_XY_IPA.msigdb_CC$Regional=unlist(lapply(lapply(tanzi_knownAD_XY_IPA.msigdb_CC$genes_HGNC,function(x)strsplit(x = x,split=",")[[1]]),function(y)paste(intersect(y,tanzi_region_genes),collapse = ",")))
fwrite(tanzi_knownAD_XY_IPA.msigdb_CC,"/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/IPA-Output/amended.FINAL.TABLEXY.interactions/tanzi_knownAD_XY_IPA.msigdb_CC.txt",sep = "\t",col.names = T,row.names = F,quote = F)

tanzi_knownAD_XY_IPA.msigdb_MF=enricher(gene = mapIds2(IDs = V(tanzi_knownAD_XY_IPA.subgraph)$name,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2],pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c5.MF)@result
tanzi_knownAD_XY_IPA.msigdb_MF$GeneRatio=paste(" ",tanzi_knownAD_XY_IPA.msigdb_MF$GeneRatio,sep = "")
tanzi_knownAD_XY_IPA.msigdb_MF$BgRatio=paste(" ",tanzi_knownAD_XY_IPA.msigdb_MF$BgRatio,sep = "")
tanzi_knownAD_XY_IPA.msigdb_MF$genes_HGNC=unlist(lapply(tanzi_knownAD_XY_IPA.msigdb_MF$geneID,function(y)paste(mapIds2(IDs = strsplit(x = y,split = "/")[[1]],IDFrom = "ENTREZID",IDTo = "SYMBOL")[[1]][,2],collapse = ",")))
tanzi_knownAD_XY_IPA.msigdb_MF$knownAD=unlist(lapply(lapply(tanzi_knownAD_XY_IPA.msigdb_MF$genes_HGNC,function(x)strsplit(x = x,split=",")[[1]]),function(y)paste(intersect(y,known_AD_genes.cons),collapse = ",")))
tanzi_knownAD_XY_IPA.msigdb_MF$SAGs=unlist(lapply(lapply(tanzi_knownAD_XY_IPA.msigdb_MF$genes_HGNC,function(x)strsplit(x = x,split=",")[[1]]),function(y)paste(intersect(y,tanzi_ranked_genes.union),collapse = ",")))
tanzi_knownAD_XY_IPA.msigdb_MF$XY=unlist(lapply(lapply(tanzi_knownAD_XY_IPA.msigdb_MF$genes_HGNC,function(x)strsplit(x = x,split=",")[[1]]),function(y)paste(intersect(y,tableXY_GREAT_genes),collapse = ",")))
tanzi_knownAD_XY_IPA.msigdb_MF$Regional=unlist(lapply(lapply(tanzi_knownAD_XY_IPA.msigdb_MF$genes_HGNC,function(x)strsplit(x = x,split=",")[[1]]),function(y)paste(intersect(y,tanzi_region_genes),collapse = ",")))
fwrite(tanzi_knownAD_XY_IPA.msigdb_MF,"/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/IPA-Output/amended.FINAL.TABLEXY.interactions/tanzi_knownAD_XY_IPA.msigdb_MF.txt",sep = "\t",col.names = T,row.names = F,quote = F)

tanzi_knownAD_XY_IPA.H=enricher(gene = mapIds2(IDs = V(tanzi_knownAD_XY_IPA.subgraph)$name,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2],pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_h)@result
tanzi_knownAD_XY_IPA.H$GeneRatio=paste(" ",tanzi_knownAD_XY_IPA.H$GeneRatio,sep = "")
tanzi_knownAD_XY_IPA.H$BgRatio=paste(" ",tanzi_knownAD_XY_IPA.H$BgRatio,sep = "")
tanzi_knownAD_XY_IPA.H$genes_HGNC=unlist(lapply(tanzi_knownAD_XY_IPA.H$geneID,function(y)paste(mapIds2(IDs = strsplit(x = y,split = "/")[[1]],IDFrom = "ENTREZID",IDTo = "SYMBOL")[[1]][,2],collapse = ",")))
tanzi_knownAD_XY_IPA.H$knownAD=unlist(lapply(lapply(tanzi_knownAD_XY_IPA.H$genes_HGNC,function(x)strsplit(x = x,split=",")[[1]]),function(y)paste(intersect(y,known_AD_genes.cons),collapse = ",")))
tanzi_knownAD_XY_IPA.H$SAGs=unlist(lapply(lapply(tanzi_knownAD_XY_IPA.H$genes_HGNC,function(x)strsplit(x = x,split=",")[[1]]),function(y)paste(intersect(y,tanzi_ranked_genes.union),collapse = ",")))
tanzi_knownAD_XY_IPA.H$XY=unlist(lapply(lapply(tanzi_knownAD_XY_IPA.H$genes_HGNC,function(x)strsplit(x = x,split=",")[[1]]),function(y)paste(intersect(y,tableXY_GREAT_genes),collapse = ",")))
tanzi_knownAD_XY_IPA.H$Regional=unlist(lapply(lapply(tanzi_knownAD_XY_IPA.H$genes_HGNC,function(x)strsplit(x = x,split=",")[[1]]),function(y)paste(intersect(y,tanzi_region_genes),collapse = ",")))
fwrite(tanzi_knownAD_XY_IPA.H,"/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/IPA-Output/amended.FINAL.TABLEXY.interactions/tanzi_knownAD_XY_IPA_msigdb_H.txt",sep = "\t",col.names = T,row.names = F,quote = F)

g1=igraph.to.graphNEL(tanzi_knownAD_XY_IPA.subgraph)
#g1=initEdgeAttribute(graph = g,attribute.name = "weight",attribute.type = "numeric",default.value = 0.0)
# g1=initNodeAttribute(graph = g,attribute.name = "significance",attribute.type = "char",default.value = V(tanzi_fav_genes.pSubnet)$name)
# g1=initNodeAttribute(graph = g,attribute.name = "significance",attribute.type = "numeric",default.value = as.numeric(V(tanzi_fav_genes.pSubnet)$significance))
# g1=initNodeAttribute(graph = g,attribute.name = "score",attribute.type = "numeric",default.value = as.numeric(V(tanzi_fav_genes.pSubnet)$score))
#g1=initNodeAttribute(graph = g,attribute.name = "p",attribute.type = "numeric",default.value = as.numeric(V(tanzi_fav_genes.pSubnet)$priority))

g1=initNodeAttribute(graph = g1,attribute.name = "label",attribute.type = "char",default.value = V(tanzi_knownAD_XY_IPA.subgraph)$name)
cw=CytoscapeWindow(title = "tanzi_fav_genes_1stN_subgraph",graph = g,overwriteWindow = T)

tanzi_knownAD_XY_IPA.KEGG=enrichr(genes = V(tanzi_knownAD_XY_IPA.subgraph)$name,databases = kegg_dbs)[[1]]%>%filter(Adjusted.P.value<0.1)
tanzi_knownAD_XY_IPA.KEGG$knownAD=unlist(lapply(lapply(tanzi_knownAD_XY_IPA.KEGG$Genes,function(x)strsplit(x = x,split=";")[[1]]),function(y)paste(intersect(y,known_AD_genes.cons),collapse = ",")))
tanzi_knownAD_XY_IPA.KEGG$SAGs=unlist(lapply(lapply(tanzi_knownAD_XY_IPA.KEGG$Genes,function(x)strsplit(x = x,split=";")[[1]]),function(y)paste(intersect(y,tanzi_ranked_genes.union),collapse = ",")))
tanzi_knownAD_XY_IPA.KEGG$XY=unlist(lapply(lapply(tanzi_knownAD_XY_IPA.KEGG$Genes,function(x)strsplit(x = x,split=";")[[1]]),function(y)paste(intersect(y,tableXY_GREAT_genes),collapse = ",")))
tanzi_knownAD_XY_IPA.KEGG$Regional=unlist(lapply(lapply(tanzi_knownAD_XY_IPA.KEGG$Genes,function(x)strsplit(x = x,split=";")[[1]]),function(y)paste(intersect(y,tanzi_region_genes),collapse = ",")))
tanzi_knownAD_XY_IPA.KEGG$Pathway=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = tanzi_knownAD_XY_IPA.KEGG$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[1]]))
tanzi_knownAD_XY_IPA.KEGG$Overlap=paste(" ",tanzi_knownAD_XY_IPA.KEGG$Overlap,sep = "")
tanzi_knownAD_XY_IPA.KEGG$KEGGID=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = tanzi_knownAD_XY_IPA.KEGG$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[2]]))


tanzi_knownAD_XY_IPA.Panther=enrichr(genes = V(tanzi_knownAD_XY_IPA.subgraph)$name,databases = panther_dbs)[[1]]%>%filter(Adjusted.P.value<0.1)
tanzi_knownAD_XY_IPA.Panther$knownAD=unlist(lapply(lapply(tanzi_knownAD_XY_IPA.Panther$Genes,function(x)strsplit(x = x,split=";")[[1]]),function(y)paste(intersect(y,known_AD_genes.cons),collapse = ",")))
tanzi_knownAD_XY_IPA.Panther$SAGs=unlist(lapply(lapply(tanzi_knownAD_XY_IPA.Panther$Genes,function(x)strsplit(x = x,split=";")[[1]]),function(y)paste(intersect(y,tanzi_ranked_genes.union),collapse = ",")))
tanzi_knownAD_XY_IPA.Panther$XY=unlist(lapply(lapply(tanzi_knownAD_XY_IPA.Panther$Genes,function(x)strsplit(x = x,split=";")[[1]]),function(y)paste(intersect(y,tableXY_GREAT_genes),collapse = ",")))
tanzi_knownAD_XY_IPA.Panther$Regional=unlist(lapply(lapply(tanzi_knownAD_XY_IPA.Panther$Genes,function(x)strsplit(x = x,split=";")[[1]]),function(y)paste(intersect(y,tanzi_region_genes),collapse = ",")))
tanzi_knownAD_XY_IPA.Panther$Pathway=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = tanzi_knownAD_XY_IPA.Panther$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[1]]))
tanzi_knownAD_XY_IPA.Panther$Overlap=paste(" ",tanzi_knownAD_XY_IPA.Panther$Overlap,sep = "")
tanzi_knownAD_XY_IPA.Panther$PantherID=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = tanzi_knownAD_XY_IPA.Panther$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[2]]))

tanzi_knownAD_XY_IPA.Biocarta=enrichr(genes = V(tanzi_knownAD_XY_IPA.subgraph)$name,databases = biocarta_dbs)[[1]]%>%filter(Adjusted.P.value<0.1)
tanzi_knownAD_XY_IPA.Biocarta$knownAD=unlist(lapply(lapply(tanzi_knownAD_XY_IPA.Biocarta$Genes,function(x)strsplit(x = x,split=";")[[1]]),function(y)paste(intersect(y,known_AD_genes.cons),collapse = ",")))
tanzi_knownAD_XY_IPA.Biocarta$SAGs=unlist(lapply(lapply(tanzi_knownAD_XY_IPA.Biocarta$Genes,function(x)strsplit(x = x,split=";")[[1]]),function(y)paste(intersect(y,tanzi_ranked_genes.union),collapse = ",")))
tanzi_knownAD_XY_IPA.Biocarta$XY=unlist(lapply(lapply(tanzi_knownAD_XY_IPA.Biocarta$Genes,function(x)strsplit(x = x,split=";")[[1]]),function(y)paste(intersect(y,tableXY_GREAT_genes),collapse = ",")))
tanzi_knownAD_XY_IPA.Biocarta$Regional=unlist(lapply(lapply(tanzi_knownAD_XY_IPA.Biocarta$Genes,function(x)strsplit(x = x,split=";")[[1]]),function(y)paste(intersect(y,tanzi_region_genes),collapse = ",")))
tanzi_knownAD_XY_IPA.Biocarta$Pathway=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = tanzi_knownAD_XY_IPA.Biocarta$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[1]]))
tanzi_knownAD_XY_IPA.Biocarta$Overlap=paste(" ",tanzi_knownAD_XY_IPA.Biocarta$Overlap,sep = "")
tanzi_knownAD_XY_IPA.Biocarta$BiocartaID=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = tanzi_knownAD_XY_IPA.Biocarta$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[3]]))

write.table(tanzi_knownAD_XY_IPA.KEGG%>%dplyr::select(Pathway,KEGGID,Overlap,P.value,Adjusted.P.value,Genes,knownAD,SAGs,Regional,XY),"/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/IPA-Output/amended.FINAL.TABLEXY.interactions/tanzi_knownAD_XY_IPA_enrichedKEGG_Enrichr.txt",sep = "\t",col.names = T,row.names = F)
write.table(tanzi_knownAD_XY_IPA.Biocarta%>%dplyr::select(Pathway,BiocartaID,Overlap,P.value,Adjusted.P.value,Genes,knownAD,SAGs,Regional,XY),"/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/IPA-Output/amended.FINAL.TABLEXY.interactions/tanzi_knownAD_XY_IPA_enrichedBiocarta_Enrichr.txt",sep = "\t",col.names = T,row.names = F)
write.table(tanzi_knownAD_XY_IPA.Panther%>%dplyr::select(Pathway,PantherID,Overlap,P.value,Adjusted.P.value,Genes,knownAD,SAGs,Regional,XY),"/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/IPA-Output/amended.FINAL.TABLEXY.interactions/tanzi_knownAD_XY_IPA_enrichedPanther_Enrichr.txt",sep = "\t",col.names = T,row.names = F)
