library(org.Hs.eg.db)
library(igraph)
library(data.table)
#library(centiserve)
library(sets)
#library(gProfileR)
library(dplyr)
library(gdata)
library(Pi)
library(limma)
#library(RDAVIDWebService)
#library(FGNet)
library(enrichR)
library(tibble)

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
setwd("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/")
dbs <- listEnrichrDbs()
kegg_dbs=grep(pattern = "KEGG",x = dbs$libraryName,value = T)[3]
biocarta_dbs=grep(pattern = "BioCarta|Panther",x = dbs$libraryName,value = T)[4]
panther_dbs=grep(pattern = "BioCarta|Panther",x = dbs$libraryName,value = T)[5]

#Read curated lists from WGS-AD
known_AD_genes=scan("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Gene_Lists/Known_AD_genes_2018_Final.txt",what = "char",sep = "\n")
known_AD_genes.LD=scan("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Gene_Lists/AD_Genes_LD.txt",what = "char",sep = "\n")
known_AD_genes.cons=scan("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Gene_Lists/Known_AD_genes_ConservativeList.txt",what = "char",sep = "\n")
names(known_AD_genes.cons)=mapIds2(IDs = known_AD_genes.cons,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2]
tanzi_ranked_genes=read.table("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Gene_Lists/All_genes_ranked.txt",header = T,sep = "\t",stringsAsFactors = F)
tanzi_ranked_genes.great=read.table("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Gene_Lists/Ranked_genes_corrected_allele_top_4000_SNPs_genes_from_great.txt",header = T,sep = "\t",stringsAsFactors = F)
tanzi_ranked_regions.great=read.table("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Gene_Lists/Region_genes_top_1000.txt",header = T,sep = "\t",stringsAsFactors = F)
tanzi_ranked_regions.great$Region_Rank=as.numeric(gsub(pattern = "Region_Rank_",replacement = "",x = tanzi_ranked_regions.great$Region_Rank))
tanzi_ranked_genes.great$SNP_Rank=as.numeric(gsub(pattern = "SNP_rank_",replacement = "",x = tanzi_ranked_genes.great$SNP_Rank))
tanzi_variant_neighbours=read.table("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/table4.region_based.rare.TOP_and_NOMINAL_with_sv_p005.txt",header = T,sep = "\t",stringsAsFactors = F)
tanzi_fav_genes=scan("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Gene_Lists/Tanzi_FavGenes_2018.txt",sep = "\n",what = "char")
tableXY_top_regional=scan("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Gene_Lists/TableXY_genes_GREAT_SNP_region.txt",what = "char",sep="\n")
tableXY_top_SNP=scan("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Gene_Lists/tableXY_top_SNPs.txt",what = "char",sep="\n")

tableXY_genes_great_snp_region=read.table("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Gene_Lists/TableXY_genes_GREAT_SNP_region.txt",sep = "\t",header = T,as.is = T)%>%select(c(gene1,gene2))
tableXY_genes_great_snp_region.vector=grep(pattern = "\\.",x = as.vector(as.matrix(tableXY_genes_great_snp_region[,c(1:2)])),value = T,invert = T)
tableXY_GREAT_genes=scan("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Gene_Lists/Genes_TableXY.txt",what="char",sep = "\t")
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
write.table(data.frame(as_edgelist(hs_ppi),stringsAsFactors = F),"/Users/sandeepamberkar/Work/Data/PPI-Data/iref14_Human_UP_noDup_table_HGNC_forMCL.txt",sep = "\t",col.names = F,row.names=F,quote=F)

tanzi_ranked_genes.AD=unique(tanzi_ranked_genes.great%>%mutate(SNP_Rank=as.integer(gsub(pattern = "SNP_rank_",replacement = "",x = SNP_Rank)))%>%filter(Corrected_association=="Case"&SNP_Rank<=1000)%>%pull(Gene))
tanzi_ranked_genes.Control=unique(tanzi_ranked_genes.great%>%mutate(SNP_Rank=as.integer(gsub(pattern = "SNP_rank_",replacement = "",x = SNP_Rank)))%>%filter(Corrected_association=="Control"&SNP_Rank<=1000)%>%pull(Gene))
tanzi_region_genes=tanzi_ranked_regions.great%>%filter(Region_Rank<=800)%>%pull(Gene_GREAT)%>%unique
tanzi_random_ranked=sample(x = toTable(org.Hs.egSYMBOL)[,2],size = length(union(tanzi_ranked_genes.AD,tanzi_ranked_genes.Control)))


tanzi_ranked_genes.Seed_df=data.frame(V(hs_ppi)$name,rep(0,length(V(hs_ppi)$name)),stringsAsFactors=F)
tanzi_ranked_genes.Seed_df[grep(pattern=paste(union(tanzi_ranked_genes.Control,tanzi_ranked_genes.AD),collapse="|"),x=tanzi_ranked_genes.Seed_df[,1]),2]=1
tanzi_ranked_genes.Seed_df[grep(pattern=paste(union(tanzi_ranked_genes.Control,tanzi_ranked_genes.AD),collapse="|"),x=tanzi_ranked_genes.Seed_df[,1]),2]=1
colnames(tanzi_ranked_genes.Seed_df)=c()
tanzi_ranked_genes.pData=xPierGenes(data=tanzi_ranked_genes.Seed_df,network.customised=hs_ppi,restart=0.9,verbose=T)


tanzi_ranked_genes.pSubnet=xPierSubnet(pNode=tanzi_ranked_genes.pData,network.customised=hs_ppi,subnet.significance=0.05,verbose=T)
tanzi_ranked_genes_AD.pSubnet=xPierSubnet(pNode=tanzi_ranked_genes.pData,network.customised=hs_ppi,subnet.significance=0.05,verbose=T)
tanzi_ranked_genes_pSubnet.KEGG=enrichr(genes = V(tanzi_ranked_genes.pSubnet)$name,databases = kegg_dbs)[[1]]%>%filter(Adjusted.P.value<0.1)
tanzi_ranked_genes_pSubnet.KEGG$KEGGID=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = tanzi_ranked_genes_pSubnet.KEGG$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[2]]))
tanzi_ranked_genes_pSubnet.KEGG$Pathway=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = tanzi_ranked_genes_pSubnet.KEGG$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[1]]))
tanzi_ranked_genes_pSubnet.KEGG$SNPgene=unlist(lapply(lapply(tanzi_ranked_genes_pSubnet.KEGG$Genes,function(x)strsplit(x = x,split=";")[[1]]),function(y)paste(intersect(y,tanzi_ranked_genes.union),collapse = ";")))
tanzi_ranked_genes_pSubnet.KEGG$knownAD=unlist(lapply(lapply(tanzi_ranked_genes_pSubnet.KEGG$Genes,function(x)strsplit(x = x,split=";")[[1]]),function(y)paste(intersect(y,known_AD_genes.cons),collapse = ";")))
tanzi_ranked_genes_pSubnet.KEGG$Overlap=paste(" ",tanzi_ranked_genes_pSubnet.KEGG$Overlap,sep = "")
fwrite(tanzi_ranked_genes_pSubnet.KEGG%>%dplyr::select(Pathway,KEGGID,Overlap,P.value,Adjusted.P.value,Genes,SNPgene,knownAD),"/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Tanzi_Ranked_SNPgenes_pSubnet_KEGG.txt",sep="\t",col.names = T,row.names = T)



known_ranked_AD_subgraph=induced.subgraph(graph = hs_ppi,vids = which(V(hs_ppi)$name%in%union(tanzi_ranked_genes.AD,known_AD_genes.cons)))
known_random_ranked_AD_subgraph=induced.subgraph(graph = hs_ppi,vids = which(V(hs_ppi)$name%in%tanzi_random_ranked))
known_region_AD_subgraph=induced.subgraph(graph = hs_ppi,vids = which(V(hs_ppi)$name%in%union(tanzi_ranked_regions.great$Gene_GREAT,known_AD_genes.cons)))
known_ranked_region_AD_subgraph=induced.subgraph(graph = hs_ppi,vids = which(V(hs_ppi)$name%in%union(union(tanzi_ranked_regions.great$Gene_GREAT,known_AD_genes.cons),tanzi_ranked_genes.AD)))

# V(known_ranked_AD_subgraph)$colour="#BDBDBD"
# V(known_ranked_AD_subgraph)$colour[which(V(known_ranked_AD_subgraph)$name%in%known_AD_genes.cons)]="#388E3C"
# V(known_ranked_AD_subgraph)$colour[which(V(known_ranked_AD_subgraph)$name%in%tanzi_ranked_genes.AD)]="#E64A19"
known_ranked_AD_subgraph2=delete_vertices(graph = known_ranked_AD_subgraph,v = which((igraph::degree(known_ranked_AD_subgraph)>2)==F))
known_random_ranked_AD_subgraph2=delete_vertices(graph = known_random_ranked_AD_subgraph,v = which((igraph::degree(known_random_ranked_AD_subgraph)>2)==F))
known_region_AD_subgraph2=delete_vertices(graph = known_region_AD_subgraph,v = which((igraph::degree(known_region_AD_subgraph)>2)==F))
known_ranked_region_AD_subgraph2=delete_vertices(graph = known_ranked_region_AD_subgraph,v = which((igraph::degree(known_ranked_region_AD_subgraph)>2)==F))
write.graph(known_ranked_AD_subgraph2,"KnownAD_Top1k_GREAT_ADvariants_PPI_Subnet_pairwiseInteractions.gml",format = "gml")
write.graph(known_region_AD_subgraph2,"KnownAD_Top1k_Regional_ADvariants_PPI_Subnet_pairwiseInteractions.gml",format = "gml")
write.graph(known_ranked_region_AD_subgraph2,"KnownAD_Top1k_Regional_GREAT_ADvariants_PPI_Subnet_pairwiseInteractions.gml",format = "gml")

# Functional analysis of knownAD-rankedAD subnetwork
kegg_genes_by_pathway.df=toTable(org.Hs.egPATH)
kegg_genes_by_pathway.df$path_id=paste("hsa",kegg_genes_by_pathway.df$path_id,sep = "")
kegg_genes_by_pathway.list=lapply(unique(kegg_genes_by_pathway.df$path_id),function(x)kegg_genes_by_pathway.df%>%filter(path_id==x))

kegg_pathways_by_gene=lapply(kegg_pathways,function(x)kegg_genes_by_pathway.df%>%dplyr::filter(path_id==x))
names(kegg_pathways_by_gene)=unlist(lapply(kegg_pathways_by_gene,function(x)unique(x[,2])))
kegg_pathways_by_gene=lapply(kegg_pathways_by_gene,function(x){x<-x[,-2];x})
kegg_pathways_by_gene=lapply(kegg_pathways_by_gene,function(x)mapIds2(IDs = x,IDFrom = "ENTREZID",IDTo = "SYMBOL")[[1]][,2])
known_genes_kegg_pathways_by_gene=kegg_pathways_by_gene[which(names(kegg_pathways_by_gene)%in%names(known_AD_genes.cons))]
names(known_genes_kegg_pathways_by_gene)=mapIds2(IDs = names(known_genes_kegg_pathways_by_gene),IDFrom = "ENTREZID",IDTo = "SYMBOL")[[1]][,2]


plcg2_kegg_pathways=c("Inositol phosphate metabolism","Metabolic pathways","EGFR tyrosine kinase inhibitor resistance","ErbB signaling pathway","Ras signaling pathway","Calcium signaling pathway","NF-kappa B signaling pathway","HIF-1 signaling pathway","Phosphatidylinositol signaling system","Phospholipase D signaling pathway","Axon guidance","VEGF signaling pathway","Osteoclast differentiation","Platelet activation","Natural killer cell mediated cytotoxicity","B cell receptor signaling pathway","Fc epsilon RI signaling pathway","Fc gamma R-mediated phagocytosis","Leukocyte transendothelial migration","Neurotrophin signaling pathway","Inflammatory mediator regulation of TRP channels","Thyroid hormone signaling pathway","AGE-RAGE signaling pathway in diabetic complications","Vibrio cholerae infection","Epithelial cell signaling in Helicobacter pylori infection","Kaposi's sarcoma-associated herpesvirus infection","Epstein-Barr virus infection","Pathways in cancer","Proteoglycans in cancer","MicroRNAs in cancer","Glioma","Non-small cell lung cancer","Hepatocellular carcinoma")
knownAD_cons.KEGG=enrichr(genes = known_AD_genes.cons,databases = kegg_dbs)[[1]]%>%filter(Adjusted.P.value<0.1)
knownAD_cons.KEGG$Pathway=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = knownAD_cons.KEGG$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[1]]))
knownAD_cons.KEGG_PLCG2_pathways=lapply(knownAD_cons.KEGG$Genes,function(y)grep(pattern = "PLCG2",x = strsplit(x = y,split = ";")[[1]]))
names(knownAD_cons.KEGG_PLCG2_pathways)=knownAD_cons.KEGG$Term
knownAD_cons.Biocarta=enrichr(genes = known_AD_genes.cons,databases = biocarta_dbs)[[1]]%>%filter(Adjusted.P.value<0.1)
knownAD_cons.Panther=enrichr(genes = known_AD_genes.cons,databases = panther_dbs)[[1]]%>%filter(Adjusted.P.value<0.1)

knownAD_genes_LD.KEGG=enrichr(genes = known_AD_genes.LD,databases = kegg_dbs)[[1]]%>%filter(Adjusted.P.value<0.1)
knownAD_genes_LD.KEGG$Pathway=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = knownAD_genes_LD.KEGG$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[1]]))


known_ranked_AD_subgraph2.KEGG=enrichr(genes =V(known_ranked_AD_subgraph2)$name,databases = kegg_dbs)[[1]]%>%filter(Adjusted.P.value<0.1)
known_random_ranked_AD_subgraph2.KEGG=enrichr(genes =V(known_random_ranked_AD_subgraph2)$name,databases = kegg_dbs)[[1]]%>%filter(Adjusted.P.value<0.1)
known_ranked_AD_subgraph2.BioCarta=enrichr(genes =V(known_ranked_AD_subgraph2)$name,databases = biocarta_dbs)[[1]]%>%filter(Adjusted.P.value<0.1)
known_ranked_AD_subgraph2.Panther=enrichr(genes =V(known_ranked_AD_subgraph2)$name,databases = panther_dbs)[[1]]%>%filter(Adjusted.P.value<0.1)

known_ranked_AD_subgraph2.KEGG$KEGGID=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = known_ranked_AD_subgraph2.KEGG$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[2]]))
known_random_ranked_AD_subgraph2.KEGG$KEGGID=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = known_random_ranked_AD_subgraph2.KEGG$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[2]]))
known_ranked_AD_subgraph2.BioCarta$BiocartaID=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = known_ranked_AD_subgraph2.BioCarta$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[3]]))
known_ranked_AD_subgraph2.Panther$PantherID=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = known_ranked_AD_subgraph2.Panther$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[2]]))

known_ranked_AD_subgraph2.KEGG$Pathway=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = known_ranked_AD_subgraph2.KEGG$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[1]]))
known_random_ranked_AD_subgraph2.KEGG$Pathway=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = known_random_ranked_AD_subgraph2.KEGG$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[1]]))
known_ranked_AD_subgraph2.BioCarta$Pathway=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = known_ranked_AD_subgraph2.BioCarta$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[1]]))
known_ranked_AD_subgraph2.Panther$Pathway=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = known_ranked_AD_subgraph2.Panther$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[1]]))

known_ranked_AD_subgraph2.KEGG$knownAD=unlist(lapply(lapply(known_ranked_AD_subgraph2.KEGG$Genes,function(x)strsplit(x = x,split=";")[[1]]),function(y)paste(intersect(y,known_AD_genes.cons),collapse = ";")))
known_random_ranked_AD_subgraph2.KEGG$knownAD=unlist(lapply(lapply(known_random_ranked_AD_subgraph2.KEGG$Genes,function(x)strsplit(x = x,split=";")[[1]]),function(y)paste(intersect(y,known_AD_genes.cons),collapse = ";")))
known_ranked_AD_subgraph2.BioCarta$knownAD=unlist(lapply(lapply(known_ranked_AD_subgraph2.BioCarta$Genes,function(x)strsplit(x = x,split=";")[[1]]),function(y)paste(intersect(y,known_AD_genes.cons),collapse = ";")))
known_ranked_AD_subgraph2.Panther$knownAD=unlist(lapply(lapply(known_ranked_AD_subgraph2.Panther$Genes,function(x)strsplit(x = x,split=";")[[1]]),function(y)paste(intersect(y,known_AD_genes.cons),collapse = ";")))


write.table(known_ranked_AD_subgraph2.KEGG%>%select(Pathway,KEGGID,Overlap,P.value,Adjusted.P.value,Genes,knownAD),"knownAD_rankedAD_subgraph_enrichedKEGG_Enrichr.txt",sep = "\t",col.names = T,row.names = F)
write.table(known_ranked_AD_subgraph2.BioCarta%>%select(Pathway,BiocartaID,Overlap,P.value,Adjusted.P.value,Genes,knownAD),"knownAD_rankedAD_subgraph_enrichedBiocarta_Enrichr.txt",sep = "\t",col.names = T,row.names = F)
write.table(known_ranked_AD_subgraph2.Panther%>%select(Pathway,PantherID,Overlap,P.value,Adjusted.P.value,Genes,knownAD),"knownAD_rankedAD_subgraph_enrichedPanther_Enrichr.txt",sep = "\t",col.names = T,row.names = F)


# Functional analysis of knownAD-regionAD subnetwork
known_region_AD_subgraph2.KEGG=enrichr(genes =V(known_region_AD_subgraph2)$name,databases = kegg_dbs)[[1]]%>%filter(Adjusted.P.value<0.1)
known_region_AD_subgraph2.BioCarta=enrichr(genes =V(known_region_AD_subgraph2)$name,databases = biocarta_dbs)[[1]]%>%filter(Adjusted.P.value<0.1)
known_region_AD_subgraph2.Panther=enrichr(genes =V(known_region_AD_subgraph2)$name,databases = panther_dbs)[[1]]%>%filter(Adjusted.P.value<0.1)

known_region_AD_subgraph2.KEGG$KEGGID=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = known_region_AD_subgraph2.KEGG$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[2]]))
known_region_AD_subgraph2.BioCarta$BiocartaID=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = known_region_AD_subgraph2.BioCarta$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[3]]))
known_region_AD_subgraph2.Panther$PantherID=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = known_region_AD_subgraph2.Panther$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[2]]))

known_region_AD_subgraph2.KEGG$Pathway=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = known_region_AD_subgraph2.KEGG$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[1]]))
known_region_AD_subgraph2.BioCarta$Pathway=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = known_region_AD_subgraph2.BioCarta$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[1]]))
known_region_AD_subgraph2.Panther$Pathway=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = known_region_AD_subgraph2.Panther$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[1]]))

known_region_AD_subgraph2.KEGG$knownAD=unlist(lapply(lapply(known_region_AD_subgraph2.KEGG$Genes,function(x)strsplit(x = x,split=";")[[1]]),function(y)paste(intersect(y,known_AD_genes.cons),collapse = ";")))
known_region_AD_subgraph2.BioCarta$knownAD=unlist(lapply(lapply(known_region_AD_subgraph2.BioCarta$Genes,function(x)strsplit(x = x,split=";")[[1]]),function(y)paste(intersect(y,known_AD_genes.cons),collapse = ";")))
known_region_AD_subgraph2.Panther$knownAD=unlist(lapply(lapply(known_region_AD_subgraph2.Panther$Genes,function(x)strsplit(x = x,split=";")[[1]]),function(y)paste(intersect(y,known_AD_genes.cons),collapse = ";")))

known_region_AD_subgraph2.KEGG$Pathway=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = known_region_AD_subgraph2.KEGG$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[1]]))
known_region_AD_subgraph2.BioCarta$Pathway=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = known_region_AD_subgraph2.BioCarta$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[1]]))
known_region_AD_subgraph2.Panther$Pathway=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = known_region_AD_subgraph2.Panther$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[1]]))

known_region_AD_subgraph2.KEGG$knownAD=unlist(lapply(lapply(known_region_AD_subgraph2.KEGG$Genes,function(x)strsplit(x = x,split=";")[[1]]),function(y)paste(intersect(y,known_AD_genes.cons),collapse = ";")))
known_region_AD_subgraph2.BioCarta$knownAD=unlist(lapply(lapply(known_region_AD_subgraph2.BioCarta$Genes,function(x)strsplit(x = x,split=";")[[1]]),function(y)paste(intersect(y,known_AD_genes.cons),collapse = ";")))
known_region_AD_subgraph2.Panther$knownAD=unlist(lapply(lapply(known_region_AD_subgraph2.Panther$Genes,function(x)strsplit(x = x,split=";")[[1]]),function(y)paste(intersect(y,known_AD_genes.cons),collapse = ";")))

write.table(known_region_AD_subgraph2.KEGG%>%select(Pathway,KEGGID,Overlap,P.value,Adjusted.P.value,Genes,knownAD),"knownAD_regionAD_subgraph_enrichedKEGG_Enrichr.txt",sep = "\t",col.names = T,row.names = F)
write.table(known_region_AD_subgraph2.BioCarta%>%select(Pathway,BiocartaID,Overlap,P.value,Adjusted.P.value,Genes,knownAD),"knownAD_regionAD_subgraph_enrichedBiocarta_Enrichr.txt",sep = "\t",col.names = T,row.names = F)
write.table(known_region_AD_subgraph2.Panther%>%select(Pathway,PantherID,Overlap,P.value,Adjusted.P.value,Genes,knownAD),"knownAD_regionAD_subgraph_enrichedPanther_Enrichr.txt",sep = "\t",col.names = T,row.names = F)

# Functional analysis of knownAD-regionAD-rankedAD subnetwork
known_ranked_region_AD_subgraph2.KEGG=enrichr(genes =V(known_ranked_region_AD_subgraph2)$name,databases = kegg_dbs)[[1]]%>%filter(Adjusted.P.value<0.1)
known_ranked_region_AD_subgraph2.BioCarta=enrichr(genes =V(known_ranked_region_AD_subgraph2)$name,databases = biocarta_dbs)[[1]]%>%filter(Adjusted.P.value<0.1)
known_ranked_region_AD_subgraph2.Panther=enrichr(genes =V(known_ranked_region_AD_subgraph2)$name,databases = panther_dbs)[[1]]%>%filter(Adjusted.P.value<0.1)

known_ranked_region_AD_subgraph2.KEGG$KEGGID=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = known_ranked_region_AD_subgraph2.KEGG$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[2]]))
known_ranked_region_AD_subgraph2.BioCarta$BiocartaID=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = known_ranked_region_AD_subgraph2.BioCarta$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[3]]))
known_ranked_region_AD_subgraph2.Panther$PantherID=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = known_ranked_region_AD_subgraph2.Panther$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[2]]))

known_ranked_region_AD_subgraph2.KEGG$Pathway=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = known_ranked_region_AD_subgraph2.KEGG$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[1]]))
known_ranked_region_AD_subgraph2.BioCarta$Pathway=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = known_ranked_region_AD_subgraph2.BioCarta$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[1]]))
known_ranked_region_AD_subgraph2.Panther$Pathway=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = known_ranked_region_AD_subgraph2.Panther$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[1]]))

known_ranked_region_AD_subgraph2.KEGG$knownAD=unlist(lapply(lapply(known_ranked_region_AD_subgraph2.KEGG$Genes,function(x)strsplit(x = x,split=";")[[1]]),function(y)paste(intersect(y,known_AD_genes.cons),collapse = ";")))
known_ranked_region_AD_subgraph2.BioCarta$knownAD=unlist(lapply(lapply(known_ranked_region_AD_subgraph2.BioCarta$Genes,function(x)strsplit(x = x,split=";")[[1]]),function(y)paste(intersect(y,known_AD_genes.cons),collapse = ";")))
known_ranked_region_AD_subgraph2.Panther$knownAD=unlist(lapply(lapply(known_ranked_region_AD_subgraph2.Panther$Genes,function(x)strsplit(x = x,split=";")[[1]]),function(y)paste(intersect(y,known_AD_genes.cons),collapse = ";")))


write.table(known_ranked_region_AD_subgraph2.KEGG%>%select(Pathway,KEGGID,Overlap,P.value,Adjusted.P.value,Genes,knownAD),"knownAD_rankedAD_regionAD_subgraph_enrichedKEGG_Enrichr.txt",sep = "\t",col.names = T,row.names = F)
write.table(known_ranked_region_AD_subgraph2.BioCarta%>%select(Pathway,BiocartaID,Overlap,P.value,Adjusted.P.value,Genes,knownAD),"knownAD_rankedAD_regionAD_subgraph_enrichedBiocarta_Enrichr.txt",sep = "\t",col.names = T,row.names = F)
write.table(known_ranked_region_AD_subgraph2.Panther%>%select(Pathway,PantherID,Overlap,P.value,Adjusted.P.value,Genes,knownAD),"knownAD_rankedAD_regionAD_subgraph_enrichedPanther_Enrichr.txt",sep = "\t",col.names = T,row.names = F)


#List wise enrichment
tanzi_ranked_genes.KEGG=enrichr(genes = tanzi_ranked_genes.union,databases = kegg_dbs)[[1]]%>%filter(Adjusted.P.value<0.1)

kegg2016=read.delim2("/Users/sandeepamberkar/Work/Data/KEGG_2016.txt",header = F,as.is = T)
kegg2016_list=foreach(i=1:dim(kegg2016)[1])%do%{
  unname(unlist(kegg2016[i,which(kegg2016[i,]!="")][,-1]))
}
names(kegg2016_list)=kegg2016[,1]
kegg2016_pathway_member_genes=unique(unlist(kegg2016_list))
kegg2016_pathway_member_genes.knownAD=intersect(known_AD_genes.cons,kegg2016_pathway_member_genes)
tanzi_ranked_genes_KEGG.knownAD_mapped=intersect(names(lapply(kegg2016_list,function(x)length(intersect(x,kegg2016_pathway_member_genes.knownAD)))[unlist(lapply(kegg2016_list,function(x)length(intersect(x,kegg2016_pathway_member_genes.knownAD))))>0]),tanzi_ranked_genes.KEGG$Term)
tanzi_ranked_genes.KEGG$knownAD_in_KEGG=NA
tanzi_ranked_genes.KEGG$knownAD_in_KEGG[which(tanzi_ranked_genes.KEGG$Term%in%tanzi_ranked_genes_KEGG.knownAD_mapped)]=unname(unlist(lapply(kegg2016_list[tanzi_ranked_genes.KEGG[tanzi_ranked_genes.KEGG$Term%in%tanzi_ranked_genes_KEGG.knownAD_mapped,1]],function(x)paste(intersect(x,known_AD_genes.cons),collapse = ","))))
tanzi_ranked_genes_KEGG.FreqGenes=data.frame(Genes=names(sort(table(unlist(lapply(tanzi_ranked_genes.KEGG$Genes[which(is.na(tanzi_ranked_genes.KEGG$knownAD_in_KEGG)==F)],function(y)strsplit(x = y,split = ";")[[1]]))),decreasing = T)),Freq=sort(table(unlist(lapply(tanzi_ranked_genes.KEGG$Genes[which(is.na(tanzi_ranked_genes.KEGG$knownAD_in_KEGG)==F)],function(y)strsplit(x = y,split = ";")[[1]]))),decreasing = T),stringsAsFactors = F)
tanzi_ranked_genes.Biocarta=enrichr(genes = union(tanzi_ranked_genes.AD,tanzi_ranked_genes.Control),databases = biocarta_dbs)[[1]]%>%filter(Adjusted.P.value<0.1)
tanzi_ranked_genes.Panther=enrichr(genes = union(tanzi_ranked_genes.AD,tanzi_ranked_genes.Control),databases = panther_dbs)[[1]]%>%filter(Adjusted.P.value<0.1)

tanzi_ranked_genes.KEGG$KEGGID=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = tanzi_ranked_genes.KEGG$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[2]]))
tanzi_ranked_genes.Biocarta$BiocartaID=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = tanzi_ranked_genes.Biocarta$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[3]]))
tanzi_ranked_genes.Panther$PantherID=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = tanzi_ranked_genes.Panther$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[2]]))

tanzi_ranked_genes.KEGG$Pathway=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = tanzi_ranked_genes.KEGG$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[1]]))
tanzi_ranked_genes.Biocarta$Pathway=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = tanzi_ranked_genes.Biocarta$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[1]]))
tanzi_ranked_genes.Panther$Pathway=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = tanzi_ranked_genes.Panther$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[1]]))

tanzi_ranked_genes.KEGG$knownAD=unlist(lapply(lapply(tanzi_ranked_genes.KEGG$Genes,function(x)strsplit(x = x,split=";")[[1]]),function(y)paste(intersect(y,known_AD_genes.cons),collapse = ";")))
tanzi_ranked_genes.Biocarta$knownAD=unlist(lapply(lapply(tanzi_ranked_genes.Biocarta$Genes,function(x)strsplit(x = x,split=";")[[1]]),function(y)paste(intersect(y,known_AD_genes.cons),collapse = ";")))
tanzi_ranked_genes.Panther$knownAD=unlist(lapply(lapply(tanzi_ranked_genes.Panther$Genes,function(x)strsplit(x = x,split=";")[[1]]),function(y)paste(intersect(y,known_AD_genes.cons),collapse = ";")))

tanzi_ranked_genes.KEGG$XY_genes=unlist(lapply(lapply(tanzi_ranked_genes.KEGG$Genes,function(x)strsplit(x = x,split=";")[[1]]),function(y)paste(intersect(y,tableXY_GREAT_genes),collapse = ";")))

write.table(tanzi_ranked_genes.KEGG%>%select(Pathway,KEGGID,Overlap,P.value,Adjusted.P.value,Genes),"tanzi_ranked_SNP_KEGG.txt",sep = "\t",col.names = T,row.names = F)
write.table(tanzi_ranked_genes.Biocarta%>%select(Pathway,BiocartaID,Overlap,P.value,Adjusted.P.value,Genes),"tanzi_ranked_SNP_Biocarta.txt",sep = "\t",col.names = T,row.names = F)
write.table(tanzi_ranked_genes.Panther%>%select(Pathway,PantherID,Overlap,P.value,Adjusted.P.value,Genes),"tanzi_ranked_SNP_Panther.txt",sep = "\t",col.names = T,row.names = F)

tanzi_region_genes.KEGG=enrichr(genes = tanzi_region_genes,databases = kegg_dbs)[[1]]%>%filter(Adjusted.P.value<0.1)
tanzi_region_genes.Biocarta=enrichr(genes = tanzi_region_genes,databases = biocarta_dbs)[[1]]%>%filter(Adjusted.P.value<0.1)
tanzi_region_genes.Panther=enrichr(genes = tanzi_region_genes,databases = panther_dbs)[[1]]%>%filter(Adjusted.P.value<0.1)

tanzi_region_genes.KEGG$KEGGID=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = tanzi_region_genes.KEGG$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[2]]))
tanzi_region_genes.Biocarta$BiocartaID=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = tanzi_region_genes.Biocarta$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[3]]))
tanzi_region_genes.Panther$PantherID=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = tanzi_region_genes.Panther$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[2]]))

tanzi_region_genes.KEGG$Pathway=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = tanzi_region_genes.KEGG$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[1]]))
tanzi_region_genes.Biocarta$Pathway=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = tanzi_region_genes.Biocarta$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[1]]))
tanzi_region_genes.Panther$Pathway=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = tanzi_region_genes.Panther$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[1]]))

write.table(tanzi_region_genes.KEGG%>%select(Pathway,KEGGID,Overlap,P.value,Adjusted.P.value,Genes),"tanzi_region_SNP_KEGG.txt",sep = "\t",col.names = T,row.names = F)
write.table(tanzi_region_genes.Biocarta%>%select(Pathway,BiocartaID,Overlap,P.value,Adjusted.P.value,Genes),"tanzi_region_SNP_Biocarta.txt",sep = "\t",col.names = T,row.names = F)
write.table(tanzi_region_genes.Panther%>%select(Pathway,PantherID,Overlap,P.value,Adjusted.P.value,Genes),"tanzi_region_SNP_Panther.txt",sep = "\t",col.names = T,row.names = F)


names(known_ranked_AD_subgraph2.KEGG_genes)=known_ranked_AD_subgraph2.KEGG$Description
known_ranked_AD_subgraph2.KEGG_genes_compareMatrix=matrix(NA,nrow = 12,ncol = 12)
for(l in 1:12){
  known_ranked_AD_subgraph2.KEGG_genes_compareMatrix[l,]=unlist(lapply(known_ranked_AD_subgraph2.KEGG_genes,function(x)length(intersect(x,known_ranked_AD_subgraph2.KEGG_genes[[l]]))))
}
colnames(known_ranked_AD_subgraph2.KEGG_genes_compareMatrix)=rownames(known_ranked_AD_subgraph2.KEGG_genes_compareMatrix)=known_ranked_AD_subgraph2.KEGG$Description
V(known_ranked_AD_subgraph2)$colour="#BDBDBD"
V(known_ranked_AD_subgraph2)$border=1.0
subnet_colours=sample(x = unique(col_vector),size = length(known_ranked_AD_subgraph2.KEGG_genes),replace = F)
for(c in 1:length(subnet_colours)){
  V(known_ranked_AD_subgraph2)$colour[which(V(known_ranked_AD_subgraph2)$name%in%known_ranked_AD_subgraph2.KEGG_genes[[c]])]=subnet_colours[c]
}

thy_KEGG_genes.expr=rep(NA,length(known_ranked_region_AD_subgraph2.KEGG_genes$`Thyroid hormone signaling pathway`))
names(thy_KEGG_genes.expr)=mapIds2(IDs = known_ranked_region_AD_subgraph2.KEGG_genes$`Thyroid hormone signaling pathway`,IDFrom = "SYMBOL",IDTo = "ENSEMBL")[[1]][,2]
names(thy_KEGG_genes.expr)=known_ranked_region_AD_subgraph2.KEGG_genes$`Thyroid hormone signaling pathway`
thy_KEGG_genes.expr[which(names(thy_KEGG_genes.expr)%in%tanzi_ranked_genes.AD)]=-1
thy_KEGG_genes.expr[which(names(thy_KEGG_genes.expr)%in%tanzi_ranked_regions.great$Gene_GREAT)]=1
plotKegg(keggIDs = "hsa04919",geneExpr = thy_KEGG_genes.expr,geneIDtype = "ENSEMBL")


pld_KEGG_genes.expr=c(rep(1,length(strsplit(tanzi_ranked_genes.KEGG[7,9],split = ";")[[1]])),rep(-1,length(strsplit(knownAD_cons.KEGG[26,9],split = ";")[[1]])))
names(pld_KEGG_genes.expr)=c(strsplit(tanzi_ranked_genes.KEGG[7,9],split = ";")[[1]],strsplit(knownAD_cons.KEGG[26,9],split = ";")[[1]])
plotKegg(keggIDs = "hsa04072",geneExpr = pld_KEGG_genes.expr)

bcell_KEGG_genes.expr=c(rep(1,length(strsplit(tanzi_ranked_genes.KEGG[19,9],split = ";")[[1]])),rep(-1,length(strsplit(knownAD_genes_LD.KEGG[21,9],split = ";")[[1]])))
names(bcell_KEGG_genes.expr)=c(strsplit(tanzi_ranked_genes.KEGG[19,9],split = ";")[[1]],strsplit(knownAD_genes_LD.KEGG[21,9],split = ";")[[1]])
plotKegg(keggIDs = "hsa04662",geneExpr = bcell_KEGG_genes.expr)
#Write 3D network layout file





geneID_KEGG_map=toTable(org.Hs.egPATH2EG)
geneID_KEGG_map.list=lapply(lapply(unique(geneID_KEGG_map$path_id),function(x)geneID_KEGG_map%>%filter(path_id==x)%>%pull(gene_id)),function(y)mapIds2(IDs = y,IDFrom = "ENTREZID",IDTo = "SYMBOL")[[1]][,2])
names(geneID_KEGG_map.list)=unique(geneID_KEGG_map$path_id)

V(tanzi_ranked_genes_AD.pSubnet)$colour="#8B8989"
V(tanzi_ranked_genes_AD.pSubnet)$colour[which(V(tanzi_ranked_genes_AD.pSubnet)$name%in%tanzi_ranked_genes.AD)]="#EE4000"
write.graph(tanzi_ranked_genes_AD.pSubnet,"/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Tanzi_RankedGenes_AD_p9_pSubnet.gml",format = "gml")
V(tanzi_ranked_genes_Control.pSubnet)$colour="#8B8989"
V(tanzi_ranked_genes_Control.pSubnet)$colour[which(V(tanzi_ranked_genes_Control.pSubnet)$name%in%tanzi_ranked_genes.Control)]="#EE4000"
write.graph(tanzi_ranked_genes_Control.pSubnet,"/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Tanzi_RankedGenes_Control_p9_pSubnet.gml",format = "gml")

data.frame(gprofiler(query=V(tanzi_ranked_genes_AD.pSubnet)$name,organism = "hsapiens",significant = T,src_filter = "KEGG",correction_method = "bonferroni",max_p_value = 0.05))
write.table(data.frame(gprofiler(query=V(tanzi_ranked_genes_Control.pSubnet)$name,organism = "hsapiens",significant = T,src_filter = "KEGG",correction_method = "bonferroni",max_p_value = 0.05)),"Tanzi_Control_Great_Subnet_p9_KEGG.txt",sep = "\t",col.names = T,row.names = F)

# Distance based comparisons
random_variants.sp=shortest.paths(graph=hs_ppi,v=V(hs_ppi)$name[which(V(hs_ppi)$name%in%known_AD_genes.cons)],to = V(hs_ppi)$name[which(V(hs_ppi)$name%in%sample(x = V(hs_ppi)$name,size = length(union(tanzi_ranked_genes.AD,tanzi_ranked_genes.Control)),replace = T))])
random_variants.sp=upperTriangle(random_variants.sp)[upperTriangle(random_variants.sp)!="Inf"]
known_AD.ShortestPaths=shortest.paths(graph = hs_ppi,v = which(V(hs_ppi)$name %in% known_AD_genes),to = which(V(hs_ppi)$name %in% known_AD_genes.cons))
known_AD.ShortestPaths=upperTriangle(known_AD.ShortestPaths)[upperTriangle(known_AD.ShortestPaths)!="Inf"]
known_AD_ranked_AD.ShortestPaths=shortest.paths(graph = hs_ppi,v = which(V(hs_ppi)$name %in% known_AD_genes.cons),to = which(V(hs_ppi)$name %in% tanzi_ranked_genes.AD))
known_AD_region_AD.ShortestPaths=shortest.paths(graph = hs_ppi,v = which(V(hs_ppi)$name %in% known_AD_genes.cons),to = which(V(hs_ppi)$name %in% tanzi_ranked_regions.great$Gene_GREAT))
known_AD_ranked_onlyAD.ShortestPaths=shortest.paths(graph = hs_ppi,v = which(V(hs_ppi)$name %in% known_AD_genes),to = which(V(hs_ppi)$name %in% setdiff(tanzi_ranked_genes.AD,tanzi_ranked_genes.Control)))
known_AD_ranked_Control.ShortestPaths=shortest.paths(graph = hs_ppi,v = which(V(hs_ppi)$name %in% known_AD_genes),to = which(V(hs_ppi)$name %in% tanzi_ranked_genes.Control))
known_AD_ranked_onlyControl.ShortestPaths=shortest.paths(graph = hs_ppi,v = which(V(hs_ppi)$name %in% known_AD_genes),to = which(V(hs_ppi)$name %in% setdiff(tanzi_ranked_genes.Control,tanzi_ranked_genes.AD)))
#known_AD_ranked_AD.ShortestPaths=upperTriangle(known_AD_ranked_AD.ShortestPaths)[upperTriangle(known_AD_ranked_AD.ShortestPaths)!="Inf"]
known_AD_ranked_AD.ShortestPaths=upperTriangle(known_AD_ranked_AD.ShortestPaths)[upperTriangle(known_AD_ranked_AD.ShortestPaths)!="Inf"]
known_AD_region_AD.ShortestPaths=upperTriangle(known_AD_region_AD.ShortestPaths)[upperTriangle(known_AD_region_AD.ShortestPaths)!="Inf"]
known_AD_ranked_Control.ShortestPaths=upperTriangle(known_AD_ranked_Control.ShortestPaths)[upperTriangle(known_AD_ranked_Control.ShortestPaths)!="Inf"]
known_AD_ranked_onlyControl.ShortestPaths=upperTriangle(known_AD_ranked_onlyControl.ShortestPaths)[upperTriangle(known_AD_ranked_onlyControl.ShortestPaths)!="Inf"]
hs_ppi.ShortestPaths=shortest.paths(graph = hs_ppi,v = V(hs_ppi)$name,to = V(hs_ppi)$name)
hs_ppi.ShortestPaths=upperTriangle(hs_ppi.ShortestPaths)[upperTriangle(hs_ppi.ShortestPaths)!="Inf"]

# Compare AD variants in expression space
setwd("/Users/sandeepamberkar/Work/Data/MSMM-MSBB-HBTRC-Data/MSBB_Array19/Normalised_Data/")
msbb_array19.files=list.files(pattern = "*.tsv",full.names = T)
msbb_array19.files=msbb_array19.files[-c(18:19)]
msbb_array19=mapply(FUN = fread,msbb_array19.files,MoreArgs = list(header=T,sep="\t",data.table=F,showProgress=T))
names(msbb_array19)=c("AC","CN","DLPFC","FP","HP","IFG","ITG","MTG","OVC","PHG","PCC","PCG","PFC","PTMN","SPL","STG","TP")
msbb_array19.covariates=read.delim2("/Users/sandeepamberkar/Work/Data/MSMM-MSBB-HBTRC-Data/MSBB_Array19/Normalised_Data/AMP-AD_MSBB_MSSM_covariates_mRNA_AffymetrixU133AB.txt",header = T,as.is = T)%>%
  mutate(Age=replace(Age,Age=="89+",90))%>%
  mutate(pH=as.numeric(pH))%>%
  mutate(CDR=as.numeric(CDR))%>%
  mutate(PLQ_Mn=as.numeric(PLQ_Mn))%>%
  mutate(BrainBank=paste("X",BrainBank,sep=""))


msbb_array19.Control=msbb_array19.covariates%>%filter(CDR<=1&Braak<=3&NP1<=1)
msbb_array19.Control$BrainBank=gsub(pattern = "X",replacement = "",x = msbb_array19.Control$BrainBank)
msbb_array19.AD=msbb_array19.covariates%>%filter(CDR>2&Braak>=4&NP1>=2)
msbb_array19.AD$BrainBank=gsub(pattern = "X",replacement = "",x = msbb_array19.AD$BrainBank)

msbb_array19.2=lapply(msbb_array19[-18],function(x){rownames(x)<-x$ID;x})
multiID.u133a=msbb_array19.2[[1]]$Gene.Symbol[grep(pattern="///",msbb_array19.2[[1]]$Gene.Symbol)]
msbb_array19.2.agg=lapply(msbb_array19.2[-18],function(y)aggregate(x=y[,-c(1:4)],by=list(Symbol=y$Gene.Symbol),mean))
msbb_array19.2.agg2=lapply(msbb_array19.2.agg,function(x){rownames(x)<- x$Symbol;x <- x[-grep(pattern = "///",x = rownames(x)),-1];x})
msbb_array19.2.agg2=readRDS("/Users/sandeepamberkar/Work/Data/MSMM-MSBB-HBTRC-Data/MSBB_Array19/Normalised_Data/msbb_array19_Agg.RDS")
names(msbb_array19.2.agg2)=names(msbb_array19)
select_brainRegions=names(msbb_array19)[c(5,7:8,10,16)]
msbb_array19.AllSamples=do.call("cbind",msbb_array19.2.agg2)
saveRDS(msbb_array19.AllSamples,"MSBB_array19_AllSamples.RDS")
msbb_array19.AllSamples=readRDS("/Users/sandeepamberkar/Work/Data/MSMM-MSBB-HBTRC-Data/MSBB_Array19/Normalised_Data/MSBB_array19_AllSamples.RDS")
msbb_array19.AllSamples_Control=msbb_array19.AllSamples[,grep(pattern = paste(msbb_array19.Control$BrainBank,collapse = "|"),x = colnames(msbb_array19.AllSamples),value = T)]
msbb_array19.AllSamples_AD=msbb_array19.AllSamples[,grep(pattern = paste(msbb_array19.AD$BrainBank,collapse = "|"),x = colnames(msbb_array19.AllSamples),value = T)]

msbb_array19.AllSamples_select_brainRegion=data.frame(cbind(msbb_array19.AllSamples_Control,msbb_array19.AllSamples_AD),stringsAsFactors = F)
msbb_array19.AllSamples_select_brainRegion.Control=msbb_array19.AllSamples_Control[,grep(pattern = paste(select_brainRegions,collapse = "|"),x = colnames(msbb_array19.AllSamples_Control),value = T)]
msbb_array19.AllSamples_select_brainRegion.AD=msbb_array19.AllSamples_AD[,grep(pattern = paste(select_brainRegions,collapse = "|"),x = colnames(msbb_array19.AllSamples_AD),value = T)]
annotation_df=data.frame(SampleType=factor(c(rep("Control",dim(msbb_array19.AllSamples_Control)[2]),rep("AD",dim(msbb_array19.AllSamples_AD)[2]))),BrainRegion=unlist(lapply(strsplit(x = colnames(msbb_array19.AllSamples_select_brainRegion),split = "\\."),`[[`,1)),stringsAsFactors = F)
rownames(annotation_df)=colnames(msbb_array19.AllSamples_select_brainRegion)
pheatmap(msbb_array19.AllSamples_select_brainRegion[rownames(msbb_array19.AllSamples_select_brainRegion)%in%rownames(deg_df),],annotation = annotation_df)

group=factor(c(rep("Control",dim(msbb_array19.AllSamples_select_brainRegion.Control)[2]),rep("AD",dim(msbb_array19.AllSamples_select_brainRegion.AD)[2])))
design_df=model.matrix(~0+group)
colnames(design_df)=levels(group)
contrasts_matrix=makeContrasts(Control-AD,levels = design_df)
fit=lmFit(object = cbind(msbb_array19.AllSamples_select_brainRegion.Control,msbb_array19.AllSamples_select_brainRegion.AD),design = design_df)
fit2=contrasts.fit(fit,contrasts_matrix)
fit2=eBayes(fit2,trend = T)
deg_df=topTable(fit2,coef = 1,sort.by = "logFC",adjust.method = "BH",p.value = 0.05,number = dim(fit$coefficients)[1])
msbb_array_selectRegions_DEG.df=deg_df%>%rownames_to_column("Gene")
tanzi_DEV=intersect(rownames(deg_df),union(tanzi_ranked_genes.AD,tanzi_ranked_genes.Control))
msbb_array_selectRegions_DEG.KEGG=enrichr(genes = msbb_array_selectRegions_DEG.df$Gene,databases = kegg_dbs)[[1]]%>%filter(Adjusted.P.value<0.1)
msbb_array_selectRegions_DEG.Panther=enrichr(genes = msbb_array_selectRegions_DEG.df$Gene,databases = panther_dbs)[[1]]%>%filter(Adjusted.P.value<=0.1)
fwrite(deg_df[rownames(deg_df)%in%tanzi_ranked_genes.union,],"/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Diff.Expressed.Variants_AD_GreatVariants_Table.txt",sep="\t",col.names = T,row.names = T)
fwrite(deg_df[which(rownames(deg_df)%in%known_AD_genes),],"/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Diff.Expressed_knownAD_Table.txt",sep="\t",col.names = T,row.names = T)

annotation_df=data.frame(SampleType=factor(c(rep("Control",dim(msbb_array19.AllSamples_select_brainRegion.Control)[2]),rep("AD",dim(msbb_array19.AllSamples_select_brainRegion.AD)[2]))),stringsAsFactors = F)
rownames(annotation_df)=c(msbb_array19.AllSamples_select_brainRegion.Control,msbb_array19.AllSamples_select_brainRegion.AD)
annotation_row_df=data.frame(VariantType=factor(c(rep("RankAD",length(intersect(tanzi_ranked_genes.union,rownames(deg_df)))))),stringsAsFactors = F)
                                                  
                    
rownames(annotation_row_df)=c(intersect(tanzi_ranked_genes.union,rownames(deg_df)),
                                                                        paste(intersect(intersect(tanzi_ranked_genes.AD,rownames(deg_df)),intersect(tanzi_ranked_regions.great$Gene_GREAT,rownames(deg_df))),".1",sep = ""))
annotation_row_df$var_counts=unlist(lapply(unlist(lapply(strsplit(rownames(annotation_row_df),split = "\\."),`[[`,1)),function(x)dim(tanzi_ranked_genes.great[tanzi_ranked_genes.great$Gene==x,])[1]))
#Plot all DEVs
pheatmap(msbb_array19.AllSamples_select_brainRegion[rownames(msbb_array19.AllSamples_select_brainRegion)%in%unlist(lapply(strsplit(rownames(annotation_row_df),split = "\\."),`[[`,1)),],annotation = annotation_df,annotation_row = annotation_row_df,annotation_names_row = T,fontsize_row = 8,fontsize_col = 8,main = "WGS-AD,differentially expressed variants across 5 brain regions")


#Plot only Control DEVs
pheatmap(msbb_array19.AllSamples_select_brainRegion[rownames(msbb_array19.AllSamples_select_brainRegion)%in%unlist(lapply(strsplit(rownames(annotation_row_df[annotation_row_df$VariantType=="Control",]),split = "\\."),`[[`,1)),],annotation = annotation_df,annotation_row = annotation_row_df,annotation_names_row = T,fontsize_row = 7,fontsize_col = 7,main = "WGS-AD,differentially expressed variants across 5 brain regions",cellwidth = 6,cellheight = 8,cutree_cols = 2,cutree_rows = 3)

#Plot only AD DEVs
pheatmap(msbb_array19.AllSamples_select_brainRegion[rownames(msbb_array19.AllSamples_select_brainRegion)%in%unlist(lapply(strsplit(rownames(annotation_row_df[annotation_row_df$VariantType=="AD",]),split = "\\."),`[[`,1)),],annotation = annotation_df,annotation_row = annotation_row_df,annotation_names_row = T,fontsize_row = 6,fontsize_col = 6,main = "WGS-AD,differentially expressed variants across 5 brain regions",cellwidth = 5,cellheight = 5)

#Plot only Common DEVs
pheatmap(msbb_array19.AllSamples_select_brainRegion[rownames(msbb_array19.AllSamples_select_brainRegion)%in%unlist(lapply(strsplit(rownames(annotation_row_df[annotation_row_df$VariantType=="Common",]),split = "\\."),`[[`,1)),],annotation = annotation_df,annotation_row = annotation_row_df,annotation_names_row = T,fontsize_row = 8,fontsize_col = 8,main = "WGS-AD,differentially expressed variants across 5 brain regions",cellwidth = 6,cellheight = 10,cutree_cols = 2,cutree_rows = 3)

# Plot heatmaps
pheatmap(msbb_array19.AllSamples_select_brainRegion[rownames(msbb_array19.AllSamples_select_brainRegion)%in%unlist(lapply(strsplit(rownames(annotation_row_df[annotation_row_df$VariantType=="RegionAD",]),split = "\\."),`[[`,1)),],annotation = annotation_df,annotation_row = annotation_row_df,annotation_names_row = T,fontsize_row = 8,fontsize_col = 8,main = "WGS-AD,differentially expressed variants across 5 brain regions")

pheatmap(msbb_array19.AllSamples_select_brainRegion[rownames(msbb_array19.AllSamples_select_brainRegion)%in%unlist(lapply(strsplit(rownames(annotation_row_df[annotation_row_df$VariantType=="RegionAD",]),split = "\\."),`[[`,1)),],annotation = annotation_df,annotation_row = annotation_row_df,annotation_names_row = T,fontsize_row = 6,fontsize_col = 6,main = "WGS-AD,differentially expressed variants across 5 brain regions")

pheatmap(msbb_array19.AllSamples_select_brainRegion[rownames(msbb_array19.AllSamples_select_brainRegion)%in%unlist(lapply(strsplit(rownames(annotation_row_df[annotation_row_df$VariantType=="RegionAD",]),split = "\\."),`[[`,1)),],annotation = annotation_df,annotation_row = annotation_row_df,annotation_names_row = T,fontsize_row = 4,fontsize_col = 4,main = "WGS-AD,differentially expressed variants across 5 brain regions")

#Compare ranked variant genes with published datasets
setwd("/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/")
pub_datasets=vector(mode = "list",length = 7)
tmp1=read.xls("Hippocampus_Proteome_AD_Hondius_etal/mmc2.xlsx",sheet = 1,header=T)
tmp2=read.xls("Hippocampus_Proteome_AD_Hondius_etal/mmc2.xlsx",sheet = 2,header=T)
#Up/Down regulated proteins stored at 1st and 2nd element of the list respectively, from Hondius et al.
pub_datasets[[1]][[1]]=droplevels(tmp1$Gene..leading.protein.)
pub_datasets[[1]][[2]]=droplevels(tmp2$Gene..leading.protein.)
names(pub_datasets[[1]])=c("Increased_AD_Proteins","Decresed_AD_Proteins")
#AD associated proteins from the Drummond et al. proteomics dataset

#GW-hydroxymethylation from Zhao et al. (David Bennet), differentially hydroxymethylated regions with NP and NFT respectively
tmp3=read.table("GW-DNA-hydroxymethylation-Bennet/1-s2.0-S155252601633059X-mmc4_NP.txt",sep = "\t",header = T,as.is = T)
tmp4=read.table("GW-DNA-hydroxymethylation-Bennet/1-s2.0-S155252601633059X-mmc5_NFT.txt",sep = "\t",header = T,as.is = T)
#GW-5hmC from Bernstein et al., differentially hydroxymethylated genes, concurrently diff. expressed
tmp4.2=read.xls("5hMC_Tau_AD/ddw109_Supp/Supplemental Table 10_Bernstein et al_R1.xlsx",header=F)
pub_datasets[[2]][[1]]=tmp3$Neartest.gene
pub_datasets[[2]][[2]]=tmp4$Nearest.gene
names(pub_datasets[[2]])=c("Diff_NP","Diff_NFT")
#Plasma-proteomics dataset from Jaeger et al., plasma proteins strongly correlated with MMSE
tmp5=read.table("PlasmaProteomics_AD_Jaeger/Jaeger_PlasmaProt_MMSEcorrelation.txt",sep = "\t",header = T,as.is = T)
pub_datasets[[3]][[1]]=tmp5%>%filter(q.value<=0.1)%>%pull(HUPO.Gene.Symbol)
names(pub_datasets[[3]])=c("PlasmaProt_MMSE")
#Zhang TYROBP signature paper, eQTL modules
tmp6=read.xls("/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/Zhang_TYROBP/mmc1.xlsx",sheet = 4,as.is=T)
tmp7=read.xls("/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/Zhang_TYROBP/mmc1.xlsx",sheet = 2,as.is=T)
lapply(unique(zhang_tyrobp_modules$Module),function(x)zhang_tyrobp_modules%>%filter(Module==x)%>%pull(Gene_Symbol))
#Initialise 2nd level list to include all cis-eQTL modules
pub_datasets[[4]]=vector(mode = "list",length = length(unique(tmp6$Module)))
names(pub_datasets[[4]])=unique(tmp6$Module)
pub_datasets[[4]]=lapply(unique(tmp6$Module), function(x)tmp6%>%filter(Module==x)%>%pull(Gene_Symbol))
names(pub_datasets[[4]])=unique(tmp6$Module)
pub_datasets[[4]]=pub_datasets[[4]][-40]

zhang_module.size=lapply(pub_datasets[[4]],length)
names(zhang_module.size)= names(pub_datasets[[4]])
zhang_module.exp_SNPs=lapply(zhang_module.size,function(x)round(length(tanzi_ranked_genes.union)/20805*x,digits = 1))
zhang_module.SNP_overlap=lapply(pub_datasets[[4]],intersect,tanzi_ranked_genes.union)
#zhang_module.SNP_overlap2=zhang_module.SNP_overlap[lapply(zhang_module.SNP_overlap,length)>2]
fisher_results.zhang_modules=vector(mode = "list",length = length(zhang_module.SNP_overlap))
for(m in 1:length(zhang_module.SNP_overlap)){
  fisher_results.zhang_modules[[m]]=fisher.test(matrix(c(length(zhang_module.SNP_overlap[[m]]),
                                                length(union(tanzi_ranked_genes.AD,tanzi_ranked_genes.Control))-length(zhang_module.SNP_overlap[[m]]),
                                                zhang_module.size[[m]]-length(zhang_module.SNP_overlap[[m]]),
                                                20805-(length(tanzi_ranked_genes.union)-length(zhang_module.SNP_overlap[[m]]))-
                                                zhang_module.size[[m]]-length(zhang_module.SNP_overlap[[m]])),nrow = 2))
  
}

res.zhang_modules=foreach(i=1:length(zhang_module.SNP_overlap))%do%{
  
  
  data.frame(Variants_in_Module=unlist(lapply(zhang_module.SNP_overlap[i],paste,collapse=",")),
             ModuleSize=unlist(zhang_module.size[[i]]),
             Expected_Variants=unname(unlist(zhang_module.exp_SNPs[i])),
             Overlap_Size=unlist(lapply(zhang_module.SNP_overlap[i],length)),
             Fisher_statistic=unname(unlist(lapply(fisher_results.zhang_modules[i],function(s)s$estimate))),
             pval=unname(unlist(lapply(fisher_results.zhang_modules[i],function(p)p$p.value))),
             stringsAsFactors = F)
}
names(res.zhang_modules)=names(zhang_module.SNP_overlap)
res.zhang_modules_df=data.frame(rbindlist(res.zhang_modules),stringsAsFactors = F)
rownames(res.zhang_modules_df)=names(zhang_module.SNP_overlap)
res.zhang_modules_df$adj.p=p.adjust(p = res.zhang_modules_df$pval,method = "fdr")
res.zhang_modules_df=rownames_to_column(res.zhang_modules_df,"ModuleName")

res.zhang_modules_df=res.zhang_modules_df%>%filter(adj.p<=0.1)
res.zhang_modules_df=res.zhang_modules_df%>%mutate(Representation=if_else(Overlap_Size>Expected_Variants,true = "Over",false = "Under"))%>%filter(Representation=="Over")

res.zhang_modules_intermod.subgraph=Reduce(graph.union,lapply(pub_datasets[[4]][res.zhang_modules_df$ModuleName],function(x)induced.subgraph(graph = hs_ppi,vids = V(hs_ppi)$name[which(V(hs_ppi)$name%in%x)])))
Reduce(graph.union,lapply(lapply(pub_datasets[[4]][res.zhang_modules_df$ModuleName],function(x)neighborhood(graph = hs_ppi,order = 1,nodes = V(hs_ppi)$name[which(V(hs_ppi)$name%in%x)],mindist = 1)[[1]]$name),function(y)induced.subgraph(graph = hs_ppi,vids = V(hs_ppi)$name[which(V(hs_ppi)$name%in%y)])))

V(res.zhang_modules_intermod.subgraph)$colour="#FAFAFA"
V(res.zhang_modules_intermod.subgraph)$colour[which(V(res.zhang_modules_intermod.subgraph)$name%in%pub_datasets[[4]][res.zhang_modules_df$ModuleName][[1]])]="#dfcadf"
V(res.zhang_modules_intermod.subgraph)$colour[which(V(res.zhang_modules_intermod.subgraph)$name%in%pub_datasets[[4]][res.zhang_modules_df$ModuleName][[2]])]="#ff8bff"
V(res.zhang_modules_intermod.subgraph)$colour[which(V(res.zhang_modules_intermod.subgraph)$name%in%pub_datasets[[4]][res.zhang_modules_df$ModuleName][[3]])]="#f2eaa0"
V(res.zhang_modules_intermod.subgraph)$colour[which(V(res.zhang_modules_intermod.subgraph)$name%in%pub_datasets[[4]][res.zhang_modules_df$ModuleName][[4]])]="#c5c5c5"
write.graph(delete_vertices(graph = res.zhang_modules_intermod.subgraph,v = names(which(degree(res.zhang_modules_intermod.subgraph)<=1))),"Zhang_SNPenriched_InterMod_Subnet.gml",format="gml")
fwrite(res.zhang_modules_df,"/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Zhang_CoexpModules_RankAD_SNPgene_enriched.txt",sep="\t",col.names = T,row.names = F)


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
zhang_celltype_ADgenes.KEGG=lapply(zhang_celltype_ADgenes.list,function(x)enrichr(genes = x,databases = kegg_dbs)[[1]]%>%filter(Adjusted.P.value<0.1))
zhang_celltype_tanzi_intercell.subgraph=Reduce(graph.union,lapply(zhang_celltype_ADgenes.list[res.tanzi_celltype.df$CellType],function(x)induced.subgraph(graph = hs_ppi,vids = V(hs_ppi)$name[which(V(hs_ppi)$name%in%x)])))
zhang_celltype_tanzi_intercell.subgraph=delete_vertices(graph = zhang_celltype_tanzi_intercell.subgraph,v = names(which(degree(zhang_celltype_tanzi_intercell.subgraph)<=1)))


V(zhang_celltype_tanzi_intercell.subgraph)$colour="#FAFAFA"
V(zhang_celltype_tanzi_intercell.subgraph)$colour[which(V(zhang_celltype_tanzi_intercell.subgraph)$name%in%zhang_celltype_ADgenes.list$Astrocytes)]="#FFF176"
V(zhang_celltype_tanzi_intercell.subgraph)$colour[which(V(zhang_celltype_tanzi_intercell.subgraph)$name%in%zhang_celltype_ADgenes.list$Microglia)]="#AED581"
V(zhang_celltype_tanzi_intercell.subgraph)$colour[which(V(zhang_celltype_tanzi_intercell.subgraph)$name%in%zhang_celltype_ADgenes.list$Neurons)]="#BA68C8"
write.graph(zhang_celltype_tanzi_intercell.subgraph,"Zhang_CellType_Tanzi_Intercell_Subnet.gml",format = "gml")


#trinitya subnetwork
tanzi_trinity.seeds=c(unname(unlist(Filter(f = function(x)length(x)>0,x = lapply(zhang_celltype_ADgenes.list,intersect,tanzi_DEV)))),known_AD_genes)
tanzi_trinity_seeds.KEGG=enrichr(genes = tanzi_trinity.seeds,databases = kegg_dbs)[[1]]%>%filter(Adjusted.P.value<0.1)%>%dplyr::select(c(Term,Genes))

write(unname(unlist(Filter(f = function(x)length(x)>0,x = lapply(zhang_celltype_ADgenes.list,intersect,tanzi_DEV)))),"./Gene_Lists/Tanzi_trinitya_Seed_genes.txt",sep = "\n")
tanzi_trinity.Seed_df=data.frame(V(hs_ppi)$name,rep(0,length(V(hs_ppi)$name)),stringsAsFactors=F)
tanzi_trinity.Seed_df[grep(pattern=paste(tanzi_trinity.seeds,collapse="|"),x=tanzi_trinity.Seed_df[,1]),2]=1

colnames(tanzi_trinity.Seed_df)=c()
tanzi_trinity.pData=xPierGenes(data=tanzi_trinity.Seed_df,network.customised=hs_ppi,restart=0.9,verbose=T)
tanzi_trinity.pSubnet=xPierSubnet(pNode=tanzi_trinity.pData,network.customised=hs_ppi,subnet.significance=0.01,verbose=T)

V(tanzi_trinity.pSubnet)$colour="#FAFAFA"
V(tanzi_trinity.pSubnet)$colour[which(V(tanzi_trinity.pSubnet)$name%in%tanzi_trinity.seeds)]="#A9DFBF"
V(tanzi_trinity.pSubnet)$colour[which(V(tanzi_trinity.pSubnet)$name%in%known_AD_genes)]="#F5B7B1"
write.graph(tanzi_trinity.pSubnet,"Tanzi_trinitya_Subgraph.gml",format = "gml")


tanzi_trinity_pSubnet.KEGG=enrichr(genes = V(tanzi_trinity.pSubnet)$name,databases = kegg_dbs)[[1]]%>%filter(Adjusted.P.value<0.1)
tanzi_trinity_pSubnet.KEGG$KEGGID=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = tanzi_trinity_pSubnet.KEGG$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[2]]))
tanzi_trinity_pSubnet.KEGG$Pathway=unlist(gsub(pattern = "_Homo sapiens_",replacement = "_",x = tanzi_trinity_pSubnet.KEGG$Term)%>%strsplit(split = "_")%>%lapply(function(x)x[[1]]))
tanzi_trinity_pSubnet.KEGG$knownAD=unlist(lapply(lapply(tanzi_trinity_pSubnet.KEGG$Genes,function(x)strsplit(x = x,split=";")[[1]]),function(y)paste(intersect(y,known_AD_genes),collapse = ";")))
tanzi_trinity_pSubnet.KEGG$DEG=unlist(lapply(lapply(tanzi_trinity_pSubnet.KEGG$Genes,function(x)strsplit(x = x,split=";")[[1]]),function(y)paste(intersect(y,rownames(deg_df)),collapse = ";")))
tanzi_trinity_pSubnet.KEGG$CTM=unlist(lapply(lapply(tanzi_trinity_pSubnet.KEGG$Genes,function(x)strsplit(x = x,split=";")[[1]]),function(y)paste(intersect(y,unname(unlist(zhang_celltype_ADgenes.list[c(1,3:4)]))),collapse = ";")))
write.table(tanzi_trinity_pSubnet.KEGG%>%select(c(Pathway,Genes,knownAD,DEG,CTM,P.value,Adjusted.P.value)),"Trinity_pSubnet_KEGG.txt",sep = "\t",col.names = T,row.names = F)

glut_synapse_KEGG_genes.expr=c(rep(1,length(strsplit(tanzi_ranked_genes.KEGG[7,9],split = ";")[[1]])),rep(-1,length(strsplit(knownAD_cons.KEGG[26,9],split = ";")[[1]])))
names(pld_KEGG_genes.expr)=c(strsplit(tanzi_ranked_genes.KEGG[7,9],split = ";")[[1]],strsplit(knownAD_cons.KEGG[26,9],split = ";")[[1]])
plotKegg(keggIDs = "hsa04072",geneExpr = pld_KEGG_genes.expr)


#PCXN analysis
pcxn_network2 <- function (object) {
  edges <- c()
  edge_thickness <- c()
  edge_color <- c()
  for (i in 1:nrow(object@data)) {
    edges <- c(edges, object@data[i, 1])
    edges <- c(edges, object@data[i, 2])
    edge_thickness <- c(edge_thickness, abs(as.numeric(object@data[i, 
                                                                   3]) * 10))
    if (as.numeric(object@data[i, 3]) < 0) 
      edge_color <- c(edge_color, "blue")
    else edge_color <- c(edge_color, "red")
  }
  graph <- igraph::make_graph(edges, directed = FALSE)
  if (object@type == "explore") {
    V(graph)$color <- "#ffaf1a"
    V(graph)[object@geneset_groups$query_geneset]$color[1] <- "#66b2ff"
  }
  if (object@type == "analyze") {
    V(graph)$color <- "white"
    for (i in 1:length(object@geneset_groups$top_correlated_genesets)) {
      tryCatch({
        cgs <- object@geneset_groups$top_correlated_genesets[i]
        V(graph)[cgs]$color <- "#ffaf1a"
      }, error = function(e) {
        print(paste("No edges containing ", cgs, " pass the correlation and p-value filters therefore\n                            the node is not included in the network", 
                    sep = ""))
      })
    }
    for (j in 1:length(object@geneset_groups$phenotype_0_genesets)) {
      tryCatch({
        cgs1 <- object@geneset_groups$phenotype_0_genesets[j]
        V(graph)[cgs1]$color <- "#66b2ff"
      }, error = function(e) {
        print(paste("No edges containing ", cgs1, " pass the correlation and p-value filters therefore\n                            the node is not included in the network", 
                    sep = ""))
      })
    }
    for (h in 1:length(object@geneset_groups$phenotype_1_genesets)) {
      tryCatch({
        cgs2 <- object@geneset_groups$phenotype_1_genesets[h]
        V(graph)[cgs2]$color <- "#19c67e"
      }, error = function(e) {
        print(paste("No edges containing ", cgs2, " pass the correlation and p-value filters therefore\n                            the node is not included in the network", 
                    sep = ""))
      })
    }
  }
  E(graph)$colour=edge_color
  E(graph)$width=edge_thickness
  #plot <- igraph::tkplot(graph, layout = igraph::layout_in_circle, 
  #                       width = NULL, height = NULL, edge.width = edge_thickness, 
   #                      edge.color = edge_color)
  #res=list(CorrGraph=graph,CorrPlot=plot)
  tna
}

tanzi_region_SNPs.KEGG_Enrichr=scan("Work/Collaborations/Tanzi_WGS/KEGG_Regional_SNPs_EnrichrOutput_pathwayList.txt",sep = "\n",what = "char")
tanzi_region_SNPs.KEGG_EnrichrTable=read.table("Work/Collaborations/Tanzi_WGS/KEGG_Ranked_SNPs_EnrichrOutput.txt",sep = "\t",header = T,as.is = T)
tanzi_ranked_SNPs.KEGG_Enrichr=scan("Work/Collaborations/Tanzi_WGS/KEGG_Ranked_SNPs_EnrichrOutput_pathwayList.txt",sep = "\n",what = "char")
pcxn_net=pcxn_network2(pcxn_analyze(collection = "pathprint",phenotype_0_genesets = grep(pattern = paste(tanzi_region_SNPs.KEGG_Enrichr,collapse = "|"),x = names(pathprint.Hs.gs),value = T),phenotype_1_genesets = grep(pattern = paste(tanzi_ranked_SNPs.KEGG_Enrichr,collapse = "|"),x = names(pathprint.Hs.gs),value = T)))
pcxn_out=pcxn_analyze(collection = "pathprint",phenotype_0_genesets = grep(pattern = paste(tanzi_region_SNPs.KEGG_Enrichr,collapse = "|"),x = names(pathprint.Hs.gs),value = T),phenotype_1_genesets = grep(pattern = paste(tanzi_ranked_SNPs.KEGG_Enrichr,collapse = "|"),x = names(pathprint.Hs.gs),value = T))
pcxn_out.df=data.frame(pcxn_out@data,stringsAsFactors = F)
pcxn_out.df_select=pcxn_out.df[union(which(pcxn_out.df$Pathway.A%in%union(pcxn_out@geneset_groups$phenotype_1_genesets,pcxn_out@geneset_groups$phenotype_0_genesets)),which(pcxn_out.df$Pathway.B%in%union(pcxn_out@geneset_groups$phenotype_1_genesets,pcxn_out@geneset_groups$phenotype_0_genesets))),]

pcxn_common_pathways=intersect(pcxn_out@geneset_groups$phenotype_0_genesets,pcxn_out@geneset_groups$phenotype_1_genesets)
pcxn_out.df_select_rankedSNP=pcxn_out.df[union(which(pcxn_out.df$Pathway.A%in%pcxn_out@geneset_groups$phenotype_1_genesets),which(pcxn_out.df$Pathway.B%in%pcxn_out@geneset_groups$phenotype_1_genesets)),]
pcxn_out.df_select_regionSNP=pcxn_out.df[union(which(pcxn_out.df$Pathway.A%in%pcxn_out@geneset_groups$phenotype_0_genesets),which(pcxn_out.df$Pathway.B%in%pcxn_out@geneset_groups$phenotype_0_genesets)),]
pcxn_common_pathways.select_commonPathways=pcxn_out.df[union(which(pcxn_out.df$Pathway.A%in%pcxn_common_pathways),which(pcxn_out.df$Pathway.B%in%pcxn_common_pathways)),]
write.table(pcxn_common_pathways.select_commonPathways,"Work/Collaborations/Tanzi_WGS/PCxN_Common_PathwayCorrTable.txt",sep = "\t",col.names = T,row.names = F,quote = F)
write.table(pcxn_out.df_select,"Work/Collaborations/Tanzi_WGS/PCxN_Ranked_Regional_AD_PathwayCorrTable.txt",sep = "\t",col.names = T,row.names = F,quote = F)
write.table(pcxn_out.df_select_rankedSNP,"Work/Collaborations/Tanzi_WGS/PCxN_RankedSNP_AD_PathwayCorrTable.txt",sep = "\t",col.names = T,row.names = F,quote = F)
write.table(pcxn_out.df_select_regionSNP,"Work/Collaborations/Tanzi_WGS/PCxN_RegionalSNP_AD_PathwayCorrTable.txt",sep = "\t",col.names = T,row.names = F,quote = F)
write.graph(pcxn_net,"Work/Collaborations/Tanzi_WGS/PCxN_Ranked_Regional_AD.gml",format="gml")
write(grep(pattern = paste(tanzi_region_SNPs.KEGG_Enrichr,collapse = "|"),x = names(pathprint.Hs.gs),value = T),"Work/Collaborations/Tanzi_WGS/KEGG_Regional_SNPs_EnrichrOutput_pathwayList_PCXNMapped.txt",sep = "\n")
write(grep(pattern = paste(tanzi_ranked_SNPs.KEGG_Enrichr,collapse = "|"),x = names(pathprint.Hs.gs),value = T),"Work/Collaborations/Tanzi_WGS/KEGG_Ranked_SNPs_EnrichrOutput_pathwayList_PCXNMapped.txt",sep = "\n")
