library(org.Hs.eg.db)
library(dplyr)
library(magrittr)
library(data.table)
library(parallel)
#library(clusterProfiler)
library(doParallel)
library(diffcoexp)
library(DCGL)
library(enrichR)

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

known_AD_genes=scan("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Gene_Lists/Known_AD_genes_V2.txt",what = "char",sep = "\n")
known_AD_genes.cons=scan("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Gene_Lists/Known_AD_genes_ConservativeList.txt",what = "char",sep = "\n")
tanzi_ranked_genes=read.table("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Gene_Lists/All_genes_ranked.txt",header = T,sep = "\t",stringsAsFactors = F)
tanzi_ranked_genes.great=read.table("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Gene_Lists/Ranked_genes_corrected_allele_top_4000_SNPs_genes_from_great.txt",header = T,sep = "\t",stringsAsFactors = F)
tanzi_ranked_genes.great$SNP_Rank=gsub(pattern = "SNP_rank_",replacement = "",x = tanzi_ranked_genes.great$SNP_Rank)
tanzi_ranked_genes.AD=unique(tanzi_ranked_genes.great%>%mutate(SNP_Rank=as.integer(gsub(pattern = "SNP_rank_",replacement = "",x = SNP_Rank)))%>%filter(Corrected_association=="Case"&SNP_Rank<=1000)%>%pull(Gene))
tanzi_ranked_genes.Control=unique(tanzi_ranked_genes.great%>%mutate(SNP_Rank=as.integer(gsub(pattern = "SNP_rank_",replacement = "",x = SNP_Rank)))%>%filter(Corrected_association=="Control"&SNP_Rank<=1000)%>%pull(Gene))
tanzi_ranked_genes.union=union(tanzi_ranked_genes.AD,tanzi_ranked_genes.Control)


setwd("/shared/hidelab2/user/md4zsa/Work/Data/MSMM_Array/MSBB_Array19/Normalised_Data")
msbb_array19.files=list.files(pattern = "*.tsv",full.names = T)
msbb_array19.files=msbb_array19.files[-c(18:19)]
msbb_array19=mapply(FUN = fread,msbb_array19.files,MoreArgs = list(header=T,sep="\t",data.table=F,showProgress=T))
names(msbb_array19)=c("AC","CN","DLPFC","FP","HP","IFG","ITG","MTG","OVC","PHG","PCC","PCG","PFC","PTMN","SPL","STG","TP")
msbb_array19.covariates=read.delim2("AMP-AD_MSBB_MSSM_covariates_mRNA_AffymetrixU133AB.txt",header = T,as.is = T)%>%
  mutate(Age=replace(Age,Age=="89+",90))%>%
  mutate(pH=as.numeric(pH))%>%
  mutate(CDR=as.numeric(CDR))%>%
  mutate(PLQ_Mn=as.numeric(PLQ_Mn))%>%
  mutate(BrainBank=paste("X",BrainBank,sep=""))
msbb_array19.2=lapply(msbb_array19[-18],function(x){rownames(x)<-x$ID;x})
#msbb_array19.3=lapply(msbb_array19.2,function(x)x[,-c(1:4)])
multiID.u133a=msbb_array19.2[[1]]$ENTREZ_GENE_ID[grep(pattern="///",msbb_array19.2[[1]]$ENTREZ_GENE_ID)]
msbb_array19.2.agg=lapply(msbb_array19.2[-18],function(y)aggregate(x=y[,-c(1:4)],by=list(EntrezID=y$ENTREZ_GENE_ID),mean))
msbb_array19.2.agg2=lapply(msbb_array19.2.agg,function(x){rownames(x)<- x$EntrezID;x <- x[-grep(pattern = "///",x = rownames(x)),];x})
msbb_array19.2.agg2=readRDS("msbb_array19_Agg.RDS")
msbb_array19.covariates$SampleType="OTHER"

msbb_array19.covariates$SampleType[(msbb_array19.covariates$CDR<=0.5&msbb_array19.covariates$NP1<=1&msbb_array19.covariates$Braak<=3)]="CONTROL"
msbb_array19.covariates$SampleType[(msbb_array19.covariates$CDR>=1&msbb_array19.covariates$NP1>=2&msbb_array19.covariates$Braak>=4)]="AD"

AD_sample.vector=paste(gsub(pattern = "X",replacement = "",msbb_array19.covariates$BrainBank[msbb_array19.covariates$SampleType=="AD"]),collapse = "|")
Control_sample.vector=paste(gsub(pattern = "X",replacement = "",msbb_array19.covariates$BrainBank[msbb_array19.covariates$SampleType=="CONTROL"]),collapse = "|")

#regnet_tf2target.HGNC=fread("/shared/hidelab2/user/md4zsa/Work/Data/TF_Databases/RegNetwork_human_regulators2.txt",header = T,sep = "\t",showProgress = T,data.table = F)%>%dplyr::filter(evidence=="Experimental")%>%dplyr::select(c(regulator_symbol,target_symbol))
regnet_tf2target_unedit=fread("/Users/sandeepamberkar/Work/Data/TF_Databases/RegNetwork_human_regulators2.txt",header = T,sep = "\t",showProgress = T,data.table = F)%>%filter(evidence=="Experimental")
regnet_tf2target=fread("/Users/sandeepamberkar/Work/Data/TF_Databases/RegNetwork_human_regulators2.txt",header = T,sep = "\t",showProgress = T,data.table = F)%>%filter(evidence=="Experimental")%>%dplyr::select(c(regulator_id,target_id))
regnet_tf2target.HGNC=fread("/Users/sandeepamberkar/Work/Data/TF_Databases/RegNetwork_human_regulators2.txt",header = T,sep = "\t",showProgress = T,data.table = F)%>%filter(evidence=="Experimental")%>%dplyr::select(c(regulator_symbol,target_symbol))
#Regroup brain regions by lobes
# lobe_bm_area.map=matrix(nrow=17,ncol=2)
# lobe_bm_area.map[1,]=c("Frontal_Lobe","FP")
# lobe_bm_area.map[2,]=c("Frontal_Lobe","AC")
# lobe_bm_area.map[3,]=c("Frontal_Lobe","PFC")
# lobe_bm_area.map[4,]=c("Occipetal_Lobe","OVC")
# lobe_bm_area.map[5,]=c("Temporal_Lobe","ITG")
# lobe_bm_area.map[6,]=c("Temporal_Lobe","MTG")
# lobe_bm_area.map[7,]=c("Temporal_Lobe","STG")
# lobe_bm_area.map[8,]=c("Parietal_Lobe","PCC")
# lobe_bm_area.map[9,]=c("Temporal_Lobe","PHG")
# lobe_bm_area.map[10,]=c("Temporal_Lobe","TP")
# lobe_bm_area.map[11,]=c("Frontal_Lobe","PCG")
# lobe_bm_area.map[12,]=c("Frontal_Lobe","IFG")
# lobe_bm_area.map[13,]=c("Frontal_Lobe","DLPFC")
# lobe_bm_area.map[14,]=c("Parietal_Lobe","SPL")
# lobe_bm_area.map[16,]=c("Temporal_Lobe","HP")
# lobe_bm_area.map[15,]=c("Dorsal_striatum","PTMN")
# lobe_bm_area.map[17,]=c("Dorsal_striatum","CN")
# lobe_bm_area.map=data.frame(Lobe=lobe_bm_area.map[,1],BrainRegion=lobe_bm_area.map[,2],stringsAsFactors = F)
# 
# 
# msbb_array.byLobe=msbb_array.byLobe2=vector(mode = "list",length = length(unique(lobe_bm_area.map$Lobe)))
# names(msbb_array.byLobe)=names(msbb_array.byLobe2)=unique(lobe_bm_area.map$Lobe)
# msbb_array.byLobe=foreach(i=1:length(names(msbb_array.byLobe)))%dopar%{
#   array_byLobe=do.call("cbind",msbb_array19.2.agg2[which(names(msbb_array19.2.agg2)%in%lobe_bm_area.map$BrainRegion[lobe_bm_area.map$Lobe==names(msbb_array.byLobe)[i]])])  
# } 

msbb_array_region.diffcoexp=vector(mode = "list",length = 17)
names(msbb_array_region.diffcoexp)=names(msbb_array19)
#msbb_array_lobe.diffcoexp=vector(mode = "list",length = 5)
for(i in 6:17){
  msbb_array_region.diffcoexp[[i]]=diffcoexp(exprs.1 = as.data.frame(msbb_array19.2.agg2[[i]])[1:500,grep(pattern = Control_sample.vector,x = colnames(as.data.frame(msbb_array19.2.agg2[[i]])))],
                 exprs.2 = as.data.frame(msbb_array19.2.agg2[[i]])[1:500,grep(pattern = AD_sample.vector,x = colnames(as.data.frame(msbb_array19.2.agg2[[i]])))],
                 rth=0.6, qth=0.2, r.diffth=0.1, q.diffth=0.1)
  saveRDS(msbb_array_region.diffcoexp[[i]],paste(names(msbb_array_region.diffcoexp)[i],"diffcoexp.RDS",sep="_"))
}

setwd("./diffcoexp_results/")
diffcoexp_results_files=list.files(path = ".",pattern = "diffcoexp")
for(f in 1:17){
  msbb_array_region.diffcoexp[[f]]=readRDS(diffcoexp_results_files[f])  
}

names(msbb_array_region.diffcoexp)=names(msbb_array19)
msbb_array.DCG=lapply(msbb_array_region.diffcoexp,function(y)y$DCGs)
msbb_array.DCL=lapply(msbb_array_region.diffcoexp,function(y)y$DCLs)
msbb_array.DRsort=vector(mode = "list",length = 17)
names(msbb_array.DRsort)=names(msbb_array19)
msbb_array.DRsort=mapply(FUN = function(a,b)DRsort(DCGs = a,DCLs = b,tf2target = regnet_tf2target,expGenes = rownames(msbb_array19.2.agg2$AC)),msbb_array.DCG,msbb_array.DCL)

msbb_array.DCG=Filter(f = function(y)dim(y)[1]>0,x = msbb_array.DCG)
msbb_array.DCL=lapply(msbb_array.DCL,function(x)x%>%mutate(Gene.1.HGNC=unname(mapIds(x = org.Hs.eg.db,keys = x$Gene.1,keytype = "ENTREZID",column = "SYMBOL",multiVals = "first")),
                                                             Gene.2.HGNC=unname(mapIds(x = org.Hs.eg.db,keys = x$Gene.2,keytype = "ENTREZID",column = "SYMBOL",multiVals = "first")),
                                                             Gene.Pair.HGNC=paste(Gene.1.HGNC,Gene.2.HGNC,sep = ","))%>%select(c(Gene.1,Gene.2,Gene.1.HGNC,Gene.2.HGNC,Gene.Pair.HGNC,cor.1,cor.2,p.1,p.2,p.diffcor,q.1,q.2,q.diffcor,cor.diff,type)))
msbb_array.DCG=lapply(msbb_array.DCG,function(x)x%>%mutate(Gene.HGNC=unname(mapIds(x = org.Hs.eg.db,keytype = "ENTREZID",column = "SYMBOL",keys = Gene,multiVals = "first")))%>%dplyr::select(Gene,Gene.HGNC,CO.links,DC.links,DCL.same,DCL.diff,DCL.switch,p,q))
msbb_array_DCG.filtered=lapply(msbb_array.DCG,function(y)y%>%pull(Gene))
msbb_array_DCG.filtered2=lapply(msbb_array.DCG,function(y)y%>%pull(Gene.HGNC))
msbb_array.DCL=lapply(msbb_array_region.diffcoexp,function(y)y$DCLs)
msbb_array.DCL=lapply(msbb_array.DCL,function(x)x%>%mutate(Gene.1.HGNC=unname(mapIds(x = org.Hs.eg.db,keys = x$Gene.1,keytype = "ENTREZID",column = "SYMBOL",multiVals = "first")),
                                                           Gene.2.HGNC=unname(mapIds(x = org.Hs.eg.db,keys = x$Gene.2,keytype = "ENTREZID",column = "SYMBOL",multiVals = "first")),
                                                           Gene.Pair.HGNC=paste(Gene.1.HGNC,Gene.2.HGNC,sep = ","))%>%dplyr::select(c(Gene.1,Gene.2,Gene.1.HGNC,Gene.2.HGNC,Gene.Pair.HGNC,cor.1,cor.2,p.1,p.2,p.diffcor,q.1,q.2,q.diffcor,cor.diff,type)))

msbb_array.DCL_types=list()
msbb_array.DCL_types$same_signed=lapply(msbb_array.DCL,function(y)y%>%filter(type=="same signed"))
msbb_array.DCL_types$diff_signed=lapply(msbb_array.DCL,function(y)y%>%filter(type=="diff signed"))
msbb_array.DCL_types$switched_opposites=lapply(msbb_array.DCL,function(y)y%>%filter(type=="switched opposites"))

msbb_array.DCL_types$same_signed=Filter(f = function(x)dim(x)[1]>2,msbb_array.DCL_types$same_signed)
msbb_array.DCL_types$diff_signed=Filter(f = function(x)dim(x)[1]>2,msbb_array.DCL_types$diff_signed)
msbb_array.DCL_types$switched_opposites=Filter(f = function(x)dim(x)[1]>2,msbb_array.DCL_types$switched_opposites)

msbb_array_top10_DCL.same_signed=sort(table(unlist(lapply(msbb_array.DCL_types$same_signed,function(x)x$Gene.Pair))),decreasing = T)[1:10]
names(msbb_array_top10_DCL.same_signed)=unlist(lapply(lapply(names(msbb_array_top10_DCL.same_signed),strsplit,split=","),function(x)paste(mapIds2(IDs = x[[1]],IDFrom = "ENTREZID",IDTo = "SYMBOL")[[1]][,2],collapse = ",")))
msbb_array_top10_DCL.diff_signed=sort(table(unlist(lapply(msbb_array.DCL_types$diff_signed,function(x)x$Gene.Pair))),decreasing = T)[1:10]
names(msbb_array_top10_DCL.diff_signed)=unlist(lapply(lapply(names(msbb_array_top10_DCL.diff_signed),strsplit,split=","),function(x)paste(mapIds2(IDs = x[[1]],IDFrom = "ENTREZID",IDTo = "SYMBOL")[[1]][,2],collapse = ",")))
msbb_array_top10_DCL.switched_opposites=sort(table(unlist(lapply(msbb_array.DCL_types$switched_opposites,function(x)x$Gene.Pair))),decreasing = T)[1:10]
names(msbb_array_top10_DCL.switched_opposites)=unlist(lapply(lapply(names(msbb_array_top10_DCL.switched_opposites),strsplit,split=","),function(x)paste(mapIds2(IDs = x[[1]],IDFrom = "ENTREZID",IDTo = "SYMBOL")[[1]][,2],collapse = ",")))

msbb_array_total_DCGs=unlist(lapply(msbb_array_DCG.filtered,length))
msbb_array_total_DCLs=unlist(lapply(msbb_array.DCL,function(x)dim(x)[1]))
dat1=data.frame(DCLs=msbb_array_total_DCLs,DCGs=msbb_array_total_DCGs,Region=Region)
dat1m=melt(dat1[,c('Region','DCLs','DCGs')],id.vars=1)
ggplot(dat1m,aes(x = Region,y = value)) + geom_bar(aes(fill = variable),stat = "identity",position = "dodge")+labs(title="MSBB array - DCGs and DCLs")+scale_y_log10()+theme(title = element_text(face = "bold",size = 15,colour = "black"),axis.title.x = element_text(face = "bold",size = 15,colour = "black"),axis.text.x = element_text(size = 15,colour = "black"),legend.text = element_text(size = 15))
# msbb_array.DCL_types_KEGG=list()
# msbb_array.DCL_types_KEGG$same_signed=lapply(msbb_array.DCL_types$same_signed,function(y)enrichr(genes = mapIds2(IDs = V(graph.data.frame(d = y[,c(1:2)],directed = F))$name,IDFrom = "ENTREZID",IDTo = "SYMBOL")[[1]][,2],databases = dbs_kegg2016)$KEGG_2016%>%filter(Adjusted.P.value<=0.1)%>%pull(Term))
# msbb_array.DCL_types_KEGG$diff_signed=lapply(msbb_array.DCL_types$diff_signed,function(y)enrichr(genes = mapIds2(IDs = V(graph.data.frame(d = y[,c(1:2)],directed = F))$name,IDFrom = "ENTREZID",IDTo = "SYMBOL")[[1]][,2],databases = dbs_kegg2016)$KEGG_2016%>%filter(Adjusted.P.value<=0.1)%>%pull(Term))
# msbb_array.DCL_types_KEGG$switched_opposites=lapply(msbb_array.DCL_types$switched_opposites,function(y)enrichr(genes = mapIds2(IDs = V(graph.data.frame(d = y[,c(1:2)],directed = F))$name,IDFrom = "ENTREZID",IDTo = "SYMBOL")[[1]][,2],databases = dbs_kegg2016)$KEGG_2016%>%filter(Adjusted.P.value<=0.1)%>%pull(Term))

saveRDS(msbb_array.DCG,"../msbb_array_DCG.RDS")
saveRDS(msbb_array.DCL,"../msbb_array_DCL.RDS")
# names(msbb_array.diffcoexp)=names(msbb_array.byLobe)
msbb_array.DRsort=vector(mode = "list",length = 17)
names(msbb_array.DRsort)=names(msbb_array19)


msbb_array.DRsort=Filter(f = function(x)length(x)>0,x = msbb_array.DRsort)

msbb_array.DRL=lapply(msbb_array.DRsort,function(x)x$DRLs%>%dplyr::filter(q.diffcor<=0.1&!is.na(common.TF)&!grepl(pattern = ";",x=common.TF))%>%dplyr::select(c(common.TF,Gene.1,Gene.2)))
msbb_array.DRL=lapply(msbb_array.DRL,function(x)x%>%mutate(common.TF.HGNC=mapIds2(IDs = as.character(common.TF),IDFrom = "ENTREZID",IDTo = "SYMBOL")[[1]][,2],Gene.1.HGNC=mapIds2(IDs = as.character(Gene.1),IDFrom = "ENTREZID",IDTo = "SYMBOL")[[1]][,2],Gene.2.HGNC=mapIds2(IDs = as.character(Gene.2),IDFrom = "ENTREZID",IDTo = "SYMBOL")[[1]][,2])%>%dplyr::select(c(common.TF.HGNC,Gene.1.HGNC,Gene.2.HGNC)))
msbb_array.DRGs=lapply(msbb_array.DRsort,function(x)mapIds2(IDs = x$DRGs%>%filter(q<0.1&DCGisTF=="TRUE")%>%pull(DCG)%>%droplevels%>%levels,IDFrom = "ENTREZID",IDTo = "SYMBOL")[[1]][,2])
msbb_array.DRGs_geneID=lapply(msbb_array.DRsort,function(x)x$DRGs%>%filter(q<0.1&DCGisTF=="TRUE")%>%pull(DCG)%>%droplevels%>%levels)


msbb_array_DRG_downstream_DCLs=vector(mode = "list",length = 10)
names(msbb_array_DRG_downstream_DCLs)=names(msbb_array.DRGs_geneID)
for(i in 1:10){
  msbb_array_DRG_downstream_DCLs[[i]]=msbb_array.DCL[names(msbb_array.DRGs_geneID)][[i]][which((msbb_array.DCL[names(msbb_array.DRGs_geneID)][[i]]$Gene.1%in%msbb_array.DRGs_geneID[[i]])|(msbb_array.DCL[names(msbb_array.DRGs_geneID)][[i]]$Gene.2%in%msbb_array.DRGs_geneID[[i]])),c('Gene.1','Gene.2')]
}
msbb_array_DRG_downstream_DCLs=lapply(msbb_array_DRG_downstream_DCLs,function(x)cbind.data.frame(Gene.1=names(mapIds(x = org.Hs.eg.db,keytype = "ENTREZID",column = "SYMBOL",keys = x$Gene.1,multiVals = "first")),
                                                                                                 Gene.1.HGNC=mapIds(x = org.Hs.eg.db,keytype = "ENTREZID",column = "SYMBOL",keys = x$Gene.1,multiVals = "first"),
                                                                                                 Gene.2=names(mapIds(x = org.Hs.eg.db,keytype = "ENTREZID",column = "SYMBOL",keys = x$Gene.2,multiVals = "first")),
                                                                                                 Gene.2.HGNC=mapIds(x = org.Hs.eg.db,keytype = "ENTREZID",column = "SYMBOL",keys = x$Gene.2,multiVals = "first"))%>%dplyr::select(c(Gene.1,Gene.2,Gene.1.HGNC,Gene.2.HGNC)))


msbb_array.common_DRGs=Reduce(intersect,lapply(msbb_array.DRsort,function(x)levels(unique(x$DCG2TF$TF))))
# mapIds2(IDs = as.character(x$DRGs$Upstream_TFofDCG[grep(pattern = ";",x$DRGs$Upstream_TFofDCG,invert = T)]),IDFrom = "ENTREZID",IDTo = "SYMBOL")[[1]][,2])
msbb_array_Upstream_TF_DCG_CO_links.df=lapply(msbb_array.DRsort,function(x)data.frame(cbind(Upstream_TFofDCG=grep(pattern = ";",x = x$DRGs$Upstream_TFofDCG[which(is.na(x$DRGs$Upstream_TFofDCG)==F)],invert = T,value = T),DC.links=x$DRGs$DC.links[grep(pattern = ";",x = x$DRGs$Upstream_TFofDCG[which(is.na(x$DRGs$Upstream_TFofDCG)==F)],invert = T)]),stringsAsFactors = F))

msbb_array.DCG2TF=lapply(msbb_array.DRsort,function(x)x$DCG2TF%>%mutate(DCG.HGNC=mapIds2(IDs = as.character(DCG),IDFrom = "ENTREZID",IDTo = "SYMBOL")[[1]][,2],
                                                                        TF.HGNC=mapIds2(IDs = as.character(TF),IDFrom = "ENTREZID",IDTo = "SYMBOL")[[1]][,2]))
msbb_array_TF_bridged.DCLs=lapply(msbb_array.DRsort,function(x)x$TF_bridged_DCL%>%filter(q.diffcor<=0.1))
msbb_array_TF_bridged.DCLs=lapply(msbb_array_TF_bridged.DCLs,function(x)x%>%mutate(common.TF.HGNC=mapIds2(IDs = as.character(common.TF),IDFrom = "ENTREZID",IDTo = "SYMBOL")[[1]][,2],
                                                                                   Gene.1.HGNC=mapIds2(IDs = as.character(Gene.1),IDFrom = "ENTREZID",IDTo = "SYMBOL")[[1]][,2],
                                                                                   Gene.2.HGNC=mapIds2(IDs = as.character(Gene.2),IDFrom = "ENTREZID",IDTo = "SYMBOL")[[1]][,2],
                                                                                   DCG.HGNC=unlist(lapply(strsplit(x=DCG,split = ";"),function(y)paste(mapIds2(IDs = gsub(pattern = " ",replacement = "",x = y),IDFrom = "ENTREZID",IDTo = "SYMBOL")[[1]][,2],collapse = ";"))))%>%select(c(common.TFinexprs,common.TFisDCG,DCG.HGNC,common.TF.HGNC,Gene.1.HGNC,Gene.2.HGNC)))
                                 


dbs <- listEnrichrDbs()
dbs_kegg2016=grep(pattern = "KEGG",x = dbs$libraryName,value = T)[3]
dbs_biocarta2016=grep(pattern = "BioCarta|Panther",x = dbs$libraryName,value = T)[4]
dbs_panther2016=grep(pattern = "BioCarta|Panther",x = dbs$libraryName,value = T)[5]
msbb_array_DCG.KEGG=lapply(lapply(lapply(msbb_array_DCG.filtered,function(x)mapIds2(IDs = x,IDFrom = "ENTREZID",IDTo = "SYMBOL")[[1]][,2]),function(y)enrichr(genes = y,kegg_dbs)),function(y)y$KEGG_2016%>%dplyr::filter(Adjusted.P.value<=0.1))
msbb_array_DCG.KEGG_List=lapply(Filter(f = function(x)dim(x)[1]>0,x = msbb_array_DCG.KEGG),function(y)y$Term)
msbb_array_DRG.KEGG=lapply(lapply(msbb_array_AD.DRGs,function(y)enrichr(genes = y,kegg_dbs)),function(y)y$KEGG_2016%>%filter(Adjusted.P.value<=0.1))
msbb_array_DRG.KEGG_List=lapply(Filter(f = function(x)dim(x)[1]>0,x = msbb_array_DRG.KEGG),function(y)y$Term)
msbb_DRsort_commonDRGs_downstreamDCG.KEGG=lapply(msbb_array.DRsort,function(x)enrichr(genes = mapIds2(IDs = subset(x$DCG2TF,TF%in%msbb_array.common_DRGs)%>%pull(DCG)%>%droplevels%>%levels,IDFrom = "ENTREZID",IDTo = "SYMBOL")[[1]][,2],databases = kegg_dbs)[[1]]%>%filter(Adjusted.P.value<=0.1)%>%pull(Term)%>%sort)
msbb_DCG_DRG_common_pathways.list=mapply(FUN = function(a,b)intersect(a,b),msbb_array_DCG.KEGG_List[intersect(names(msbb_array_DCG.KEGG_List),names(msbb_array_DRG.KEGG_List))],msbb_array_DRG.KEGG_List[intersect(names(msbb_array_DCG.KEGG_List),names(msbb_array_DRG.KEGG_List))])

common_regions_DCG_DRG=intersect(names(msbb_array_DCG.KEGG_List),names(msbb_array_DRG.KEGG_List))
msbb_array_DCG.KEGG_List2=msbb_array_DCG.KEGG_List[common_regions_DCG_DRG]
msbb_array_DRG.KEGG_List2=msbb_array_DRG.KEGG_List[common_regions_DCG_DRG]


msbb_array.DEG=vector(mode = "list",length = 17)
names(msbb_array.DEG)=names(msbb_array19.2.agg2)
for(t1 in 1:17){
  
  control_data=msbb_array19.2.agg2[[t1]][,grep(pattern = paste('\\b',gsub(pattern = "X",replacement = "",x = msbb_array19.covariates$BrainBank[which(msbb_array19.covariates$SampleType=="CONTROL")]),'\\b',sep="",collapse = '|'),x = colnames(msbb_array19.2.agg2[[t1]]))]
  case_data=msbb_array19.2.agg2[[t1]][,grep(pattern = paste('\\b',gsub(pattern = "X",replacement = "",x = msbb_array19.covariates$BrainBank[which(msbb_array19.covariates$SampleType=="AD")]),'\\b',sep="",collapse = '|'),x = colnames(msbb_array19.2.agg2[[t1]]))]
  group=factor(c(rep("CONTROL",length=length(colnames(control_data))),(rep("AD",length=length(colnames(case_data))))))
  design_df=model.matrix(~0+group)
  colnames(design_df)=levels(group)
  contrasts_matrix=limma::makeContrasts(CONTROL-AD,levels = design_df)
  fit=limma::lmFit(object = cbind(control_data,case_data),design = design_df)
  fit2=limma::contrasts.fit(fit,contrasts_matrix)
  fit2=limma::eBayes(fit2)
  fit2=limma::eBayes(fit2,trend = T)
  msbb_array.DEG[[t1]]=limma::topTable(fit2,coef = 1,sort.by = "logFC",adjust.method = "BH",p.value = 0.1,number = dim(control_data)[1])%>%rownames_to_column("Gene")
 
  if(dim(msbb_array.DEG[[t1]])[1]==0){
    msbb_array.DEG[[t1]]=0
    }
      else
        msbb_array.DEG[[t1]]=limma::topTable(fit2,coef = 1,sort.by = "logFC",adjust.method = "BH",p.value = 0.1,number = dim(control_data)[1])
}

msbb_array.DEG=lapply(Filter(f = function(x)length(x)>1,x = msbb_array.DEG),function(y)mapIds2(IDs = rownames(y),IDFrom = "ENTREZID",IDTo = "SYMBOL")[[1]][,2])


msbb_array_DRrank.TDD=mcmapply(FUN=function(a,b)DRrank(DCGs = a,DCLs = b,tf2target = regnet_tf2target,expGenes = rownames(msbb_array19.2.agg2$AC),rank.method = "TDD",Nperm = 500),msbb_array.DCG,msbb_array.DCL,mc.cores = 8)

msbb_array_DRrank.TDD=msbb_array_DRrank.TED=vector(mode = "list",length = 17)
names(msbb_array_DRrank.TDD)=names(msbb_array_DRrank.TED)=names(msbb_array19)

for(tdd in 1:17){
  msbb_array_DRrank.TED[[tdd]]=DRrank(DCGs = msbb_array.DCG[[tdd]],DCLs = msbb_array.DCL[[tdd]],tf2target = regnet_tf2target,expGenes = rownames(msbb_array19.2.agg2$AC),rank.method = "TED",Nperm = 1000)
  msbb_array_DRrank.TDD[[tdd]]=DRrank(DCGs = msbb_array.DCG[[tdd]],DCLs = msbb_array.DCL[[tdd]],tf2target = regnet_tf2target,expGenes = rownames(msbb_array19.2.agg2$AC),rank.method = "TDD",Nperm = 1000)
  saveRDS(msbb_array_DRrank.TED[[tdd]],paste(names(msbb_array_DRrank.TED)[tdd],"DRrank_TED.RDS",sep="_"))
  saveRDS(msbb_array_DRrank.TDD[[tdd]],paste(names(msbb_array_DRrank.TDD)[tdd],"DRrank_TDD.RDS",sep="_"))
}
setwd("../")
msbb_ted_files=list.files(path = "./diffcoexp_results/",pattern = "TED",full.names = T)
msbb_tdd_files=list.files(path = "./diffcoexp_results/",pattern = "TDD",full.names = T)
for(tdd in 1:17){
  msbb_array_DRrank.TED[[tdd]]=readRDS(msbb_ted_files[[tdd]])
  msbb_array_DRrank.TDD[[tdd]]=readRDS(msbb_tdd_files[[tdd]])
}
msbb_array_DRrank.TED=Filter(f = function(x)dim(x)[1]>0,x=lapply(msbb_array_DRrank.TED,function(x)x%>%rownames_to_column("Gene")%>%mutate(Gene.HGNC=mapIds2(IDs = as.character(Gene),IDFrom = "ENTREZID",IDTo = "SYMBOL")[[1]][,2])%>%filter(p.adj<=0.2)%>%select(c(Gene.HGNC,score,p,p.adj))))

#merge.data.frame(x = msbb_array.DCG2TF$FP[msbb_array.DCG2TF$FP$TF.HGNC%in%msbb_array_DRrank.TED$FP$Gene.HGNC,],y = msbb_array.DCG$FP[which(msbb_array.DCG$FP$Gene.HGNC%in%regnet_tf2target.HGNC[which(regnet_tf2target.HGNC$regulator_symbol%in%msbb_array_DRrank.TED$FP$Gene.HGNC),2]),],by.x = "DCG.HGNC",by.y = "Gene.HGNC")
msbb_array_DRrank.TDD=Filter(f = function(x)dim(x)[1]>0,x = lapply(msbb_array_DRrank.TDD,function(x)x%>%rownames_to_column("Gene")%>%mutate(Gene.HGNC=mapIds2(IDs = as.character(Gene),IDFrom = "ENTREZID",IDTo = "SYMBOL")[[1]][,2])%>%filter(p.adj<=0.2)%>%select(c(Gene.HGNC,score,p,p.adj))))

#Validation with transcriptional signature Ciryam et al.

ms_subprot_AD=vector(mode = "list",length = 2)
ms_subprot_AD[[1]]=read.table("/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/TF_signature_AD_metastable_proteome/ms_subproteome_AD_Downreg.tsv",sep = "\t",header = T,as.is = T)%>%pull(To)%>%as.character
ms_subprot_AD[[2]]=read.table("/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/TF_signature_AD_metastable_proteome/ms_subproteome_AD_Upreg.tsv",sep = "\t",header = T,as.is = T)%>%pull(To)%>%as.character
names(ms_subprot_AD)=c("DownReg","UpReg")

dhmc_bennet.NP=read.table("/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/GW-DNA-hydroxymethylation-Bennet/1-s2.0-S155252601633059X-mmc4_NP.txt",sep = "\t",header = T,as.is = T)
dhmc_bennet.NFT=read.table("/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/GW-DNA-hydroxymethylation-Bennet/1-s2.0-S155252601633059X-mmc5_NFT.txt",sep = "\t",header = T,as.is = T)
d5mc_bernstein.Dhml=scan("/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/5hMC_Tau_AD/ddw109_Supp/Supplemental Table 6_Bernstein et al_R1.txt",sep = "\n",what = "char")

#Brain cell type markers Zhang et al.
zhang_celltype_ADgenes=read.xls('/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/Zhang_19BrainRegions_Paper/Zhang_BrainCelltype_Markers.xlsx',skip=1,sheet=3,header=T,as.is=T)
zhang_celltype_ADgenes.list=vector(mode = 'list',length = 5)
names(zhang_celltype_ADgenes.list)=sort(unique(zhang_celltype_ADgenes$Cell.type))
zhang_celltype_ADgenes.list$Astrocytes=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_ADgenes.list)[1],zhang_celltype_ADgenes$Cell.type)]
zhang_celltype_ADgenes.list$Endothelial=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_ADgenes.list)[2],zhang_celltype_ADgenes$Cell.type)]
zhang_celltype_ADgenes.list$Microglia=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_ADgenes.list)[3],zhang_celltype_ADgenes$Cell.type)]
zhang_celltype_ADgenes.list$Neurons=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_ADgenes.list)[4],zhang_celltype_ADgenes$Cell.type)]
zhang_celltype_ADgenes.list$Oligodendrocytes=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_ADgenes.list)[5],zhang_celltype_ADgenes$Cell.type)]


msbb_array_DCG_celltype.Astrocytes=Filter(f = function(x)length(x)>0,x = lapply(lapply(msbb_array_DCG.filtered,function(x)mapIds2(IDs = x,IDFrom = "ENTREZID",IDTo = "SYMBOL")[[1]][,2]),intersect,zhang_celltype_ADgenes.list$Astrocytes))
msbb_array_DCG_celltype.Endothelial=Filter(f = function(x)length(x)>0,x = lapply(lapply(msbb_array_DCG.filtered,function(x)mapIds2(IDs = x,IDFrom = "ENTREZID",IDTo = "SYMBOL")[[1]][,2]),intersect,zhang_celltype_ADgenes.list$Endothelial))
msbb_array_DCG_celltype.Microglia=Filter(f = function(x)length(x)>0,x = lapply(lapply(msbb_array_DCG.filtered,function(x)mapIds2(IDs = x,IDFrom = "ENTREZID",IDTo = "SYMBOL")[[1]][,2]),intersect,zhang_celltype_ADgenes.list$Microglia))
msbb_array_DCG_celltype.Neurons=Filter(f = function(x)length(x)>0,x = lapply(lapply(msbb_array_DCG.filtered,function(x)mapIds2(IDs = x,IDFrom = "ENTREZID",IDTo = "SYMBOL")[[1]][,2]),intersect,zhang_celltype_ADgenes.list$Neurons))
msbb_array_DCG_celltype.Oligodendrocytes=Filter(f = function(x)length(x)>0,x = lapply(lapply(msbb_array_DCG.filtered,function(x)mapIds2(IDs = x,IDFrom = "ENTREZID",IDTo = "SYMBOL")[[1]][,2]),intersect,zhang_celltype_ADgenes.list$Oligodendrocytes))



#Plot TF bridged DCLs

for(i in 1:17){
  DRplot(DCGs = msbb_array.DCG[[i]],DCLs = msbb_array.DCL[[i]],tf2target = regnet_tf2target,expGenes = rownames(msbb_array19.2.agg2$AC),type = "TF_bridged_DCL",vsize = 10,figname = paste(names(msbb_array.DCG)[i],"TF_bridged_DCL.pdf",sep = "_"))
}

#Reproduce results in Berchtold dataset

gse11882_series_matrix=fread("/Users/sandeepamberkar/Work/Data/AD_GSE11882/GSE11882_series_matrix.txt",sep = "\t",header = T,stringsAsFactors = F,showProgress = T,data.table = F)
gpl570_annot=fread("/Users/sandeepamberkar/Work/Data/AD_GSE11882/GPL570.annot",skip = 27,sep = "\t",header = T,stringsAsFactors = F,showProgress = T,data.table = F)
gse11882_series_matrix$Symbol=gpl570_annot$`Gene symbol`
gse11882_series_matrix.agg=aggregate(x = gse11882_series_matrix[,-c(1,175)],by=list(symbol=gse11882_series_matrix$Symbol),mean)
rownames(gse11882_series_matrix.agg)=gse11882_series_matrix.agg$symbol

gse11882_soft=getGEO(filename = "/Users/sandeepamberkar/Work/Data/AD_GSE11882/GSE11882_family.soft.gz",destdir = "/Users/sandeepamberkar/Work/Data/AD_GSE11882/GSE11882_family.soft")


gse48350_soft=getGEO(filename = "/Users/sandeepamberkar/Work/Data/AD_GSE48350/GSE48350_family.soft.gz",destdir = "/Users/sandeepamberkar/Work/Data/AD_GSE48350/GSE48350_family.soft")


#Read TF-Target interactions, preprocess data
