library(diffcoexp)
library(org.Hs.eg.db)
library(dplyr)
library(magrittr)
library(enrichR)
library(data.table)
library(foreach)
library(doParallel)
library(DCGL)
library(clusterProfiler)

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

tanzi_ranked_genes=read.table("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Gene_Lists/All_genes_ranked.txt",header = T,sep = "\t",stringsAsFactors = F)
tanzi_ranked_genes.great=read.table("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Gene_Lists/Ranked_genes_corrected_allele_top_4000_SNPs_genes_from_great.txt",header = T,sep = "\t",stringsAsFactors = F)
tanzi_ranked_genes.great$SNP_Rank=gsub(pattern = "SNP_rank_",replacement = "",x = tanzi_ranked_genes.great$SNP_Rank)
tanzi_ranked_genes.AD=unique(tanzi_ranked_genes.great%>%mutate(SNP_Rank=as.integer(gsub(pattern = "SNP_rank_",replacement = "",x = SNP_Rank)))%>%filter(Corrected_association=="Case"&SNP_Rank<=1000)%>%pull(Gene))
tanzi_ranked_genes.Control=unique(tanzi_ranked_genes.great%>%mutate(SNP_Rank=as.integer(gsub(pattern = "SNP_rank_",replacement = "",x = SNP_Rank)))%>%filter(Corrected_association=="Control"&SNP_Rank<=1000)%>%pull(Gene))
tanzi_ranked_genes.union=union(tanzi_ranked_genes.AD,tanzi_ranked_genes.Control)

setwd("/Users/sandeepamberkar/Work/Data/MSMM-MSBB-HBTRC-Data/MSBB_Array19/GSE84422/")
msbb_gse84422_GPL96_97_samplesToAnalyse.exprs=readRDS("msbb_gse84422_GPL96_97_samplesToAnalyse_exprs.RDS")
names(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs)=gsub(pattern = " ",replacement = "_",x = names(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs))
msbb_gse84422_GPL570_samplesToAnalyse.exprs=readRDS("msbb_gse84422_GPL570_samplesToAnalyse_exprs.RDS")
names(msbb_gse84422_GPL570_samplesToAnalyse.exprs)=gsub(pattern = " ",replacement = "_",x = names(msbb_gse84422_GPL570_samplesToAnalyse.exprs))
regnet_tf2target.HGNC=fread("/shared/hidelab2/user/md4zsa/Work/Data/TF_Databases/RegNetwork_human_regulators2.txt",header = T,sep = "\t",showProgress = T,data.table = F)%>%dplyr::filter(evidence=="Experimental")%>%dplyr::select(c(regulator_symbol,target_symbol))
regnet_tf2target.HGNC=fread("/Users/sandeepamberkar/Work/Data/TF_Databases/RegNetwork_human_regulators2.txt",header = T,sep = "\t",showProgress = T,data.table = F)%>%filter(evidence=="Experimental")%>%dplyr::select(c(regulator_symbol,target_symbol))

msbb_gse84422_diffcoexp_results_files=list.files(path = ".",pattern = "diffcoexp")%>%grep(pattern = ".RDS",value = T)%>%sort
msbb_gse84422_diffcoexp_results=vector(mode = "list",length = length(msbb_gse84422_diffcoexp_results_files))
names(msbb_gse84422_diffcoexp_results)=gsub(pattern = " ",unlist(lapply(lapply(msbb_gse84422_diffcoexp_results_files,function(y)strsplit(x = y,split = "_")[[1]]),`[[`,1)),replacement = "_")
for(f in 1:19){
  msbb_gse84422_diffcoexp_results[[f]]=readRDS(msbb_gse84422_diffcoexp_results_files[f])  
}
msbb_gse84422.DCGs=lapply(msbb_gse84422_diffcoexp_results,function(x)x$DCGs)
msbb_gse84422.DCGs_list=lapply(msbb_gse84422.DCGs,function(x)x%>%dplyr::filter(q<=0.05)%>%pull(Gene))
msbb_gse84422.DCGs_list_Entrez=lapply(msbb_gse84422.DCGs_list,function(y)unname(mapIds(x = org.Hs.eg.db,keys = y,keytype = "SYMBOL",column = "ENTREZID")))
msbb_gse84422.DCLs=lapply(msbb_gse84422_diffcoexp_results,function(x)x$DCLs%>%filter(q.diffcor<=0.05)%>%rownames_to_column("GenePair"))
msbb_gse84422_total_DCGs=unlist(lapply(msbb_gse84422.DCGs_list,length))
msbb_gse84422_total_DCLs=unlist(lapply(msbb_gse84422.DCLs,function(x)dim(x)[1]))
dat1=data.frame(DCLs=msbb_gse84422_total_DCLs,DCGs=msbb_gse84422_total_DCGs,Region=names(msbb_gse84422_total_DCLs))
dat1m=melt(dat1[,c('Region','DCLs','DCGs')],id.vars=1)
ggplot(dat1m,aes(x = Region,y = value)) + geom_bar(aes(fill = variable),stat = "identity",position = "dodge")+
  labs(title="MSBB array - DCGs and DCLs")+
  theme(title = element_text(face = "bold",size = 15,colour = "black"),axis.title.x = element_text(face = "bold",size = 15,colour = "black"),axis.text.x = element_text(size = 10,colour = "black",angle = 90),legend.text = element_text(size = 15))

#Read DRsort results
msbb_gse84422.DRsort=vector(mode = "list",length = 19)
names(msbb_gse84422.DRsort)=names(msbb_gse84422_diffcoexp_results)
msbb_gse84422.DRsort$Amygdala=DRsort(DCGs = msbb_gse84422.DCGs$Amygdala,DCLs = msbb_gse84422.DCLs$Amygdala,tf2target = regnet_tf2target.HGNC,expGenes = rownames(msbb_gse84422_GPL570_samplesToAnalyse.exprs$Amygdala))

msbb_gse84422.DRsort$Anterior_Cingulate=DRsort(DCGs = msbb_gse84422.DCGs$Anterior_Cingulate,DCLs = msbb_gse84422.DCLs$Anterior_Cingulate,tf2target = regnet_tf2target.HGNC,expGenes = rownames(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs$Anterior_Cingulate))

msbb_gse84422.DRsort$Caudate_Nucleus=DRsort(DCGs = msbb_gse84422.DCGs$Caudate_Nucleus,DCLs = msbb_gse84422.DCLs$Caudate_Nucleus,tf2target = regnet_tf2target.HGNC,expGenes = rownames(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs$Caudate_Nucleus))

msbb_gse84422.DRsort$Dorsolateral_Prefrontal_Cortex=DRsort(DCGs = msbb_gse84422.DCGs$Dorsolateral_Prefrontal_Cortex,DCLs = msbb_gse84422.DCLs$Dorsolateral_Prefrontal_Cortex,tf2target = regnet_tf2target.HGNC,expGenes = rownames(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs$Dorsolateral_Prefrontal_Cortex))

msbb_gse84422.DRsort$Frontal_Pole=DRsort(DCGs = msbb_gse84422.DCGs$Frontal_Pole,DCLs = msbb_gse84422.DCLs$Frontal_Pole,tf2target = regnet_tf2target.HGNC,expGenes = rownames(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs$Frontal_Pole))

msbb_gse84422.DRsort$Hippocampus=DRsort(DCGs = msbb_gse84422.DCGs$Hippocampus,DCLs = msbb_gse84422.DCLs$Hippocampus,tf2target = regnet_tf2target.HGNC,expGenes = rownames(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs$Hippocampus))

msbb_gse84422.DRsort$Inferior_Frontal_Gyrus=DRsort(DCGs = msbb_gse84422.DCGs$Inferior_Frontal_Gyrus,DCLs = msbb_gse84422.DCLs$Inferior_Frontal_Gyrus,tf2target = regnet_tf2target.HGNC,expGenes = rownames(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs$Inferior_Frontal_Gyrus))

msbb_gse84422.DRsort$Inferior_Temporal_Gyrus=DRsort(DCGs = msbb_gse84422.DCGs$Inferior_Temporal_Gyrus,DCLs = msbb_gse84422.DCLs$Inferior_Temporal_Gyrus,tf2target = regnet_tf2target.HGNC,expGenes = rownames(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs$Inferior_Temporal_Gyrus))

msbb_gse84422.DRsort$Middle_Temporal_Gyrus=DRsort(DCGs = msbb_gse84422.DCGs$Middle_Temporal_Gyrus,DCLs = msbb_gse84422.DCLs$Middle_Temporal_Gyrus,tf2target = regnet_tf2target.HGNC,expGenes = rownames(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs$Middle_Temporal_Gyrus))

msbb_gse84422.DRsort$Nucleus_Accumbens=DRsort(DCGs = msbb_gse84422.DCGs$Nucleus_Accumbens,DCLs = msbb_gse84422.DCLs$Nucleus_Accumbens,tf2target = regnet_tf2target.HGNC,expGenes = rownames(msbb_gse84422_GPL570_samplesToAnalyse.exprs$Nucleus_Accumbens))

msbb_gse84422.DRsort$Occipital_Visual_Cortex=DRsort(DCGs = msbb_gse84422.DCGs$Occipital_Visual_Cortex,DCLs = msbb_gse84422.DCLs$Occipital_Visual_Cortex,tf2target = regnet_tf2target.HGNC,expGenes = rownames(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs$Occipital_Visual_Cortex))

msbb_gse84422.DRsort$Parahippocampal_Gyrus=DRsort(DCGs = msbb_gse84422.DCGs$Parahippocampal_Gyrus,DCLs = msbb_gse84422.DCLs$Parahippocampal_Gyrus,tf2target = regnet_tf2target.HGNC,expGenes = rownames(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs$Parahippocampal_Gyrus))

msbb_gse84422.DRsort$Posterior_Cingulate_Cortex=DRsort(DCGs = msbb_gse84422.DCGs$Posterior_Cingulate_Cortex,DCLs = msbb_gse84422.DCLs$Posterior_Cingulate_Cortex,tf2target = regnet_tf2target.HGNC,expGenes = rownames(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs$Posterior_Cingulate_Cortex))

msbb_gse84422.DRsort$Precentral_Gyrus=DRsort(DCGs = msbb_gse84422.DCGs$Precentral_Gyrus,DCLs = msbb_gse84422.DCLs$Precentral_Gyrus,tf2target = regnet_tf2target.HGNC,expGenes = rownames(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs$Precentral_Gyrus))

msbb_gse84422.DRsort$Prefrontal_Cortex=DRsort(DCGs = msbb_gse84422.DCGs$Prefrontal_Cortex,DCLs = msbb_gse84422.DCLs$Prefrontal_Cortex,tf2target = regnet_tf2target.HGNC,expGenes = rownames(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs$Prefrontal_Cortex))

msbb_gse84422.DRsort$Putamen=DRsort(DCGs = msbb_gse84422.DCGs$Putamen,DCLs = msbb_gse84422.DCLs$Putamen,tf2target = regnet_tf2target.HGNC,expGenes = rownames(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs$Putamen))

msbb_gse84422.DRsort$Superior_Parietal_Lobule=DRsort(DCGs = msbb_gse84422.DCGs$Superior_Parietal_Lobule,DCLs = msbb_gse84422.DCLs$Superior_Parietal_Lobule,tf2target = regnet_tf2target.HGNC,expGenes = rownames(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs$Superior_Parietal_Lobule))

msbb_gse84422.DRsort$Superior_Temporal_Gyrus=DRsort(DCGs = msbb_gse84422.DCGs$Superior_Temporal_Gyrus,DCLs = msbb_gse84422.DCLs$Superior_Temporal_Gyrus,tf2target = regnet_tf2target.HGNC,expGenes = rownames(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs$Superior_Temporal_Gyrus))

msbb_gse84422.DRsort$Temporal_Pole=DRsort(DCGs = msbb_gse84422.DCGs$Temporal_Pole,DCLs = msbb_gse84422.DCLs$Temporal_Pole,tf2target = regnet_tf2target.HGNC,expGenes = rownames(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs$Temporal_Pole))


# for(i in c(1,10)){msbb_gse84422.DRsort$Amygdala=DRsort(DCGs = msbb_gse84422.DCGs$Amygdala,DCLs = msbb_gse84422.DCLs$Amygdala,tf2target = regnet_tf2target.HGNC,expGenes = rownames(msbb_gse84422_GPL570_samplesToAnalyse.exprs$Amygdala))
#   msbb_gse84422.DRsort[[i]]=DRsort(DCGs = msbb_gse84422.DCGs[[i]],DCLs = msbb_gse84422.DCLs[[i]],tf2target = regnet_tf2target.HGNC,expGenes = rownames(msbb_gse84422_GPL570_samplesToAnalyse.exprs$Amygdala))
# }
# for(i in c(2:9,11:19)){
#   msbb_gse84422.DRsort[[i]]=DRsort(DCGs = msbb_gse84422.DCGs[[i]],DCLs = msbb_gse84422.DCLs[[i]],tf2target = regnet_tf2target.HGNC,expGenes = rownames(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs$Frontal_Pole))
# }

# msbb_gse84422.DRsort=readRDS("msbb_gse84422_DRsort.RDS")
msbb_gse84422.DRGs=lapply(msbb_gse84422.DRsort,function(x)x$DRGs%>%filter(DCGisTF=="TRUE"&q<=0.05))
msbb_gse84422.DRGs_list=lapply(msbb_gse84422.DRsort,function(x)x$DRGs%>%filter(DCGisTF=="TRUE"&q<=0.05)%>%pull(DCG)%>%droplevels%>%levels)
msbb_gse84422.DRGs_list_Entrez=lapply(msbb_gse84422.DRGs_list,function(y)unname(mapIds(x = org.Hs.eg.db,keys = y,keytype = "SYMBOL",column = "ENTREZID")))
msbb_gse84422.DRLs=lapply(msbb_gse84422.DRsort,function(x)x$DRLs%>%filter(q.diffcor<=0.05))
msbb_gse84422.TF_bridged_DCL=Filter(f = function(x)dim(x)[1]>0,x = lapply(msbb_gse84422.DRsort,function(x)x$TF_bridged_DCL%>%filter(q.diffcor<=0.05)))
msbb_gse84422.DCG2TF=lapply(msbb_gse84422.DRsort,function(x)x$DCG2TF)
msbb_gse84422_bridge_DCL_by_TF=vector(mode = "list",length = length(names(msbb_gse84422.TF_bridged_DCL)))
names(msbb_gse84422_bridge_DCL_by_TF)=names(msbb_gse84422.TF_bridged_DCL)
for(i in 1:length(names(msbb_gse84422_bridge_DCL_by_TF))){
  msbb_gse84422_bridge_DCL_by_TF[[i]]=vector(mode = "list",length = length(lapply(unique(msbb_gse84422.TF_bridged_DCL[[i]]$common.TF),function(x)msbb_gse84422.TF_bridged_DCL[[i]]%>%filter(common.TF==x)%>%select(c(Gene.1,Gene.2)))))
  names(msbb_gse84422_bridge_DCL_by_TF[[i]])=unique(msbb_gse84422.TF_bridged_DCL[[i]]$common.TF)
  msbb_gse84422_bridge_DCL_by_TF[[i]][1:length(names(msbb_gse84422_bridge_DCL_by_TF[[i]]))]=lapply(unique(msbb_gse84422.TF_bridged_DCL[[i]]$common.TF),function(x)msbb_gse84422.TF_bridged_DCL[[i]]%>%filter(common.TF==x)%>%select(c(Gene.1,Gene.2)))
}


msbb_array_DRrank.TDD=msbb_array_DRrank.TED=vector(mode = "list",length = 19)
names(msbb_array_DRrank.TDD)=names(msbb_array_DRrank.TED)=names(msbb_gse84422_diffcoexp_results)

msbb_array_DRrank.TDD[c(1,10)]=mcmapply(FUN=function(a,b)DRrank(DCGs = a,DCLs = b,tf2target = regnet_tf2target.HGNC,expGenes = rownames(msbb_gse84422_GPL570_samplesToAnalyse.exprs$Amygdala),rank.method = "TDD",Nperm = 500),msbb_gse84422.DCGs[c(1,10)],msbb_gse84422.DCLs[c(1,10)],mc.cores = mc)

#Functional enrichments
msigdb_c7=read.gmt(gmtfile = "/Users/sandeepamberkar/Work/Software/MSigDB_Collections/c7.all.v6.1.symbols.gmt")
msigdb_c2=read.gmt(gmtfile = "/Users/sandeepamberkar/Work/Software/MSigDB_Collections/c2.all.v6.1.symbols.gmt")
msigdb_h=read.gmt(gmtfile = "/Users/sandeepamberkar/Work/Software/MSigDB_Collections/h.all.v6.1.symbols.gmt")
msigdb_c5.BP=read.gmt(gmtfile = "/Users/sandeepamberkar/Work/Software/MSigDB_Collections/c5.bp.v6.1.symbols.gmt")
msigdb_c5.CC=read.gmt(gmtfile = "/Users/sandeepamberkar/Work/Software/MSigDB_Collections/c5.cc.v6.1.symbols.gmt")
msigdb_c5.MF=read.gmt(gmtfile = "/Users/sandeepamberkar/Work/Software/MSigDB_Collections/c5.mf.v6.1.symbols.gmt")

msbb_gse84422_DCGs.c7=lapply(msbb_gse84422.DCGs_list,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.01,pAdjustMethod = "BH",TERM2GENE = msigdb_c7),stringsAsFactors = F))
msbb_gse84422_DCGs.c2=lapply(msbb_gse84422.DCGs_list,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.01,pAdjustMethod = "BH",TERM2GENE = msigdb_c2),stringsAsFactors = F))
msbb_gse84422_DCGs.h=lapply(msbb_gse84422.DCGs_list,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.01,pAdjustMethod = "BH",TERM2GENE = msigdb_h),stringsAsFactors = F))
msbb_gse84422_DCGs.KEGG=vector(mode = "list",length = 19)
names(msbb_gse84422_DCGs.KEGG)=names(msbb_gse84422_diffcoexp_results)
msbb_gse84422_DCGs.KEGG=lapply(msbb_gse84422.DCGs_list_Entrez,function(x)data.frame(enrichKEGG(gene = x,organism = "hsa",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1),stringsAsFactors = F))
msbb_gse84422_DCGs.KEGG=Filter(f = function(x)dim(x)[1]>0,x = msbb_gse84422_DCGs.KEGG)
msbb_gse84422_DRGs.KEGG=lapply(msbb_gse84422.DRGs_list_Entrez,function(x)data.frame(enrichKEGG(gene = x,organism = "hsa",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1),stringsAsFactors = F))
msbb_gse84422_DRGs.KEGG=Filter(f = function(x)dim(x)[1]>0,x = msbb_gse84422_DRGs.KEGG)
msbb_gse84422_DCGs_KEGG.list=lapply(msbb_gse84422_DCGs.KEGG,function(x)x$Description)
msbb_gse84422_DRGs_KEGG.list=lapply(msbb_gse84422_DRGs.KEGG,function(x)x$Description)

msbb_gse84422_DCG_DRG_common_pathways.list=mapply(FUN = function(a,b)intersect(a,b),msbb_gse84422_DCGs_KEGG.list[intersect(names(msbb_gse84422_DCGs_KEGG.list),names(msbb_gse84422_DRGs_KEGG.list))],msbb_gse84422_DRGs_KEGG.list[intersect(names(msbb_gse84422_DCGs_KEGG.list),names(msbb_gse84422_DRGs_KEGG.list))])


#Genetic basis for AD
msbb_gse84422_tanzi_SAGs.overlap=lapply(msbb_gse84422.DCGs_list,intersect,tanzi_ranked_genes.union)
dat2=data.frame(DCG_SAG_overlap=unlist(lapply(msbb_gse84422_tanzi_SAGs.overlap,length)),Total_DCGs=msbb_gse84422_total_DCGs,Region=names(msbb_gse84422_total_DCGs))
dat2m=melt(dat1[,c('Region','DCG_SAG_overlap','Total_DCGs')],id.vars=1)
ggplot(dat1m,aes(x = Region,y = value)) + geom_bar(aes(fill = variable),stat = "identity",position = "dodge")+
  labs(title="MSBB - Overlap with Tanzi SAGs")+scale_y_log10()+
  theme(title = element_text(face = "bold",size = 12,colour = "black"),axis.title.x = element_text(face = "bold",size = 15,colour = "black"),axis.text.x = element_text(size = 10,colour = "black",angle = 90),legend.text = element_text(size = 15))


msbb_array_DRrank.TDD=msbb_array_DRrank.TED=vector(mode = "list",length = 19)
names(msbb_array_DRrank.TDD)=names(msbb_array_DRrank.TED)=names(msbb_gse84422_diffcoexp_results)
msbb_array_DRrank_TED.files=list.files(pattern = "TED",full.names = T)
msbb_array_DRrank_TDD.files=list.files(pattern = "TDD",full.names = T)
for(i in 1:19){
  msbb_array_DRrank.TED[[i]]=readRDS(msbb_array_DRrank_TED.files[[i]])
  msbb_array_DRrank.TDD[[i]]=readRDS(msbb_array_DRrank_TDD.files[[i]])
}
msbb_array_DRrank.TED=lapply(msbb_array_DRrank.TED,function(x)x%>%filter(p<=0.05)%>%mutate(Rank=with_order(order_by = score,fun = row_number,x = -score))%>%arrange(desc(-Rank))%>%dplyr::select(c(Gene,Rank,score,p,p.adj)))
msbb_array_DRrank.TED=Filter(f = function(x)dim(x)[1]>0,x = msbb_array_DRrank.TED)
msbb_array_DRrank.TDD=lapply(msbb_array_DRrank.TDD,function(x)x%>%rownames_to_column("Gene")%>%filter(p<=0.05)%>%mutate(Rank=with_order(order_by = score,fun = row_number,x = -score))%>%arrange(desc(-Rank))%>%dplyr::select(c(Gene,Rank,score,p,p.adj)))
msbb_array_DRrank.TDD=Filter(f = function(x)dim(x)[1]>0,x = msbb_array_DRrank.TDD)
#NOTCH subunits appear as causal
msbb_gse84422_NOTCH_downstream.DCGs=lapply(msbb_gse84422.DCG2TF[names(Filter(f = function(y)length(y)>0,x = lapply(msbb_array_DRrank.TED,function(x)grep(pattern = "NOTCH",x = x$Gene))))],function(x)x[grep(pattern = "NOTCH",x = x$TF),])
msbb_gse84422_PIAS_downstream.DCGs=lapply(msbb_gse84422.DCG2TF[names(Filter(f = function(y)length(y)>0,x = lapply(msbb_array_DRrank.TED,function(x)grep(pattern = "PIAS",x = x$Gene))))],function(x)x[grep(pattern = "PIAS",x = x$TF),])
##########################################################################################################################################################
# Validation with genetic, epigenomic and cell-type markers
##########################################################################################################################################################
msbb_exp_DCGs=lapply(lapply(msbb_gse84422.DCGs_list,length),function(x)round(length(tanzi_ranked_genes.union)/21982*x,digits = 1))
msbb_gse84422_tanzi_SAGs.overlap=lapply(msbb_gse84422.DCGs_list,function(x)length(intersect(x,tanzi_ranked_genes.union)))
msbb_gse84422_tanzi_SAGs.fisher_results=vector(mode = "list",length = 19)
for(m in 1:length(msbb_gse84422.DCGs_list)){
  msbb_gse84422_tanzi_SAGs.fisher_results[[m]]=fisher.test(matrix(c(length(msbb_gse84422_tanzi_SAGs.overlap[[m]]),
                                                              length(tanzi_ranked_genes.union)-length(msbb_gse84422_tanzi_SAGs.overlap[[m]]),
                                                              lapply(msbb_gse84422.DCGs_list,length)[[m]]-length(msbb_gse84422_tanzi_SAGs.overlap[[m]]),
                                                              21982-(length(tanzi_ranked_genes.union)-length(msbb_gse84422_tanzi_SAGs.overlap[[m]]))-
                                                              lapply(msbb_gse84422.DCGs_list,length)[[m]]-length(msbb_gse84422_tanzi_SAGs.overlap[[m]])),nrow = 2))
}

res.SAG_DCG=foreach(i=1:length(msbb_gse84422_tanzi_SAGs.overlap))%do%{
  
  
  data.frame(SAGs_in_DCGs=unlist(lapply(msbb_gse84422_tanzi_SAGs.overlap[i],paste,collapse=",")),
             DCGs=unlist(lapply(msbb_gse84422.DCGs_list,length)[[i]]),
             Expected_SAGs_in_DCGs=unname(unlist(msbb_exp_DCGs[i])),
             Fisher_statistic=unname(unlist(lapply(msbb_gse84422_tanzi_SAGs.fisher_results[i],function(s)s$estimate))),
             pval=unname(unlist(lapply(msbb_gse84422_tanzi_SAGs.fisher_results[i],function(p)p$p.value))),
             stringsAsFactors = F)
}
names(res.SAG_DCG)=names(msbb_gse84422_tanzi_SAGs.overlap)
res.SAG_DCG_df=data.frame(rbindlist(res.SAG_DCG),stringsAsFactors = F)
rownames(res.SAG_DCG_df)=names(msbb_gse84422_tanzi_SAGs.overlap)
res.SAG_DCG_df$adj.p=p.adjust(p = res.SAG_DCG_df$pval,method = "fdr")
res.SAG_DCG_df$SAGs_in_DCGs=as.numeric(res.SAG_DCG_df$SAGs_in_DCGs)
res.SAG_DCG_df=res.SAG_DCG_df%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(SAGs_in_DCGs>Expected_SAGs_in_DCGs,true = "Over",false = "Under"))
fwrite(res.SAG_DCG_df,"MSBB_SAG_DCG_EnrichmentSummary.txt",sep = "\t",col.names = T,row.names = F,quote = F)

dhmc_bennet.NP=read.table("/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/GW-DNA-hydroxymethylation-Bennet/1-s2.0-S155252601633059X-mmc4_NP.txt",sep = "\t",header = T,as.is = T)
dhmc_bennet.NFT=read.table("/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/GW-DNA-hydroxymethylation-Bennet/1-s2.0-S155252601633059X-mmc5_NFT.txt",sep = "\t",header = T,as.is = T)
d5mc_bernstein.Dhml=scan("/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/5hMC_Tau_AD/ddw109_Supp/Supplemental Table 6_Bernstein et al_R1.txt",sep = "\n",what = "char")

msbb_gse84422_DCG_d5mc.overlap=lapply(msbb_gse84422.DCGs_list,function(x)length(intersect(x,d5mc_bernstein.Dhml)))
dhmc_bennet_NFT.genes=dhmc_bennet.NFT%>%filter(q.value_<=0.1)%>%pull(Nearest.gene)
dhmc_bennet_NP.genes=dhmc_bennet.NP%>%filter(q.vlaue_<=0.1)%>%pull(Neartest.gene)
msbb_gse84422_DCG_dhmc_NFT.overlap=lapply(msbb_gse84422.DCGs_list,function(x)length(intersect(x,dhmc_bennet_NFT.genes)))
msbb_gse84422_DCG_dhmc_NP.overlap=lapply(msbb_gse84422.DCGs_list,function(x)length(intersect(x,dhmc_bennet_NP.genes)))
msbb_exp_dhmc_NP_DCGs=lapply(lapply(msbb_gse84422.DCGs_list,length),function(x)round(length(dhmc_bennet_NP.genes)/21982*x,digits = 1))

msbb_gse84422_dhmc_NP.fisher_results=msbb_gse84422_dhmc_NFT.fisher_results=msbb_gse84422_d5mc.fisher_results=vector(mode = "list",length = length(msbb_gse84422.DCGs_list))
names(msbb_gse84422_dhmc_NP.fisher_results)=names(msbb_gse84422_dhmc_NFT.fisher_results)=names(msbb_gse84422_d5mc.fisher_results)=names(msbb_gse84422.DCGs_list)
for(m in 1:length(msbb_gse84422.DCGs_list)){
  msbb_gse84422_dhmc_NP.fisher_results[[m]]=fisher.test(matrix(c(length(msbb_gse84422_tanzi_SAGs.overlap[[m]]),
                                                                    length(dhmc_bennet_NP.genes)-length(msbb_gse84422_DCG_dhmc_NP.overlap[[m]]),
                                                                    lapply(msbb_gse84422.DCGs_list,length)[[m]]-length(msbb_gse84422_DCG_dhmc_NP.overlap[[m]]),
                                                                    21982-(length(dhmc_bennet_NP.genes)-length(msbb_gse84422_DCG_dhmc_NP.overlap[[m]]))-
                                                                    lapply(msbb_gse84422.DCGs_list,length)[[m]]-length(msbb_gse84422_DCG_dhmc_NP.overlap[[m]])),nrow = 2))
}

res.dhmc_NP_DCG=foreach(i=1:length(msbb_gse84422_tanzi_SAGs.overlap))%do%{
  
  
  data.frame(NPgenes_in_DCGs=unlist(lapply(msbb_gse84422_DCG_dhmc_NP.overlap[i],paste,collapse=",")),
             DCGs=unlist(lapply(msbb_gse84422.DCGs_list,length)[[i]]),
             Expected_dhmc_NPgenes_in_DCGs=unname(unlist(msbb_exp_DCGs[i])),
             Fisher_statistic=unname(unlist(lapply(msbb_gse84422_dhmc_NP.fisher_results[i],function(s)s$estimate))),
             pval=unname(unlist(lapply(msbb_gse84422_dhmc_NP.fisher_results[i],function(p)p$p.value))),
             stringsAsFactors = F)
}


names(res.dhmc_NP_DCG)=names(msbb_gse84422_DCG_dhmc_NP.overlap)
res.dhmc_NP_DCG_df=data.frame(rbindlist(res.dhmc_NP_DCG),stringsAsFactors = F)
rownames(res.dhmc_NP_DCG_df)=names(msbb_gse84422_DCG_dhmc_NP.overlap)
res.dhmc_NP_DCG_df$adj.p=p.adjust(p = res.dhmc_NP_DCG_df$pval,method = "fdr")
res.dhmc_NP_DCG_df$NPgenes_in_DCGs=as.numeric(res.dhmc_NP_DCG_df$NPgenes_in_DCGs)
res.dhmc_NP_DCG_df=res.dhmc_NP_DCG_df%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(NPgenes_in_DCGs>Expected_dhmc_NPgenes_in_DCGs,true = "Over",false = "Under"))
fwrite(res.dhmc_NP_DCG_df,"MSBB_DhMR_NP_DCG_EnrichmentSummary.txt",sep = "\t",col.names = T,row.names = F,quote = F)

for(m in 1:length(msbb_gse84422.DCGs_list)){
  msbb_gse84422_dhmc_NFT.fisher_results[[m]]=fisher.test(matrix(c(length(msbb_gse84422_tanzi_SAGs.overlap[[m]]),
                                                                  length(dhmc_bennet_NFT.genes)-length(msbb_gse84422_DCG_dhmc_NFT.overlap[[m]]),
                                                                  lapply(msbb_gse84422.DCGs_list,length)[[m]]-length(msbb_gse84422_DCG_dhmc_NFT.overlap[[m]]),
                                                                  21982-(length(dhmc_bennet_NFT.genes)-length(msbb_gse84422_DCG_dhmc_NFT.overlap[[m]]))-
                                                                    lapply(msbb_gse84422.DCGs_list,length)[[m]]-length(msbb_gse84422_DCG_dhmc_NFT.overlap[[m]])),nrow = 2))
}

res.dhmc_NFT_DCG=foreach(i=1:length(msbb_gse84422_tanzi_SAGs.overlap))%do%{
  
  
  data.frame(NFTgenes_in_DCGs=unlist(lapply(msbb_gse84422_DCG_dhmc_NFT.overlap[i],paste,collapse=",")),
             DCGs=unlist(lapply(msbb_gse84422.DCGs_list,length)[[i]]),
             Expected_dhmc_NFTgenes_in_DCGs=unname(unlist(msbb_exp_DCGs[i])),
             Fisher_statistic=unname(unlist(lapply(msbb_gse84422_dhmc_NFT.fisher_results[i],function(s)s$estimate))),
             pval=unname(unlist(lapply(msbb_gse84422_dhmc_NFT.fisher_results[i],function(p)p$p.value))),
             stringsAsFactors = F)
}


names(res.dhmc_NFT_DCG)=names(msbb_gse84422_DCG_dhmc_NFT.overlap)
res.dhmc_NFT_DCG_df=data.frame(rbindlist(res.dhmc_NFT_DCG),stringsAsFactors = F)
rownames(res.dhmc_NFT_DCG_df)=names(msbb_gse84422_DCG_dhmc_NFT.overlap)
res.dhmc_NFT_DCG_df$adj.p=p.adjust(p = res.dhmc_NFT_DCG_df$pval,method = "fdr")
res.dhmc_NFT_DCG_df$NFTgenes_in_DCGs=as.numeric(res.dhmc_NFT_DCG_df$NFTgenes_in_DCGs)
res.dhmc_NFT_DCG_df=res.dhmc_NFT_DCG_df%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(NFTgenes_in_DCGs>Expected_dhmc_NFTgenes_in_DCGs,true = "Over",false = "Under"))
fwrite(res.dhmc_NFT_DCG_df,"MSBB_DhMR_NFT_DCG_EnrichmentSummary.txt",sep = "\t",col.names = T,row.names = F,quote = F)
#Brain cell type markers Zhang et al.
zhang_celltype_ADgenes=read.xls('/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/Zhang_19BrainRegions_Paper/Zhang_BrainCelltype_Markers.xlsx',skip=1,sheet=3,header=T,as.is=T)
zhang_celltype_ADgenes.list=vector(mode = 'list',length = 5)
names(zhang_celltype_ADgenes.list)=sort(unique(zhang_celltype_ADgenes$Cell.type))
zhang_celltype_ADgenes.list$Astrocytes=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_ADgenes.list)[1],zhang_celltype_ADgenes$Cell.type)]
zhang_celltype_ADgenes.list$Endothelial=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_ADgenes.list)[2],zhang_celltype_ADgenes$Cell.type)]
zhang_celltype_ADgenes.list$Microglia=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_ADgenes.list)[3],zhang_celltype_ADgenes$Cell.type)]
zhang_celltype_ADgenes.list$Neurons=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_ADgenes.list)[4],zhang_celltype_ADgenes$Cell.type)]
zhang_celltype_ADgenes.list$Oligodendrocytes=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_ADgenes.list)[5],zhang_celltype_ADgenes$Cell.type)]

msbb_gse84422_DCG_Astrocytes.overlap=lapply(msbb_gse84422.DCGs_list,intersect,zhang_celltype_ADgenes.list$Astrocytes)
msbb_gse84422_DCG_Endothelial.overlap=lapply(msbb_gse84422.DCGs_list,intersect,zhang_celltype_ADgenes.list$Endothelial)
msbb_gse84422_DCG_Microglia.overlap=lapply(msbb_gse84422.DCGs_list,intersect,zhang_celltype_ADgenes.list$Microglia)
msbb_gse84422_DCG_Neurons.overlap=lapply(msbb_gse84422.DCGs_list,intersect,zhang_celltype_ADgenes.list$Neurons)
msbb_gse84422_DCG_Oligodendrocytes.overlap=lapply(msbb_gse84422.DCGs_list,intersect,zhang_celltype_ADgenes.list$Oligodendrocytes)

msbb_gse84422.exp_Astrocytes=lapply(lapply(msbb_gse84422.DCGs_list,length),function(x)round(length(zhang_celltype_ADgenes.list$Astrocytes)/21982*x,digits = 1))
msbb_gse84422.exp_Endothelial=lapply(lapply(msbb_gse84422.DCGs_list,length),function(x)round(length(zhang_celltype_ADgenes.list$Endothelial)/21982*x,digits = 1))
msbb_gse84422.exp_Microglia=lapply(lapply(msbb_gse84422.DCGs_list,length),function(x)round(length(zhang_celltype_ADgenes.list$Microglia)/21982*x,digits = 1))
msbb_gse84422.exp_Neurons=lapply(lapply(msbb_gse84422.DCGs_list,length),function(x)round(length(zhang_celltype_ADgenes.list$Neurons)/21982*x,digits = 1))
msbb_gse84422.exp_Oligodendrocytes=lapply(lapply(msbb_gse84422.DCGs_list,length),function(x)round(length(zhang_celltype_ADgenes.list$Oligodendrocytes)/21982*x,digits = 1))

msbb_gse84422_Astrocytes.fisher_results=msbb_gse84422_Endothelial.fisher_results=msbb_gse84422_Microglia.fisher_results=msbb_gse84422_Neurons.fisher_results=msbb_gse84422_Oligodendrocytes.fisher_results=vector(mode = "list",length = length(msbb_gse84422.DCGs_list))
names(msbb_gse84422_Astrocytes.fisher_results)=names(msbb_gse84422_Endothelial.fisher_results)=names(msbb_gse84422_Microglia.fisher_results)=names(msbb_gse84422_Neurons.fisher_results)=names(msbb_gse84422_Oligodendrocytes.fisher_results)=names(msbb_gse84422_total_DCGs)

for(m in 1:length(msbb_gse84422_Astrocytes.fisher_results)){
  msbb_gse84422_Astrocytes.fisher_results[[m]]=fisher.test(matrix(c(length(msbb_gse84422.exp_Astrocytes[[m]]),
                                                                    length(zhang_celltype_ADgenes.list$Astrocytes)-length(msbb_gse84422_DCG_Astrocytes.overlap[[m]]),
                                                                    lapply(msbb_gse84422.DCGs_list,length)[[m]]-length(msbb_gse84422_DCG_Astrocytes.overlap[[m]]),
                                                                    21982-(length(zhang_celltype_ADgenes.list$Astrocytes)-length(msbb_gse84422_DCG_Astrocytes.overlap[[m]]))-
                                                                      lapply(msbb_gse84422.DCGs_list,length)[[m]]-length(msbb_gse84422_DCG_Astrocytes.overlap[[m]])),nrow = 2))
}
for(m in 1:length(msbb_gse84422_Endothelial.fisher_results)){
  msbb_gse84422_Endothelial.fisher_results[[m]]=fisher.test(matrix(c(length(msbb_gse84422.exp_Endothelial[[m]]),
                                                                     length(zhang_celltype_ADgenes.list$Endothelial)-length(msbb_gse84422_DCG_Endothelial.overlap[[m]]),
                                                                     lapply(msbb_gse84422.DCGs_list,length)[[m]]-length(msbb_gse84422_DCG_Endothelial.overlap[[m]]),
                                                                     21982-(length(zhang_celltype_ADgenes.list$Endothelial)-length(msbb_gse84422_DCG_Endothelial.overlap[[m]]))-
                                                                       lapply(msbb_gse84422.DCGs_list,length)[[m]]-length(msbb_gse84422_DCG_Endothelial.overlap[[m]])),nrow = 2))
}
for(m in 1:length(msbb_gse84422_Microglia.fisher_results)){
  msbb_gse84422_Microglia.fisher_results[[m]]=fisher.test(matrix(c(length(msbb_gse84422.exp_Microglia[[m]]),
                                                                   length(zhang_celltype_ADgenes.list$Microglia)-length(msbb_gse84422_DCG_Microglia.overlap[[m]]),
                                                                   lapply(msbb_gse84422.DCGs_list,length)[[m]]-length(msbb_gse84422_DCG_Microglia.overlap[[m]]),
                                                                   21982-(length(zhang_celltype_ADgenes.list$Microglia)-length(msbb_gse84422_DCG_Microglia.overlap[[m]]))-
                                                                     lapply(msbb_gse84422.DCGs_list,length)[[m]]-length(msbb_gse84422_DCG_Microglia.overlap[[m]])),nrow = 2))
}
for(m in 1:length(msbb_gse84422_Neurons.fisher_results)){
  msbb_gse84422_Neurons.fisher_results[[m]]=fisher.test(matrix(c(length(msbb_gse84422.exp_Neurons[[m]]),
                                                                 length(zhang_celltype_ADgenes.list$Neurons)-length(msbb_gse84422_DCG_Neurons.overlap[[m]]),
                                                                 lapply(msbb_gse84422.DCGs_list,length)[[m]]-length(msbb_gse84422_DCG_Neurons.overlap[[m]]),
                                                                 21982-(length(zhang_celltype_ADgenes.list$Neurons)-length(msbb_gse84422_DCG_Neurons.overlap[[m]]))-
                                                                   lapply(msbb_gse84422.DCGs_list,length)[[m]]-length(msbb_gse84422_DCG_Neurons.overlap[[m]])),nrow = 2))
}
for(m in 1:length(msbb_gse84422_Oligodendrocytes.fisher_results)){
  msbb_gse84422_Oligodendrocytes.fisher_results[[m]]=fisher.test(matrix(c(length(msbb_gse84422.exp_Oligodendrocytes[[m]]),
                                                                          length(zhang_celltype_ADgenes.list$Oligodendrocytes)-length(msbb_gse84422_DCG_Oligodendrocytes.overlap[[m]]),
                                                                          lapply(msbb_gse84422.DCGs_list,length)[[m]]-length(msbb_gse84422_DCG_Oligodendrocytes.overlap[[m]]),
                                                                          21982-(length(zhang_celltype_ADgenes.list$Oligodendrocytes)-length(msbb_gse84422_DCG_Oligodendrocytes.overlap[[m]]))-
                                                                            lapply(msbb_gse84422.DCGs_list,length)[[m]]-length(msbb_gse84422_DCG_Oligodendrocytes.overlap[[m]])),nrow = 2))
}

res.DCG_Astrocytes=foreach(i=1:length(msbb_gse84422_DCG_Astrocytes.overlap))%do%{
  
  
  data.frame(AstrocyteMarkers_in_DCGs=unlist(lapply(msbb_gse84422_DCG_Astrocytes.overlap[i],paste,collapse=",")),
             Nr_AstrocyteMarkers_in_DCGs=unlist(lapply(lapply(msbb_gse84422_DCG_Astrocytes.overlap,length)[i],paste,collapse=",")),
             DCGs=unlist(lapply(msbb_gse84422.DCGs_list,length)[[i]]),
             Expected_AstrocyteMarkers_in_DCGs=unname(unlist(msbb_gse84422.exp_Astrocytes[i])),
             Fisher_statistic=unname(unlist(lapply(msbb_gse84422_Astrocytes.fisher_results[i],function(s)s$estimate))),
             pval=unname(unlist(lapply(msbb_gse84422_Astrocytes.fisher_results[i],function(p)p$p.value))),
             stringsAsFactors = F)
}
names(res.DCG_Astrocytes)=names(msbb_gse84422_DCG_Astrocytes.overlap)
res.DCG_Astrocytes_df=data.frame(rbindlist(res.DCG_Astrocytes),stringsAsFactors = F)
rownames(res.DCG_Astrocytes_df)=names(msbb_gse84422_DCG_Astrocytes.overlap)
res.DCG_Astrocytes_df$adj.p=p.adjust(p = res.DCG_Astrocytes_df$pval,method = "fdr")
res.DCG_Astrocytes_df$Expected_AstrocyteMarkers_in_DCGs=as.numeric(res.DCG_Astrocytes_df$Expected_AstrocyteMarkers_in_DCGs)
res.DCG_Astrocytes_df$Nr_AstrocyteMarkers_in_DCGs=as.numeric(res.DCG_Astrocytes_df$Nr_AstrocyteMarkers_in_DCGs)
res.DCG_Astrocytes_df=res.DCG_Astrocytes_df%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(Nr_AstrocyteMarkers_in_DCGs>Expected_AstrocyteMarkers_in_DCGs,true = "Over",false = "Under"))


res.DCG_Endothelial=foreach(i=1:length(msbb_gse84422_DCG_Endothelial.overlap))%do%{
  
  
  data.frame(EndothelialMarkers_in_DCGs=unlist(lapply(msbb_gse84422_DCG_Endothelial.overlap[i],paste,collapse=",")),
             Nr_EndothelialMarkers_in_DCGs=unlist(lapply(lapply(msbb_gse84422_DCG_Endothelial.overlap,length)[i],paste,collapse=",")),
             DCGs=unlist(lapply(msbb_gse84422.DCGs_list,length)[[i]]),
             Expected_EndothelialMarkers_in_DCGs=unname(unlist(msbb_gse84422.exp_Endothelial[i])),
             Fisher_statistic=unname(unlist(lapply(msbb_gse84422_Endothelial.fisher_results[i],function(s)s$estimate))),
             pval=unname(unlist(lapply(msbb_gse84422_Endothelial.fisher_results[i],function(p)p$p.value))),
             stringsAsFactors = F)
}
names(res.DCG_Endothelial)=names(msbb_gse84422_DCG_Endothelial.overlap)
res.DCG_Endothelial_df=data.frame(rbindlist(res.DCG_Endothelial),stringsAsFactors = F)
rownames(res.DCG_Endothelial_df)=names(msbb_gse84422_DCG_Endothelial.overlap)
res.DCG_Endothelial_df$adj.p=p.adjust(p = res.DCG_Endothelial_df$pval,method = "fdr")
res.DCG_Endothelial_df$Expected_EndothelialMarkers_in_DCGs=as.numeric(res.DCG_Endothelial_df$Expected_EndothelialMarkers_in_DCGs)
res.DCG_Endothelial_df$Nr_EndothelialMarkers_in_DCGs=as.numeric(res.DCG_Endothelial_df$Nr_EndothelialMarkers_in_DCGs)
res.DCG_Endothelial_df=res.DCG_Endothelial_df%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(Nr_EndothelialMarkers_in_DCGs>Expected_EndothelialMarkers_in_DCGs,true = "Over",false = "Under"))

res.DCG_Microglia=foreach(i=1:length(msbb_gse84422_DCG_Microglia.overlap))%do%{
  
  
  data.frame(MicrogliaMarkers_in_DCGs=unlist(lapply(msbb_gse84422_DCG_Microglia.overlap[i],paste,collapse=",")),
             Nr_MicrogliaMarkers_in_DCGs=unlist(lapply(lapply(msbb_gse84422_DCG_Microglia.overlap,length)[i],paste,collapse=",")),
             DCGs=unlist(lapply(msbb_gse84422.DCGs_list,length)[[i]]),
             Expected_MicrogliaMarkers_in_DCGs=unname(unlist(msbb_gse84422.exp_Microglia[i])),
             Fisher_statistic=unname(unlist(lapply(msbb_gse84422_Microglia.fisher_results[i],function(s)s$estimate))),
             pval=unname(unlist(lapply(msbb_gse84422_Microglia.fisher_results[i],function(p)p$p.value))),
             stringsAsFactors = F)
}
names(res.DCG_Microglia)=names(msbb_gse84422_DCG_Microglia.overlap)
res.DCG_Microglia_df=data.frame(rbindlist(res.DCG_Microglia),stringsAsFactors = F)
rownames(res.DCG_Microglia_df)=names(msbb_gse84422_DCG_Microglia.overlap)
res.DCG_Microglia_df$adj.p=p.adjust(p = res.DCG_Microglia_df$pval,method = "fdr")
res.DCG_Microglia_df$Expected_MicrogliaMarkers_in_DCGs=as.numeric(res.DCG_Microglia_df$Expected_MicrogliaMarkers_in_DCGs)
res.DCG_Microglia_df$Nr_MicrogliaMarkers_in_DCGs=as.numeric(res.DCG_Microglia_df$Nr_MicrogliaMarkers_in_DCGs)
res.DCG_Microglia_df=res.DCG_Microglia_df%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(Nr_MicrogliaMarkers_in_DCGs>Expected_MicrogliaMarkers_in_DCGs,true = "Over",false = "Under"))

res.DCG_Neurons=foreach(i=1:length(msbb_gse84422_DCG_Neurons.overlap))%do%{
  
  
  data.frame(NeuronsMarkers_in_DCGs=unlist(lapply(msbb_gse84422_DCG_Neurons.overlap[i],paste,collapse=",")),
             Nr_NeuronsMarkers_in_DCGs=unlist(lapply(lapply(msbb_gse84422_DCG_Neurons.overlap,length)[i],paste,collapse=",")),
             DCGs=unlist(lapply(msbb_gse84422.DCGs_list,length)[[i]]),
             Expected_NeuronsMarkers_in_DCGs=unname(unlist(msbb_gse84422.exp_Neurons[i])),
             Fisher_statistic=unname(unlist(lapply(msbb_gse84422_Neurons.fisher_results[i],function(s)s$estimate))),
             pval=unname(unlist(lapply(msbb_gse84422_Neurons.fisher_results[i],function(p)p$p.value))),
             stringsAsFactors = F)
}
names(res.DCG_Neurons)=names(msbb_gse84422_DCG_Neurons.overlap)
res.DCG_Neurons_df=data.frame(rbindlist(res.DCG_Neurons),stringsAsFactors = F)
rownames(res.DCG_Neurons_df)=names(msbb_gse84422_DCG_Neurons.overlap)
res.DCG_Neurons_df$adj.p=p.adjust(p = res.DCG_Neurons_df$pval,method = "fdr")
res.DCG_Neurons_df$Expected_NeuronsMarkers_in_DCGs=as.numeric(res.DCG_Neurons_df$Expected_NeuronsMarkers_in_DCGs)
res.DCG_Neurons_df$Nr_NeuronsMarkers_in_DCGs=as.numeric(res.DCG_Neurons_df$Nr_NeuronsMarkers_in_DCGs)
res.DCG_Neurons_df=res.DCG_Neurons_df%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(Nr_NeuronsMarkers_in_DCGs>Expected_NeuronsMarkers_in_DCGs,true = "Over",false = "Under"))

res.DCG_Oligodendrocytes=foreach(i=1:length(msbb_gse84422_DCG_Oligodendrocytes.overlap))%do%{
  
  
  data.frame(OligodendrocytesMarkers_in_DCGs=unlist(lapply(msbb_gse84422_DCG_Oligodendrocytes.overlap[i],paste,collapse=",")),
             Nr_OligodendrocytesMarkers_in_DCGs=unlist(lapply(lapply(msbb_gse84422_DCG_Oligodendrocytes.overlap,length)[i],paste,collapse=",")),
             DCGs=unlist(lapply(msbb_gse84422.DCGs_list,length)[[i]]),
             Expected_OligodendrocytesMarkers_in_DCGs=unname(unlist(msbb_gse84422.exp_Oligodendrocytes[i])),
             Fisher_statistic=unname(unlist(lapply(msbb_gse84422_Oligodendrocytes.fisher_results[i],function(s)s$estimate))),
             pval=unname(unlist(lapply(msbb_gse84422_Oligodendrocytes.fisher_results[i],function(p)p$p.value))),
             stringsAsFactors = F)
}
names(res.DCG_Oligodendrocytes)=names(msbb_gse84422_DCG_Oligodendrocytes.overlap)
res.DCG_Oligodendrocytes_df=data.frame(rbindlist(res.DCG_Oligodendrocytes),stringsAsFactors = F)
rownames(res.DCG_Oligodendrocytes_df)=names(msbb_gse84422_DCG_Oligodendrocytes.overlap)
res.DCG_Oligodendrocytes_df$adj.p=p.adjust(p = res.DCG_Oligodendrocytes_df$pval,method = "fdr")
res.DCG_Oligodendrocytes_df$Expected_OligodendrocytesMarkers_in_DCGs=as.numeric(res.DCG_Oligodendrocytes_df$Expected_OligodendrocytesMarkers_in_DCGs)
res.DCG_Oligodendrocytes_df$Nr_OligodendrocytesMarkers_in_DCGs=as.numeric(res.DCG_Oligodendrocytes_df$Nr_OligodendrocytesMarkers_in_DCGs)
res.DCG_Oligodendrocytes_df=res.DCG_Oligodendrocytes_df%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(Nr_OligodendrocytesMarkers_in_DCGs>Expected_OligodendrocytesMarkers_in_DCGs,true = "Over",false = "Under"))

write.table(res.DCG_Astrocytes_df,"DCG_Astrocyte_Enrichment.txt",sep = "\t",col.names = T,row.names = F,quote = F)
write.table(res.DCG_Endothelial_df,"DCG_Endothelial_Enrichment.txt",sep = "\t",col.names = T,row.names = F,quote = F)
write.table(res.DCG_Microglia_df,"DCG_Microglia_Enrichment.txt",sep = "\t",col.names = T,row.names = F,quote = F)
write.table(res.DCG_Neurons_df,"DCG_Neurons_Enrichment.txt",sep = "\t",col.names = T,row.names = F,quote = F)
write.table(res.DCG_Oligodendrocytes_df,"DCG_Oligodendrocytes_Enrichment.txt",sep = "\t",col.names = T,row.names = F,quote = F)

#Validation with Berchtold dataset

gse48350_diffcoexp_files=list.files(path = "/Users/sandeepamberkar/Work/Data/AD_GSE48350/",pattern = "*diffcoexp.RDS",full.names = T)

gse48350_brain_regions.diffcoexp=vector(mode = "list",length = 4)
names(gse48350_brain_regions.diffcoexp)=gsub(pattern = "./|diffcoexp.RDS",replacement = "",x = gse48350_diffcoexp_files)%>%gsub(pattern = "/UsersandeepamberkaWorDatAD_GSE4835/",replacement = "")
for(f in 1:4){
  gse48350_brain_regions.diffcoexp[[f]]=readRDS(gse48350_diffcoexp_files[f])  
}
gse48350.DCGs=lapply(gse48350_brain_regions.diffcoexp,function(x)x$DCGs)
gse48350.DCLs=lapply(gse48350_brain_regions.diffcoexp,function(x)x$DCLs)

