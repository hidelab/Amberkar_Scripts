library(dplyr)
library(magrittr)
library(data.table)
library(DCGL)
library(enrichR)
library(limma)
library(tibble)
library(igraph)
library(org.Hs.eg.db)
library(circlize)

dbs <- listEnrichrDbs()
kegg_dbs=grep(pattern = "KEGG",x = dbs$libraryName,value = T)[3]

jaccard=function(A,B){
  jc=set_cardinality(intersect(A,B))/set_cardinality(union(A,B))
  return(jc)
}

earlyAD_diffcoexp_files=list.files(path = "/Users/sandeepamberkar/Work/Data/MSMM-MSBB-HBTRC-Data/MSBB_Array19/GSE84422/Lobe_Analysis/Dyscorrelation_Landscape_Results/EarlyAD_diffcoexp/",pattern = "earlyAD_diffcoexp",full.names = T)
earlyAD_samples=readRDS("/Users/sandeepamberkar/Work/Data/MSMM-MSBB-HBTRC-Data/MSBB_Array19/GSE84422/Lobe_Analysis/Dyscorrelation_Landscape_Results/EarlyAD_diffcoexp/msbb_gse84422_GPL96_97_earlyAD_samplesToAnalyse.RDS")
earlyAD_samples.exprs=readRDS("/Users/sandeepamberkar/Work/Data/MSMM-MSBB-HBTRC-Data/MSBB_Array19/GSE84422/Lobe_Analysis/Dyscorrelation_Landscape_Results/EarlyAD_diffcoexp/msbb_gse84422_GPL96_97_earlyAD_samplesToAnalyse_exprs.RDS")
setwd("/Users/sandeepamberkar/Work/Data/MSMM-MSBB-HBTRC-Data/MSBB_Array19/GSE84422/Lobe_Analysis/Dyscorrelation_Landscape_Results/LateAD_diffcoexp/")
lateAD_diffcoexp_files=list.files(path = "/Users/sandeepamberkar/Work/Data/MSMM-MSBB-HBTRC-Data/MSBB_Array19/GSE84422/Lobe_Analysis/Dyscorrelation_Landscape_Results/LateAD_diffcoexp/",pattern = "lateAD_diffcoexp",full.names = T)
lateAD_samples=readRDS("/Users/sandeepamberkar/Work/Data/MSMM-MSBB-HBTRC-Data/MSBB_Array19/GSE84422/Lobe_Analysis/Dyscorrelation_Landscape_Results/LateAD_diffcoexp/msbb_gse84422_GPL96_97_lateAD_samplesToAnalyse.RDS")
lateAD_samples.exprs=readRDS("/Users/sandeepamberkar/Work/Data/MSMM-MSBB-HBTRC-Data/MSBB_Array19/GSE84422/Lobe_Analysis/Dyscorrelation_Landscape_Results/LateAD_diffcoexp/msbb_gse84422_GPL96_97_lateAD_samplesToAnalyse_exprs.RDS")

regnet_tf2target_data=fread("/Users/sandeepamberkar/Work/Data/TF_Databases/RegNetwork_human_regulators2.txt",header = T,sep = "\t",showProgress = T,data.table = F)
regnet_tf2target.HGNC=fread("/Users/sandeepamberkar/Work/Data/TF_Databases/RegNetwork_human_regulators2.txt",header = T,sep = "\t",showProgress = T,data.table = F)%>%dplyr::select(c(regulator_symbol,target_symbol))

#DEG analyses to compare with DCGs
earlyAD_DEGs=lateAD_DEGs=vector(mode = "list",length = length(earlyAD_samples.exprs))
names(earlyAD_DEGs)=names(lateAD_DEGs)=names(earlyAD_samples.exprs)

for(i in 1:length(earlyAD_DEGs)){
  c_exprs=earlyAD_samples.exprs[[i]][,earlyAD_samples[[i]]$SampleType=="CONTROL"]
  d_exprs=earlyAD_samples.exprs[[i]][,earlyAD_samples[[i]]$SampleType=="AD"]
  group=factor(c(rep("CONTROL",dim(c_exprs)[2]),rep("AD",dim(d_exprs)[2])))
  design_df=model.matrix(~0+group)
  colnames(design_df)=levels(group)
  contrasts_matrix=makeContrasts(CONTROL-AD,levels = design_df)
  fit=lmFit(object = cbind.data.frame(c_exprs,d_exprs),design = design_df)
  fit2=contrasts.fit(fit,contrasts_matrix)
  fit2=eBayes(fit2,trend = T)
  earlyAD_DEGs[[i]]=topTable(fit2,coef = 1,number = 19530,p.value = 0.1,adjust.method = "BH")%>%rownames_to_column("Gene")
}
for(i in 1:length(lateAD_DEGs)){
  c_exprs=lateAD_samples.exprs[[i]][,lateAD_samples[[i]]$SampleType=="CONTROL"]
  d_exprs=lateAD_samples.exprs[[i]][,lateAD_samples[[i]]$SampleType=="AD"]
  group=factor(c(rep("CONTROL",dim(c_exprs)[2]),rep("AD",dim(d_exprs)[2])))
  design_df=model.matrix(~0+group)
  colnames(design_df)=levels(group)
  contrasts_matrix=makeContrasts(CONTROL-AD,levels = design_df)
  fit=lmFit(object = cbind.data.frame(c_exprs,d_exprs),design = design_df)
  fit2=contrasts.fit(fit,contrasts_matrix)
  fit2=eBayes(fit2,trend = T)
  lateAD_DEGs[[i]]=topTable(fit2,coef = 1,number = 19530,p.value = 0.1,adjust.method = "BH")%>%rownames_to_column("Gene")
}


earlyAD_diffcoexp_results=lateAD_diffcoexp_results=vector(mode = "list",length = length(earlyAD_diffcoexp_files))
for(i in 1:length(earlyAD_diffcoexp_files)){
  earlyAD_diffcoexp_results[[i]]=readRDS(earlyAD_diffcoexp_files[[i]])
  lateAD_diffcoexp_results[[i]]=readRDS(lateAD_diffcoexp_files[[i]])
}
names(earlyAD_diffcoexp_results)=names(lateAD_diffcoexp_results)=names(earlyAD_samples.exprs)

random_earlyAD_DCGs=replicate(n = 1000,expr = sample(x = rownames(lateAD_samples.exprs$Frontal_Pole),size = length(earlyAD_diffcoexp_results$Frontal_Pole$DCGs$Gene),replace = F),simplify = T)
random_lateAD_DCGs=replicate(n = 1000,expr = sample(x = rownames(lateAD_samples.exprs$Frontal_Pole),size = length(lateAD_diffcoexp_results$Frontal_Pole$DCGs$Gene),replace = F),simplify = T)
common_DCGs=intersect(earlyAD_diffcoexp_results$Frontal_Pole$DCGs$Gene,lateAD_diffcoexp_results$Frontal_Pole$DCGs$Gene)

res=foreach(i=1:1000)%do%{
  length(intersect(random_earlyAD_DCGs[,i],random_lateAD_DCGs[,i]))
}

earlyAD_DCGs=earlyAD_DCLs=earlyAD_DCGs.list=earlyAD_DCLs.list=earlyAD_DRsort=earlyAD_DRGs=earlyAD_DRGs.list=earlyAD_DRLs=earlyAD_DRLs.list=earlyAD_bridged_DCLs=vector(mode = "list",length = length(earlyAD_diffcoexp_results))
names(earlyAD_DCGs)=names(earlyAD_DCLs)=names(earlyAD_DCGs.list)=names(earlyAD_DCLs.list)=names(earlyAD_DRsort)=names(earlyAD_DRGs)=names(earlyAD_DRGs.list)=names(earlyAD_DRLs)=names(earlyAD_DRLs.list)=names(earlyAD_bridged_DCLs)

earlyAD_DCGs=lapply(earlyAD_diffcoexp_results,function(x)x$DCGs)
earlyAD_DCGs.filtered=lapply(earlyAD_diffcoexp_results,function(x)x$DCGs%>%dplyr::filter(q<=0.1))
earlyAD_DCLs=lapply(earlyAD_diffcoexp_results,function(x)x$DCLs)
earlyAD_DCLs.filtered=lapply(earlyAD_diffcoexp_results,function(x)x$DCLs%>%dplyr::filter(q.diffcor<=0.1))
earlyAD_DCLs.filtered=earlyAD_DCLs.filtered[-12]
earlyAD_DCGs.list=lapply(earlyAD_DCGs.filtered,function(x)x%>%pull(Gene))
earlyAD_DCGs.list=earlyAD_DCGs.list[-12]
earlyAD_DCLs.subgraph=lapply(earlyAD_DCLs.filtered,function(x)graph.data.frame(d = x[,c('Gene.1','Gene.2')],directed = F))

earlyAD_DCGs.Entrez_list=lapply(earlyAD_DCGs.list,function(x)unname(mapIds(x = org.Hs.eg.db,keys = x,column = "ENTREZID",keytype = "SYMBOL")))

# earlyAD_RS_DCGs.list=earlyAD.RS_DCGs=vector(mode = "list",length = length(earlyAD_samples.exprs[-12]))
# names(earlyAD_RS_DCGs.list)=names(earlyAD.RS_DCGs)=names(earlyAD_samples.exprs)[-12]

# for(i in 1:length(earlyAD_DCGs.list)){
#   earlyAD_RS_DCGs.list[[i]]=setdiff(earlyAD_DCGs.list[[i]],(Reduce(union,lapply(earlyAD_DCGs.list,function(x)intersect(x,earlyAD_DCGs.list[[i]]))[-i])))
#   earlyAD.RS_DCGs[[i]]=earlyAD_DCGs.filtered[[i]][earlyAD_RS_DCGs.list[[i]]%in%earlyAD_DCGs.filtered[[i]]$Gene,]
# }
# earlyAD_RS_DCGs.Entrez_list=lapply(earlyAD_RS_DCGs.list,function(x)unname(mapIds(x = org.Hs.eg.db,keys = x,column = "ENTREZID",keytype = "SYMBOL")))
# earlyAD_RS_DCGs.KEGG=lapply(earlyAD_RS_DCGs.list,function(x)enrichr(genes = x,databases = kegg_dbs)[[1]]%>%dplyr::filter(Adjusted.P.value<=0.1))
# earlyAD_RS_DCGs.KEGG2=lapply(earlyAD_RS_DCGs.Entrez_list,function(x)data.frame(enrichKEGG(gene = x,organism = "hsa",pvalueCutoff = 0.1,pAdjustMethod = "BH")))

#Functional enrichments in MSigdb collections
msigdb_c7=read.gmt(gmtfile = "/Users/sandeepamberkar/Work/Software/MSigDB_Collections/c7.all.v6.1.symbols.gmt")
msigdb_c2=read.gmt(gmtfile = "/Users/sandeepamberkar/Work/Software/MSigDB_Collections/c2.all.v6.1.symbols.gmt")
msigdb_h=read.gmt(gmtfile = "/Users/sandeepamberkar/Work/Software/MSigDB_Collections/h.all.v6.1.symbols.gmt")
msigdb_biocarta=read.gmt(gmtfile = "/Users/sandeepamberkar/Work/Software/MSigDB_Collections/c2.cp.biocarta.v6.1.symbols.gmt")
msigdb_c5.BP=read.gmt(gmtfile = "/Users/sandeepamberkar/Work/Software/MSigDB_Collections/c5.bp.v6.1.symbols.gmt")
msigdb_c5.CC=read.gmt(gmtfile = "/Users/sandeepamberkar/Work/Software/MSigDB_Collections/c5.cc.v6.1.symbols.gmt")
msigdb_c5.MF=read.gmt(gmtfile = "/Users/sandeepamberkar/Work/Software/MSigDB_Collections/c5.mf.v6.1.symbols.gmt")

# earlyAD_RS_DCGs.c7=lapply(earlyAD_RS_DCGs.list,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c7),stringsAsFactors = F))
# earlyAD_RS_DCGs.c2=lapply(earlyAD_RS_DCGs.list,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c2),stringsAsFactors = F))
# earlyAD_RS_DCGs.h=lapply(earlyAD_RS_DCGs.list,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_h),stringsAsFactors = F))
# earlyAD_RS_DCGs.biocarta=lapply(earlyAD_RS_DCGs.list,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_biocarta),stringsAsFactors = F))
# earlyAD_RS_DCGs.BP=lapply(earlyAD_RS_DCGs.list,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c5.BP),stringsAsFactors = F))
# earlyAD_RS_DCGs.CC=lapply(earlyAD_RS_DCGs.list,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c5.CC),stringsAsFactors = F))
# earlyAD_RS_DCGs.MF=lapply(earlyAD_RS_DCGs.list,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c5.MF),stringsAsFactors = F))

earlyAD_DCGs.c7=lapply(earlyAD_DCGs.list,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c7),stringsAsFactors = F))
earlyAD_DCGs.c2=lapply(earlyAD_DCGs.list,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c2),stringsAsFactors = F))
earlyAD_DCGs.h=lapply(earlyAD_DCGs.list,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_h),stringsAsFactors = F))
earlyAD_DCGs.biocarta=lapply(earlyAD_DCGs.list,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_biocarta),stringsAsFactors = F))
earlyAD_DCGs.BP=lapply(earlyAD_DCGs.list,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c5.BP),stringsAsFactors = F))
earlyAD_DCGs.CC=lapply(earlyAD_DCGs.list,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c5.CC),stringsAsFactors = F))
earlyAD_DCGs.MF=lapply(earlyAD_DCGs.list,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c5.MF),stringsAsFactors = F))




for(t in 1:length(earlyAD_samples.exprs)){
  earlyAD_DRsort[[t]]=DRsort(DCGs = earlyAD_DCGs[[t]],DCLs = earlyAD_DCLs[[t]],tf2target = regnet_tf2target.HGNC,expGenes = rownames(earlyAD_samples.exprs[[t]]))
}

names(earlyAD_DRsort)=names(earlyAD_samples.exprs)

earlyAD_DRGs=lapply(earlyAD_DRsort,function(x)x$DRGs%>%dplyr::filter(q<=0.1&DCGisTF=="TRUE"))
earlyAD_DRGs.list=lapply(earlyAD_DRGs,function(x)x$DCG%>%droplevels)

# earlyAD_RS_DRGs=mapply(FUN = function(a,b)a[which(a$DCG%in%b),],earlyAD_DRGs,earlyAD_RS_DCGs.list,SIMPLIFY = F)
# earlyAD_RS_DRGs.list=lapply(earlyAD_RS_DRGs,function(x)x%>%pull(1)%>%droplevels%>%levels)
# earlyAD_RS_DRGs.KEGG=lapply(earlyAD_RS_DRGs.list,function(x)enrichr(genes = x,databases = kegg_dbs)[[1]]%>%dplyr::filter(Adjusted.P.value<=0.1))
# 

earlyAD_DRLs=lapply(earlyAD_DRsort,function(x)x$DRLs%>%dplyr::filter(q.diffcor<=0.1))
earlyAD_DRLs.subgraph=lapply(earlyAD_DRLs,function(x)graph.data.frame(d = x[,c('Gene.1','Gene.2')],directed = F))
earlyAD_TF_bridged_DCL=lapply(earlyAD_DRsort,function(x)x$TF_bridged_DCL%>%filter(q.diffcor<=0.1)%>%dplyr::select(c(common.TFisDCG,common.TF,'Gene.1','Gene.2',DCG,'q.diffcor',type)))
earlyAD_TF_bridged_DCL.commonTFs_list=lapply(earlyAD_TF_bridged_DCL,function(x)x%>%pull('common.TF')%>%unique)
#Build graphs for TF-bridged DCLs, taking the top10 TFs

earlyAD_TF_bridged_DCL.top10_TFs=sort(table(unlist(lapply(earlyAD_TF_bridged_DCL,function(x)x%>%pull('common.TF')%>%unique))),decreasing = T)[1:10]

for(g in 1:10){
  filt_TF_bridged_DCL=Filter(f = function(y)dim(y)[1]>2,x = lapply(earlyAD_TF_bridged_DCL,function(x)x[x$common.TF%in%names(earlyAD_TF_bridged_DCL.top10_TFs[g]),c(2:5)]))
  filt_TF_bridged_DCL.subgraph=lapply(filt_TF_bridged_DCL,function(x)graph.data.frame(d = x[,c(2:3)],directed = F))
  #filt_TF_bridged_DCL.subgraph=lapply(filt_TF_bridged_DCL.subgraph,function(x){V(x)$borderWidth=5;x})
  #filt_TF_bridged_DCL.subgraph=lapply(filt_TF_bridged_DCL.subgraph,function(x){V(x)$borderWidth[V(x)$name%in%filt_TF_bridged_DCL[[g]]$DCG]=15;x})
  
  write.graph(graph =Reduce(graph.union,filt_TF_bridged_DCL.subgraph),file = paste("earlyAD",names(earlyAD_TF_bridged_DCL.top10_TFs[g]),"bridgedDCL.gml",sep = "_"),format = "gml")  
  for(l in 1:length(filt_TF_bridged_DCL)){
    write(filt_TF_bridged_DCL[[l]]$DCG,paste("earlyAD",names(filt_TF_bridged_DCL)[l],names(earlyAD_TF_bridged_DCL.top10_TFs[g]),"bridgedDCL_regulated_DCGs.txt",sep = "_"),sep = "\n")
  }
}





earlyAD_DCG.CompMatrix=matrix(NA,ncol=length(earlyAD_DCGs.list),nrow=length(earlyAD_DCGs.list))
rownames(earlyAD_DCG.CompMatrix)=colnames(earlyAD_DCG.CompMatrix)=c("FP","OVC","ITG","MTG","STG","PCC","AC","PHG","TP","IFG","DLPFC","SPL","HIPP")
for(i in 1:length(earlyAD_DCGs.list)){
  earlyAD_DCG.CompMatrix[i,]=unlist(lapply(earlyAD_DCGs.list,function(x)length(intersect(x,earlyAD_DCGs.list[[i]]))))
}
diag(earlyAD_DCG.CompMatrix)=0
chordDiagramFromMatrix(earlyAD_DCG.CompMatrix[-12,-12])




earlyAD_DRrank_TDD.files=earlyAD_DRrank_TED.files=earlyAD_DRrank.TDD=earlyAD_DRrank.TED=vector(mode = "list",length = length(earlyAD_DCGs))
names(earlyAD_DRrank_TDD.files)=names(earlyAD_DRrank_TED.files)=names(earlyAD_DRrank.TDD)=names(earlyAD_DRrank.TED)=names(earlyAD_DCGs)

earlyAD_DRrank_TDD.files=list.files(path = "EarlyAD_TDD",pattern = "TED",full.names = T)
earlyAD_DRrank_TED.files=list.files(path = "EarlyAD_TED",pattern = "TED_new_pbinom",full.names = T)
for(t in 1:length(earlyAD_DRrank_TED.files)){
  #earlyAD_DRrank.TDD[[t]]=readRDS(earlyAD_DRrank_TDD.files[[t]])
  earlyAD_DRrank.TED[[t]]=readRDS(earlyAD_DRrank_TED.files[[t]])
}
earlyAD_sigTF_bridgedDCLs=lapply(earlyAD_DRrank.TED,function(x)x%>%filter(`p.adjusted`<=0.1))
earlyAD_sigTF_bridgedDCLs=Filter(f = function(x)dim(x)[1]>0,x=earlyAD_sigTF_bridgedDCLs)
earlyAD_sigTF_bridgedDCLs.list=lapply(earlyAD_sigTF_bridgedDCLs,function(x)x%>%pull(regulator_gene_symbol)%>%droplevels%>%levels)


#In silico validation with genetic and epigenomic data

known_AD_genes=scan("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Gene_Lists/Known_AD_genes_2018_Final.txt",sep = "\n",what = "char")
janssen_GWGWAS_MAGMA_hits=scan("/Users/sandeepamberkar/Work/Data/Janssen-AD/Janssen-GWGWAS-MAGMA-SNPgenes.txt",sep = "\n",what = "char")
janssen_meta_hits=scan("/Users/sandeepamberkar/Work/Data/Janssen-AD/Janssen-Meta-Annotated-SNPgenes.txt",sep = "\n",what = "char")
earlyAD_janssenAD.overlap=lapply(earlyAD_RS_DCGs.list,intersect,janssen_GWGWAS_MAGMA_hits)

earlyAD_exp_DCGs=lapply(lapply(earlyAD_RS_DCGs.list,length),function(x)round(length(janssen_meta_hits)/19530*x,digits = 1))
earlyAD_janssenAD.overlap=lapply(earlyAD_RS_DCGs.list,function(x)intersect(x,janssen_meta_hits))
earlyAD_janssenAD.fisher_results=vector(mode = "list",length = 13)
for(m in 1:length(earlyAD_RS_DCGs.list)){
  earlyAD_janssenAD.fisher_results[[m]]=fisher.test(matrix(c(length(earlyAD_janssenAD.overlap[[m]]),
                                                             length(known_AD_genes)-length(earlyAD_janssenAD.overlap[[m]]),
                                                             lapply(earlyAD_RS_DCGs.list,length)[[m]]-length(earlyAD_janssenAD.overlap[[m]]),
                                                             19530-(length(known_AD_genes)-length(earlyAD_janssenAD.overlap[[m]]))-
                                                             lapply(earlyAD_RS_DCGs.list,length)[[m]]-length(earlyAD_janssenAD.overlap[[m]])),nrow = 2))
  
}
res.janssen_SAG_earlyAD_DCG=foreach(i=1:length(earlyAD_janssenAD.overlap))%do%{
  
  data.frame(SAGs_in_earlyAD_DCGs=unlist(lapply(earlyAD_janssenAD.overlap[i],length)),
             earlyAD_DCGs=unlist(lapply(earlyAD_RS_DCGs.list,length)[[i]]),
             Expected_SAGs_in_earlyAD_DCGs=unname(unlist(earlyAD_exp_DCGs[i])),
             Fisher_statistic=unname(unlist(lapply(earlyAD_janssenAD.fisher_results[i],function(s)s$estimate))),
             pval=unname(unlist(lapply(earlyAD_janssenAD.fisher_results[i],function(p)p$p.value))),
             stringsAsFactors = F)
}
names(res.janssen_SAG_earlyAD_DCG)=names(earlyAD_janssenAD.overlap)
res.SAG_earlyAD_DCG_df=data.frame(rbindlist(res.janssen_SAG_earlyAD_DCG),stringsAsFactors = F)
rownames(res.SAG_earlyAD_DCG_df)=names(earlyAD_janssenAD.overlap)
res.SAG_earlyAD_DCG_df$adj.p=p.adjust(p = res.SAG_earlyAD_DCG_df$pval,method = "fdr")
res.SAG_earlyAD_DCG_df$SAGs_in_earlyAD_DCGs=as.numeric(res.SAG_earlyAD_DCG_df$SAGs_in_earlyAD_DCGs)
res.SAG_earlyAD_DCG_df=res.SAG_earlyAD_DCG_df%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(SAGs_in_earlyAD_DCGs>Expected_SAGs_in_earlyAD_DCGs,true = "Over",false = "Under"))  
fwrite(res.SAG_earlyAD_DCG_df,"EarlyAD_SAG_earlyAD_DCG_EnrichmentSummary.txt",sep = "\t",col.names = T,row.names = F,quote = F)

dhmc_bennet.NP=read.table("/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/GW-DNA-hydroxymethylation-Bennet/1-s2.0-S155252601633059X-mmc4_NP.txt",sep = "\t",header = T,as.is = T)
dhmc_bennet.NFT=read.table("/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/GW-DNA-hydroxymethylation-Bennet/1-s2.0-S155252601633059X-mmc5_NFT.txt",sep = "\t",header = T,as.is = T)
d5mc_bernstein.Dhml=scan("/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/5hMC_Tau_AD/ddw109_Supp/Supplemental Table 6_Bernstein et al_R1.txt",sep = "\n",what = "char")

earlyAD_DCG_d5mc.overlap=lapply(earlyAD_DCGs.list,function(x)length(intersect(x,d5mc_bernstein.Dhml)))
dhmc_bennet_NFT.genes=dhmc_bennet.NFT%>%dplyr::filter(q.value_<=0.1)%>%pull(Nearest.gene)
dhmc_bennet_NP.genes=dhmc_bennet.NP%>%dplyr::filter(q.vlaue_<=0.1)%>%pull(Neartest.gene)
earlyAD_DCG_dhmc_NFT.overlap=lapply(earlyAD_DCGs.list,function(x)intersect(x,dhmc_bennet_NFT.genes))
earlyAD_DCG_dhmc_NP.overlap=lapply(earlyAD_DCGs.list,function(x)intersect(x,dhmc_bennet_NP.genes))
earlyAD_exp_dhmc_NP_DCGs=lapply(lapply(earlyAD_DCGs.list,length),function(x)round(length(dhmc_bennet_NP.genes)/19530*x,digits = 1))
earlyAD_exp_dhmc_NFT_DCGs=lapply(lapply(earlyAD_DCGs.list,length),function(x)round(length(dhmc_bennet_NFT.genes)/19530*x,digits = 1))

earlyAD_dhmc_NP.fisher_results=earlyAD_dhmc_NFT.fisher_results=earlyAD_d5mc.fisher_results=vector(mode = "list",length = length(earlyAD_DCGs.list))
names(earlyAD_dhmc_NP.fisher_results)=names(earlyAD_dhmc_NFT.fisher_results)=names(earlyAD_d5mc.fisher_results)=names(earlyAD_DCGs.list)
for(m in 1:length(earlyAD_RS_DCGs.list)){
  earlyAD_dhmc_NP.fisher_results[[m]]=fisher.test(matrix(c(length(earlyAD_DCG_dhmc_NP.overlap[[m]]),
                                                           length(dhmc_bennet_NP.genes)-length(earlyAD_DCG_dhmc_NP.overlap[[m]]),
                                                           lapply(earlyAD_DCGs.list,length)[[m]]-length(earlyAD_DCG_dhmc_NP.overlap[[m]]),
                                                           19530-(length(dhmc_bennet_NP.genes)-length(earlyAD_DCG_dhmc_NP.overlap[[m]]))-
                                                           lapply(earlyAD_DCGs.list,length)[[m]]-length(earlyAD_DCG_dhmc_NP.overlap[[m]])),nrow = 2))
}

res.dhmc_NP_DCG=foreach(i=1:length(earlyAD_DCG_dhmc_NP.overlap))%do%{
  
  
  data.frame(NPgenes_in_DCGs=unlist(lapply(earlyAD_DCG_dhmc_NP.overlap[i],length)),
             DCGs=unlist(lapply(earlyAD_DCGs.list,length)[[i]]),
             Expected_dhmc_NPgenes_in_DCGs=unname(unlist(earlyAD_exp_dhmc_NP_DCGs[i])),
             Fisher_statistic=unname(unlist(lapply(earlyAD_dhmc_NP.fisher_results[i],function(s)s$estimate))),
             pval=unname(unlist(lapply(earlyAD_dhmc_NP.fisher_results[i],function(p)p$p.value))),
             stringsAsFactors = F)
}


names(res.dhmc_NP_DCG)=names(earlyAD_DCG_dhmc_NP.overlap)
res.dhmc_NP_DCG_df=data.frame(rbindlist(res.dhmc_NP_DCG),stringsAsFactors = F)
rownames(res.dhmc_NP_DCG_df)=names(earlyAD_DCG_dhmc_NP.overlap)
res.dhmc_NP_DCG_df$adj.p=p.adjust(p = res.dhmc_NP_DCG_df$pval,method = "fdr")
res.dhmc_NP_DCG_df$NPgenes_in_DCGs=as.numeric(res.dhmc_NP_DCG_df$NPgenes_in_DCGs)
res.dhmc_NP_DCG_df=res.dhmc_NP_DCG_df%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(NPgenes_in_DCGs>Expected_dhmc_NPgenes_in_DCGs,true = "Over",false = "Under"))
fwrite(res.dhmc_NP_DCG_df,"EarlyAD_DhMR_NP_earlyAD_DCG_EnrichmentSummary.txt",sep = "\t",col.names = T,row.names = F,quote = F)


earlyAD_DCG_dhmc_NFT.overlap=lapply(earlyAD_DCGs.list,function(x)intersect(x,dhmc_bennet_NFT.genes))
for(m in 1:length(res.dhmc_NP_DCG_df)){
  earlyAD_dhmc_NFT.fisher_results[[m]]=fisher.test(matrix(c(length(earlyAD_DCG_dhmc_NFT.overlap[[m]]),
                                                            length(dhmc_bennet_NFT.genes)-length(earlyAD_DCG_dhmc_NFT.overlap[[m]]),
                                                            lapply(earlyAD_DCGs.list,length)[[m]]-length(earlyAD_DCG_dhmc_NFT.overlap[[m]]),
                                                            19530-(length(dhmc_bennet_NFT.genes)-length(earlyAD_DCG_dhmc_NFT.overlap[[m]]))-
                                                            lapply(earlyAD_DCGs.list,length)[[m]]-length(earlyAD_DCG_dhmc_NFT.overlap[[m]])),nrow = 2))
}

res.dhmc_NFT_DCG=foreach(i=1:length(earlyAD_DCG_dhmc_NFT.overlap))%do%{
  
  
  data.frame(NFTgenes_in_DCGs=unlist(lapply(earlyAD_DCG_dhmc_NFT.overlap[i],length)),
             DCGs=unlist(lapply(earlyAD_DCGs.list,length)[[i]]),
             Expected_dhmc_NFTgenes_in_DCGs=unname(unlist(earlyAD_exp_dhmc_NFT_DCGs[i])),
             Fisher_statistic=unname(unlist(lapply(earlyAD_dhmc_NFT.fisher_results[i],function(s)s$estimate))),
             pval=unname(unlist(lapply(earlyAD_dhmc_NFT.fisher_results[i],function(p)p$p.value))),
             stringsAsFactors = F)
}


names(res.dhmc_NFT_DCG)=names(earlyAD_DCG_dhmc_NFT.overlap)
res.dhmc_NFT_DCG_df=data.frame(rbindlist(res.dhmc_NFT_DCG),stringsAsFactors = F)
rownames(res.dhmc_NFT_DCG_df)=names(earlyAD_DCG_dhmc_NFT.overlap)
res.dhmc_NFT_DCG_df$adj.p=p.adjust(p = res.dhmc_NFT_DCG_df$pval,method = "fdr")
res.dhmc_NFT_DCG_df$NFTgenes_in_DCGs=as.numeric(res.dhmc_NFT_DCG_df$NFTgenes_in_DCGs)
res.dhmc_NFT_DCG_df=res.dhmc_NFT_DCG_df%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(NFTgenes_in_DCGs>Expected_dhmc_NFTgenes_in_DCGs,true = "Over",false = "Under"))
fwrite(res.dhmc_NFT_DCG_df,"EarlyAD_DhMR_NFT_earlyAD_DCG_EnrichmentSummary.txt",sep = "\t",col.names = T,row.names = F,quote = F)


#STG GW-Methylation dataset
stg_gw_methylation=fread("/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/STG-GW-Methylation-Watson/13073_2015_258_MOESM6_ESM.txt",sep = "\t",header = T,data.table = F)
stg_gw_methylation_dmr.genes=stg_gw_methylation%>%dplyr::filter(`FDR_DMR.q.value`<=0.1)%>%pull(`Closest.Gene`)%>%unique
earlyAD_DCG_stg_methylation_dmr.overlap=lapply(earlyAD_DCGs.list,function(x)intersect(x,stg_gw_methylation_dmr.genes))
earlyAD_exp_stg_methylation_dmr_DCGs=lapply(lapply(earlyAD_DCGs.list,length),function(x)round(length(stg_gw_methylation_dmr.genes)/19530*x,digits = 1))

earlAD_stg_methylation.fisher_results=vector(mode = "list",length = length(earlyAD_DCGs.list))
names(earlAD_stg_methylation.fisher_results)=names(earlyAD_DCGs.list)
for(m in 1:length(earlyAD_DCGs.list)){
  earlAD_stg_methylation.fisher_results[[m]]=fisher.test(matrix(c(length(earlyAD_DCG_stg_methylation_dmr.overlap[[m]]),
                                                            length(stg_gw_methylation_dmr.genes)-length(earlyAD_DCG_stg_methylation_dmr.overlap[[m]]),
                                                            lapply(earlyAD_DCGs.list,length)[[m]]-length(earlyAD_DCG_stg_methylation_dmr.overlap[[m]]),
                                                            19530-(length(stg_gw_methylation_dmr.genes)-length(earlyAD_DCG_stg_methylation_dmr.overlap[[m]]))-
                                                            lapply(earlyAD_DCGs.list,length)[[m]]-length(earlyAD_DCG_stg_methylation_dmr.overlap[[m]])),nrow = 2))
}

res.earlyAD_stg_methylation=foreach(i=1:length(earlyAD_DCG_stg_methylation_dmr.overlap))%do%{
  
  
  data.frame(Methylated_genes_in_DCGs=unlist(lapply(earlyAD_DCG_stg_methylation_dmr.overlap[i],length)),
             DCGs=unlist(lapply(earlyAD_DCGs.list,length)[[i]]),
             Expected_STG_methylation_in_DCGs=unname(unlist(earlyAD_exp_stg_methylation_dmr_DCGs[i])),
             Fisher_statistic=unname(unlist(lapply(earlAD_stg_methylation.fisher_results[i],function(s)s$estimate))),
             pval=unname(unlist(lapply(earlAD_stg_methylation.fisher_results[i],function(p)p$p.value))),
             stringsAsFactors = F)
}


names(res.earlyAD_stg_methylation)=names(earlyAD_DCG_stg_methylation_dmr.overlap)
res_earlyAD_stg_methylation.df=data.frame(rbindlist(res.earlyAD_stg_methylation),stringsAsFactors = F)
rownames(res_earlyAD_stg_methylation.df)=names(earlyAD_DCG_stg_methylation_dmr.overlap)
res_earlyAD_stg_methylation.df$adj.p=p.adjust(p = res_earlyAD_stg_methylation.df$pval,method = "fdr")
res_earlyAD_stg_methylation.df$Methylated_genes_in_DCGs=as.numeric(res_earlyAD_stg_methylation.df$Methylated_genes_in_DCGs)
res_earlyAD_stg_methylation.df=res_earlyAD_stg_methylation.df%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(Methylated_genes_in_DCGs>Expected_STG_methylation_in_DCGs,true = "Over",false = "Under"))






#Brain cell type markers Zhang et al.
zhang_celltype_ADgenes=read.xls('/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/Zhang_19BrainRegions_Paper/13073_2016_355_MOESM1_ESM.xlsx',skip=1,sheet=3,header=T,as.is=T)
zhang_celltype_ADgenes.list=vector(mode = 'list',length = 5)
names(zhang_celltype_ADgenes.list)=c("ast","end","mic","neu","oli")
zhang_celltype_ADgenes.list$ast=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = "Astrocytes",zhang_celltype_ADgenes$Cell.type)]
zhang_celltype_ADgenes.list$end=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = "Endothelial",zhang_celltype_ADgenes$Cell.type)]
zhang_celltype_ADgenes.list$mic=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = "Microglia",zhang_celltype_ADgenes$Cell.type)]
zhang_celltype_ADgenes.list$neu=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = "Neurons",zhang_celltype_ADgenes$Cell.type)]
zhang_celltype_ADgenes.list$oli=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = "Oligodendrocytes",zhang_celltype_ADgenes$Cell.type)]


earlyAD_DCG_Astrocytes.overlap=lapply(earlyAD_DCGs.list,intersect,zhang_celltype_ADgenes.list$ast)
earlyAD_DCG_Endothelial.overlap=lapply(earlyAD_DCGs.list,intersect,zhang_celltype_ADgenes.list$end)
earlyAD_DCG_Microglia.overlap=lapply(earlyAD_DCGs.list,intersect,zhang_celltype_ADgenes.list$mic)
earlyAD_DCG_Neurons.overlap=lapply(earlyAD_DCGs.list,intersect,zhang_celltype_ADgenes.list$neu)
earlyAD_DCG_Oligodendrocytes.overlap=lapply(earlyAD_DCGs.list,intersect,zhang_celltype_ADgenes.list$oli)


earlyAD.exp_Astrocytes=lapply(lapply(earlyAD_DCGs.list,length),function(x)round(length(zhang_celltype_ADgenes.list$ast)/19530*x,digits = 1))
earlyAD.exp_Endothelial=lapply(lapply(earlyAD_DCGs.list,length),function(x)round(length(zhang_celltype_ADgenes.list$end)/19530*x,digits = 1))
earlyAD.exp_Microglia=lapply(lapply(earlyAD_DCGs.list,length),function(x)round(length(zhang_celltype_ADgenes.list$mic)/19530*x,digits = 1))
earlyAD.exp_Neurons=lapply(lapply(earlyAD_DCGs.list,length),function(x)round(length(zhang_celltype_ADgenes.list$neu)/19530*x,digits = 1))
earlyAD.exp_Oligodendrocytes=lapply(lapply(earlyAD_DCGs.list,length),function(x)round(length(zhang_celltype_ADgenes.list$oli)/19530*x,digits = 1))


earlyAD_Astrocytes.fisher_results=earlyAD_Endothelial.fisher_results=earlyAD_Microglia.fisher_results=earlyAD_Neurons.fisher_results=earlyAD_Oligodendrocytes.fisher_results=vector(mode = "list",length = length(earlyAD_DCGs.list))
names(earlyAD_Astrocytes.fisher_results)=names(earlyAD_Endothelial.fisher_results)=names(earlyAD_Microglia.fisher_results)=names(earlyAD_Neurons.fisher_results)=names(earlyAD_Oligodendrocytes.fisher_results)=names(earlyAD_DCGs.list)

for(m in 1:length(earlyAD_Astrocytes.fisher_results)){
  earlyAD_Astrocytes.fisher_results[[m]]=fisher.test(matrix(c(length(earlyAD_DCG_Astrocytes.overlap[[m]]),
                                                              length(zhang_celltype_ADgenes.list$ast)-length(earlyAD_DCG_Astrocytes.overlap[[m]]),
                                                              lapply(earlyAD_DCGs.list,length)[[m]]-length(earlyAD_DCG_Astrocytes.overlap[[m]]),
                                                              19530-(length(zhang_celltype_ADgenes.list$ast)-length(earlyAD_DCG_Astrocytes.overlap[[m]]))-
                                                                lapply(earlyAD_DCGs.list,length)[[m]]-length(earlyAD_DCG_Astrocytes.overlap[[m]])),nrow = 2))
}
for(m in 1:length(earlyAD_Endothelial.fisher_results)){
  earlyAD_Endothelial.fisher_results[[m]]=fisher.test(matrix(c(length(earlyAD_DCG_Endothelial.overlap[[m]]),
                                                               length(zhang_celltype_ADgenes.list$end)-length(earlyAD_DCG_Endothelial.overlap[[m]]),
                                                               lapply(earlyAD_DCGs.list,length)[[m]]-length(earlyAD_DCG_Endothelial.overlap[[m]]),
                                                               19530-(length(zhang_celltype_ADgenes.list$end)-length(earlyAD_DCG_Endothelial.overlap[[m]]))-
                                                                 lapply(earlyAD_DCGs.list,length)[[m]]-length(earlyAD_DCG_Endothelial.overlap[[m]])),nrow = 2))
}
for(m in 1:length(earlyAD_Microglia.fisher_results)){
  earlyAD_Microglia.fisher_results[[m]]=fisher.test(matrix(c(length(earlyAD_DCG_Microglia.overlap[[m]]),
                                                             length(zhang_celltype_ADgenes.list$mic)-length(earlyAD_DCG_Microglia.overlap[[m]]),
                                                             lapply(earlyAD_DCGs.list,length)[[m]]-length(earlyAD_DCG_Microglia.overlap[[m]]),
                                                             19530-(length(zhang_celltype_ADgenes.list$mic)-length(earlyAD_DCG_Microglia.overlap[[m]]))-
                                                               lapply(earlyAD_DCGs.list,length)[[m]]-length(earlyAD_DCG_Microglia.overlap[[m]])),nrow = 2))
}
for(m in 1:length(earlyAD_Neurons.fisher_results)){
  earlyAD_Neurons.fisher_results[[m]]=fisher.test(matrix(c(length(earlyAD_DCG_Neurons.overlap[[m]]),
                                                           length(zhang_celltype_ADgenes.list$neu)-length(earlyAD_DCG_Neurons.overlap[[m]]),
                                                           lapply(earlyAD_DCGs.list,length)[[m]]-length(earlyAD_DCG_Neurons.overlap[[m]]),
                                                           19530-(length(zhang_celltype_ADgenes.list$neu)-length(earlyAD_DCG_Neurons.overlap[[m]]))-
                                                             lapply(earlyAD_DCGs.list,length)[[m]]-length(earlyAD_DCG_Neurons.overlap[[m]])),nrow = 2))
}
for(m in 1:length(earlyAD_Oligodendrocytes.fisher_results)){
  earlyAD_Oligodendrocytes.fisher_results[[m]]=fisher.test(matrix(c(length(earlyAD_DCG_Oligodendrocytes.overlap[[m]]),
                                                                    length(zhang_celltype_ADgenes.list$oli)-length(earlyAD_DCG_Oligodendrocytes.overlap[[m]]),
                                                                    lapply(earlyAD_DCGs.list,length)[[m]]-length(earlyAD_DCG_Oligodendrocytes.overlap[[m]]),
                                                                    19530-(length(zhang_celltype_ADgenes.list$oli)-length(earlyAD_DCG_Oligodendrocytes.overlap[[m]]))-
                                                                    lapply(earlyAD_DCGs.list,length)[[m]]-length(earlyAD_DCG_Oligodendrocytes.overlap[[m]])),nrow = 2))
}




res.earlyAD_DCG_Astrocytes=foreach(i=1:length(earlyAD_DCG_Astrocytes.overlap))%do%{
  
  
  data.frame(AstrocyteMarkers_in_DCGs=unlist(lapply(earlyAD_DCG_Astrocytes.overlap[i],paste,collapse=",")),
             Nr_AstrocyteMarkers_in_DCGs=unlist(lapply(lapply(earlyAD_DCG_Astrocytes.overlap,length)[i],paste,collapse=",")),
             DCGs=unlist(lapply(earlyAD_DCGs.list,length)[[i]]),
             Expected_AstrocyteMarkers_in_DCGs=unname(unlist(earlyAD.exp_Astrocytes[i])),
             Fisher_statistic=unname(unlist(lapply(earlyAD_Astrocytes.fisher_results[i],function(s)s$estimate))),
             pval=unname(unlist(lapply(earlyAD_Astrocytes.fisher_results[i],function(p)p$p.value))),
             stringsAsFactors = F)
}
names(res.earlyAD_DCG_Astrocytes)=names(earlyAD_DCG_Astrocytes.overlap)
res.earlyAD_DCG_Astrocytes_df=data.frame(rbindlist(res.earlyAD_DCG_Astrocytes),stringsAsFactors = F)
rownames(res.earlyAD_DCG_Astrocytes_df)=names(earlyAD_DCG_Astrocytes.overlap)
res.earlyAD_DCG_Astrocytes_df$adj.p=p.adjust(p = res.earlyAD_DCG_Astrocytes_df$pval,method = "fdr")
res.earlyAD_DCG_Astrocytes_df$Expected_AstrocyteMarkers_in_DCGs=as.numeric(res.earlyAD_DCG_Astrocytes_df$Expected_AstrocyteMarkers_in_DCGs)
res.earlyAD_DCG_Astrocytes_df$Nr_AstrocyteMarkers_in_DCGs=as.numeric(res.earlyAD_DCG_Astrocytes_df$Nr_AstrocyteMarkers_in_DCGs)
res.earlyAD_DCG_Astrocytes_df=res.earlyAD_DCG_Astrocytes_df%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(Nr_AstrocyteMarkers_in_DCGs>Expected_AstrocyteMarkers_in_DCGs,true = "Over",false = "Under"))%>%filter(Representation=="Over")


res.earlyAD_DCG_Endothelial=foreach(i=1:length(earlyAD_DCG_Endothelial.overlap))%do%{
  
  
  data.frame(EndothelialMarkers_in_DCGs=unlist(lapply(earlyAD_DCG_Endothelial.overlap[i],paste,collapse=",")),
             Nr_EndothelialMarkers_in_DCGs=unlist(lapply(lapply(earlyAD_DCG_Endothelial.overlap,length)[i],paste,collapse=",")),
             DCGs=unlist(lapply(earlyAD_DCGs.list,length)[[i]]),
             Expected_EndothelialMarkers_in_DCGs=unname(unlist(earlyAD.exp_Endothelial[i])),
             Fisher_statistic=unname(unlist(lapply(earlyAD_Endothelial.fisher_results[i],function(s)s$estimate))),
             pval=unname(unlist(lapply(earlyAD_Endothelial.fisher_results[i],function(p)p$p.value))),
             stringsAsFactors = F)
}
names(res.earlyAD_DCG_Endothelial)=names(earlyAD_DCG_Endothelial.overlap)
res.earlyAD_DCG_Endothelial_df=data.frame(rbindlist(res.earlyAD_DCG_Endothelial),stringsAsFactors = F)
rownames(res.earlyAD_DCG_Endothelial_df)=names(earlyAD_DCG_Endothelial.overlap)
res.earlyAD_DCG_Endothelial_df$adj.p=p.adjust(p = res.earlyAD_DCG_Endothelial_df$pval,method = "fdr")
res.earlyAD_DCG_Endothelial_df$Expected_EndothelialMarkers_in_DCGs=as.numeric(res.earlyAD_DCG_Endothelial_df$Expected_EndothelialMarkers_in_DCGs)
res.earlyAD_DCG_Endothelial_df$Nr_EndothelialMarkers_in_DCGs=as.numeric(res.earlyAD_DCG_Endothelial_df$Nr_EndothelialMarkers_in_DCGs)
res.earlyAD_DCG_Endothelial_df=res.earlyAD_DCG_Endothelial_df%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(Nr_EndothelialMarkers_in_DCGs>Expected_EndothelialMarkers_in_DCGs,true = "Over",false = "Under"))%>%filter(Representation=="Over")

res.earlyAD_DCG_Microglia=foreach(i=1:length(earlyAD_DCG_Microglia.overlap))%do%{
  
  
  data.frame(MicrogliaMarkers_in_DCGs=unlist(lapply(earlyAD_DCG_Microglia.overlap[i],paste,collapse=",")),
             Nr_MicrogliaMarkers_in_DCGs=unlist(lapply(lapply(earlyAD_DCG_Microglia.overlap,length)[i],paste,collapse=",")),
             DCGs=unlist(lapply(earlyAD_DCGs.list,length)[[i]]),
             Expected_MicrogliaMarkers_in_DCGs=unname(unlist(earlyAD.exp_Microglia[i])),
             Fisher_statistic=unname(unlist(lapply(earlyAD_Microglia.fisher_results[i],function(s)s$estimate))),
             pval=unname(unlist(lapply(earlyAD_Microglia.fisher_results[i],function(p)p$p.value))),
             stringsAsFactors = F)
}
names(res.earlyAD_DCG_Microglia)=names(earlyAD_DCG_Microglia.overlap)
res.earlyAD_DCG_Microglia_df=data.frame(rbindlist(res.earlyAD_DCG_Microglia),stringsAsFactors = F)
rownames(res.earlyAD_DCG_Microglia_df)=names(earlyAD_DCG_Microglia.overlap)
res.earlyAD_DCG_Microglia_df$adj.p=p.adjust(p = res.earlyAD_DCG_Microglia_df$pval,method = "fdr")
res.earlyAD_DCG_Microglia_df$Expected_MicrogliaMarkers_in_DCGs=as.numeric(res.earlyAD_DCG_Microglia_df$Expected_MicrogliaMarkers_in_DCGs)
res.earlyAD_DCG_Microglia_df$Nr_MicrogliaMarkers_in_DCGs=as.numeric(res.earlyAD_DCG_Microglia_df$Nr_MicrogliaMarkers_in_DCGs)
res.earlyAD_DCG_Microglia_df=res.earlyAD_DCG_Microglia_df%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(Nr_MicrogliaMarkers_in_DCGs>Expected_MicrogliaMarkers_in_DCGs,true = "Over",false = "Under"))%>%filter(Representation=="Over")

res.earlyAD_DCG_Neurons=foreach(i=1:length(earlyAD_DCG_Neurons.overlap))%do%{
  
  
  data.frame(NeuronsMarkers_in_DCGs=unlist(lapply(earlyAD_DCG_Neurons.overlap[i],paste,collapse=",")),
             Nr_NeuronsMarkers_in_DCGs=unlist(lapply(lapply(earlyAD_DCG_Neurons.overlap,length)[i],paste,collapse=",")),
             DCGs=unlist(lapply(earlyAD_DCGs.list,length)[[i]]),
             Expected_NeuronsMarkers_in_DCGs=unname(unlist(earlyAD.exp_Neurons[i])),
             Fisher_statistic=unname(unlist(lapply(earlyAD_Neurons.fisher_results[i],function(s)s$estimate))),
             pval=unname(unlist(lapply(earlyAD_Neurons.fisher_results[i],function(p)p$p.value))),
             stringsAsFactors = F)
}
names(res.earlyAD_DCG_Neurons)=names(earlyAD_DCG_Neurons.overlap)
res.earlyAD_DCG_Neurons_df=data.frame(rbindlist(res.earlyAD_DCG_Neurons),stringsAsFactors = F)
rownames(res.earlyAD_DCG_Neurons_df)=names(earlyAD_DCG_Neurons.overlap)
res.earlyAD_DCG_Neurons_df$adj.p=p.adjust(p = res.earlyAD_DCG_Neurons_df$pval,method = "fdr")
res.earlyAD_DCG_Neurons_df$Expected_NeuronsMarkers_in_DCGs=as.numeric(res.earlyAD_DCG_Neurons_df$Expected_NeuronsMarkers_in_DCGs)
res.earlyAD_DCG_Neurons_df$Nr_NeuronsMarkers_in_DCGs=as.numeric(res.earlyAD_DCG_Neurons_df$Nr_NeuronsMarkers_in_DCGs)
res.earlyAD_DCG_Neurons_df=res.earlyAD_DCG_Neurons_df%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(Nr_NeuronsMarkers_in_DCGs>Expected_NeuronsMarkers_in_DCGs,true = "Over",false = "Under"))%>%filter(Representation=="Over")

res.earlyAD_DCG_Oligodendrocytes=foreach(i=1:length(earlyAD_DCG_Oligodendrocytes.overlap))%do%{
  
  
  data.frame(OligodendrocytesMarkers_in_DCGs=unlist(lapply(earlyAD_DCG_Oligodendrocytes.overlap[i],paste,collapse=",")),
             Nr_OligodendrocytesMarkers_in_DCGs=unlist(lapply(lapply(earlyAD_DCG_Oligodendrocytes.overlap,length)[i],paste,collapse=",")),
             DCGs=unlist(lapply(earlyAD_DCGs.list,length)[[i]]),
             Expected_OligodendrocytesMarkers_in_DCGs=unname(unlist(earlyAD.exp_Oligodendrocytes[i])),
             Fisher_statistic=unname(unlist(lapply(earlyAD_Oligodendrocytes.fisher_results[i],function(s)s$estimate))),
             pval=unname(unlist(lapply(earlyAD_Oligodendrocytes.fisher_results[i],function(p)p$p.value))),
             stringsAsFactors = F)
}
names(res.earlyAD_DCG_Oligodendrocytes)=names(earlyAD_DCG_Oligodendrocytes.overlap)
res.earlyAD_DCG_Oligodendrocytes_df=data.frame(rbindlist(res.earlyAD_DCG_Oligodendrocytes),stringsAsFactors = F)
rownames(res.earlyAD_DCG_Oligodendrocytes_df)=names(earlyAD_DCG_Oligodendrocytes.overlap)
res.earlyAD_DCG_Oligodendrocytes_df$adj.p=p.adjust(p = res.earlyAD_DCG_Oligodendrocytes_df$pval,method = "fdr")
res.earlyAD_DCG_Oligodendrocytes_df$Expected_OligodendrocytesMarkers_in_DCGs=as.numeric(res.earlyAD_DCG_Oligodendrocytes_df$Expected_OligodendrocytesMarkers_in_DCGs)
res.earlyAD_DCG_Oligodendrocytes_df$Nr_OligodendrocytesMarkers_in_DCGs=as.numeric(res.earlyAD_DCG_Oligodendrocytes_df$Nr_OligodendrocytesMarkers_in_DCGs)
res.earlyAD_DCG_Oligodendrocytes_df=res.earlyAD_DCG_Oligodendrocytes_df%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(Nr_OligodendrocytesMarkers_in_DCGs>Expected_OligodendrocytesMarkers_in_DCGs,true = "Over",false = "Under"))%>%filter(Representation=="Over")



#Redhead et al. viral network drivers, loss-gain genes
redhead_nnet_driver_genes.df=read.xls("/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/Redhead_AD_HerpesVirus_Dataset/mmc2.xlsx",sheet = 3,header=T)
redhead_nnet_driver_genes.list=earlyAD_DCG_nnet_driver_genes.overlap=earlyAD_DCG_nnet_driver_genes.expected=earlyAD_DCG_nnet_driver_genes.fisher=vector(mode = "list",length = 3)
names(redhead_nnet_driver_genes.list)=names(earlyAD_DCG_nnet_driver_genes.overlap)=names(earlyAD_DCG_nnet_driver_genes.expected)=names(earlyAD_DCG_nnet_driver_genes.fisher)=unique(redhead_nnet_driver_genes.df$driver_type)

for(i in 1:3){
  redhead_nnet_driver_genes.list[[i]]=redhead_nnet_driver_genes.df%>%filter(driver_type==names(redhead_nnet_driver_genes.list)[i])%>%pull(driver_symbol)%>%droplevels%>%levels
  earlyAD_DCG_nnet_driver_genes.overlap[[i]]=lapply(earlyAD_DCGs.list,intersect,redhead_nnet_driver_genes.list[[i]])  
  earlyAD_DCG_nnet_driver_genes.expected[[i]]=lapply(lapply(earlyAD_DCGs.list,length),function(x)round(length(redhead_nnet_driver_genes.list[[i]])/19530*x,digits = 1))
  earlyAD_DCG_nnet_driver_genes.fisher[[i]]=vector(mode = "list",length = 12)
  names(earlyAD_DCG_nnet_driver_genes.fisher[[i]])=names(earlyAD_DCGs.list)
  for(m in 1:length(earlyAD_DCGs.list)){
    earlyAD_DCG_nnet_driver_genes.fisher[[i]][[m]]=fisher.test(matrix(c(length(earlyAD_DCG_nnet_driver_genes.overlap[[i]][[m]]),
                                                                        length(redhead_nnet_driver_genes.list[[i]])-length(earlyAD_DCG_nnet_driver_genes.overlap[[i]][[m]]),
                                                                        lapply(earlyAD_DCGs.list,length)[[m]]-length(earlyAD_DCG_nnet_driver_genes.overlap[[i]][[m]]),
                                                                        19530-(length(redhead_nnet_driver_genes.list[[i]])-length(earlyAD_DCG_nnet_driver_genes.overlap[[i]][[m]]))-
                                                                        lapply(earlyAD_DCGs.list,length)[[m]]-length(earlyAD_DCG_nnet_driver_genes.overlap[[i]][[m]])),nrow = 2))}
}
  

  earlyAD_DCG_nnet_driver_genes.fisher[[i]]=vector(mode = "list",length = length(earlyAD_DCGs.list))
  earlyAD_DCG_nnet_driver_genes.fisher[[i]]=







