library(org.Hs.eg.db)
library(igraph)
library(data.table)
library(centiserve)
library(sets)
library(clusterProfiler)
library(dplyr)
library(gdata)
library(foreach)
library(doParallel)

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
setwd("/Users/sandeepamberkar/Work/Data/AMP-AD_RNAseq_ReSeq/Normalised_covariate_corrected_NoResiduals/")

ConsensusDEGs=vector(mode = "list",length = 7)
ConsensusDEGs_df=vector(mode = "list",length = 7)
names(ConsensusDEGs)=names(ConsensusDEGs_df)=c("DLPFC","FP","IFG","PHG","STG","CER","TCX")
ConsensusDEGs_df$DLPFC=fread("../DEG_Analyses/ROSMAP_DLPFC_DiffExpression.tsv",header = T,stringsAsFactors = F,data.table = F) %>%filter(Comparison=="cogdx1-cogdx4")
ConsensusDEGs_df$CER=fread("../DEG_Analyses/MAYO_CBE_TCX_DiffExpression.tsv",header = T,stringsAsFactors = F,data.table = F) %>%filter(BrainRegion=="CER"&Comparison=="CONTROL-AD")
ConsensusDEGs_df$TCX=fread("../DEG_Analyses/MAYO_CBE_TCX_DiffExpression.tsv",header = T,stringsAsFactors = F,data.table = F) %>%filter(BrainRegion=="TCX"&Comparison=="CONTROL-AD")
ConsensusDEGs_df$FP=fread("../DEG_Analyses/MSSM_FP_STG_PHG_IFG_DiffExpression.tsv",header = T,stringsAsFactors = F,data.table = F)%>%filter(BrainRegion=="FP"&Comparison=="CONTROL-AD")
ConsensusDEGs_df$IFG=fread("../DEG_Analyses/MSSM_FP_STG_PHG_IFG_DiffExpression.tsv",header = T,stringsAsFactors = F,data.table = F)%>%filter(BrainRegion=="IFG"&Comparison=="CONTROL-AD")
ConsensusDEGs_df$PHG=fread("../DEG_Analyses/MSSM_FP_STG_PHG_IFG_DiffExpression.tsv",header = T,stringsAsFactors = F,data.table = F)%>%filter(BrainRegion=="PHG"&Comparison=="CONTROL-AD")
ConsensusDEGs_df$STG=fread("../DEG_Analyses/MSSM_FP_STG_PHG_IFG_DiffExpression.tsv",header = T,stringsAsFactors = F,data.table = F)%>%filter(BrainRegion=="STG"&Comparison=="CONTROL-AD")
ConsensusDEGs=lapply(ConsensusDEGs_df,function(x)x%>%filter(Direction=="UP"|Direction=="DOWN")%>%filter(hgnc_symbol!="")%>%select(hgnc_symbol))
tanzi_damaging_variants=scan("../../../Collaborations/Tanzi_WGS/Damaging_exonic_GeneSymbol.txt",sep = "\n",what = "char")
tanzi_protective_variants=scan("../../../Collaborations/Tanzi_WGS/Protective_exonic_GeneSymbol.txt",sep = "\n",what = "char")

metaAD_DCLs=metaAD_DCGs=metaAD_DRLs=metaAD_DRGs=metaAD_DRsort_results=vector(mode = "list",length=7)
names(metaAD_DCLs)=names(metaAD_DCGs)=names(metaAD_DRLs)=names(metaAD_DRGs)=names(metaAD_DRsort_results)=c("Mayo_TCX","Mayo_CER","ROSMAP_DLPFC","MSBB_FP","MSBB_IFG","MSBB_PHG","MSBB_STG")

metaAD_DRsort_results$Mayo_TCX=readRDS("MAYO/DCGL_results/TCX_DiffCoexp_DRsort_res.RDS")
metaAD_DRsort_results$Mayo_CER=readRDS("MAYO/DCGL_results/CER_DRsort_noLFC.RDS")
metaAD_DRsort_results$ROSMAP_DLPFC=readRDS("ROSMAP/DCGL_results/ROSMAP_DiffCoexp_DRsort_noLFC.res.RDS")
metaAD_DRsort_results$MSBB_FP=readRDS("MSMM/DCGL_results/MSBB_FP_DRsort_noLFC.res.RDS")
metaAD_DRsort_results$MSBB_IFG=readRDS("MSMM/DCGL_results/MSBB_IFG_DRsort_noLFC.res.RDS")
metaAD_DRsort_results$MSBB_PHG=readRDS("MSMM/DCGL_results/MSBB_PHG_DRsort_noLFC.res.RDS")
metaAD_DRsort_results$MSBB_STG=readRDS("MSMM/DCGL_results/MSBB_STG_DRsort_noLFC.res.RDS")

metaAD_DCLs=lapply(metaAD_DRsort_results,function(x)x$DCLs)
metaAD_DCLs.graph=lapply(metaAD_DCLs,function(x)data.frame(x)%>%filter(q.diffcor<0.05)%>%select(c(Gene.1,Gene.2))%>%graph.data.frame(directed = F))
metaAD_DCGs=lapply(metaAD_DRsort_results,function(x)data.frame(x$DCGs,stringsAsFactors = F))
metaAD_DCGs.filtered=lapply(lapply(metaAD_DCGs,function(x)x%>%filter(q<0.05)%>%select(DCG)),`[[`,1)
metaAD_DRLs=lapply(metaAD_DRsort_results,function(x)x$DRLs)
metaAD_DRGs=lapply(metaAD_DRsort_results,function(x)x$DRGs)
metaAD_TF_bridged_DCLs=lapply(metaAD_DRsort_results,function(x)x$TF_bridged_DCL)


ampAD_AggModules.data=fread("../AggregateModules/AMP-AD_AggModules.tsv",sep = "\t",header = T,stringsAsFactors = F,showProgress = T,data.table = F)
tmp=lapply(unique(ampAD_AggModules.data$brainRegion),function(x)ampAD_AggModules.data%>%filter(brainRegion==x)%>%select(c(Module,external_gene_name)))
ampAD_AggModules_data_brainRegion=ampAD_AggModules_data_brainRegion_OlpGenes=vector(mode="list",length=7)
names(ampAD_AggModules_data_brainRegion)=names(ampAD_AggModules_data_brainRegion_OlpGenes)=unique(ampAD_AggModules.data$brainRegion)

for(r in 1:length(ampAD_AggModules_data_brainRegion)){
  ampAD_AggModules_data_brainRegion[[r]]=lapply(lapply(unique(tmp[[r]]$Module),function(x)tmp[[r]]%>%filter(Module==x)%>%select(external_gene_name)),`[[`,1)
  names(ampAD_AggModules_data_brainRegion[[r]])=unique(tmp[[r]]$Module)
}

comp_ConsensusModules_damaging_DLPFC.res=foreach(i=1:length(ConsensusUniqueModules$DLPFC))%do%{
  jcd=lapply(ConsensusUniqueModules$DLPFC[i],jaccard,tanzi_damaging_variants)
  module_size=unlist(lapply(ConsensusUniqueModules$DLPFC[i],length))
  module_match=lapply(lapply(ConsensusUniqueModules$DLPFC[i],intersect,tanzi_damaging_variants),length)>1
  overlapping_genes=lapply(lapply(ConsensusUniqueModules$DLPFC[i],intersect,tanzi_damaging_variants)[module_match],function(y)mapIds2(IDs = y,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2])
  olp_KEGG=lapply(overlapping_genes,function(k)data.frame(enrichKEGG(gene = k,organism = "hsa",pvalueCutoff = 0.05,pAdjustMethod = "BH")))
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T]=lapply(olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T],`[[`,2)
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==F]="No enriched KEGG pathways found!"
  pval_vector=unlist(lapply(overlapping_genes[i],function(x)phyper(length(intersect(x,tanzi_damaging_variants)),length(x),19000,length(tanzi_damaging_variants),lower.tail = F)))
  data.frame(ModuleName=names(ConsensusUniqueModules$DLPFC[i][module_match]),
             Jaccard=unlist(jcd[module_match]),
             ModuleSize=module_size[module_match],
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes[module_match],length)),
             Nr_Of_OlpKEGG=unlist(lapply(olp_KEGG,paste,sep=";",collapse=";")),
             p.overlap=pval_vector[module_match],stringsAsFactors = F)
  
}
comp_ConsensusModules_damaging_DLPFC.df=data.frame(rbindlist(comp_ConsensusModules_damaging_DLPFC.res[lapply(comp_ConsensusModules_damaging_DLPFC.res,function(x)dim(x)[1])>0]),stringsAsFactors = F)%>%arrange(desc(-p.overlap))
write.table(comp_ConsensusModules_damaging_DLPFC.df,"../../../Collaborations/Tanzi_WGS/ConsensusModules_damaging_DLPFC_comparison.txt",col.names = T,row.names = F,sep = "\t",quote = F)

comp_ConsensusModules_damaging_TCX.res=foreach(i=1:length(ConsensusUniqueModules$TCX))%do%{
  jcd=lapply(ConsensusUniqueModules$TCX[i],jaccard,tanzi_damaging_variants)
  module_size=unlist(lapply(ConsensusUniqueModules$TCX[i],length))
  module_match=lapply(lapply(ConsensusUniqueModules$TCX[i],intersect,tanzi_damaging_variants),length)>1
  overlapping_genes=lapply(lapply(ConsensusUniqueModules$TCX[i],intersect,tanzi_damaging_variants)[module_match],function(y)mapIds2(IDs = y,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2])
  olp_KEGG=lapply(overlapping_genes,function(k)data.frame(enrichKEGG(gene = k,organism = "hsa",pvalueCutoff = 0.05,pAdjustMethod = "BH")))
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T]=lapply(olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T],`[[`,2)
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==F]="No enriched KEGG pathways found!"
  pval_vector=unlist(lapply(overlapping_genes[i],function(x)phyper(length(intersect(x,tanzi_damaging_variants)),length(x),19000,length(tanzi_damaging_variants),lower.tail = F)))
  data.frame(ModuleName=names(ConsensusUniqueModules$TCX[i][module_match]),
            Jaccard=unlist(jcd[module_match]),
             ModuleSize=module_size[module_match],
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes[module_match],length)),
             Nr_Of_OlpKEGG=unlist(lapply(olp_KEGG,paste,sep=";",collapse=";")),
             p.overlap=pval_vector[module_match],stringsAsFactors = F)
  
}
comp_ConsensusModules_damaging_TCX.df=data.frame(rbindlist(comp_ConsensusModules_damaging_TCX.res[lapply(comp_ConsensusModules_damaging_TCX.res,function(x)dim(x)[1])>0]),stringsAsFactors = F)%>%arrange(desc(-p.overlap))
write.table(comp_ConsensusModules_damaging_DLPFC.df,"../../../Collaborations/Tanzi_WGS/ConsensusModules_damaging_DLPFC_comparison.txt",col.names = T,row.names = F,sep = "\t",quote = F)

comp_ConsensusModules_damaging_CER.res=foreach(i=1:length(ConsensusUniqueModules$CER))%do%{
  jcd=lapply(ConsensusUniqueModules$CER[i],jaccard,tanzi_damaging_variants)
  module_size=unlist(lapply(ConsensusUniqueModules$CER[i],length))
  module_match=lapply(lapply(ConsensusUniqueModules$CER[i],intersect,tanzi_damaging_variants),length)>1
  overlapping_genes=lapply(lapply(ConsensusUniqueModules$CER[i],intersect,tanzi_damaging_variants)[module_match],function(y)mapIds2(IDs = y,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2])
  olp_KEGG=lapply(overlapping_genes,function(k)data.frame(enrichKEGG(gene = k,organism = "hsa",pvalueCutoff = 0.05,pAdjustMethod = "BH")))
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T]=lapply(olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T],`[[`,2)
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==F]="No enriched KEGG pathways found!"
  pval_vector=unlist(lapply(overlapping_genes[i],function(x)phyper(length(intersect(x,tanzi_damaging_variants)),length(x),19000,length(tanzi_damaging_variants),lower.tail = F)))
  data.frame(ModuleName=names(ConsensusUniqueModules$CER[i][module_match]),
             Jaccard=unlist(jcd[module_match]),
             ModuleSize=module_size[module_match],
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes[module_match],length)),
             Nr_Of_OlpKEGG=unlist(lapply(olp_KEGG,paste,sep=";",collapse=";")),
             p.overlap=pval_vector[module_match],stringsAsFactors = F)
  
}
comp_ConsensusModules_damaging_CER.df=data.frame(rbindlist(comp_ConsensusModules_damaging_CER.res[lapply(comp_ConsensusModules_damaging_CER.res,function(x)dim(x)[1])>0]),stringsAsFactors = F)%>%arrange(desc(-p.overlap))
write.table(comp_ConsensusModules_damaging_DLPFC.df,"../../../Collaborations/Tanzi_WGS/ConsensusModules_damaging_DLPFC_comparison.txt",col.names = T,row.names = F,sep = "\t",quote = F)

comp_ConsensusModules_damaging_FP.res=foreach(i=1:length(ConsensusUniqueModules$FP))%do%{
  jcd=lapply(ConsensusUniqueModules$FP[i],jaccard,tanzi_damaging_variants)
  module_size=unlist(lapply(ConsensusUniqueModules$FP[i],length))
  module_match=lapply(lapply(ConsensusUniqueModules$FP[i],intersect,tanzi_damaging_variants),length)>1
  overlapping_genes=lapply(lapply(ConsensusUniqueModules$FP[i],intersect,tanzi_damaging_variants)[module_match],function(y)mapIds2(IDs = y,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2])
  olp_KEGG=lapply(overlapping_genes,function(k)data.frame(enrichKEGG(gene = k,organism = "hsa",pvalueCutoff = 0.05,pAdjustMethod = "BH")))
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T]=lapply(olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T],`[[`,2)
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==F]="No enriched KEGG pathways found!"
  pval_vector=unlist(lapply(overlapping_genes[i],function(x)phyper(length(intersect(x,tanzi_damaging_variants)),length(x),19000,length(tanzi_damaging_variants),lower.tail = F)))
  data.frame(ModuleName=names(ConsensusUniqueModules$FP[i][module_match]),
             Jaccard=unlist(jcd[module_match]),
             ModuleSize=module_size[module_match],
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes[module_match],length)),
             Nr_Of_OlpKEGG=unlist(lapply(olp_KEGG,paste,sep=";",collapse=";")),
             p.overlap=pval_vector[module_match],stringsAsFactors = F)
  
}
comp_ConsensusModules_damaging_FP.df=data.frame(rbindlist(comp_ConsensusModules_damaging_FP.res[lapply(comp_ConsensusModules_damaging_FP.res,function(x)dim(x)[1])>0]),stringsAsFactors = F)%>%arrange(desc(-p.overlap))
write.table(comp_ConsensusModules_damaging_FP.df,"../../../Collaborations/Tanzi_WGS/ConsensusModules_damaging_FP_comparison.txt",col.names = T,row.names = F,sep = "\t",quote = F)

comp_ConsensusModules_damaging_IFG.res=foreach(i=1:length(ConsensusUniqueModules$IFG))%do%{
  jcd=lapply(ConsensusUniqueModules$IFG[i],jaccard,tanzi_damaging_variants)
  module_size=unlist(lapply(ConsensusUniqueModules$IFG[i],length))
  module_match=lapply(lapply(ConsensusUniqueModules$IFG[i],intersect,tanzi_damaging_variants),length)>1
  overlapping_genes=lapply(lapply(ConsensusUniqueModules$IFG[i],intersect,tanzi_damaging_variants)[module_match],function(y)mapIds2(IDs = y,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2])
  olp_KEGG=lapply(overlapping_genes,function(k)data.frame(enrichKEGG(gene = k,organism = "hsa",pvalueCutoff = 0.05,pAdjustMethod = "BH")))
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T]=lapply(olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T],`[[`,2)
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==F]="No enriched KEGG pathways found!"
  pval_vector=unlist(lapply(overlapping_genes[i],function(x)phyper(length(intersect(x,tanzi_damaging_variants)),length(x),19000,length(tanzi_damaging_variants),lower.tail = F)))
  data.frame(ModuleName=names(ConsensusUniqueModules$IFG[i][module_match]),
             Jaccard=unlist(jcd[module_match]),
             ModuleSize=module_size[module_match],
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes[module_match],length)),
             Nr_Of_OlpKEGG=unlist(lapply(olp_KEGG,paste,sep=";",collapse=";")),
             p.overlap=pval_vector[module_match],stringsAsFactors = F)
  
}
comp_ConsensusModules_damaging_IFG.df=data.frame(rbindlist(comp_ConsensusModules_damaging_IFG.res[lapply(comp_ConsensusModules_damaging_IFG.res,function(x)dim(x)[1])>0]),stringsAsFactors = F)%>%arrange(desc(-p.overlap))
write.table(comp_ConsensusModules_damaging_IFG.df,"../../../Collaborations/Tanzi_WGS/ConsensusModules_damaging_IFG_comparison.txt",col.names = T,row.names = F,sep = "\t",quote = F)

comp_ConsensusModules_damaging_PHG.res=foreach(i=1:length(ConsensusUniqueModules$PHG))%do%{
  jcd=lapply(ConsensusUniqueModules$PHG[i],jaccard,tanzi_damaging_variants)
  module_size=unlist(lapply(ConsensusUniqueModules$PHG[i],length))
  module_match=lapply(lapply(ConsensusUniqueModules$PHG[i],intersect,tanzi_damaging_variants),length)>1
  overlapping_genes=lapply(lapply(ConsensusUniqueModules$PHG[i],intersect,tanzi_damaging_variants)[module_match],function(y)mapIds2(IDs = y,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2])
  olp_KEGG=lapply(overlapping_genes,function(k)data.frame(enrichKEGG(gene = k,organism = "hsa",pvalueCutoff = 0.05,pAdjustMethod = "BH")))
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T]=lapply(olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T],`[[`,2)
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==F]="No enriched KEGG pathways found!"
  pval_vector=unlist(lapply(overlapping_genes[i],function(x)phyper(length(intersect(x,tanzi_damaging_variants)),length(x),19000,length(tanzi_damaging_variants),lower.tail = F)))
  data.frame(ModuleName=names(ConsensusUniqueModules$PHG[i][module_match]),
             Jaccard=unlist(jcd[module_match]),
             ModuleSize=module_size[module_match],
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes[module_match],length)),
             Nr_Of_OlpKEGG=unlist(lapply(olp_KEGG,paste,sep=";",collapse=";")),
             p.overlap=pval_vector[module_match],stringsAsFactors = F)
  
}
comp_ConsensusModules_damaging_PHG.df=data.frame(rbindlist(comp_ConsensusModules_damaging_PHG.res[lapply(comp_ConsensusModules_damaging_PHG.res,function(x)dim(x)[1])>0]),stringsAsFactors = F)%>%arrange(desc(-p.overlap))
write.table(comp_ConsensusModules_damaging_PHG.df,"../../../Collaborations/Tanzi_WGS/ConsensusModules_damaging_PHG_comparison.txt",col.names = T,row.names = F,sep = "\t",quote = F)

comp_ConsensusModules_damaging_STG.res=foreach(i=1:length(ConsensusUniqueModules$STG))%do%{
  jcd=lapply(ConsensusUniqueModules$STG[i],jaccard,tanzi_damaging_variants)
  module_size=unlist(lapply(ConsensusUniqueModules$STG[i],length))
  module_match=lapply(lapply(ConsensusUniqueModules$STG[i],intersect,tanzi_damaging_variants),length)>1
  overlapping_genes=lapply(lapply(ConsensusUniqueModules$STG[i],intersect,tanzi_damaging_variants)[module_match],function(y)mapIds2(IDs = y,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2])
  olp_KEGG=lapply(overlapping_genes,function(k)data.frame(enrichKEGG(gene = k,organism = "hsa",pvalueCutoff = 0.05,pAdjustMethod = "BH")))
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T]=lapply(olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T],`[[`,2)
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==F]="No enriched KEGG pathways found!"
  pval_vector=unlist(lapply(overlapping_genes[i],function(x)phyper(length(intersect(x,tanzi_damaging_variants)),length(x),19000,length(tanzi_damaging_variants),lower.tail = F)))
  data.frame(ModuleName=names(ConsensusUniqueModules$STG[i][module_match]),
             Jaccard=unlist(jcd[module_match]),
             ModuleSize=module_size[module_match],
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes[module_match],length)),
             Nr_Of_OlpKEGG=unlist(lapply(olp_KEGG,paste,sep=";",collapse=";")),
             p.overlap=pval_vector[module_match],stringsAsFactors = F)
  
}
comp_ConsensusModules_damaging_STG.df=data.frame(rbindlist(comp_ConsensusModules_damaging_STG.res[lapply(comp_ConsensusModules_damaging_STG.res,function(x)dim(x)[1])>0]),stringsAsFactors = F)%>%arrange(desc(-p.overlap))
write.table(comp_ConsensusModules_damaging_STG.df,"../../../Collaborations/Tanzi_WGS/ConsensusModules_damaging_STG_comparison.txt",col.names = T,row.names = F,sep = "\t",quote = F)

#Compare with Aggregate modules
comp_AggModules_damaging_DLPFC.res=foreach(i=1:length(ampAD_AggModules_data_brainRegion$DLPFC))%do%{
  jcd=lapply(ampAD_AggModules_data_brainRegion$DLPFC[i],jaccard,tanzi_damaging_variants)
  module_size=unlist(lapply(ampAD_AggModules_data_brainRegion$DLPFC[i],length))
  module_match=lapply(lapply(ampAD_AggModules_data_brainRegion$DLPFC[i],intersect,tanzi_damaging_variants),length)>1
  overlapping_genes=lapply(lapply(ampAD_AggModules_data_brainRegion$DLPFC[i],intersect,tanzi_damaging_variants)[module_match],function(y)mapIds2(IDs = y,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2])
  olp_KEGG=lapply(overlapping_genes,function(k)data.frame(enrichKEGG(gene = k,organism = "hsa",pvalueCutoff = 0.05,pAdjustMethod = "BH")))
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T]=lapply(olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T],`[[`,2)
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==F]="No enriched KEGG pathways found!"
  pval_vector=unlist(lapply(overlapping_genes[i],function(x)phyper(length(intersect(x,tanzi_damaging_variants)),length(x),19000,length(tanzi_damaging_variants),lower.tail = F)))
  data.frame(ModuleName=names(ampAD_AggModules_data_brainRegion$DLPFC[i][module_match]),
             Jaccard=unlist(jcd[module_match]),
             ModuleSize=module_size[module_match],
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes[module_match],length)),
             Nr_Of_OlpKEGG=unlist(lapply(olp_KEGG,paste,sep=";",collapse=";")),
             p.overlap=pval_vector[module_match],stringsAsFactors = F)
  
}
comp_AggModules_damaging_DLPFC.df=data.frame(rbindlist(comp_AggModules_damaging_DLPFC.res[lapply(comp_AggModules_damaging_DLPFC.res,function(x)dim(x)[1])>0]),stringsAsFactors = F)%>%arrange(desc(-p.overlap))
write.table(comp_AggModules_damaging_DLPFC.df,"../../../Collaborations/Tanzi_WGS/AggModules_damaging_DLPFC_comparison.txt",col.names = T,row.names = F,sep = "\t",quote = F)

comp_AggModules_damaging_CER.res=foreach(i=1:length(ampAD_AggModules_data_brainRegion$CBE))%do%{
  jcd=lapply(ampAD_AggModules_data_brainRegion$CBE[i],jaccard,tanzi_damaging_variants)
  module_size=unlist(lapply(ampAD_AggModules_data_brainRegion$CBE[i],length))
  module_match=lapply(lapply(ampAD_AggModules_data_brainRegion$CBE[i],intersect,tanzi_damaging_variants),length)>1
  overlapping_genes=lapply(lapply(ampAD_AggModules_data_brainRegion$CBE[i],intersect,tanzi_damaging_variants)[module_match],function(y)mapIds2(IDs = y,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2])
  olp_KEGG=lapply(overlapping_genes,function(k)data.frame(enrichKEGG(gene = k,organism = "hsa",pvalueCutoff = 0.05,pAdjustMethod = "BH")))
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T]=lapply(olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T],`[[`,2)
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==F]="No enriched KEGG pathways found!"
  pval_vector=unlist(lapply(overlapping_genes[i],function(x)phyper(length(intersect(x,tanzi_damaging_variants)),length(x),19000,length(tanzi_damaging_variants),lower.tail = F)))
  data.frame(ModuleName=names(ampAD_AggModules_data_brainRegion$CBE[i][module_match]),
             Jaccard=unlist(jcd[module_match]),
             ModuleSize=module_size[module_match],
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes[module_match],length)),
             Nr_Of_OlpKEGG=unlist(lapply(olp_KEGG,paste,sep=";",collapse=";")),
             p.overlap=pval_vector[module_match],stringsAsFactors = F)
  
}
comp_AggModules_damaging_CER.df=data.frame(rbindlist(comp_AggModules_damaging_CER.res[lapply(comp_AggModules_damaging_CER.res,function(x)dim(x)[1])>0]),stringsAsFactors = F)%>%arrange(desc(-p.overlap))
write.table(comp_AggModules_damaging_CER.df,"../../../Collaborations/Tanzi_WGS/AggModules_damaging_CER_comparison.txt",col.names = T,row.names = F,sep = "\t",quote = F)

comp_AggModules_damaging_TCX.res=foreach(i=1:length(ampAD_AggModules_data_brainRegion$TCX))%do%{
  jcd=lapply(ampAD_AggModules_data_brainRegion$TCX[i],jaccard,tanzi_damaging_variants)
  module_size=unlist(lapply(ampAD_AggModules_data_brainRegion$TCX[i],length))
  module_match=lapply(lapply(ampAD_AggModules_data_brainRegion$TCX[i],intersect,tanzi_damaging_variants),length)>1
  overlapping_genes=lapply(lapply(ampAD_AggModules_data_brainRegion$TCX[i],intersect,tanzi_damaging_variants)[module_match],function(y)mapIds2(IDs = y,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2])
  olp_KEGG=lapply(overlapping_genes,function(k)data.frame(enrichKEGG(gene = k,organism = "hsa",pvalueCutoff = 0.05,pAdjustMethod = "BH")))
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T]=lapply(olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T],`[[`,2)
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==F]="No enriched KEGG pathways found!"
  pval_vector=unlist(lapply(overlapping_genes[i],function(x)phyper(length(intersect(x,tanzi_damaging_variants)),length(x),19000,length(tanzi_damaging_variants),lower.tail = F)))
  data.frame(ModuleName=names(ampAD_AggModules_data_brainRegion$TCX[i][module_match]),
             Jaccard=unlist(jcd[module_match]),
             ModuleSize=module_size[module_match],
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes[module_match],length)),
             Nr_Of_OlpKEGG=unlist(lapply(olp_KEGG,paste,sep=";",collapse=";")),
             p.overlap=pval_vector[module_match],stringsAsFactors = F)
  
}
comp_AggModules_damaging_TCX.df=data.frame(rbindlist(comp_AggModules_damaging_TCX.res[lapply(comp_AggModules_damaging_TCX.res,function(x)dim(x)[1])>0]),stringsAsFactors = F)%>%arrange(desc(-p.overlap))
write.table(comp_AggModules_damaging_TCX.df,"../../../Collaborations/Tanzi_WGS/AggModules_damaging_TCX_comparison.txt",col.names = T,row.names = F,sep = "\t",quote = F)

comp_AggModules_damaging_FP.res=foreach(i=1:length(ampAD_AggModules_data_brainRegion$FP))%do%{
  jcd=lapply(ampAD_AggModules_data_brainRegion$FP[i],jaccard,tanzi_damaging_variants)
  module_size=unlist(lapply(ampAD_AggModules_data_brainRegion$FP[i],length))
  module_match=lapply(lapply(ampAD_AggModules_data_brainRegion$FP[i],intersect,tanzi_damaging_variants),length)>1
  overlapping_genes=lapply(lapply(ampAD_AggModules_data_brainRegion$FP[i],intersect,tanzi_damaging_variants)[module_match],function(y)mapIds2(IDs = y,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2])
  olp_KEGG=lapply(overlapping_genes,function(k)data.frame(enrichKEGG(gene = k,organism = "hsa",pvalueCutoff = 0.05,pAdjustMethod = "BH")))
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T]=lapply(olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T],`[[`,2)
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==F]="No enriched KEGG pathways found!"
  pval_vector=unlist(lapply(overlapping_genes[i],function(x)phyper(length(intersect(x,tanzi_damaging_variants)),length(x),19000,length(tanzi_damaging_variants),lower.tail = F)))
  data.frame(ModuleName=names(ampAD_AggModules_data_brainRegion$FP[i][module_match]),
             Jaccard=unlist(jcd[module_match]),
             ModuleSize=module_size[module_match],
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes[module_match],length)),
             Nr_Of_OlpKEGG=unlist(lapply(olp_KEGG,paste,sep=";",collapse=";")),
             p.overlap=pval_vector[module_match],stringsAsFactors = F)
  
}
comp_AggModules_damaging_FP.df=data.frame(rbindlist(comp_AggModules_damaging_FP.res[lapply(comp_AggModules_damaging_FP.res,function(x)dim(x)[1])>0]),stringsAsFactors = F)%>%arrange(desc(-p.overlap))
write.table(comp_AggModules_damaging_FP.df,"../../../Collaborations/Tanzi_WGS/AggModules_damaging_FP_comparison.txt",col.names = T,row.names = F,sep = "\t",quote = F)

comp_AggModules_damaging_IFG.res=foreach(i=1:length(ampAD_AggModules_data_brainRegion$IFG))%do%{
  jcd=lapply(ampAD_AggModules_data_brainRegion$IFG[i],jaccard,tanzi_damaging_variants)
  module_size=unlist(lapply(ampAD_AggModules_data_brainRegion$IFG[i],length))
  module_match=lapply(lapply(ampAD_AggModules_data_brainRegion$IFG[i],intersect,tanzi_damaging_variants),length)>1
  overlapping_genes=lapply(lapply(ampAD_AggModules_data_brainRegion$IFG[i],intersect,tanzi_damaging_variants)[module_match],function(y)mapIds2(IDs = y,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2])
  olp_KEGG=lapply(overlapping_genes,function(k)data.frame(enrichKEGG(gene = k,organism = "hsa",pvalueCutoff = 0.05,pAdjustMethod = "BH")))
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T]=lapply(olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T],`[[`,2)
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==F]="No enriched KEGG pathways found!"
  pval_vector=unlist(lapply(overlapping_genes[i],function(x)phyper(length(intersect(x,tanzi_damaging_variants)),length(x),19000,length(tanzi_damaging_variants),lower.tail = F)))
  data.frame(ModuleName=names(ampAD_AggModules_data_brainRegion$IFG[i][module_match]),
             Jaccard=unlist(jcd[module_match]),
             ModuleSize=module_size[module_match],
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes[module_match],length)),
             Nr_Of_OlpKEGG=unlist(lapply(olp_KEGG,paste,sep=";",collapse=";")),
             p.overlap=pval_vector[module_match],stringsAsFactors = F)
  
}
comp_AggModules_damaging_IFG.df=data.frame(rbindlist(comp_AggModules_damaging_IFG.res[lapply(comp_AggModules_damaging_IFG.res,function(x)dim(x)[1])>0]),stringsAsFactors = F)%>%arrange(desc(-p.overlap))
write.table(comp_AggModules_damaging_IFG.df,"../../../Collaborations/Tanzi_WGS/AggModules_damaging_IFG_comparison.txt",col.names = T,row.names = F,sep = "\t",quote = F)

comp_AggModules_damaging_PHG.res=foreach(i=1:length(ampAD_AggModules_data_brainRegion$PHG))%do%{
  jcd=lapply(ampAD_AggModules_data_brainRegion$PHG[i],jaccard,tanzi_damaging_variants)
  module_size=unlist(lapply(ampAD_AggModules_data_brainRegion$PHG[i],length))
  module_match=lapply(lapply(ampAD_AggModules_data_brainRegion$PHG[i],intersect,tanzi_damaging_variants),length)>1
  overlapping_genes=lapply(lapply(ampAD_AggModules_data_brainRegion$PHG[i],intersect,tanzi_damaging_variants)[module_match],function(y)mapIds2(IDs = y,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2])
  olp_KEGG=lapply(overlapping_genes,function(k)data.frame(enrichKEGG(gene = k,organism = "hsa",pvalueCutoff = 0.05,pAdjustMethod = "BH")))
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T]=lapply(olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T],`[[`,2)
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==F]="No enriched KEGG pathways found!"
  pval_vector=unlist(lapply(overlapping_genes[i],function(x)phyper(length(intersect(x,tanzi_damaging_variants)),length(x),19000,length(tanzi_damaging_variants),lower.tail = F)))
  data.frame(ModuleName=names(ampAD_AggModules_data_brainRegion$PHG[i][module_match]),
             Jaccard=unlist(jcd[module_match]),
             ModuleSize=module_size[module_match],
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes[module_match],length)),
             Nr_Of_OlpKEGG=unlist(lapply(olp_KEGG,paste,sep=";",collapse=";")),
             p.overlap=pval_vector[module_match],stringsAsFactors = F)
  
}
comp_AggModules_damaging_PHG.df=data.frame(rbindlist(comp_AggModules_damaging_PHG.res[lapply(comp_AggModules_damaging_PHG.res,function(x)dim(x)[1])>0]),stringsAsFactors = F)%>%arrange(desc(-p.overlap))
write.table(comp_AggModules_damaging_PHG.df,"../../../Collaborations/Tanzi_WGS/AggModules_damaging_PHG_comparison.txt",col.names = T,row.names = F,sep = "\t",quote = F)

comp_AggModules_damaging_STG.res=foreach(i=1:length(ampAD_AggModules_data_brainRegion$STG))%do%{
  jcd=lapply(ampAD_AggModules_data_brainRegion$STG[i],jaccard,tanzi_damaging_variants)
  module_size=unlist(lapply(ampAD_AggModules_data_brainRegion$STG[i],length))
  module_match=lapply(lapply(ampAD_AggModules_data_brainRegion$STG[i],intersect,tanzi_damaging_variants),length)>1
  overlapping_genes=lapply(lapply(ampAD_AggModules_data_brainRegion$STG[i],intersect,tanzi_damaging_variants)[module_match],function(y)mapIds2(IDs = y,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2])
  olp_KEGG=lapply(overlapping_genes,function(k)data.frame(enrichKEGG(gene = k,organism = "hsa",pvalueCutoff = 0.05,pAdjustMethod = "BH")))
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T]=lapply(olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T],`[[`,2)
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==F]="No enriched KEGG pathways found!"
  pval_vector=unlist(lapply(overlapping_genes[i],function(x)phyper(length(intersect(x,tanzi_damaging_variants)),length(x),19000,length(tanzi_damaging_variants),lower.tail = F)))
  data.frame(ModuleName=names(ampAD_AggModules_data_brainRegion$STG[i][module_match]),
             Jaccard=unlist(jcd[module_match]),
             ModuleSize=module_size[module_match],
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes[module_match],length)),
             Nr_Of_OlpKEGG=unlist(lapply(olp_KEGG,paste,sep=";",collapse=";")),
             p.overlap=pval_vector[module_match],stringsAsFactors = F)
  
}
comp_AggModules_damaging_STG.df=data.frame(rbindlist(comp_AggModules_damaging_STG.res[lapply(comp_AggModules_damaging_STG.res,function(x)dim(x)[1])>0]),stringsAsFactors = F)%>%arrange(desc(-p.overlap))
write.table(comp_AggModules_damaging_STG.df,"../../../Collaborations/Tanzi_WGS/AggModules_damaging_STG_comparison.txt",col.names = T,row.names = F,sep = "\t",quote = F)

########################################################################################################################
# All comparisons for protective variants #
########################################################################################################################
comp_ConsensusModules_protective_DLPFC.res=foreach(i=1:length(ConsensusUniqueModules$DLPFC))%do%{
  jcd=lapply(ConsensusUniqueModules$DLPFC[i],jaccard,tanzi_protective_variants)
  module_size=unlist(lapply(ConsensusUniqueModules$DLPFC[i],length))
  module_match=lapply(lapply(ConsensusUniqueModules$DLPFC[i],intersect,tanzi_protective_variants),length)>1
  overlapping_genes=lapply(lapply(ConsensusUniqueModules$DLPFC[i],intersect,tanzi_protective_variants)[module_match],function(y)mapIds2(IDs = y,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2])
  olp_KEGG=lapply(overlapping_genes,function(k)data.frame(enrichKEGG(gene = k,organism = "hsa",pvalueCutoff = 0.05,pAdjustMethod = "BH")))
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T]=lapply(olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T],`[[`,2)
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==F]="No enriched KEGG pathways found!"
  pval_vector=unlist(lapply(overlapping_genes,function(x)phyper(length(intersect(x,tanzi_protective_variants)),length(x),19000,length(tanzi_protective_variants),lower.tail = F)))
  data.frame(ModuleName=names(ConsensusUniqueModules$DLPFC[i][module_match]),
             Jaccard=unlist(jcd[module_match]),
             ModuleSize=module_size[module_match],
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes[module_match],length)),
             Nr_Of_OlpKEGG=unlist(lapply(olp_KEGG,paste,sep=";",collapse=";")),
             p.overlap=pval_vector[module_match],stringsAsFactors = F)
  
}
comp_ConsensusModules_protective_DLPFC.df=data.frame(rbindlist(comp_ConsensusModules_protective_DLPFC.res[lapply(comp_ConsensusModules_protective_DLPFC.res,function(x)dim(x)[1])>0]),stringsAsFactors = F)%>%arrange(desc(-p.overlap))
write.table(comp_ConsensusModules_protective_DLPFC.df,"../../../Collaborations/Tanzi_WGS/ConsensusModules_protective_DLPFC_comparison.txt",col.names = T,row.names = F,sep = "\t",quote = F)

comp_ConsensusModules_protective_TCX.res=foreach(i=1:length(ConsensusUniqueModules$TCX))%do%{
  jcd=lapply(ConsensusUniqueModules$TCX[i],jaccard,tanzi_protective_variants)
  module_size=unlist(lapply(ConsensusUniqueModules$TCX[i],length))
  module_match=lapply(lapply(ConsensusUniqueModules$TCX[i],intersect,tanzi_protective_variants),length)>1
  overlapping_genes=lapply(lapply(ConsensusUniqueModules$TCX[i],intersect,tanzi_protective_variants)[module_match],function(y)mapIds2(IDs = y,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2])
  olp_KEGG=lapply(overlapping_genes,function(k)data.frame(enrichKEGG(gene = k,organism = "hsa",pvalueCutoff = 0.05,pAdjustMethod = "BH")))
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T]=lapply(olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T],`[[`,2)
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==F]="No enriched KEGG pathways found!"
  pval_vector=unlist(lapply(overlapping_genes,function(x)phyper(length(intersect(x,tanzi_protective_variants)),length(x),19000,length(tanzi_protective_variants),lower.tail = F)))
  data.frame(ModuleName=names(ConsensusUniqueModules$TCX[i][module_match]),
             Jaccard=unlist(jcd[module_match]),
             ModuleSize=module_size[module_match],
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes[module_match],length)),
             Nr_Of_OlpKEGG=unlist(lapply(olp_KEGG,paste,sep=";",collapse=";")),
             p.overlap=pval_vector[module_match],stringsAsFactors = F)
  
}
comp_ConsensusModules_protective_TCX.df=data.frame(rbindlist(comp_ConsensusModules_protective_TCX.res[lapply(comp_ConsensusModules_protective_TCX.res,function(x)dim(x)[1])>0]),stringsAsFactors = F)%>%arrange(desc(-p.overlap))
write.table(comp_ConsensusModules_protective_DLPFC.df,"../../../Collaborations/Tanzi_WGS/ConsensusModules_protective_DLPFC_comparison.txt",col.names = T,row.names = F,sep = "\t",quote = F)

comp_ConsensusModules_protective_CER.res=foreach(i=1:length(ConsensusUniqueModules$CER))%do%{
  jcd=lapply(ConsensusUniqueModules$CER[i],jaccard,tanzi_protective_variants)
  module_size=unlist(lapply(ConsensusUniqueModules$CER[i],length))
  module_match=lapply(lapply(ConsensusUniqueModules$CER[i],intersect,tanzi_protective_variants),length)>1
  overlapping_genes=lapply(lapply(ConsensusUniqueModules$CER[i],intersect,tanzi_protective_variants)[module_match],function(y)mapIds2(IDs = y,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2])
  olp_KEGG=lapply(overlapping_genes,function(k)data.frame(enrichKEGG(gene = k,organism = "hsa",pvalueCutoff = 0.05,pAdjustMethod = "BH")))
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T]=lapply(olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T],`[[`,2)
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==F]="No enriched KEGG pathways found!"
  pval_vector=unlist(lapply(overlapping_genes,function(x)phyper(length(intersect(x,tanzi_protective_variants)),length(x),19000,length(tanzi_protective_variants),lower.tail = F)))
  data.frame(ModuleName=names(ConsensusUniqueModules$CER[i][module_match]),
             Jaccard=unlist(jcd[module_match]),
             ModuleSize=module_size[module_match],
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes[module_match],length)),
             Nr_Of_OlpKEGG=unlist(lapply(olp_KEGG,paste,sep=";",collapse=";")),
             p.overlap=pval_vector[module_match],stringsAsFactors = F)
  
}
comp_ConsensusModules_protective_CER.df=data.frame(rbindlist(comp_ConsensusModules_protective_CER.res[lapply(comp_ConsensusModules_protective_CER.res,function(x)dim(x)[1])>0]),stringsAsFactors = F)%>%arrange(desc(-p.overlap))
write.table(comp_ConsensusModules_protective_DLPFC.df,"../../../Collaborations/Tanzi_WGS/ConsensusModules_protective_DLPFC_comparison.txt",col.names = T,row.names = F,sep = "\t",quote = F)

comp_ConsensusModules_protective_FP.res=foreach(i=1:length(ConsensusUniqueModules$FP))%do%{
  jcd=lapply(ConsensusUniqueModules$FP[i],jaccard,tanzi_protective_variants)
  module_size=unlist(lapply(ConsensusUniqueModules$FP[i],length))
  module_match=lapply(lapply(ConsensusUniqueModules$FP[i],intersect,tanzi_protective_variants),length)>1
  overlapping_genes=lapply(lapply(ConsensusUniqueModules$FP[i],intersect,tanzi_protective_variants)[module_match],function(y)mapIds2(IDs = y,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2])
  olp_KEGG=lapply(overlapping_genes,function(k)data.frame(enrichKEGG(gene = k,organism = "hsa",pvalueCutoff = 0.05,pAdjustMethod = "BH")))
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T]=lapply(olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T],`[[`,2)
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==F]="No enriched KEGG pathways found!"
  pval_vector=unlist(lapply(overlapping_genes,function(x)phyper(length(intersect(x,tanzi_protective_variants)),length(x),19000,length(tanzi_protective_variants),lower.tail = F)))
  data.frame(ModuleName=names(ConsensusUniqueModules$FP[i][module_match]),
             Jaccard=unlist(jcd[module_match]),
             ModuleSize=module_size[module_match],
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes[module_match],length)),
             Nr_Of_OlpKEGG=unlist(lapply(olp_KEGG,paste,sep=";",collapse=";")),
             p.overlap=pval_vector[module_match],stringsAsFactors = F)
  
}
comp_ConsensusModules_protective_FP.df=data.frame(rbindlist(comp_ConsensusModules_protective_FP.res[lapply(comp_ConsensusModules_protective_FP.res,function(x)dim(x)[1])>0]),stringsAsFactors = F)%>%arrange(desc(-p.overlap))
write.table(comp_ConsensusModules_protective_FP.df,"../../../Collaborations/Tanzi_WGS/ConsensusModules_protective_FP_comparison.txt",col.names = T,row.names = F,sep = "\t",quote = F)

comp_ConsensusModules_protective_IFG.res=foreach(i=1:length(ConsensusUniqueModules$IFG))%do%{
  jcd=lapply(ConsensusUniqueModules$IFG[i],jaccard,tanzi_protective_variants)
  module_size=unlist(lapply(ConsensusUniqueModules$IFG[i],length))
  module_match=lapply(lapply(ConsensusUniqueModules$IFG[i],intersect,tanzi_protective_variants),length)>1
  overlapping_genes=lapply(lapply(ConsensusUniqueModules$IFG[i],intersect,tanzi_protective_variants)[module_match],function(y)mapIds2(IDs = y,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2])
  olp_KEGG=lapply(overlapping_genes,function(k)data.frame(enrichKEGG(gene = k,organism = "hsa",pvalueCutoff = 0.05,pAdjustMethod = "BH")))
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T]=lapply(olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T],`[[`,2)
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==F]="No enriched KEGG pathways found!"
  pval_vector=unlist(lapply(overlapping_genes,function(x)phyper(length(intersect(x,tanzi_protective_variants)),length(x),19000,length(tanzi_protective_variants),lower.tail = F)))
  data.frame(ModuleName=names(ConsensusUniqueModules$IFG[i][module_match]),
             Jaccard=unlist(jcd[module_match]),
             ModuleSize=module_size[module_match],
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes[module_match],length)),
             Nr_Of_OlpKEGG=unlist(lapply(olp_KEGG,paste,sep=";",collapse=";")),
             p.overlap=pval_vector[module_match],stringsAsFactors = F)
  
}
comp_ConsensusModules_protective_IFG.df=data.frame(rbindlist(comp_ConsensusModules_protective_IFG.res[lapply(comp_ConsensusModules_protective_IFG.res,function(x)dim(x)[1])>0]),stringsAsFactors = F)%>%arrange(desc(-p.overlap))
write.table(comp_ConsensusModules_protective_IFG.df,"../../../Collaborations/Tanzi_WGS/ConsensusModules_protective_IFG_comparison.txt",col.names = T,row.names = F,sep = "\t",quote = F)

comp_ConsensusModules_protective_PHG.res=foreach(i=1:length(ConsensusUniqueModules$PHG))%do%{
  jcd=lapply(ConsensusUniqueModules$PHG[i],jaccard,tanzi_protective_variants)
  module_size=unlist(lapply(ConsensusUniqueModules$PHG[i],length))
  module_match=lapply(lapply(ConsensusUniqueModules$PHG[i],intersect,tanzi_protective_variants),length)>1
  overlapping_genes=lapply(lapply(ConsensusUniqueModules$PHG[i],intersect,tanzi_protective_variants)[module_match],function(y)mapIds2(IDs = y,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2])
  olp_KEGG=lapply(overlapping_genes,function(k)data.frame(enrichKEGG(gene = k,organism = "hsa",pvalueCutoff = 0.05,pAdjustMethod = "BH")))
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T]=lapply(olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T],`[[`,2)
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==F]="No enriched KEGG pathways found!"
  pval_vector=unlist(lapply(overlapping_genes,function(x)phyper(length(intersect(x,tanzi_protective_variants)),length(x),19000,length(tanzi_protective_variants),lower.tail = F)))
  data.frame(ModuleName=names(ConsensusUniqueModules$PHG[i][module_match]),
             Jaccard=unlist(jcd[module_match]),
             ModuleSize=module_size[module_match],
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes[module_match],length)),
             Nr_Of_OlpKEGG=unlist(lapply(olp_KEGG,paste,sep=";",collapse=";")),
             p.overlap=pval_vector[module_match],stringsAsFactors = F)
  
}
comp_ConsensusModules_protective_PHG.df=data.frame(rbindlist(comp_ConsensusModules_protective_PHG.res[lapply(comp_ConsensusModules_protective_PHG.res,function(x)dim(x)[1])>0]),stringsAsFactors = F)%>%arrange(desc(-p.overlap))
write.table(comp_ConsensusModules_protective_PHG.df,"../../../Collaborations/Tanzi_WGS/ConsensusModules_protective_PHG_comparison.txt",col.names = T,row.names = F,sep = "\t",quote = F)

comp_ConsensusModules_protective_STG.res=foreach(i=1:length(ConsensusUniqueModules$STG))%do%{
  jcd=lapply(ConsensusUniqueModules$STG[i],jaccard,tanzi_protective_variants)
  module_size=unlist(lapply(ConsensusUniqueModules$STG[i],length))
  module_match=lapply(lapply(ConsensusUniqueModules$STG[i],intersect,tanzi_protective_variants),length)>1
  overlapping_genes=lapply(lapply(ConsensusUniqueModules$STG[i],intersect,tanzi_protective_variants)[module_match],function(y)mapIds2(IDs = y,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2])
  olp_KEGG=lapply(overlapping_genes,function(k)data.frame(enrichKEGG(gene = k,organism = "hsa",pvalueCutoff = 0.05,pAdjustMethod = "BH")))
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T]=lapply(olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T],`[[`,2)
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==F]="No enriched KEGG pathways found!"
  pval_vector=unlist(lapply(overlapping_genes,function(x)phyper(length(intersect(x,tanzi_protective_variants)),length(x),19000,length(tanzi_protective_variants),lower.tail = F)))
  data.frame(ModuleName=names(ConsensusUniqueModules$STG[i][module_match]),
             Jaccard=unlist(jcd[module_match]),
             ModuleSize=module_size[module_match],
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes[module_match],length)),
             Nr_Of_OlpKEGG=unlist(lapply(olp_KEGG,paste,sep=";",collapse=";")),
             p.overlap=pval_vector[module_match],stringsAsFactors = F)
  
}
comp_ConsensusModules_protective_STG.df=data.frame(rbindlist(comp_ConsensusModules_protective_STG.res[lapply(comp_ConsensusModules_protective_STG.res,function(x)dim(x)[1])>0]),stringsAsFactors = F)%>%arrange(desc(-p.overlap))
write.table(comp_ConsensusModules_protective_STG.df,"../../../Collaborations/Tanzi_WGS/ConsensusModules_protective_STG_comparison.txt",col.names = T,row.names = F,sep = "\t",quote = F)

#Compare with Aggregate modules
comp_AggModules_protective_DLPFC.res=foreach(i=1:length(ampAD_AggModules_data_brainRegion$DLPFC))%do%{
  jcd=lapply(ampAD_AggModules_data_brainRegion$DLPFC[i],jaccard,tanzi_protective_variants)
  module_size=unlist(lapply(ampAD_AggModules_data_brainRegion$DLPFC[i],length))
  module_match=lapply(lapply(ampAD_AggModules_data_brainRegion$DLPFC[i],intersect,tanzi_protective_variants),length)>1
  overlapping_genes=lapply(lapply(ampAD_AggModules_data_brainRegion$DLPFC[i],intersect,tanzi_protective_variants)[module_match],function(y)mapIds2(IDs = y,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2])
  olp_KEGG=lapply(overlapping_genes,function(k)data.frame(enrichKEGG(gene = k,organism = "hsa",pvalueCutoff = 0.05,pAdjustMethod = "BH")))
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T]=lapply(olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T],`[[`,2)
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==F]="No enriched KEGG pathways found!"
  pval_vector=unlist(lapply(overlapping_genes,function(x)phyper(length(intersect(x,tanzi_protective_variants)),length(x),19000,length(tanzi_protective_variants),lower.tail = F)))
  data.frame(ModuleName=names(ampAD_AggModules_data_brainRegion$DLPFC[i][module_match]),
             Jaccard=unlist(jcd[module_match]),
             ModuleSize=module_size[module_match],
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes[module_match],length)),
             Nr_Of_OlpKEGG=unlist(lapply(olp_KEGG,paste,sep=";",collapse=";")),
             p.overlap=pval_vector[module_match],stringsAsFactors = F)
  
}
comp_AggModules_protective_DLPFC.df=data.frame(rbindlist(comp_AggModules_protective_DLPFC.res[lapply(comp_AggModules_protective_DLPFC.res,function(x)dim(x)[1])>0]),stringsAsFactors = F)%>%arrange(desc(-p.overlap))
write.table(comp_AggModules_protective_DLPFC.df,"../../../Collaborations/Tanzi_WGS/AggModules_protective_DLPFC_comparison.txt",col.names = T,row.names = F,sep = "\t",quote = F)

comp_AggModules_protective_CER.res=foreach(i=1:length(ampAD_AggModules_data_brainRegion$CBE))%do%{
  jcd=lapply(ampAD_AggModules_data_brainRegion$CBE[i],jaccard,tanzi_protective_variants)
  module_size=unlist(lapply(ampAD_AggModules_data_brainRegion$CBE[i],length))
  module_match=lapply(lapply(ampAD_AggModules_data_brainRegion$CBE[i],intersect,tanzi_protective_variants),length)>1
  overlapping_genes=lapply(lapply(ampAD_AggModules_data_brainRegion$CBE[i],intersect,tanzi_protective_variants)[module_match],function(y)mapIds2(IDs = y,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2])
  olp_KEGG=lapply(overlapping_genes,function(k)data.frame(enrichKEGG(gene = k,organism = "hsa",pvalueCutoff = 0.05,pAdjustMethod = "BH")))
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T]=lapply(olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T],`[[`,2)
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==F]="No enriched KEGG pathways found!"
  pval_vector=unlist(lapply(overlapping_genes,function(x)phyper(length(intersect(x,tanzi_protective_variants)),length(x),19000,length(tanzi_protective_variants),lower.tail = F)))
  data.frame(ModuleName=names(ampAD_AggModules_data_brainRegion$CBE[i][module_match]),
             Jaccard=unlist(jcd[module_match]),
             ModuleSize=module_size[module_match],
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes[module_match],length)),
             Nr_Of_OlpKEGG=unlist(lapply(olp_KEGG,paste,sep=";",collapse=";")),
             p.overlap=pval_vector[module_match],stringsAsFactors = F)
  
}
comp_AggModules_protective_CER.df=data.frame(rbindlist(comp_AggModules_protective_CER.res[lapply(comp_AggModules_protective_CER.res,function(x)dim(x)[1])>0]),stringsAsFactors = F)%>%arrange(desc(-p.overlap))
write.table(comp_AggModules_protective_CER.df,"../../../Collaborations/Tanzi_WGS/AggModules_protective_CER_comparison.txt",col.names = T,row.names = F,sep = "\t",quote = F)

comp_AggModules_protective_TCX.res=foreach(i=1:length(ampAD_AggModules_data_brainRegion$TCX))%do%{
  jcd=lapply(ampAD_AggModules_data_brainRegion$TCX[i],jaccard,tanzi_protective_variants)
  module_size=unlist(lapply(ampAD_AggModules_data_brainRegion$TCX[i],length))
  module_match=lapply(lapply(ampAD_AggModules_data_brainRegion$TCX[i],intersect,tanzi_protective_variants),length)>1
  overlapping_genes=lapply(lapply(ampAD_AggModules_data_brainRegion$TCX[i],intersect,tanzi_protective_variants)[module_match],function(y)mapIds2(IDs = y,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2])
  olp_KEGG=lapply(overlapping_genes,function(k)data.frame(enrichKEGG(gene = k,organism = "hsa",pvalueCutoff = 0.05,pAdjustMethod = "BH")))
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T]=lapply(olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T],`[[`,2)
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==F]="No enriched KEGG pathways found!"
  pval_vector=unlist(lapply(overlapping_genes,function(x)phyper(length(intersect(x,tanzi_protective_variants)),length(x),19000,length(tanzi_protective_variants),lower.tail = F)))
  data.frame(ModuleName=names(ampAD_AggModules_data_brainRegion$TCX[i][module_match]),
             Jaccard=unlist(jcd[module_match]),
             ModuleSize=module_size[module_match],
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes[module_match],length)),
             Nr_Of_OlpKEGG=unlist(lapply(olp_KEGG,paste,sep=";",collapse=";")),
             p.overlap=pval_vector[module_match],stringsAsFactors = F)
  
}
comp_AggModules_protective_TCX.df=data.frame(rbindlist(comp_AggModules_protective_TCX.res[lapply(comp_AggModules_protective_TCX.res,function(x)dim(x)[1])>0]),stringsAsFactors = F)%>%arrange(desc(-p.overlap))
write.table(comp_AggModules_protective_TCX.df,"../../../Collaborations/Tanzi_WGS/AggModules_protective_TCX_comparison.txt",col.names = T,row.names = F,sep = "\t",quote = F)

comp_AggModules_protective_FP.res=foreach(i=1:length(ampAD_AggModules_data_brainRegion$FP))%do%{
  jcd=lapply(ampAD_AggModules_data_brainRegion$FP[i],jaccard,tanzi_protective_variants)
  module_size=unlist(lapply(ampAD_AggModules_data_brainRegion$FP[i],length))
  module_match=lapply(lapply(ampAD_AggModules_data_brainRegion$FP[i],intersect,tanzi_protective_variants),length)>1
  overlapping_genes=lapply(lapply(ampAD_AggModules_data_brainRegion$FP[i],intersect,tanzi_protective_variants)[module_match],function(y)mapIds2(IDs = y,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2])
  olp_KEGG=lapply(overlapping_genes,function(k)data.frame(enrichKEGG(gene = k,organism = "hsa",pvalueCutoff = 0.05,pAdjustMethod = "BH")))
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T]=lapply(olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T],`[[`,2)
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==F]="No enriched KEGG pathways found!"
  pval_vector=unlist(lapply(overlapping_genes,function(x)phyper(length(intersect(x,tanzi_protective_variants)),length(x),19000,length(tanzi_protective_variants),lower.tail = F)))
  data.frame(ModuleName=names(ampAD_AggModules_data_brainRegion$FP[i][module_match]),
             Jaccard=unlist(jcd[module_match]),
             ModuleSize=module_size[module_match],
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes[module_match],length)),
             Nr_Of_OlpKEGG=unlist(lapply(olp_KEGG,paste,sep=";",collapse=";")),
             p.overlap=pval_vector[module_match],stringsAsFactors = F)
  
}
comp_AggModules_protective_FP.df=data.frame(rbindlist(comp_AggModules_protective_FP.res[lapply(comp_AggModules_protective_FP.res,function(x)dim(x)[1])>0]),stringsAsFactors = F)%>%arrange(desc(-p.overlap))
write.table(comp_AggModules_protective_FP.df,"../../../Collaborations/Tanzi_WGS/AggModules_protective_FP_comparison.txt",col.names = T,row.names = F,sep = "\t",quote = F)

comp_AggModules_protective_IFG.res=foreach(i=1:length(ampAD_AggModules_data_brainRegion$IFG))%do%{
  jcd=lapply(ampAD_AggModules_data_brainRegion$IFG[i],jaccard,tanzi_protective_variants)
  module_size=unlist(lapply(ampAD_AggModules_data_brainRegion$IFG[i],length))
  module_match=lapply(lapply(ampAD_AggModules_data_brainRegion$IFG[i],intersect,tanzi_protective_variants),length)>1
  overlapping_genes=lapply(lapply(ampAD_AggModules_data_brainRegion$IFG[i],intersect,tanzi_protective_variants)[module_match],function(y)mapIds2(IDs = y,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2])
  olp_KEGG=lapply(overlapping_genes,function(k)data.frame(enrichKEGG(gene = k,organism = "hsa",pvalueCutoff = 0.05,pAdjustMethod = "BH")))
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T]=lapply(olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T],`[[`,2)
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==F]="No enriched KEGG pathways found!"
  pval_vector=unlist(lapply(overlapping_genes,function(x)phyper(length(intersect(x,tanzi_protective_variants)),length(x),19000,length(tanzi_protective_variants),lower.tail = F)))
  data.frame(ModuleName=names(ampAD_AggModules_data_brainRegion$IFG[i][module_match]),
             Jaccard=unlist(jcd[module_match]),
             ModuleSize=module_size[module_match],
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes[module_match],length)),
             Nr_Of_OlpKEGG=unlist(lapply(olp_KEGG,paste,sep=";",collapse=";")),
             p.overlap=pval_vector[module_match],stringsAsFactors = F)
  
}
comp_AggModules_protective_IFG.df=data.frame(rbindlist(comp_AggModules_protective_IFG.res[lapply(comp_AggModules_protective_IFG.res,function(x)dim(x)[1])>0]),stringsAsFactors = F)%>%arrange(desc(-p.overlap))
write.table(comp_AggModules_protective_IFG.df,"../../../Collaborations/Tanzi_WGS/AggModules_protective_IFG_comparison.txt",col.names = T,row.names = F,sep = "\t",quote = F)

comp_AggModules_protective_PHG.res=foreach(i=1:length(ampAD_AggModules_data_brainRegion$PHG))%do%{
  jcd=lapply(ampAD_AggModules_data_brainRegion$PHG[i],jaccard,tanzi_protective_variants)
  module_size=unlist(lapply(ampAD_AggModules_data_brainRegion$PHG[i],length))
  module_match=lapply(lapply(ampAD_AggModules_data_brainRegion$PHG[i],intersect,tanzi_protective_variants),length)>1
  overlapping_genes=lapply(lapply(ampAD_AggModules_data_brainRegion$PHG[i],intersect,tanzi_protective_variants),function(y)mapIds2(IDs = y,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2])
  olp_KEGG=lapply(overlapping_genes,function(k)data.frame(enrichKEGG(gene = k,organism = "hsa",pvalueCutoff = 0.05,pAdjustMethod = "BH")))
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T]=lapply(olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T],`[[`,2)
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==F]="No enriched KEGG pathways found!"
  pval_vector=unlist(lapply(overlapping_genes,function(x)phyper(length(intersect(x,tanzi_protective_variants)),length(x),19000,length(tanzi_protective_variants),lower.tail = F)))
  data.frame(ModuleName=names(ampAD_AggModules_data_brainRegion$PHG[i]),
             Jaccard=unlist(jcd),
             ModuleSize=module_size,
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes,length)),
             Nr_Of_OlpKEGG=unlist(lapply(olp_KEGG,paste,sep=";",collapse=";")),
             p.overlap=pval_vector,stringsAsFactors = F)
  
}
comp_AggModules_protective_PHG.df=data.frame(rbindlist(comp_AggModules_protective_PHG.res[lapply(comp_AggModules_protective_PHG.res,function(x)dim(x)[1])>0]),stringsAsFactors = F)%>%arrange(desc(-p.overlap))
write.table(comp_AggModules_protective_PHG.df,"../../../Collaborations/Tanzi_WGS/AggModules_protective_PHG_comparison.txt",col.names = T,row.names = F,sep = "\t",quote = F)

comp_AggModules_protective_STG.res=foreach(i=1:length(ampAD_AggModules_data_brainRegion$STG))%do%{
  jcd=lapply(ampAD_AggModules_data_brainRegion$STG[i],jaccard,tanzi_protective_variants)
  module_size=unlist(lapply(ampAD_AggModules_data_brainRegion$STG[i],length))
  module_match=lapply(lapply(ampAD_AggModules_data_brainRegion$STG[i],intersect,tanzi_protective_variants),length)>1
  overlapping_genes=lapply(lapply(ampAD_AggModules_data_brainRegion$STG[i],intersect,tanzi_protective_variants)[module_match],function(y)mapIds2(IDs = y,IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2])
  olp_KEGG=lapply(overlapping_genes,function(k)data.frame(enrichKEGG(gene = k,organism = "hsa",pvalueCutoff = 0.05,pAdjustMethod = "BH")))
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T]=lapply(olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==T],`[[`,2)
  olp_KEGG[lapply(olp_KEGG,function(d)dim(d)[1]>0)==F]="No enriched KEGG pathways found!"
  pval_vector=unlist(lapply(overlapping_genes,function(x)phyper(length(intersect(x,tanzi_protective_variants)),length(x),19000,length(tanzi_protective_variants),lower.tail = F)))
  data.frame(ModuleName=names(ampAD_AggModules_data_brainRegion$STG[i][module_match]),
             Jaccard=unlist(jcd[module_match]),
             ModuleSize=module_size[module_match],
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes[module_match],length)),
             Nr_Of_OlpKEGG=unlist(lapply(olp_KEGG,paste,sep=";",collapse=";")),
             p.overlap=pval_vector,stringsAsFactors = F)
  
}
comp_AggModules_protective_STG.df=data.frame(rbindlist(comp_AggModules_protective_STG.res[lapply(comp_AggModules_protective_STG.res,function(x)dim(x)[1])>0]),stringsAsFactors = F)%>%arrange(desc(-p.overlap))
write.table(comp_AggModules_protective_STG.df,"../../../Collaborations/Tanzi_WGS/AggModules_protective_STG_comparison.txt",col.names = T,row.names = F,sep = "\t",quote = F)
