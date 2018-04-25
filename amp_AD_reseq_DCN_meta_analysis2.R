library(org.Hs.eg.db)
library(igraph)
library(data.table)
library(centiserve)
library(sets)
library(gProfileR)
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
names(ConsensusDEGs)=names(ConsensusDEGs_df)=c("TCX","CER","DLPFC","FP","IFG","PHG","STG")

ConsensusDEGs_df$TCX=fread("../DEG_Analyses/MAYO_CBE_TCX_DiffExpression.tsv",header = T,stringsAsFactors = F,data.table = F) %>%filter(BrainRegion=="TCX"&Comparison=="CONTROL-AD")
ConsensusDEGs_df$CER=fread("../DEG_Analyses/MAYO_CBE_TCX_DiffExpression.tsv",header = T,stringsAsFactors = F,data.table = F) %>%filter(BrainRegion=="CER"&Comparison=="CONTROL-AD")
ConsensusDEGs_df$DLPFC=fread("../DEG_Analyses/ROSMAP_DLPFC_DiffExpression.tsv",header = T,stringsAsFactors = F,data.table = F) %>%filter(Comparison=="cogdx1-cogdx4")
ConsensusDEGs_df$FP=fread("../DEG_Analyses/MSSM_FP_STG_PHG_IFG_DiffExpression.tsv",header = T,stringsAsFactors = F,data.table = F)%>%filter(BrainRegion=="FP"&Comparison=="CONTROL-AD")
ConsensusDEGs_df$IFG=fread("../DEG_Analyses/MSSM_FP_STG_PHG_IFG_DiffExpression.tsv",header = T,stringsAsFactors = F,data.table = F)%>%filter(BrainRegion=="IFG"&Comparison=="CONTROL-AD")
ConsensusDEGs_df$PHG=fread("../DEG_Analyses/MSSM_FP_STG_PHG_IFG_DiffExpression.tsv",header = T,stringsAsFactors = F,data.table = F)%>%filter(BrainRegion=="PHG"&Comparison=="CONTROL-AD")
ConsensusDEGs_df$STG=fread("../DEG_Analyses/MSSM_FP_STG_PHG_IFG_DiffExpression.tsv",header = T,stringsAsFactors = F,data.table = F)%>%filter(BrainRegion=="STG"&Comparison=="CONTROL-AD")
ConsensusDEGs=lapply(ConsensusDEGs_df,function(x)x%>%filter(Direction=="UP"|Direction=="DOWN")%>%filter(hgnc_symbol!="")%>%pull(hgnc_symbol))

metaAD_DCLs=metaAD_DCGs=metaAD_DRLs=metaAD_DRGs=metaAD_DRsort_results=metaAD_DiffCoexp=vector(mode = "list",length=7)
names(metaAD_DCLs)=names(metaAD_DCGs)=names(metaAD_DRLs)=names(metaAD_DRGs)=names(metaAD_DRsort_results)=names(metaAD_DiffCoexp)=c("Mayo_TCX","Mayo_CER","ROSMAP_DLPFC","MSBB_FP","MSBB_IFG","MSBB_PHG","MSBB_STG")

metaAD_DRsort_results$Mayo_TCX=readRDS("./DCGL_results/TCX_DiffCoexp_DRsort_noLFC_bugfix2.res.RDS")
metaAD_DRsort_results$Mayo_CER=readRDS("./DCGL_results/CER_DiffCoexp_DRsort_noLFC_bugfix2.res.RDS")
metaAD_DRsort_results$ROSMAP_DLPFC=readRDS("./DCGL_results/ROSMAP_DiffCoexp_DRsort_noLFC_bugfix2.res.RDS")
metaAD_DRsort_results$MSBB_FP=readRDS("./DCGL_results/MSBB_FP_DRsort_noLFC_bugfix2.res.RDS")
metaAD_DRsort_results$MSBB_IFG=readRDS("./DCGL_results/MSBB_IFG_DRsort_noLFC_bugfix2.res.RDS")
metaAD_DRsort_results$MSBB_PHG=readRDS("./DCGL_results/MSBB_PHG_DRsort_noLFC_bugfix2.res.RDS")
metaAD_DRsort_results$MSBB_STG=readRDS("./DCGL_results/MSBB_STG_DiffCoexp_DRsort_noLFC_bugfix2.res.RDS")

metaAD_DCLs=lapply(metaAD_DRsort_results,function(x)x$DCLs)
metaAD_DCLs.graph=lapply(metaAD_DCLs,function(x)data.frame(x)%>%filter(q.diffcor<0.1)%>%select(c(Gene.1,Gene.2))%>%graph.data.frame(directed = F))
metaAD_DCGs=lapply(metaAD_DRsort_results,function(x)data.frame(x$DCGs,stringsAsFactors = F))
metaAD_DCGs.filtered=lapply(metaAD_DCGs,function(x)x%>%filter(q<0.1)%>%pull(DCG)%>%droplevels%>%levels)
metaAD_DRLs=lapply(metaAD_DRsort_results,function(x)x$DRLs%>%filter(q.diffcor<0.1))
metaAD_DRLs.graph=lapply(metaAD_DRLs,function(x)data.frame(x)%>%filter(q.diffcor<0.1)%>%select(c(Gene.1,Gene.2))%>%graph.data.frame(directed = F))
metaAD_DRGs=lapply(metaAD_DRsort_results,function(x)x$DRGs)
metaAD_DRGs.filtered=lapply(metaAD_DRsort_results,function(x)x$DRGs%>%filter(q<0.1)%>%pull(DCG)%>%droplevels%>%levels)
metaAD_TF_bridged_DCLs=lapply(metaAD_DRsort_results,function(x)x$TF_bridged_DCL)
metaAD_TF_bridged_DCLs.filtered=lapply(metaAD_TF_bridged_DCLs,function(x)x%>%filter(q.diffcor<0.1)%>%select(c(Gene.1,Gene.2)))

jaccard_DCG_DEG=mapply(FUN=function(a,b)jaccard(a,b),metaAD_DCGs.filtered,ConsensusDEGs,SIMPLIFY = F)
nr_DEGs=lapply(ConsensusDEGs,length)
nr_DCGs=lapply(metaAD_DCGs.filtered,length)
overlapping_genes=mapply(FUN=function(a,b)intersect(a,b),metaAD_DCGs.filtered,ConsensusDEGs,SIMPLIFY = F)
overlapping_genes_size=lapply(overlapping_genes,length)
olp_KEGG=lapply(lapply(lapply(overlapping_genes,gprofiler,organism = "hsapiens",significant = T,src_filter = "KEGG",correction_method = "bonferroni",max_p_value = 0.05),`[`,c(9,12)),function(x)paste(x[,1],x[,2],sep = "--",collapse = ";"))
olp_KEGG[lapply(olp_KEGG,function(x)x=="")==T]="No enriched KEGG pathways found!"
fisher_results=exp_DCG=vector(mode="list",length = length(overlapping_genes_size))
names(fisher_results)=names(exp_DCG)=names(metaAD_DCGs.filtered)
#Run Fisher's test individually by dataset as each dataset had uniquely profiled genes, after normalisation
fisher_results$ROSMAP_DLPFC=fisher.test(matrix(c(overlapping_genes_size$ROSMAP_DLPFC,
                                                 length(metaAD_DCGs.filtered$ROSMAP_DLPFC)-overlapping_genes_size$ROSMAP_DLPFC,
                                                 length(ConsensusDEGs$DLPFC)-overlapping_genes_size$ROSMAP_DLPFC,
                                                 14228-(length(metaAD_DCGs.filtered$ROSMAP_DLPFC)-overlapping_genes_size$ROSMAP_DLPFC)-(length(ConsensusDEGs$DLPFC)-overlapping_genes_size$ROSMAP_DLPFC)-overlapping_genes_size$ROSMAP_DLPFC),nrow = 2))
exp_DCG$ROSMAP_DLPFC=(length(metaAD_DCGs.filtered$ROSMAP_DLPFC)/14228)*length(ConsensusDEGs$DLPFC)

fisher_results$Mayo_TCX=fisher.test(matrix(c(overlapping_genes_size$Mayo_TCX,
                                           length(metaAD_DCGs.filtered$Mayo_TCX)-overlapping_genes_size$Mayo_TCX,
                                           length(ConsensusDEGs$TCX)-overlapping_genes_size$Mayo_TCX,
                                           15160-(length(metaAD_DCGs.filtered$Mayo_TCX)-overlapping_genes_size$Mayo_TCX)-(length(ConsensusDEGs$TCX)-overlapping_genes_size$Mayo_TCX)-overlapping_genes_size$Mayo_TCX),nrow = 2))

exp_DCG=mapply(FUN=function(a,b)(a/15160)*b,nr_DCGs,nr_DEGs,SIMPLIFY = F)

fisher_results$Mayo_CER=fisher.test(matrix(c(overlapping_genes_size$Mayo_CER,
                                             length(metaAD_DCGs.filtered$Mayo_CER)-overlapping_genes_size$Mayo_CER,
                                             length(ConsensusDEGs$CER)-overlapping_genes_size$Mayo_CER,
                                             15160-(length(metaAD_DCGs.filtered$Mayo_CER)-overlapping_genes_size$Mayo_CER)-(length(ConsensusDEGs$CER)-overlapping_genes_size$Mayo_CER)-overlapping_genes_size$Mayo_CER),nrow = 2))
exp_DCG$Mayo_CER=(length(metaAD_DCGs.filtered$Mayo_CER)/15160)*length(ConsensusDEGs$CER)
fisher_results$MSBB_FP=fisher.test(matrix(c(overlapping_genes_size$MSBB_FP,
                                            length(metaAD_DCGs.filtered$MSBB_FP)-overlapping_genes_size$MSBB_FP,
                                            length(ConsensusDEGs$FP)-overlapping_genes_size$MSBB_FP,
                                            14681-(length(metaAD_DCGs.filtered$MSBB_FP)-overlapping_genes_size$MSBB_FP)-(length(ConsensusDEGs$FP)-overlapping_genes_size$MSBB_FP)-overlapping_genes_size$MSBB_FP),nrow = 2))
exp_DCG$MSBB_FP=(length(metaAD_DCGs.filtered$MSBB_FP)/14681)*length(ConsensusDEGs$FP)

fisher_results$MSBB_IFG=fisher.test(matrix(c(overlapping_genes_size$MSBB_IFG,
                                             length(metaAD_DCGs.filtered$MSBB_IFG)-overlapping_genes_size$MSBB_IFG,
                                             length(ConsensusDEGs$IFG)-overlapping_genes_size$MSBB_IFG,
                                             14681-(length(metaAD_DCGs.filtered$MSBB_IFG)-overlapping_genes_size$MSBB_IFG)-(length(ConsensusDEGs$IFG)-overlapping_genes_size$MSBB_IFG)-overlapping_genes_size$MSBB_IFG),nrow = 2))
exp_DCG$MSBB_IFG=(length(metaAD_DCGs.filtered$MSBB_IFG)/14681)*length(ConsensusDEGs$IFG)

fisher_results$MSBB_PHG=fisher.test(matrix(c(overlapping_genes_size$MSBB_PHG,
                                             length(metaAD_DCGs.filtered$MSBB_PHG)-overlapping_genes_size$MSBB_PHG,
                                             length(ConsensusDEGs$PHG)-overlapping_genes_size$MSBB_PHG,
                                             14681-(length(metaAD_DCGs.filtered$MSBB_PHG)-overlapping_genes_size$MSBB_PHG)-(length(ConsensusDEGs$PHG)-overlapping_genes_size$MSBB_PHG)-overlapping_genes_size$MSBB_PHG),nrow = 2))
exp_DCG$MSBB_PHG=(length(metaAD_DCGs.filtered$MSBB_PHG)/14681)*length(ConsensusDEGs$PHG)

fisher_results$MSBB_STG=fisher.test(matrix(c(overlapping_genes_size$MSBB_STG,
                                             length(metaAD_DCGs.filtered$MSBB_STG)-overlapping_genes_size$MSBB_STG,
                                             length(ConsensusDEGs$STG)-overlapping_genes_size$MSBB_STG,
                                             14681-(length(metaAD_DCGs.filtered$MSBB_STG)-overlapping_genes_size$MSBB_STG)-(length(ConsensusDEGs$STG)-overlapping_genes_size$MSBB_STG)-overlapping_genes_size$MSBB_STG),nrow = 2))
exp_DCG$MSBB_STG=(length(metaAD_DCGs.filtered$MSBB_STG)/14681)*length(ConsensusDEGs$STG)



for(m in 1:length(overlapping_genes_size)){
  if(overlapping_genes_size[[m]]>exp_DCG[[m]]){
    DCG_representation[[m]]="Over"  
  }
  else {DCG_representation[[m]]="Under" }
} 
comp_DCG_DEG.df=data.frame(Jcd=unlist(jaccard_DCG_DEG),
                           Nr_DCGs=unlist(nr_DCGs),
                           Nr_DEGs=unlist(nr_DEGs),
                           Nr_OverlappedGenes=unlist(overlapping_genes_size),
                           Expected_DCGs=unlist(exp_DCG),
                           Nr_Of_OlpKEGG=unlist(olp_KEGG),
                           p.overlap=unlist(lapply(fisher_results,function(p)p$p.value)),
                           DCG_Representation=unlist(DCG_representation),
                           stringsAsFactors = F)
write.table(comp_DCG_DEG.df,"Comp_DCG_DEG_df.txt",sep="\t",col.names = T,row.names = T,quote = F)

#Compare with AMP-AD aggregate modules
ampAD_AggModules.data=fread("../AggregateModules/AMP-AD_AggModules.tsv",sep = "\t",header = T,stringsAsFactors = F,showProgress = T,data.table = F)
tmp=lapply(unique(ampAD_AggModules.data$brainRegion),function(x)ampAD_AggModules.data%>%filter(brainRegion==x)%>%select(c(Module,external_gene_name)))
ampAD_AggModules_data_brainRegion=ampAD_AggModules_data_brainRegion_OlpGenes=vector(mode="list",length=7)
names(ampAD_AggModules_data_brainRegion)=names(ampAD_AggModules_data_brainRegion_OlpGenes)=unique(ampAD_AggModules.data$brainRegion)

for(r in 1:length(ampAD_AggModules_data_brainRegion)){
  ampAD_AggModules_data_brainRegion[[r]]=lapply(lapply(unique(tmp[[r]]$Module),function(x)tmp[[r]]%>%filter(Module==x)%>%select(external_gene_name)),`[[`,1)
  names(ampAD_AggModules_data_brainRegion[[r]])=unique(tmp[[r]]$Module)
}


comp_DCG_DLPFC.res=foreach(i=1:length(ampAD_AggModules_data_brainRegion$DLPFC))%dopar%{
  jaccard_DCG=lapply(ampAD_AggModules_data_brainRegion$DLPFC[i],jaccard,metaAD_DCGs.filtered$ROSMAP_DLPFC)
  module_size=lapply(ampAD_AggModules_data_brainRegion$DLPFC[i],length)
  overlapping_genes=lapply(ampAD_AggModules_data_brainRegion$DLPFC[i],intersect,metaAD_DCGs.filtered$ROSMAP_DLPFC)
  overlapping_genes_size=lapply(overlapping_genes,length)
  olp_KEGG=lapply(lapply(lapply(overlapping_genes,gprofiler,organism = "hsapiens",significant = T,src_filter = "KEGG",correction_method = "bonferroni",max_p_value = 0.05),`[`,c(9,12)),function(x)paste(x[,1],x[,2],sep = "--",collapse = ";"))
  olp_KEGG[lapply(olp_KEGG,function(x)x=="")==T]="No enriched KEGG pathways found!"
  fisher_results=DCG_representation=vector(mode="list",length = length(overlapping_genes_size))
  for(m in 1:length(overlapping_genes_size)){
    fisher_results[[m]]=fisher.test(matrix(c(overlapping_genes_size[[m]],
                                             length(metaAD_DCGs.filtered$ROSMAP_DLPFC)-overlapping_genes_size[[m]],
                                             module_size[[m]]-overlapping_genes_size[[m]],
                                             14228-(length(metaAD_DCGs.filtered$ROSMAP_DLPFC)-overlapping_genes_size[[m]])-length(module_size[[m]]-overlapping_genes_size[[m]])-overlapping_genes_size[[m]]),nrow = 2))
  }
  
  exp_DCG=as.list((length(metaAD_DCGs.filtered$ROSMAP_DLPFC)/14228)*unlist(module_size))
  for(m in 1:length(overlapping_genes_size)){
    if(overlapping_genes_size[[m]]>exp_DCG[[m]]){
           DCG_representation[[m]]="Over"  
    }
    else {DCG_representation[[m]]="Under" }
  } 
  data.frame(Jcd=unlist(jaccard_DCG),
             ModuleSize=unlist(module_size),
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes,length)),
             Expected_DCG_in_Module=unlist(exp_DCG),
             Nr_Of_OlpKEGG=unlist(olp_KEGG),
             p.overlap=unlist(lapply(fisher_results,function(p)p$p.value)),
             DCG_Representation=unlist(DCG_representation),
             stringsAsFactors = F)
}
comp_DCG_DLPFC.df=data.frame(rbindlist(comp_DCG_DLPFC.res),stringsAsFactors = F)
rownames(comp_DCG_DLPFC.df)=names(ampAD_AggModules_data_brainRegion$DLPFC)
write.table(comp_DCG_DLPFC.df,"AggModules_DCG_DLPFC_df.txt",sep="\t",col.names = T,row.names = T,quote = F)

comp_DCG_TCX.res=foreach(i=1:length(ampAD_AggModules_data_brainRegion$TCX))%dopar%{
  jaccard_DCG=lapply(ampAD_AggModules_data_brainRegion$TCX[i],jaccard,metaAD_DCGs.filtered$Mayo_TCX)
  module_size=lapply(ampAD_AggModules_data_brainRegion$TCX[i],length)
  overlapping_genes=lapply(ampAD_AggModules_data_brainRegion$TCX[i],intersect,metaAD_DCGs.filtered$Mayo_TCX)
  overlapping_genes_size=lapply(overlapping_genes,length)
  olp_KEGG=lapply(lapply(lapply(overlapping_genes,gprofiler,organism = "hsapiens",significant = T,src_filter = "KEGG",correction_method = "bonferroni",max_p_value = 0.05),`[`,c(9,12)),function(x)paste(x[,1],x[,2],sep = "--",collapse = ";"))
  olp_KEGG[lapply(olp_KEGG,function(x)x=="")==T]="No enriched KEGG pathways found!"
  fisher_results=DCG_representation=vector(mode="list",length = length(overlapping_genes_size))
  for(m in 1:length(overlapping_genes_size)){
    fisher_results[[m]]=fisher.test(matrix(c(overlapping_genes_size[[m]],
                                             length(metaAD_DCGs.filtered$Mayo_TCX)-overlapping_genes_size[[m]],
                                             module_size[[m]]-overlapping_genes_size[[m]],
                                             15160-(length(metaAD_DCGs.filtered$Mayo_TCX)-overlapping_genes_size[[m]])-length(module_size[[m]]-overlapping_genes_size[[m]])-overlapping_genes_size[[m]]),nrow = 2))
  }
  
  exp_DCG=as.list((length(metaAD_DCGs.filtered$Mayo_TCX)/15160)*unlist(module_size))
  for(m in 1:length(overlapping_genes_size)){
    if(overlapping_genes_size[[m]]>exp_DCG[[m]]){
      DCG_representation[[m]]="Over"  
    }
    else {DCG_representation[[m]]="Under" }
  } 
  data.frame(Jcd=unlist(jaccard_DCG),
             ModuleSize=unlist(module_size),
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes,length)),
             Expected_DCG_in_Module=unlist(exp_DCG),
             Nr_Of_OlpKEGG=unlist(olp_KEGG),
             p.overlap=unlist(lapply(fisher_results,function(p)p$p.value)),
             DCG_Representation=unlist(DCG_representation),
             stringsAsFactors = F)
}
comp_DCG_TCX.df=data.frame(rbindlist(comp_DCG_TCX.res),stringsAsFactors = F)
rownames(comp_DCG_TCX.df)=names(ampAD_AggModules_data_brainRegion$TCX)
write.table(comp_DCG_TCX.df,"AggModules_DCG_TCX_df.txt",sep="\t",col.names = T,row.names = T,quote = F)

comp_DCG_CER.res=foreach(i=1:length(ampAD_AggModules_data_brainRegion$CBE))%dopar%{
  jaccard_DCG=lapply(ampAD_AggModules_data_brainRegion$CBE[i],jaccard,metaAD_DCGs.filtered$Mayo_CER)
  module_size=lapply(ampAD_AggModules_data_brainRegion$CBE[i],length)
  overlapping_genes=lapply(ampAD_AggModules_data_brainRegion$CBE[i],intersect,metaAD_DCGs.filtered$Mayo_CER)
  overlapping_genes_size=lapply(overlapping_genes,length)
  olp_KEGG=lapply(lapply(lapply(overlapping_genes,gprofiler,organism = "hsapiens",significant = T,src_filter = "KEGG",correction_method = "bonferroni",max_p_value = 0.05),`[`,c(9,12)),function(x)paste(x[,1],x[,2],sep = "--",collapse = ";"))
  olp_KEGG[lapply(olp_KEGG,function(x)x=="")==T]="No enriched KEGG pathways found!"
  fisher_results=DCG_representation=vector(mode="list",length = length(overlapping_genes_size))
  for(m in 1:length(overlapping_genes_size)){
    fisher_results[[m]]=fisher.test(matrix(c(overlapping_genes_size[[m]],
                                             length(metaAD_DCGs.filtered$Mayo_CER)-overlapping_genes_size[[m]],
                                             module_size[[m]]-overlapping_genes_size[[m]],
                                             15160-(length(metaAD_DCGs.filtered$Mayo_CER)-overlapping_genes_size[[m]])-length(module_size[[m]]-overlapping_genes_size[[m]])-overlapping_genes_size[[m]]),nrow = 2))
  }
  
  exp_DCG=as.list((length(metaAD_DCGs.filtered$Mayo_CER)/15160)*unlist(module_size))
  for(m in 1:length(overlapping_genes_size)){
    if(overlapping_genes_size[[m]]>exp_DCG[[m]]){
      DCG_representation[[m]]="Over"  
    }
    else {DCG_representation[[m]]="Under" }
  } 
  data.frame(Jcd=unlist(jaccard_DCG),
             ModuleSize=unlist(module_size),
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes,length)),
             Expected_DCG_in_Module=unlist(exp_DCG),
             Nr_Of_OlpKEGG=unlist(olp_KEGG),
             p.overlap=unlist(lapply(fisher_results,function(p)p$p.value)),
             DCG_Representation=unlist(DCG_representation),
             stringsAsFactors = F)
}
comp_DCG_CER.df=data.frame(rbindlist(comp_DCG_CER.res),stringsAsFactors = F)
rownames(comp_DCG_CER.df)=names(ampAD_AggModules_data_brainRegion$CBE)
write.table(comp_DCG_CER.df,"AggModules_DCG_CER_df.txt",sep="\t",col.names = T,row.names = T,quote = F)

comp_DCG_FP.res=foreach(i=1:length(ampAD_AggModules_data_brainRegion$FP))%dopar%{
  jaccard_DCG=lapply(ampAD_AggModules_data_brainRegion$FP[i],jaccard,metaAD_DCGs.filtered$MSBB_FP)
  module_size=lapply(ampAD_AggModules_data_brainRegion$FP[i],length)
  overlapping_genes=lapply(ampAD_AggModules_data_brainRegion$FP[i],intersect,metaAD_DCGs.filtered$MSBB_FP)
  overlapping_genes_size=lapply(overlapping_genes,length)
  olp_KEGG=lapply(lapply(lapply(overlapping_genes,gprofiler,organism = "hsapiens",significant = T,src_filter = "KEGG",correction_method = "bonferroni",max_p_value = 0.05),`[`,c(9,12)),function(x)paste(x[,1],x[,2],sep = "--",collapse = ";"))
  olp_KEGG[lapply(olp_KEGG,function(x)x=="")==T]="No enriched KEGG pathways found!"
  fisher_results=DCG_representation=vector(mode="list",length = length(overlapping_genes_size))
  for(m in 1:length(overlapping_genes_size)){
    fisher_results[[m]]=fisher.test(matrix(c(overlapping_genes_size[[m]],
                                             length(metaAD_DCGs.filtered$MSBB_FP)-overlapping_genes_size[[m]],
                                             module_size[[m]]-overlapping_genes_size[[m]],
                                             14681-(length(metaAD_DCGs.filtered$MSBB_FP)-overlapping_genes_size[[m]])-length(module_size[[m]]-overlapping_genes_size[[m]])-overlapping_genes_size[[m]]),nrow = 2))
  }
  
  exp_DCG=as.list((length(metaAD_DCGs.filtered$MSBB_FP)/14681)*unlist(module_size))
  for(m in 1:length(overlapping_genes_size)){
    if(overlapping_genes_size[[m]]>exp_DCG[[m]]){
      DCG_representation[[m]]="Over"  
    }
    else {DCG_representation[[m]]="Under" }
  } 
  data.frame(Jcd=unlist(jaccard_DCG),
             ModuleSize=unlist(module_size),
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes,length)),
             Expected_DCG_in_Module=unlist(exp_DCG),
             Nr_Of_OlpKEGG=unlist(olp_KEGG),
             p.overlap=unlist(lapply(fisher_results,function(p)p$p.value)),
             DCG_Representation=unlist(DCG_representation),
             stringsAsFactors = F)
}
comp_DCG_FP.df=data.frame(rbindlist(comp_DCG_FP.res),stringsAsFactors = F)
rownames(comp_DCG_FP.df)=names(ampAD_AggModules_data_brainRegion$FP)
write.table(comp_DCG_FP.df,"AggModules_DCG_FP_df.txt",sep="\t",col.names = T,row.names = T,quote = F)

comp_DCG_IFG.res=foreach(i=1:length(ampAD_AggModules_data_brainRegion$IFG))%dopar%{
  jaccard_DCG=lapply(ampAD_AggModules_data_brainRegion$IFG[i],jaccard,metaAD_DCGs.filtered$MSBB_IFG)
  module_size=lapply(ampAD_AggModules_data_brainRegion$IFG[i],length)
  overlapping_genes=lapply(ampAD_AggModules_data_brainRegion$IFG[i],intersect,metaAD_DCGs.filtered$MSBB_IFG)
  overlapping_genes_size=lapply(overlapping_genes,length)
  olp_KEGG=lapply(lapply(lapply(overlapping_genes,gprofiler,organism = "hsapiens",significant = T,src_filter = "KEGG",correction_method = "bonferroni",max_p_value = 0.05),`[`,c(9,12)),function(x)paste(x[,1],x[,2],sep = "--",collapse = ";"))
  olp_KEGG[lapply(olp_KEGG,function(x)x=="")==T]="No enriched KEGG pathways found!"
  fisher_results=DCG_representation=vector(mode="list",length = length(overlapping_genes_size))
  for(m in 1:length(overlapping_genes_size)){
    fisher_results[[m]]=fisher.test(matrix(c(overlapping_genes_size[[m]],
                                             length(metaAD_DCGs.filtered$MSBB_IFG)-overlapping_genes_size[[m]],
                                             module_size[[m]]-overlapping_genes_size[[m]],
                                             14681-(length(metaAD_DCGs.filtered$MSBB_IFG)-overlapping_genes_size[[m]])-length(module_size[[m]]-overlapping_genes_size[[m]])-overlapping_genes_size[[m]]),nrow = 2))
  }
  
  exp_DCG=as.list((length(metaAD_DCGs.filtered$MSBB_IFG)/14681)*unlist(module_size))
  for(m in 1:length(overlapping_genes_size)){
    if(overlapping_genes_size[[m]]>exp_DCG[[m]]){
      DCG_representation[[m]]="Over"  
    }
    else {DCG_representation[[m]]="Under" }
  } 
  data.frame(Jcd=unlist(jaccard_DCG),
             ModuleSize=unlist(module_size),
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes,length)),
             Expected_DCG_in_Module=unlist(exp_DCG),
             Nr_Of_OlpKEGG=unlist(olp_KEGG),
             p.overlap=unlist(lapply(fisher_results,function(p)p$p.value)),
             DCG_Representation=unlist(DCG_representation),
             stringsAsFactors = F)
}
comp_DCG_IFG.df=data.frame(rbindlist(comp_DCG_IFG.res),stringsAsFactors = F)
rownames(comp_DCG_IFG.df)=names(ampAD_AggModules_data_brainRegion$IFG)
write.table(comp_DCG_IFG.df,"AggModules_DCG_IFG_df.txt",sep="\t",col.names = T,row.names = T,quote = F)

comp_DCG_PHG.res=foreach(i=1:length(ampAD_AggModules_data_brainRegion$PHG))%dopar%{
  jaccard_DCG=lapply(ampAD_AggModules_data_brainRegion$PHG[i],jaccard,metaAD_DCGs.filtered$MSBB_PHG)
  module_size=lapply(ampAD_AggModules_data_brainRegion$PHG[i],length)
  overlapping_genes=lapply(ampAD_AggModules_data_brainRegion$PHG[i],intersect,metaAD_DCGs.filtered$MSBB_PHG)
  overlapping_genes_size=lapply(overlapping_genes,length)
  olp_KEGG=lapply(lapply(lapply(overlapping_genes,gprofiler,organism = "hsapiens",significant = T,src_filter = "KEGG",correction_method = "bonferroni",max_p_value = 0.05),`[`,c(9,12)),function(x)paste(x[,1],x[,2],sep = "--",collapse = ";"))
  olp_KEGG[lapply(olp_KEGG,function(x)x=="")==T]="No enriched KEGG pathways found!"
  fisher_results=DCG_representation=vector(mode="list",length = length(overlapping_genes_size))
  for(m in 1:length(overlapping_genes_size)){
    fisher_results[[m]]=fisher.test(matrix(c(overlapping_genes_size[[m]],
                                             length(metaAD_DCGs.filtered$MSBB_PHG)-overlapping_genes_size[[m]],
                                             module_size[[m]]-overlapping_genes_size[[m]],
                                             14681-(length(metaAD_DCGs.filtered$MSBB_PHG)-overlapping_genes_size[[m]])-length(module_size[[m]]-overlapping_genes_size[[m]])-overlapping_genes_size[[m]]),nrow = 2))
  }
  
  exp_DCG=as.list((length(metaAD_DCGs.filtered$MSBB_PHG)/14681)*unlist(module_size))
  for(m in 1:length(overlapping_genes_size)){
    if(overlapping_genes_size[[m]]>exp_DCG[[m]]){
      DCG_representation[[m]]="Over"  
    }
    else {DCG_representation[[m]]="Under" }
  } 
  data.frame(Jcd=unlist(jaccard_DCG),
             ModuleSize=unlist(module_size),
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes,length)),
             Expected_DCG_in_Module=unlist(exp_DCG),
             Nr_Of_OlpKEGG=unlist(olp_KEGG),
             p.overlap=unlist(lapply(fisher_results,function(p)p$p.value)),
             DCG_Representation=unlist(DCG_representation),
             stringsAsFactors = F)
}
comp_DCG_PHG.df=data.frame(rbindlist(comp_DCG_PHG.res),stringsAsFactors = F)
rownames(comp_DCG_PHG.df)=names(ampAD_AggModules_data_brainRegion$PHG)
write.table(comp_DCG_PHG.df,"AggModules_DCG_PHG_df.txt",sep="\t",col.names = T,row.names = T,quote = F)

comp_DCG_STG.res=foreach(i=1:length(ampAD_AggModules_data_brainRegion$STG))%dopar%{
  jaccard_DCG=lapply(ampAD_AggModules_data_brainRegion$STG[i],jaccard,metaAD_DCGs.filtered$MSBB_STG)
  module_size=lapply(ampAD_AggModules_data_brainRegion$STG[i],length)
  overlapping_genes=lapply(ampAD_AggModules_data_brainRegion$STG[i],intersect,metaAD_DCGs.filtered$MSBB_STG)
  overlapping_genes_size=lapply(overlapping_genes,length)
  olp_KEGG=lapply(lapply(lapply(overlapping_genes,gprofiler,organism = "hsapiens",significant = T,src_filter = "KEGG",correction_method = "bonferroni",max_p_value = 0.05),`[`,c(9,12)),function(x)paste(x[,1],x[,2],sep = "--",collapse = ";"))
  olp_KEGG[lapply(olp_KEGG,function(x)x=="")==T]="No enriched KEGG pathways found!"
  fisher_results=DCG_representation=vector(mode="list",length = length(overlapping_genes_size))
  for(m in 1:length(overlapping_genes_size)){
    fisher_results[[m]]=fisher.test(matrix(c(overlapping_genes_size[[m]],
                                             length(metaAD_DCGs.filtered$MSBB_STG)-overlapping_genes_size[[m]],
                                             module_size[[m]]-overlapping_genes_size[[m]],
                                             14681-(length(metaAD_DCGs.filtered$MSBB_STG)-overlapping_genes_size[[m]])-length(module_size[[m]]-overlapping_genes_size[[m]])-overlapping_genes_size[[m]]),nrow = 2))
  }
  
  exp_DCG=as.list((length(metaAD_DCGs.filtered$MSBB_STG)/14681)*unlist(module_size))
  for(m in 1:length(overlapping_genes_size)){
    if(overlapping_genes_size[[m]]>exp_DCG[[m]]){
      DCG_representation[[m]]="Over"  
    }
    else {DCG_representation[[m]]="Under" }
  } 
  data.frame(Jcd=unlist(jaccard_DCG),
             ModuleSize=unlist(module_size),
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes,length)),
             Expected_DCG_in_Module=unlist(exp_DCG),
             Nr_Of_OlpKEGG=unlist(olp_KEGG),
             p.overlap=unlist(lapply(fisher_results,function(p)p$p.value)),
             DCG_Representation=unlist(DCG_representation),
             stringsAsFactors = F)
}
comp_DCG_STG.df=data.frame(rbindlist(comp_DCG_STG.res),stringsAsFactors = F)
rownames(comp_DCG_STG.df)=names(ampAD_AggModules_data_brainRegion$STG)
write.table(comp_DCG_STG.df,"AggModules_DCG_STG_df.txt",sep="\t",col.names = T,row.names = T,quote = F)


AggModules_DCG_Comparison=vector(mode = "list",length = 7)
names(AggModules_DCG_Comparison)=c("DLPFC","TCX","CER","FP","IFG","PHG","STG")
AggModules_DCG_Comparison$DLPFC=comp_DCG_DLPFC.df
AggModules_DCG_Comparison$TCX=comp_DCG_TCX.df
AggModules_DCG_Comparison$CER=comp_DCG_CER.df
AggModules_DCG_Comparison$FP=comp_DCG_FP.df
AggModules_DCG_Comparison$IFG=comp_DCG_IFG.df
AggModules_DCG_Comparison$PHG=comp_DCG_PHG.df
AggModules_DCG_Comparison$STG=comp_DCG_STG.df

for(r in 1:7){
  AggModules_DCG_Comparison[[r]]=AggModules_DCG_Comparison[[r]]%>%mutate(BrainRegion=names(AggModules_DCG_Comparison)[r])
}




# Compare consensus modules to DCGs and DCLs
ConsensusModules=ConsensusUniqueModules=vector(mode = "list",length = 7)
ConsensusModules_df=vector(mode = "list",length = 7)
names(ConsensusModules)=names(ConsensusModules_df)=names(ConsensusUniqueModules)=c("DLPFC","FP","IFG","PHG","STG","CER","TCX")
ConsensusModules_df$DLPFC=read.csv("../ConsensusModules/consensusDLPFC.csv",header = T,as.is = T,stringsAsFactors = F)
ConsensusModules_df$CER=read.csv("../ConsensusModules/consensusCBE.csv",header = T,as.is = T,stringsAsFactors = F)
ConsensusModules_df$TCX=read.csv("../ConsensusModules/consensusTCX.csv",header = T,as.is = T,stringsAsFactors = F)
ConsensusModules_df$FP=read.csv("../ConsensusModules/consensusFP.csv",header = T,as.is = T,stringsAsFactors = F)
ConsensusModules_df$IFG=read.csv("../ConsensusModules/consensusIFG.csv",header = T,as.is = T,stringsAsFactors = F)
ConsensusModules_df$PHG=read.csv("../ConsensusModules/consensusPHG.csv",header = T,as.is = T,stringsAsFactors = F)
ConsensusModules_df$STG=read.csv("../ConsensusModules/consensusSTG.csv",header = T,as.is = T,stringsAsFactors = F)

ConsensusUniqueModules$DLPFC=ConsensusUniqueModules$CBE=ConsensusUniqueModules$TCX=ConsensusUniqueModules$FP=ConsensusUniqueModules$IFG=ConsensusUniqueModules$PHG=ConsensusUniqueModules$STG=list()

for(i in 1:7){
  ConsensusUniqueModules[[i]]=lapply(lapply(lapply(unique(ConsensusModules_df[[i]]$moduleLabel),function(x)ConsensusModules_df[[i]]%>%filter(moduleLabel==x)%>%select(Gene.ID)),`[[`,1),function(y)mapIds2(IDs = y,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[1]][,2])
  names(ConsensusUniqueModules[[i]])=unique(ConsensusModules_df[[i]]$moduleLabel)
}

comp2_DCG_DLPFC.res=foreach(i=1:length(ConsensusUniqueModules$DLPFC))%dopar%{
  jaccard_DCG=lapply(ConsensusUniqueModules$DLPFC[i],jaccard,metaAD_DCGs.filtered$ROSMAP_DLPFC)
  module_size=lapply(ConsensusUniqueModules$DLPFC[i],length)
  overlapping_genes=lapply(ConsensusUniqueModules$DLPFC[i],intersect,metaAD_DCGs.filtered$ROSMAP_DLPFC)
  overlapping_genes_size=lapply(overlapping_genes,length)
  olp_KEGG=lapply(lapply(lapply(overlapping_genes,gprofiler,organism = "hsapiens",significant = T,src_filter = "KEGG",correction_method = "bonferroni",max_p_value = 0.05),`[`,c(9,12)),function(x)paste(x[,1],x[,2],sep = "--",collapse = ";"))
  olp_KEGG[lapply(olp_KEGG,function(x)x=="")==T]="No enriched KEGG pathways found!"
  fisher_results=DCG_representation=vector(mode="list",length = length(overlapping_genes_size))
  for(m in 1:length(overlapping_genes_size)){
    fisher_results[[m]]=fisher.test(matrix(c(overlapping_genes_size[[m]],
                                             length(metaAD_DCGs.filtered$ROSMAP_DLPFC)-overlapping_genes_size[[m]],
                                             module_size[[m]]-overlapping_genes_size[[m]],
                                             14228-(length(metaAD_DCGs.filtered$ROSMAP_DLPFC)-overlapping_genes_size[[m]])-length(module_size[[m]]-overlapping_genes_size[[m]])-overlapping_genes_size[[m]]),nrow = 2))
  }
  
  exp_DCG=as.list((length(metaAD_DCGs.filtered$ROSMAP_DLPFC)/14228)*unlist(module_size))
  for(m in 1:length(overlapping_genes_size)){
    if(overlapping_genes_size[[m]]>exp_DCG[[m]]){
      DCG_representation[[m]]="Over"  
    }
    else {DCG_representation[[m]]="Under" }
  } 
  data.frame(Jcd=unlist(jaccard_DCG),
             ModuleSize=unlist(module_size),
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes,length)),
             Expected_DCG_in_Module=unlist(exp_DCG),
             Nr_Of_OlpKEGG=unlist(olp_KEGG),
             p.overlap=unlist(lapply(fisher_results,function(p)p$p.value)),
             DCG_Representation=unlist(DCG_representation),
             stringsAsFactors = F)
}
comp2_DCG_DLPFC.df=data.frame(rbindlist(comp2_DCG_DLPFC.res),stringsAsFactors = F)
rownames(comp2_DCG_DLPFC.df)=names(ConsensusUniqueModules$DLPFC)
write.table(comp2_DCG_DLPFC.df,"ConsensusModules_DCG_DLPFC_df.txt",sep="\t",col.names = T,row.names = T,quote = F)


comp2_DCG_TCX.res=foreach(i=1:length(ConsensusUniqueModules$TCX))%dopar%{
  jaccard_DCG=lapply(ConsensusUniqueModules$TCX[i],jaccard,metaAD_DCGs.filtered$Mayo_TCX)
  module_size=lapply(ConsensusUniqueModules$TCX[i],length)
  overlapping_genes=lapply(ConsensusUniqueModules$TCX[i],intersect,metaAD_DCGs.filtered$Mayo_TCX)
  overlapping_genes_size=lapply(overlapping_genes,length)
  olp_KEGG=lapply(lapply(lapply(overlapping_genes,gprofiler,organism = "hsapiens",significant = T,src_filter = "KEGG",correction_method = "bonferroni",max_p_value = 0.05),`[`,c(9,12)),function(x)paste(x[,1],x[,2],sep = "--",collapse = ";"))
  olp_KEGG[lapply(olp_KEGG,function(x)x=="")==T]="No enriched KEGG pathways found!"
  fisher_results=DCG_representation=vector(mode="list",length = length(overlapping_genes_size))
  for(m in 1:length(overlapping_genes_size)){
    fisher_results[[m]]=fisher.test(matrix(c(overlapping_genes_size[[m]],
                                             length(metaAD_DCGs.filtered$Mayo_TCX)-overlapping_genes_size[[m]],
                                             module_size[[m]]-overlapping_genes_size[[m]],
                                             14228-(length(metaAD_DCGs.filtered$Mayo_TCX)-overlapping_genes_size[[m]])-length(module_size[[m]]-overlapping_genes_size[[m]])-overlapping_genes_size[[m]]),nrow = 2))
  }
  
  exp_DCG=as.list((length(metaAD_DCGs.filtered$Mayo_TCX)/14228)*unlist(module_size))
  for(m in 1:length(overlapping_genes_size)){
    if(overlapping_genes_size[[m]]>exp_DCG[[m]]){
      DCG_representation[[m]]="Over"  
    }
    else {DCG_representation[[m]]="Under" }
  } 
  data.frame(Jcd=unlist(jaccard_DCG),
             ModuleSize=unlist(module_size),
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes,length)),
             Expected_DCG_in_Module=unlist(exp_DCG),
             Nr_Of_OlpKEGG=unlist(olp_KEGG),
             p.overlap=unlist(lapply(fisher_results,function(p)p$p.value)),
             DCG_Representation=unlist(DCG_representation),
             stringsAsFactors = F)
}
comp2_DCG_TCX.df=data.frame(rbindlist(comp2_DCG_TCX.res),stringsAsFactors = F)
rownames(comp2_DCG_TCX.df)=names(ConsensusUniqueModules$TCX)
write.table(comp2_DCG_TCX.df[order(comp2_DCG_TCX.df$p.overlap,decreasing = F),],"ConsensusModules_DCG_TCX_df.txt",sep="\t",col.names = T,row.names = T,quote = F)


comp2_DCG_CER.res=foreach(i=1:length(ConsensusUniqueModules$CER))%dopar%{
  jaccard_DCG=lapply(ConsensusUniqueModules$CER[i],jaccard,metaAD_DCGs.filtered$Mayo_CER)
  module_size=lapply(ConsensusUniqueModules$CER[i],length)
  overlapping_genes=lapply(ConsensusUniqueModules$CER[i],intersect,metaAD_DCGs.filtered$Mayo_CER)
  overlapping_genes_size=lapply(overlapping_genes,length)
  olp_KEGG=lapply(lapply(lapply(overlapping_genes,gprofiler,organism = "hsapiens",significant = T,src_filter = "KEGG",correction_method = "bonferroni",max_p_value = 0.05),`[`,c(9,12)),function(x)paste(x[,1],x[,2],sep = "--",collapse = ";"))
  olp_KEGG[lapply(olp_KEGG,function(x)x=="")==T]="No enriched KEGG pathways found!"
  fisher_results=DCG_representation=vector(mode="list",length = length(overlapping_genes_size))
  for(m in 1:length(overlapping_genes_size)){
    fisher_results[[m]]=fisher.test(matrix(c(overlapping_genes_size[[m]],
                                             length(metaAD_DCGs.filtered$Mayo_CER)-overlapping_genes_size[[m]],
                                             module_size[[m]]-overlapping_genes_size[[m]],
                                             14228-(length(metaAD_DCGs.filtered$Mayo_CER)-overlapping_genes_size[[m]])-length(module_size[[m]]-overlapping_genes_size[[m]])-overlapping_genes_size[[m]]),nrow = 2))
  }
  
  exp_DCG=as.list((length(metaAD_DCGs.filtered$Mayo_CER)/14228)*unlist(module_size))
  for(m in 1:length(overlapping_genes_size)){
    if(overlapping_genes_size[[m]]>exp_DCG[[m]]){
      DCG_representation[[m]]="Over"  
    }
    else {DCG_representation[[m]]="Under" }
  } 
  data.frame(Jcd=unlist(jaccard_DCG),
             ModuleSize=unlist(module_size),
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes,length)),
             Expected_DCG_in_Module=unlist(exp_DCG),
             Nr_Of_OlpKEGG=unlist(olp_KEGG),
             p.overlap=unlist(lapply(fisher_results,function(p)p$p.value)),
             DCG_Representation=unlist(DCG_representation),
             stringsAsFactors = F)
}
comp2_DCG_CER.df=data.frame(rbindlist(comp2_DCG_CER.res),stringsAsFactors = F)
rownames(comp2_DCG_CER.df)=names(ConsensusUniqueModules$CER)
write.table(comp2_DCG_CER.df[order(comp2_DCG_CER.df$p.overlap,decreasing = F),],"ConsensusModules_DCG_CER_df.txt",sep="\t",col.names = T,row.names = T,quote = F)

comp2_DCG_FP.res=foreach(i=1:length(ConsensusUniqueModules$FP))%dopar%{
  jaccard_DCG=lapply(ConsensusUniqueModules$FP[i],jaccard,metaAD_DCGs.filtered$MSBB_FP)
  module_size=lapply(ConsensusUniqueModules$FP[i],length)
  overlapping_genes=lapply(ConsensusUniqueModules$FP[i],intersect,metaAD_DCGs.filtered$MSBB_FP)
  overlapping_genes_size=lapply(overlapping_genes,length)
  olp_KEGG=lapply(lapply(lapply(overlapping_genes,gprofiler,organism = "hsapiens",significant = T,src_filter = "KEGG",correction_method = "bonferroni",max_p_value = 0.05),`[`,c(9,12)),function(x)paste(x[,1],x[,2],sep = "--",collapse = ";"))
  olp_KEGG[lapply(olp_KEGG,function(x)x=="")==T]="No enriched KEGG pathways found!"
  fisher_results=DCG_representation=vector(mode="list",length = length(overlapping_genes_size))
  for(m in 1:length(overlapping_genes_size)){
    fisher_results[[m]]=fisher.test(matrix(c(overlapping_genes_size[[m]],
                                             length(metaAD_DCGs.filtered$MSBB_FP)-overlapping_genes_size[[m]],
                                             module_size[[m]]-overlapping_genes_size[[m]],
                                             14228-(length(metaAD_DCGs.filtered$MSBB_FP)-overlapping_genes_size[[m]])-length(module_size[[m]]-overlapping_genes_size[[m]])-overlapping_genes_size[[m]]),nrow = 2))
  }
  
  exp_DCG=as.list((length(metaAD_DCGs.filtered$MSBB_FP)/14228)*unlist(module_size))
  for(m in 1:length(overlapping_genes_size)){
    if(overlapping_genes_size[[m]]>exp_DCG[[m]]){
      DCG_representation[[m]]="Over"  
    }
    else {DCG_representation[[m]]="Under" }
  } 
  data.frame(Jcd=unlist(jaccard_DCG),
             ModuleSize=unlist(module_size),
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes,length)),
             Expected_DCG_in_Module=unlist(exp_DCG),
             Nr_Of_OlpKEGG=unlist(olp_KEGG),
             p.overlap=unlist(lapply(fisher_results,function(p)p$p.value)),
             DCG_Representation=unlist(DCG_representation),
             stringsAsFactors = F)
}
comp2_DCG_FP.df=data.frame(rbindlist(comp2_DCG_FP.res),stringsAsFactors = F)
rownames(comp2_DCG_FP.df)=names(ConsensusUniqueModules$FP)
write.table(comp2_DCG_FP.df[order(comp2_DCG_FP.df$p.overlap,decreasing = F),],"ConsensusModules_DCG_FP_df.txt",sep="\t",col.names = T,row.names = T,quote = F)

comp2_DCG_IFG.res=foreach(i=1:length(ConsensusUniqueModules$IFG))%dopar%{
  jaccard_DCG=lapply(ConsensusUniqueModules$IFG[i],jaccard,metaAD_DCGs.filtered$MSBB_IFG)
  module_size=lapply(ConsensusUniqueModules$IFG[i],length)
  overlapping_genes=lapply(ConsensusUniqueModules$IFG[i],intersect,metaAD_DCGs.filtered$MSBB_IFG)
  overlapping_genes_size=lapply(overlapping_genes,length)
  olp_KEGG=lapply(lapply(lapply(overlapping_genes,gprofiler,organism = "hsapiens",significant = T,src_filter = "KEGG",correction_method = "bonferroni",max_p_value = 0.05),`[`,c(9,12)),function(x)paste(x[,1],x[,2],sep = "--",collapse = ";"))
  olp_KEGG[lapply(olp_KEGG,function(x)x=="")==T]="No enriched KEGG pathways found!"
  fisher_results=DCG_representation=vector(mode="list",length = length(overlapping_genes_size))
  for(m in 1:length(overlapping_genes_size)){
    fisher_results[[m]]=fisher.test(matrix(c(overlapping_genes_size[[m]],
                                             length(metaAD_DCGs.filtered$MSBB_IFG)-overlapping_genes_size[[m]],
                                             module_size[[m]]-overlapping_genes_size[[m]],
                                             14228-(length(metaAD_DCGs.filtered$MSBB_IFG)-overlapping_genes_size[[m]])-length(module_size[[m]]-overlapping_genes_size[[m]])-overlapping_genes_size[[m]]),nrow = 2))
  }
  
  exp_DCG=as.list((length(metaAD_DCGs.filtered$MSBB_IFG)/14228)*unlist(module_size))
  for(m in 1:length(overlapping_genes_size)){
    if(overlapping_genes_size[[m]]>exp_DCG[[m]]){
      DCG_representation[[m]]="Over"  
    }
    else {DCG_representation[[m]]="Under" }
  } 
  data.frame(Jcd=unlist(jaccard_DCG),
             ModuleSize=unlist(module_size),
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes,length)),
             Expected_DCG_in_Module=unlist(exp_DCG),
             Nr_Of_OlpKEGG=unlist(olp_KEGG),
             p.overlap=unlist(lapply(fisher_results,function(p)p$p.value)),
             DCG_Representation=unlist(DCG_representation),
             stringsAsFactors = F)
}
comp2_DCG_IFG.df=data.frame(rbindlist(comp2_DCG_IFG.res),stringsAsFactors = F)
rownames(comp2_DCG_IFG.df)=names(ConsensusUniqueModules$IFG)
write.table(comp2_DCG_IFG.df[order(comp2_DCG_IFG.df$p.overlap,decreasing = F),],"ConsensusModules_DCG_IFG_df.txt",sep="\t",col.names = T,row.names = T,quote = F)


comp2_DCG_PHG.res=foreach(i=1:length(ConsensusUniqueModules$PHG))%dopar%{
  jaccard_DCG=lapply(ConsensusUniqueModules$PHG[i],jaccard,metaAD_DCGs.filtered$MSBB_PHG)
  module_size=lapply(ConsensusUniqueModules$PHG[i],length)
  overlapping_genes=lapply(ConsensusUniqueModules$PHG[i],intersect,metaAD_DCGs.filtered$MSBB_PHG)
  overlapping_genes_size=lapply(overlapping_genes,length)
  olp_KEGG=lapply(lapply(lapply(overlapping_genes,gprofiler,organism = "hsapiens",significant = T,src_filter = "KEGG",correction_method = "bonferroni",max_p_value = 0.05),`[`,c(9,12)),function(x)paste(x[,1],x[,2],sep = "--",collapse = ";"))
  olp_KEGG[lapply(olp_KEGG,function(x)x=="")==T]="No enriched KEGG pathways found!"
  fisher_results=DCG_representation=vector(mode="list",length = length(overlapping_genes_size))
  for(m in 1:length(overlapping_genes_size)){
    fisher_results[[m]]=fisher.test(matrix(c(overlapping_genes_size[[m]],
                                             length(metaAD_DCGs.filtered$MSBB_PHG)-overlapping_genes_size[[m]],
                                             module_size[[m]]-overlapping_genes_size[[m]],
                                             14228-(length(metaAD_DCGs.filtered$MSBB_PHG)-overlapping_genes_size[[m]])-length(module_size[[m]]-overlapping_genes_size[[m]])-overlapping_genes_size[[m]]),nrow = 2))
  }
  
  exp_DCG=as.list((length(metaAD_DCGs.filtered$MSBB_PHG)/14228)*unlist(module_size))
  for(m in 1:length(overlapping_genes_size)){
    if(overlapping_genes_size[[m]]>exp_DCG[[m]]){
      DCG_representation[[m]]="Over"  
    }
    else {DCG_representation[[m]]="Under" }
  } 
  data.frame(Jcd=unlist(jaccard_DCG),
             ModuleSize=unlist(module_size),
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes,length)),
             Expected_DCG_in_Module=unlist(exp_DCG),
             Nr_Of_OlpKEGG=unlist(olp_KEGG),
             p.overlap=unlist(lapply(fisher_results,function(p)p$p.value)),
             DCG_Representation=unlist(DCG_representation),
             stringsAsFactors = F)
}
comp2_DCG_PHG.df=data.frame(rbindlist(comp2_DCG_PHG.res),stringsAsFactors = F)
rownames(comp2_DCG_PHG.df)=names(ConsensusUniqueModules$PHG)
write.table(comp2_DCG_PHG.df[order(comp2_DCG_PHG.df$p.overlap,decreasing = F),],"ConsensusModules_DCG_PHG_df.txt",sep="\t",col.names = T,row.names = T,quote = F)


comp2_DCG_STG.res=foreach(i=1:length(ConsensusUniqueModules$STG))%dopar%{
  jaccard_DCG=lapply(ConsensusUniqueModules$STG[i],jaccard,metaAD_DCGs.filtered$MSBB_STG)
  module_size=lapply(ConsensusUniqueModules$STG[i],length)
  overlapping_genes=lapply(ConsensusUniqueModules$STG[i],intersect,metaAD_DCGs.filtered$MSBB_STG)
  overlapping_genes_size=lapply(overlapping_genes,length)
  olp_KEGG=lapply(lapply(lapply(overlapping_genes,gprofiler,organism = "hsapiens",significant = T,src_filter = "KEGG",correction_method = "bonferroni",max_p_value = 0.05),`[`,c(9,12)),function(x)paste(x[,1],x[,2],sep = "--",collapse = ";"))
  olp_KEGG[lapply(olp_KEGG,function(x)x=="")==T]="No enriched KEGG pathways found!"
  fisher_results=DCG_representation=vector(mode="list",length = length(overlapping_genes_size))
  for(m in 1:length(overlapping_genes_size)){
    fisher_results[[m]]=fisher.test(matrix(c(overlapping_genes_size[[m]],
                                             length(metaAD_DCGs.filtered$MSBB_STG)-overlapping_genes_size[[m]],
                                             module_size[[m]]-overlapping_genes_size[[m]],
                                             14228-(length(metaAD_DCGs.filtered$MSBB_STG)-overlapping_genes_size[[m]])-length(module_size[[m]]-overlapping_genes_size[[m]])-overlapping_genes_size[[m]]),nrow = 2))
  }
  
  exp_DCG=as.list((length(metaAD_DCGs.filtered$MSBB_STG)/14228)*unlist(module_size))
  for(m in 1:length(overlapping_genes_size)){
    if(overlapping_genes_size[[m]]>exp_DCG[[m]]){
      DCG_representation[[m]]="Over"  
    }
    else {DCG_representation[[m]]="Under" }
  } 
  data.frame(Jcd=unlist(jaccard_DCG),
             ModuleSize=unlist(module_size),
             Nr_OverlappedGenes=unlist(lapply(overlapping_genes,length)),
             Expected_DCG_in_Module=unlist(exp_DCG),
             Nr_Of_OlpKEGG=unlist(olp_KEGG),
             p.overlap=unlist(lapply(fisher_results,function(p)p$p.value)),
             DCG_Representation=unlist(DCG_representation),
             stringsAsFactors = F)
}
comp2_DCG_STG.df=data.frame(rbindlist(comp2_DCG_STG.res),stringsAsFactors = F)
rownames(comp2_DCG_STG.df)=names(ConsensusUniqueModules$STG)
write.table(comp2_DCG_STG.df[order(comp2_DCG_STG.df$p.overlap,decreasing = F),],"ConsensusModules_DCG_STG_df.txt",sep="\t",col.names = T,row.names = T,quote = F)

