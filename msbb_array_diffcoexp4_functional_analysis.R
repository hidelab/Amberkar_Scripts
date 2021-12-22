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

lateAD_diffcoexp_files=lateAD_diffcoexp_files[-11]
lateAD_samples.exprs=lateAD_samples.exprs[-12]
lateAD_samples=lateAD_samples[-12]

regnet_tf2target_data=fread("/Users/sandeepamberkar/Work/Data/TF_Databases/RegNetwork_human_regulators2.txt",header = T,sep = "\t",showProgress = T,data.table = F)%>%dplyr::filter(evidence=="Experimental")
regnet_tf2target.HGNC=fread("/Users/sandeepamberkar/Work/Data/TF_Databases/RegNetwork_human_regulators2.txt",header = T,sep = "\t",showProgress = T,data.table = F)%>%dplyr::select(c(regulator_symbol,target_symbol))

#DEG analyses to compare with DCGs
earlyAD_DEGs=lateAD_DEGs=vector(mode = "list",length = length(earlyAD_samples.exprs))
earlyAD_samples.exprs=earlyAD_samples.exprs[names(earlyAD_diffcoexp_results)]
earlyAD_samples=earlyAD_samples[names(earlyAD_diffcoexp_results)]
lateAD_samples.exprs=lateAD_samples.exprs[names(lateAD_diffcoexp_results)]
lateAD_samples=lateAD_samples[names(lateAD_diffcoexp_results)]
names(earlyAD_DEGs)=unlist(lapply(strsplit(split = "_earlyAD",basename(earlyAD_diffcoexp_files)),`[[`,1))
names(lateAD_DEGs)=unlist(lapply(strsplit(split = "_lateAD",basename(lateAD_diffcoexp_files)),`[[`,1))

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
  earlyAD_DEGs[[i]]=topTable(fit2,coef = 1,number = 19530,p.value = 0.2,adjust.method = "BH")%>%rownames_to_column("Gene")
}
for(i in 2:length(lateAD_DEGs)){
  c_exprs=lateAD_samples.exprs[[i]][,lateAD_samples[[i]]$SampleType=="CONTROL"]
  d_exprs=lateAD_samples.exprs[[i]][,lateAD_samples[[i]]$SampleType=="AD"]
  group=factor(c(rep("CONTROL",dim(c_exprs)[2]),rep("AD",dim(d_exprs)[2])))
  design_df=model.matrix(~0+group)
  colnames(design_df)=levels(group)
  contrasts_matrix=makeContrasts(CONTROL-AD,levels = design_df)
  fit=lmFit(object = cbind.data.frame(c_exprs,d_exprs),design = design_df)
  fit2=contrasts.fit(fit,contrasts_matrix)
  fit2=eBayes(fit2,trend = T)
  lateAD_DEGs[[i]]=topTable(fit2,coef = 1,number = 19530,p.value = 0.2,adjust.method = "BH")%>%rownames_to_column("Gene")
}


earlyAD_diffcoexp_results=vector(mode = "list",length = length(earlyAD_diffcoexp_files))
lateAD_diffcoexp_results=vector(mode = "list",length = length(lateAD_diffcoexp_files))
for(i in 1:length(earlyAD_diffcoexp_files)){
  earlyAD_diffcoexp_results[[i]]=readRDS(earlyAD_diffcoexp_files[[i]])
  lateAD_diffcoexp_results[[i]]=readRDS(lateAD_diffcoexp_files[[i]])
}
names(earlyAD_diffcoexp_results)=unlist(lapply(strsplit(split = "_earlyAD",basename(earlyAD_diffcoexp_files)),`[[`,1))
names(lateAD_diffcoexp_results)=unlist(lapply(strsplit(split = "_lateAD",basename(lateAD_diffcoexp_files)),`[[`,1))

earlyAD_DCGs=lapply(earlyAD_diffcoexp_results,function(x)x$DCGs)
earlyAD_DCGs.filtered=lapply(earlyAD_diffcoexp_results,function(x)x$DCGs%>%dplyr::filter(q<=0.1))
earlyAD_DCGs.filtered=lapply(earlyAD_DCGs.filtered,function(x){x$EnsemblID <- unname(mapIds(x = EnsDb.Hsapiens.v79,keys = x$Gene,keytype = "SYMBOL",symbol = "GENEID"));x[c(1,9,2:8)]})
earlyAD_DCGs.filtered=lapply(earlyAD_DCGs.filtered,function(x)x%>%dplyr::filter(!grepl(pattern = "///|\\-|\\.",x = Gene)))
earlyAD_DCLs=lapply(earlyAD_diffcoexp_results,function(x)x$DCLs)
earlyAD_DCLs.filtered=lapply(earlyAD_diffcoexp_results,function(x)x$DCLs%>%dplyr::filter(q.diffcor<=0.1))
earlyAD_DCLs.filtered=lapply(earlyAD_DCLs.filtered,function(x)x%>%dplyr::filter(!grepl(pattern = "///",x = `Gene.1`)&!grepl(pattern = "///",x = `Gene.2`)))
earlyAD_DCLs.filtered=lapply(earlyAD_DCLs.filtered,function(x){x$Ensembl.1 <- unname(mapIds(x = EnsDb.Hsapiens.v79,keys = x$Gene.1,keytype = "SYMBOL",symbol = "GENEID"));
                                                               x$Ensembl.2 <- unname(mapIds(x = EnsDb.Hsapiens.v79,keys = x$Gene.2,keytype = "SYMBOL",symbol = "GENEID"));x[c(1,13,2,14,3:12)]})
earlyAD_DCGs.list=lapply(earlyAD_DCGs.filtered,function(x)x%>%dplyr::filter(!grepl(pattern = "///|\\-|\\.",x = `Gene`))%>%pull(Gene))
earlyAD_DCLs.subgraph=lapply(earlyAD_DCLs.filtered,function(x)graph.data.frame(d = x[,c('Gene.1','Gene.2')],directed = F))
earlyAD_DCGs.Entrez_list=lapply(earlyAD_DCGs.list,function(x)unname(mapIds(x = org.Hs.eg.db,keys = x,column = "ENTREZID",keytype = "SYMBOL")))

lateAD_DCGs=lapply(lateAD_diffcoexp_results,function(x)x$DCGs)
lateAD_DCGs.filtered=lapply(lateAD_diffcoexp_results,function(x)x$DCGs%>%dplyr::filter(q<=0.1))
lateAD_DCGs.filtered=lapply(lateAD_DCGs.filtered,function(x){x$EnsemblID <- unname(mapIds(x = EnsDb.Hsapiens.v79,keys = x$Gene,keytype = "SYMBOL",symbol = "GENEID"));x[c(1,9,2:8)]})
lateAD_DCGs.filtered=lapply(lateAD_DCGs.filtered,function(x)x%>%dplyr::filter(!grepl(pattern = "///|\\-|\\.",x = Gene)))
lateAD_DCLs=lapply(lateAD_diffcoexp_results,function(x)x$DCLs)
lateAD_DCLs.filtered=lapply(lateAD_diffcoexp_results,function(x)x$DCLs%>%dplyr::filter(q.diffcor<=0.1))
lateAD_DCLs.filtered=lapply(earlyAD_DCLs.filtered,function(x)x%>%dplyr::filter(!grepl(pattern = "///",x = `Gene.1`)&!grepl(pattern = "///",x = `Gene.2`)))
lateAD_DCLs.filtered=lapply(earlyAD_DCLs.filtered,function(x){x$Ensembl.1 <- unname(mapIds(x = EnsDb.Hsapiens.v79,keys = x$Gene.1,keytype = "SYMBOL",symbol = "GENEID"));
                                                              x$Ensembl.2 <- unname(mapIds(x = EnsDb.Hsapiens.v79,keys = x$Gene.2,keytype = "SYMBOL",symbol = "GENEID"));x[c(1,13,2,14,3:12)]})


lateAD_DCGs.list=lapply(lateAD_DCGs.filtered,function(x)x%>%dplyr::filter(!grepl(pattern = "///|\\-|\\.",x = `Gene`))%>%pull(Gene))
lateAD_DCLs.subgraph=lapply(lateAD_DCLs.filtered,function(x)graph.data.frame(d = x[,c('Gene.1','Gene.2')],directed = F))
lateAD_DCGs.Entrez_list=lapply(lateAD_DCGs.list,function(x)unname(mapIds(x = org.Hs.eg.db,keys = x,column = "ENTREZID",keytype = "SYMBOL")))



common_DCGs=Filter(f = function(x)length(x)>0,x = mapply(FUN = function(a,b)intersect(a$Gene,b$Gene),earlyAD_DCGs.filtered,lateAD_DCGs.filtered))
for(i in 1:length(common_DCGs)){
  early_common_DCGs=earlyAD_DCGs.filtered[[i]][earlyAD_DCGs.filtered[[i]]$Gene%in%common_DCGs[[i]],]%>%arrange(Gene)
  late_common_DCGs=lateAD_DCGs.filtered[[i]][lateAD_DCGs.filtered[[i]]$Gene%in%common_DCGs[[i]],]%>%arrange(Gene)
  merge_common_DCGs=merge.data.frame(early_common_DCGs,late_common_DCGs,by = 'Gene')%>%dplyr::select(c(-10))
  colnames(merge_common_DCGs)[c(3:9)]=gsub(pattern = ".x",replacement = ".earlyAD",x = colnames(merge_common_DCGs)[c(3:9)])
  colnames(merge_common_DCGs)[c(10:16)]=gsub(pattern = ".y",replacement = ".lateAD",x = colnames(merge_common_DCGs)[c(10:16)])
  colnames(merge_common_DCGs)[2]="EnsemblID"
  fwrite(merge_common_DCGs,paste(names(common_DCGs)[i],"persistent_DCGs.txt",sep = "_"),sep = "\t",col.names = T,row.names = F)
  write(common_DCGs[[i]],paste(names(common_DCGs)[i],"persistent_DCGs_list.txt",sep = "_"),sep = "\n")
}


#Aging signature

aging_signature.Up=read.xls("../../../../../HAGR_Meta_Aging_Signature/btp073_Supplementary_Data/bioinf-2008-1490-File002.xls",sheet = 2)
aging_signature.Down=read.xls("../../../../../HAGR_Meta_Aging_Signature/btp073_Supplementary_Data/bioinf-2008-1490-File002.xls",sheet = 4)
aging_signature=c(aging_signature.Up%>%dplyr::filter(q.value<=0.1)%>%pull(Symbol)%>%droplevels%>%levels,aging_signature.Down%>%dplyr::filter(q.value<=0.1)%>%pull(Symbol)%>%droplevels%>%levels)

earlyAD_DCG_aging.overlap=lapply(earlyAD_DCGs.list,intersect,aging_signature)
lateAD_DCG_aging.overlap=lapply(lateAD_DCGs.list,intersect,aging_signature)

#Correct for aging and late AD persisting DCGs
earlyAD_DCGs.uniq=Filter(f = function(x)length(x)>0,x = mapply(FUN = function(a,b)setdiff(a$Gene,b$Gene),earlyAD_DCGs.filtered,lateAD_DCGs.filtered))
earlyAD_DCGs.uniq=Filter(f = function(x)length(x)>0,x = mapply(FUN = function(a,b)setdiff(a,b),earlyAD_DCGs.uniq,earlyAD_DCG_aging.overlap))
lateAD_DCGs.uniq=Filter(f = function(x)length(x)>0,x = mapply(FUN = function(a,b)setdiff(a$Gene,b$Gene),lateAD_DCGs.filtered,earlyAD_DCGs.filtered))
lateAD_DCGs.uniq=Filter(f = function(x)length(x)>0,x = mapply(FUN = function(a,b)setdiff(a,b),lateAD_DCGs.uniq,lateAD_DCG_aging.overlap))
lateAD_DCGs.uniq=lapply(lateAD_DCGs.uniq,function(x)grep(pattern = "///|\\-|\\.",x = x,invert = T,value = T))


earlyAD_DCGs.Uniq_df=lateAD_DCGs.Uniq_df=vector(mode = "list",length = 12)
names(earlyAD_DCGs.Uniq_df)=names(lateAD_DCGs.Uniq_df)=names(earlyAD_DCGs.filtered)
for(i in 1:length(earlyAD_DCGs.Uniq_df)){
  earlyAD_DCGs.Uniq_df[[i]]=earlyAD_DCGs.filtered[[i]][earlyAD_DCGs.filtered[[i]]$Gene%in%earlyAD_DCGs.uniq[[i]],]
  lateAD_DCGs.Uniq_df[[i]]=lateAD_DCGs.filtered[[i]][lateAD_DCGs.filtered[[i]]$Gene%in%lateAD_DCGs.uniq[[i]],]
  fwrite(earlyAD_DCGs.Uniq_df[[i]],paste("./EarlyAD_diffcoexp/EarlyAD",names(earlyAD_DCGs.Uniq_df)[i],"restricted_DCGs.txt",sep = "_"),sep = "\t",col.names = T,row.names = F)
  write(earlyAD_DCGs.Uniq_df[[i]]$Gene,paste("./EarlyAD_diffcoexp/EarlyAD",names(earlyAD_DCGs.Uniq_df)[i],"restricted_DCGs_list.txt",sep = "_"),sep = "\n")
  fwrite(lateAD_DCGs.Uniq_df[[i]],paste("./LateAD_diffcoexp/LateAD",names(lateAD_DCGs.Uniq_df)[i],"restricted_DCGs.txt",sep = "_"),sep = "\t",col.names = T,row.names = F)
  write(lateAD_DCGs.Uniq_df[[i]]$Gene,paste("./LateAD_diffcoexp/LateAD",names(lateAD_DCGs.Uniq_df)[i],"restricted_DCGs_list.txt",sep = "_"),sep = "\n")
}
  
earlyAD_janssenAD.fisher=lateAD_janssenAD.fisher=




random_common_DCGs.list=vector(mode = "list",length = length(names(common_DCGs)))
names(random_common_DCGs.list)=names(common_DCGs)
for(t in names(common_DCGs)){
  random_earlyAD_DCGs=replicate(n = 1000,expr = sample(x = rownames(earlyAD_samples.exprs[[t]]),size = length(earlyAD_diffcoexp_results[[t]]$DCGs$Gene),replace = F),simplify = T)
  random_lateAD_DCGs=replicate(n = 1000,expr = sample(x = rownames(lateAD_samples.exprs[[t]]),size = length(lateAD_diffcoexp_results[[t]]$DCGs$Gene),replace = F),simplify = T)
  random_common_DCGs.list[[t]]=unlist(foreach(i=1:1000)%do%{
    length(intersect(random_earlyAD_DCGs[,i],random_lateAD_DCGs[,i]))
  })  
}

# 
# earlyAD_DCGs=earlyAD_DCLs=earlyAD_DCGs.list=earlyAD_DCLs.list=earlyAD_DRsort=earlyAD_DRGs=earlyAD_DRGs.list=earlyAD_DRLs=earlyAD_DRLs.list=earlyAD_bridged_DCLs=vector(mode = "list",length = length(earlyAD_diffcoexp_results))
# names(earlyAD_DCGs)=names(earlyAD_DCLs)=names(earlyAD_DCGs.list)=names(earlyAD_DCLs.list)=names(earlyAD_DRsort)=names(earlyAD_DRGs)=names(earlyAD_DRGs.list)=names(earlyAD_DRLs)=names(earlyAD_DRLs.list)=names(earlyAD_bridged_DCLs)

common_DCGs.fisher=vector(mode = "list",length = length(earlyAD_DCGs))
for(i in 1:length(earlyAD_DCGs)){
  common_DCGs.fisher[[i]]=fisher.test(matrix(c(length(common_DCGs[[i]]),
                       length(earlyAD_DCGs.list[[i]])-length(common_DCGs[[i]]),
                       length(lateAD_DCGs.list[[i]])-length(common_DCGs[[i]]),
                       19529-(length(earlyAD_DCGs.list[[i]])-length(common_DCGs[[i]]))-(length(lateAD_DCGs.list[[i]])-length(common_DCGs[[i]]))-length(common_DCGs[[i]])),nrow = 2))
}
names(common_DCGs.fisher)=names(earlyAD_DCGs)
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
earlyAD_DCGs.Uniq_c7=lapply(earlyAD_DCGs.uniq,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c7),stringsAsFactors = F))
earlyAD_DCGs.Uniq_c2=lapply(earlyAD_DCGs.uniq,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c2),stringsAsFactors = F))
for(i in 1:length(earlyAD_DCGs.Uniq_c2)){
  fwrite(earlyAD_DCGs.Uniq_c2[[i]][grep(pattern = "KEGG|REACTOME",x = earlyAD_DCGs.Uniq_c2[[i]]$ID),-2],file = paste("./EarlyAD_diffcoexp/",names(earlyAD_DCGs.Uniq_c2)[i],"_earlyAD_restricted_enriched_pathways.txt",sep = ""),sep = "\t",col.names = T,row.names = F)
}
earlyAD_DCGs.Uniq_h=lapply(earlyAD_DCGs.uniq,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_h),stringsAsFactors = F))
earlyAD_DCGs.Uniq_biocarta=lapply(earlyAD_DCGs.uniq,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_biocarta),stringsAsFactors = F))
earlyAD_DCGs.Uniq_BP=lapply(earlyAD_DCGs.uniq,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c5.BP),stringsAsFactors = F))
earlyAD_DCGs.Uniq_CC=lapply(earlyAD_DCGs.uniq,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c5.CC),stringsAsFactors = F))
earlyAD_DCGs.Uniq_MF=lapply(earlyAD_DCGs.uniq,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c5.MF),stringsAsFactors = F))


earlyAD_DCGs.c7=lapply(earlyAD_DCGs.list,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c7),stringsAsFactors = F))
earlyAD_DCGs.c2=lapply(earlyAD_DCGs.list,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c2),stringsAsFactors = F))
earlyAD_DCGs.h=lapply(earlyAD_DCGs.list,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_h),stringsAsFactors = F))
earlyAD_DCGs.biocarta=lapply(earlyAD_DCGs.list,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_biocarta),stringsAsFactors = F))
earlyAD_DCGs.BP=lapply(earlyAD_DCGs.list,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c5.BP),stringsAsFactors = F))
earlyAD_DCGs.CC=lapply(earlyAD_DCGs.list,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c5.CC),stringsAsFactors = F))
earlyAD_DCGs.MF=lapply(earlyAD_DCGs.list,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c5.MF),stringsAsFactors = F))

#Late DCGs functional enrichment
lateAD_DCGs.c7=lapply(lateAD_DCGs.list,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c7),stringsAsFactors = F))
lateAD_DCGs.c2=lapply(lateAD_DCGs.list,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c2),stringsAsFactors = F))

lateAD_DCGs.h=lapply(lateAD_DCGs.list,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_h),stringsAsFactors = F))
lateAD_DCGs.biocarta=lapply(lateAD_DCGs.list,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_biocarta),stringsAsFactors = F))
lateAD_DCGs.BP=lapply(lateAD_DCGs.list,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c5.BP),stringsAsFactors = F))
lateAD_DCGs.CC=lapply(lateAD_DCGs.list,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c5.CC),stringsAsFactors = F))
lateAD_DCGs.MF=lapply(lateAD_DCGs.list,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c5.MF),stringsAsFactors = F))

lateAD_DCGs.Uniq_c7=lapply(lateAD_DCGs.uniq,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c7),stringsAsFactors = F))
lateAD_DCGs.Uniq_c2=lapply(lateAD_DCGs.uniq,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c2),stringsAsFactors = F))
for(i in 1:length(lateAD_DCGs.Uniq_c2)){
  
  fwrite(lateAD_DCGs.Uniq_c2[[i]][grep(pattern = "KEGG|REACTOME",x = lateAD_DCGs.Uniq_c2[[i]]$ID),-2],file = paste("./LateAD_diffcoexp/",names(lateAD_DCGs.Uniq_c2)[i],"_lateAD_restricted_enriched_pathways.txt",sep = ""),sep = "\t",col.names = T,row.names = F)
}
lateAD_DCGs.Uniq_h=lapply(lateAD_DCGs.uniq,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_h),stringsAsFactors = F))
lateAD_DCGs.Uniq_biocarta=lapply(lateAD_DCGs.uniq,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_biocarta),stringsAsFactors = F))
lateAD_DCGs.Uniq_BP=lapply(lateAD_DCGs.uniq,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c5.BP),stringsAsFactors = F))
lateAD_DCGs.Uniq_CC=lapply(lateAD_DCGs.uniq,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c5.CC),stringsAsFactors = F))
lateAD_DCGs.Uniq_MF=lapply(lateAD_DCGs.uniq,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c5.MF),stringsAsFactors = F))

common_DCGs_c7=lapply(common_DCGs,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c7),stringsAsFactors = F))
common_DCGs_c2=lapply(common_DCGs,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c2),stringsAsFactors = F))
common_DCGs_h=lapply(common_DCGs,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_h),stringsAsFactors = F))
common_DCGs_biocarta=lapply(common_DCGs,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_biocarta),stringsAsFactors = F))
common_DCGs_BP=lapply(common_DCGs,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c5.BP),stringsAsFactors = F))
common_DCGs_CC=lapply(common_DCGs,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c5.CC),stringsAsFactors = F))
common_DCGs_MF=lapply(common_DCGs,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c5.MF),stringsAsFactors = F))

DCG_KEGG_Summary.df=data.frame(rbindlist(list(EarlyAD_DCG_KEGG=lapply(earlyAD_DCGs.Uniq_c2,function(x)grep(pattern = "KEGG|REACTOME",x = x$ID,value = T)%>%paste(collapse = "; ")),
                              Common_DCG_KEGG=lapply(common_DCGs_c2,function(x)grep(pattern = "KEGG|REACTOME",x = x$ID,value = T)%>%paste(collapse = "; ")),
                              LateAD_DCG_KEGG=lapply(lateAD_DCGs.Uniq_c2,function(x)grep(pattern = "KEGG|REACTOME",x = x$ID,value = T)%>%paste(collapse = "; ")))))
DCG_KEGG_Summary.df=data.frame(t(DCG_KEGG_Summary.df),stringsAsFactors = F)
colnames(DCG_KEGG_Summary.df)=c("EarlyAD-specific-DCG-pathways","Persistent-DCG-pathways","LateAD-specific-DCG-pathways")
WriteXLS(DCG_KEGG_Summary.df,"../DCG_Pathway_Summary.xlsx",row.names = T,col.names = T)

earlyAD_DRsort=vector(mode = "list",length = length(earlyAD_DCGs.filtered))
for(t in 1:length(earlyAD_DCGs.filtered)){
  earlyAD_DRsort[[t]]=DRsort(DCGs = earlyAD_DCGs[[names(earlyAD_DCGs.filtered)[t]]],DCLs = earlyAD_DCLs[[names(earlyAD_DCLs.filtered)[t]]],tf2target = regnet_tf2target.HGNC,expGenes = rownames(earlyAD_samples.exprs[[t]]))
}

names(earlyAD_DRsort)=names(earlyAD_DCGs.filtered)

earlyAD_DRGs=lapply(earlyAD_DRsort,function(x)x$DRGs%>%dplyr::filter(q<=0.1&DCGisTF=="TRUE"))
earlyAD_DRGs.list=lapply(earlyAD_DRGs,function(x)x$DCG%>%droplevels)
earlyAD_DRLs=lapply(earlyAD_DRsort,function(x)x$DRLs%>%dplyr::filter(q.diffcor<=0.05))
earlyAD_DRLs.subgraph=lapply(earlyAD_DRLs,function(x)graph.data.frame(d = x[,c('Gene.1','Gene.2')],directed = F))
earlyAD_TF_bridged_DCL=lapply(earlyAD_DRsort,function(x)x$TF_bridged_DCL%>%dplyr::filter(q.diffcor<=0.01))

lateAD_DRsort=vector(mode = "list",length = length(lateAD_DCGs.filtered))
for(t in 1:length(lateAD_DCGs.filtered)){
  lateAD_DRsort[[t]]=DRsort(DCGs = lateAD_DCGs[[names(lateAD_DCGs.filtered)[t]]],DCLs = lateAD_DCLs[[names(lateAD_DCLs.filtered)[t]]],tf2target = regnet_tf2target.HGNC,expGenes = rownames(lateAD_samples.exprs[[t]]))
}

names(lateAD_DRsort)=names(lateAD_DCGs.filtered)

lateAD_DRGs=lapply(lateAD_DRsort,function(x)x$DRGs%>%dplyr::filter(q<=0.1&DCGisTF=="TRUE"))
lateAD_DRGs.list=lapply(lateAD_DRGs,function(x)x$DCG%>%droplevels%>%levels)
lateAD_DRLs=lapply(lateAD_DRsort,function(x)x$DRLs%>%dplyr::filter(q.diffcor<=0.05))
lateAD_DRLs.subgraph=lapply(lateAD_DRLs,function(x)graph.data.frame(d = x[,c('Gene.1','Gene.2')],directed = F))
lateAD_TF_bridged_DCL=lapply(lateAD_DRsort,function(x)x$TF_bridged_DCL%>%dplyr::filter(q.diffcor<=0.01))


# earlyAD_DCG.CompMatrix=matrix(NA,ncol=length(earlyAD_DCGs.list),nrow=length(earlyAD_DCGs.list))
# rownames(earlyAD_DCG.CompMatrix)=colnames(earlyAD_DCG.CompMatrix)=c("FP","OVC","ITG","MTG","STG","PCC","AC","PHG","TP","IFG","DLPFC","SPL","HIPP")
# for(i in 1:length(earlyAD_DCGs.list)){
#   earlyAD_DCG.CompMatrix[i,]=unlist(lapply(earlyAD_DCGs.list,function(x)length(intersect(x,earlyAD_DCGs.list[[i]]))))
# }
# diag(earlyAD_DCG.CompMatrix)=0
# chordDiagramFromMatrix(earlyAD_DCG.CompMatrix[-12,-12])




earlyAD_DRrank.TED=vector(mode = "list",length = length(earlyAD_DCGs.filtered))
earlyAD_DRrank_allTFs.TED=vector(mode = "list",length = length(earlyAD_DCGs.filtered))
lateAD_DRrank.TED=vector(mode = "list",length = length(lateAD_DCGs.filtered))
names(earlyAD_DRrank.TED)=names(earlyAD_DCGs.filtered)
names(earlyAD_DRrank_allTFs.TED)=names(earlyAD_DCGs.filtered)
names(lateAD_DRrank.TED)=names(lateAD_DCGs.filtered)

setwd("/Users/sandeepamberkar/Work/Data/MSMM-MSBB-HBTRC-Data/MSBB_Array19/GSE84422/Lobe_Analysis/Dyscorrelation_Landscape_Results/EarlyAD_diffcoexp/")
earlyAD_DRrank_TED.files=list.files(path = "EarlyAD_TED",pattern = "TED_new_pbinom",full.names = T)
earlyAD_DRrank_allTFs_TED.files=list.files(path = "EarlyAD_TED/All_TFs",pattern = "allTFs",full.names = T)

for(t in 1:length(earlyAD_DRrank_TED.files)){
  
  earlyAD_DRrank.TED[[t]]=readRDS(earlyAD_DRrank_TED.files[[t]])
  
}
for(t in 1:length(earlyAD_DRrank_allTFs_TED.files)){
  
  earlyAD_DRrank_allTFs.TED[[t]]=readRDS(earlyAD_DRrank_allTFs_TED.files[[t]])
  
}

earlyAD_DRrank_allTFs_TED.files=list.files(path = "EarlyAD_TED/All_TFs",pattern = "allTFs",full.names = T)
setwd("/Users/sandeepamberkar/Work/Data/MSMM-MSBB-HBTRC-Data/MSBB_Array19/GSE84422/Lobe_Analysis/Dyscorrelation_Landscape_Results/LateAD_diffcoexp/")
lateAD_DRrank_TED.files=list.files(path = "LateAD_TED",pattern = "TED_new_pbinom",full.names = T)
#Remove Superior Parietal Lobule region from the TED results
lateAD_DRrank_TED.files=lateAD_DRrank_TED.files[-11]
for(t in 1:length(lateAD_DRrank_TED.files)){
  
  lateAD_DRrank.TED[[t]]=readRDS(lateAD_DRrank_TED.files[[t]])
  
}

earlyAD_sigTF_bridgedDCLs=lapply(earlyAD_DRrank.TED,function(x)x%>%filter(`p.adjusted`<=0.1))
earlyAD_sigTF_bridgedDCLs=Filter(f = function(x)dim(x)[1]>1,x=earlyAD_sigTF_bridgedDCLs)
earlyAD_sigTF_bridgedDCLs.list=lapply(earlyAD_sigTF_bridgedDCLs,function(x)x%>%pull(regulator_gene_symbol)%>%droplevels%>%levels)
common_earlyAD_sigTF_bridgedDCLs=Reduce(intersect,earlyAD_sigTF_bridgedDCLs.list)

lateAD_sigTF_bridgedDCLs=lapply(lateAD_DRrank.TED,function(x)x%>%filter(`p.adjusted`<=0.1))
lateAD_sigTF_bridgedDCLs=Filter(f = function(x)dim(x)[1]>1,x=lateAD_sigTF_bridgedDCLs)
lateAD_sigTF_bridgedDCLs.list=lapply(lateAD_sigTF_bridgedDCLs,function(x)x%>%pull(regulator_gene_symbol)%>%droplevels%>%levels)
common_lateAD_sigTF_bridgedDCLs=Reduce(intersect,lateAD_sigTF_bridgedDCLs.list)

common_earlyAD_sigTF_regulated_bridgedDCLs.df=data.frame(rbindlist(lapply(earlyAD_TF_bridged_DCL,function(x)x%>%filter(common.TF%in%common_earlyAD_sigTF_bridgedDCLs)%>%pull(common.TF)%>%table%>%as.list)))
rownames(common_earlyAD_sigTF_regulated_bridgedDCLs.df)=names(earlyAD_TF_bridged_DCL)
earlyAD_sigTF_regulated_bridgedDCLs.df=data.frame(rbindlist(lapply(earlyAD_TF_bridged_DCL,function(x)x%>%filter(common.TF%in%common_earlyAD_sigTF_bridgedDCLs)%>%pull(common.TF)%>%table%>%as.list)))
write.table(common_earlyAD_sigTF_regulated_bridgedDCLs.df,"common_earlyAD_sigTF_regulated_bridgedDCLs_df.txt",sep = "\t",col.names = T,row.names = T)
common_lateAD_sigTF_regulated_bridgedDCLs.df=data.frame(rbindlist(lapply(lateAD_TF_bridged_DCL,function(x)x%>%filter(common.TF%in%common_lateAD_sigTF_bridgedDCLs)%>%pull(common.TF)%>%table%>%as.list)))
rownames(common_lateAD_sigTF_regulated_bridgedDCLs.df)=names(lateAD_TF_bridged_DCL)
write.table(common_lateAD_sigTF_regulated_bridgedDCLs.df,"common_lateAD_sigTF_regulated_bridgedDCLs_df.txt",sep = "\t",col.names = T,row.names = T)


earlyAD_TF_bridged_DCL.TF_downstream=earlyAD_TF_bridged_DCL.TF_downstream_subgraph=earlyAD_TF_bridged_DCL.TF_downstream_subgraph_entrez=earlyAD_TF_bridged_DCL.TF_downstream_subgraph_GO=earlyAD_TF_bridged_DCL.TF_downstream_subgraph_GO_simplify=earlyAD_TF_bridged_DCL.TF_downstream_GOsummaries=vector(mode = "list",length = dim(common_earlyAD_sigTF_regulated_bridgedDCLs.df)[2])
names(earlyAD_TF_bridged_DCL.TF_downstream)=names(earlyAD_TF_bridged_DCL.TF_downstream_subgraph)=names(earlyAD_TF_bridged_DCL.TF_downstream_subgraph_entrez)=names(earlyAD_TF_bridged_DCL.TF_downstream_subgraph_GO)=names(earlyAD_TF_bridged_DCL.TF_downstream_subgraph_GO_simplify)=names(earlyAD_TF_bridged_DCL.TF_downstream_GOsummaries)=colnames(common_earlyAD_sigTF_regulated_bridgedDCLs.df)
for(l in 1:length(common_earlyAD_sigTF_regulated_bridgedDCLs.df)){
  
    earlyAD_TF_bridged_DCL.TF_downstream[[l]]=lapply(earlyAD_TF_bridged_DCL,function(x)x%>%filter(common.TF==colnames(common_earlyAD_sigTF_regulated_bridgedDCLs.df)[l])%>%dplyr::select(c(`Gene.1`,`Gene.2`)))
    earlyAD_TF_bridged_DCL.TF_downstream_subgraph[[l]]=lapply(earlyAD_TF_bridged_DCL.TF_downstream[[l]],function(x)graph.data.frame(d = x,directed = F))
    earlyAD_TF_bridged_DCL.TF_downstream_subgraph_entrez[[l]]=lapply(earlyAD_TF_bridged_DCL.TF_downstream_subgraph[[l]],function(x)unname(mapIds(x = org.Hs.eg.db,keys = V(x)$name,column = "ENTREZID",keytype = "SYMBOL")))
}

#Downstream DCLs for all sig TFs
earlyAD_TF_bridged_DCL.TF_downstream=earlyAD_TF_bridged_DCL.TF_downstream_subgraph=vector(mode = "list",length = dim(common_earlyAD_sigTF_regulated_bridgedDCLs.df)[2])


lateAD_TF_bridged_DCL.TF_downstream=lateAD_TF_bridged_DCL.TF_downstream_subgraph=lateAD_TF_bridged_DCL.TF_downstream_subgraph_entrez=lateAD_TF_bridged_DCL.TF_downstream_subgraph_GO=lateAD_TF_bridged_DCL.TF_downstream_subgraph_GO_simplify=lateAD_TF_bridged_DCL.TF_downstream_GOsummaries=vector(mode = "list",length = dim(common_lateAD_sigTF_regulated_bridgedDCLs.df)[2])
names(lateAD_TF_bridged_DCL.TF_downstream)=names(lateAD_TF_bridged_DCL.TF_downstream_subgraph)=names(lateAD_TF_bridged_DCL.TF_downstream_subgraph_entrez)=names(lateAD_TF_bridged_DCL.TF_downstream_subgraph_GO)=names(lateAD_TF_bridged_DCL.TF_downstream_subgraph_GO_simplify)=names(lateAD_TF_bridged_DCL.TF_downstream_GOsummaries)=colnames(common_lateAD_sigTF_regulated_bridgedDCLs.df)
for(l in 1:length(common_lateAD_sigTF_regulated_bridgedDCLs.df)){
  
  lateAD_TF_bridged_DCL.TF_downstream[[l]]=lapply(lateAD_TF_bridged_DCL,function(x)x%>%filter(common.TF==colnames(common_lateAD_sigTF_regulated_bridgedDCLs.df)[l])%>%dplyr::select(c(`Gene.1`,`Gene.2`)))
  lateAD_TF_bridged_DCL.TF_downstream_subgraph[[l]]=lapply(lateAD_TF_bridged_DCL.TF_downstream[[l]],function(x)graph.data.frame(d = x,directed = F))
  lateAD_TF_bridged_DCL.TF_downstream_subgraph_entrez[[l]]=lapply(lateAD_TF_bridged_DCL.TF_downstream_subgraph[[l]],function(x)unname(mapIds(x = org.Hs.eg.db,keys = V(x)$name,column = "ENTREZID",keytype = "SYMBOL")))
}

    #earlyAD_TF_bridged_DCL.TF_downstream_GOsummaries[[l]]=gosummaries(x =lapply(earlyAD_TF_bridged_DCL.TF_downstream_subgraph[[l]],function(x)V(x)$name),organism = "hsapiens",max_p_value = 0.05,go_branches = "BP",ordered_query = F)
    #plot(earlyAD_TF_bridged_DCL.TF_downstream_GOsummaries[[l]],filename=paste(names(earlyAD_TF_bridged_DCL.TF_downstream_GOsummaries)[l],"downstream_GOsummary.pdf",sep = "_"))  

for(t in 1:length(earlyAD_TF_bridged_DCL.TF_downstream_subgraph)){
  for(i in 1:length(earlyAD_TF_bridged_DCL.TF_downstream_subgraph[[t]])){
    write.graph(graph = earlyAD_TF_bridged_DCL.TF_downstream_subgraph[[t]][[i]],file = paste(names(earlyAD_TF_bridged_DCL.TF_downstream_subgraph)[t],names(earlyAD_TF_bridged_DCL.TF_downstream_subgraph[[t]])[i],"downstream_subgraph.gml",sep = "_"),format = "gml")
  }
}
  
  
  
  
gs=gosummaries(lapply(earlyAD_TF_bridged_DCL.CTCF_downstream_subgraph,function(x)V(x)$name),organism = "hsapiens",max_p_value = 0.05)
earlyAD_TF_bridged_DCL.CTCF_downstream=lapply(earlyAD_TF_bridged_DCL,function(x)x%>%filter(common.TF=="CTCF")%>%dplyr::select(c(`Gene.1`,`Gene.2`)))
earlyAD_TF_bridged_DCL.CTCF_downstream_subgraph=lapply(earlyAD_TF_bridged_DCL.CTCF_downstream,function(x)graph.data.frame(d = x,directed = F))
#In silico validation with genetic and epigenomic data

known_AD_genes=scan("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Gene_Lists/Known_AD_genes_2018_Final.txt",sep = "\n",what = "char")
janssen_GWGWAS_MAGMA_hits=scan("/Users/sandeepamberkar/Work/Data/Janssen-AD/Janssen-GWGWAS-MAGMA-SNPgenes.txt",sep = "\n",what = "char")
janssen_meta_hits=scan("/Users/sandeepamberkar/Work/Data/Janssen-AD/Janssen-Meta-Annotated-SNPgenes.txt",sep = "\n",what = "char")
janssen_meta_hits2=c(janssen_meta_hits,unique(unlist(strsplit(grep(pattern = ":",janssen_meta_hits,value = T),split = ":"))))

janssen_C2=data.frame(enricher(gene = janssen_meta_hits,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c2),stringsAsFactors = F)
janssen_C2_KEGG=janssen_C2[grep(pattern = "KEGG",x = janssen_C2$ID),]

janssen_KEGG_random=replicate(n = 1000,expr = sample(x = unique(grep(pattern = "KEGG",msigdb_c2$ont,value = T)),size = dim(janssen_C2_KEGG)[1],replace = F),simplify = T)
earlyAD_FrontalPole_KEGG_random=replicate(n = 1000,expr = sample(x = unique(grep(pattern = "KEGG",msigdb_c2$ont,value = T)),size = dim(earlyAD_DCGs.Uniq_c2$Frontal_Pole)[1],replace = F),simplify = T)


# Common enriched pathways in Jansen and early/late Uniq DCGs
janssen_earlyAD_KEGG=Filter(f = function(x)length(x)>0,x = lapply(lapply(earlyAD_DCGs.Uniq_c2,function(x)grep(pattern = "KEGG",x = x$ID,value = T)%>%gsub(pattern = "KEGG_",replacement = "")%>%gsub(pattern = "_",replacement = " ")),intersect,toupper(gsub(pattern = "_Homo sapiens_hsa[0-9]+",enrichr(genes = janssen_meta_hits,databases = kegg_dbs)[[1]]%>%dplyr::filter(`Adjusted.P.value`<=0.2)%>%pull(Term),replacement = ""))))
janssen_lateAD_KEGG=Filter(f = function(x)length(x)>0,x = lapply(lapply(lateAD_DCGs.Uniq_c2,function(x)grep(pattern = "KEGG",x = x$ID,value = T)%>%gsub(pattern = "KEGG_",replacement = "")%>%gsub(pattern = "_",replacement = " ")),intersect,toupper(gsub(pattern = "_Homo sapiens_hsa[0-9]+",enrichr(genes = janssen_meta_hits,databases = kegg_dbs)[[1]]%>%dplyr::filter(`Adjusted.P.value`<=0.1)%>%pull(Term),replacement = ""))))
earlyAD_janssenAD.overlap=lapply(earlyAD_DCGs.uniq,intersect,janssen_meta_hits)
commonDCGs_janssenAD.overlap=lapply(common_DCGs,intersect,janssen_meta_hits)
lateAD_janssenAD.overlap=lapply(lateAD_DCGs.uniq,intersect,janssen_meta_hits)
earlyAD_janssenAD.c2=lapply(Filter(f = function(x)length(x)>0,earlyAD_janssenAD.overlap),function(y)data.frame(enricher(gene =  y,pvalueCutoff = 0.05,pAdjustMethod = "BH",TERM2GENE = msigdb_c2),stringsAsFactors = F))
#earlyAD_janssenAD.c2=lapply(earlyAD_janssenAD.c2,function(x)x%>%dplyr::filter(Count>2))
earlyAD_janssenAD.c2=Filter(f = function(x)dim(x)[1]>=1,lapply(Filter(f = function(x)dim(x)[1]>1,earlyAD_janssenAD.c2),function(a)a[grep(pattern = "KEGG|REACTOME",x = a$ID),]))
commonDCGs_janssenAD.c2=lapply(Filter(f = function(x)length(x)>0,commonDCGs_janssenAD.overlap),function(y)data.frame(enricher(gene =  y,pvalueCutoff = 0.05,pAdjustMethod = "BH",TERM2GENE = msigdb_c2),stringsAsFactors = F))
lateAD_janssenAD.c2=lapply(Filter(f = function(x)length(x)>0,lateAD_janssenAD.overlap),function(y)data.frame(enricher(gene =  y,pvalueCutoff = 0.05,pAdjustMethod = "BH",TERM2GENE = msigdb_c2),stringsAsFactors = F))

tanzi_SAGs=scan("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Gene_Lists/Tanzi_Ranked_SNP_Genes.txt",sep = "\n",what="char")


janssen_DCG_pathways.df=cbind.data.frame(EarlyAD_GWAS_pathways=unlist(lapply(lateAD_janssenAD.c2,function(x)paste(grep(pattern = "REACTOME|KEGG",x = x$ID,value = T),collapse = ";"))),
                 PersistentDCGs_GWAS_pathways=unlist(lapply(commonDCGs_janssenAD.c2,function(x)paste(grep(pattern = "REACTOME|KEGG",x = x$ID,value = T),collapse = ";"))),
                 LateAD_GWAS_pathways=unlist(lapply(lateAD_janssenAD.c2,function(x)paste(grep(pattern = "REACTOME|KEGG",x = x$ID,value = T),collapse = ";"))))

earlyAD_janssenAD.exp=lapply(lapply(earlyAD_DCGs.uniq,length),function(x)round(length(janssen_meta_hits)/19529*x,digits = 1))
lateAD_janssenAD.exp=lapply(lapply(lateAD_DCGs.uniq,length),function(x)round(length(janssen_meta_hits)/19529*x,digits = 1))


res_early.janssen_DCG=foreach(i=1:length(earlyAD_janssenAD.overlap))%do%{
  
  
  data.frame(GWAS_in_DCGs=unlist(lapply(earlyAD_janssenAD.overlap[i],length)),
             DCGs=unlist(lapply(earlyAD_DCGs.uniq,length)[[i]]),
             Expected_GWAS_in_DCGs=unname(unlist(earlyAD_janssenAD.exp[i])),
             Fisher_statistic=unname(unlist(lapply(earlyAD_janssenAD.fisher[i],function(s)s$estimate))),
             pval=unname(unlist(lapply(earlyAD_janssenAD.fisher[i],function(p)p$p.value))),
             stringsAsFactors = F)
}
names(res_early.janssen_DCG)=names(earlyAD_janssenAD.overlap)
res.SAG_earlyAD_DCG_df=data.frame(rbindlist(res_early.janssen_DCG),stringsAsFactors = F)
rownames(res.SAG_earlyAD_DCG_df)=names(earlyAD_janssenAD.overlap)
res.SAG_earlyAD_DCG_df$adj.p=p.adjust(p = res.SAG_earlyAD_DCG_df$pval,method = "fdr")
res.SAG_earlyAD_DCG_df$GWAS_in_DCGs=as.numeric(res.SAG_earlyAD_DCG_df$GWAS_in_DCGs)
res.SAG_earlyAD_DCG_df=res.SAG_earlyAD_DCG_df%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(GWAS_in_DCGs>Expected_GWAS_in_DCGs,true = "Over",false = "Under"))  
fwrite(res.SAG_earlyAD_DCG_df,"EarlyAD_Janssen_DCG_EnrichmentSummary.txt",sep = "\t",col.names = T,row.names = F,quote = F)

res_late.janssen_DCG=foreach(i=1:length(lateAD_janssenAD.overlap))%do%{
  
  
  data.frame(GWAS_in_DCGs=unlist(lapply(lateAD_janssenAD.overlap[i],length)),
             DCGs=unlist(lapply(lateAD_DCGs.uniq,length)[[i]]),
             Expected_GWAS_in_DCGs=unname(unlist(lateAD_janssenAD.exp[i])),
             Fisher_statistic=unname(unlist(lapply(lateAD_janssenAD.fisher[i],function(s)s$estimate))),
             pval=unname(unlist(lapply(lateAD_janssenAD.fisher[i],function(p)p$p.value))),
             stringsAsFactors = F)
}
names(res_late.janssen_DCG)=names(lateAD_janssenAD.overlap)
res.SAG_lateAD_DCG_df=data.frame(rbindlist(res_late.janssen_DCG),stringsAsFactors = F)
rownames(res.SAG_lateAD_DCG_df)=names(lateAD_janssenAD.overlap)
res.SAG_lateAD_DCG_df$adj.p=p.adjust(p = res.SAG_lateAD_DCG_df$pval,method = "fdr")
res.SAG_lateAD_DCG_df$GWAS_in_DCGs=as.numeric(res.SAG_lateAD_DCG_df$GWAS_in_DCGs)
res.SAG_lateAD_DCG_df=res.SAG_lateAD_DCG_df%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(GWAS_in_DCGs>Expected_GWAS_in_DCGs,true = "Over",false = "Under"))  
fwrite(res.SAG_lateAD_DCG_df,"LateAD_Janssen_DCG_EnrichmentSummary.txt",sep = "\t",col.names = T,row.names = F,quote = F)


dhmc_bennet.NP=read.table("/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/GW-DNA-hydroxymethylation-Bennet/1-s2.0-S155252601633059X-mmc4_NP.txt",sep = "\t",header = T,as.is = T)
dhmc_bennet.NFT=read.table("/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/GW-DNA-hydroxymethylation-Bennet/1-s2.0-S155252601633059X-mmc5_NFT.txt",sep = "\t",header = T,as.is = T)
d5mc_bernstein.Dhml=scan("/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/5hMC_Tau_AD/ddw109_Supp/Supplemental Table 6_Bernstein et al_R1.txt",sep = "\n",what = "char")

earlyAD_DCG_d5mc.overlap=lapply(earlyAD_DCGs.list,function(x)length(intersect(x,d5mc_bernstein.Dhml)))
dhmc_bennet_NFT.genes=dhmc_bennet.NFT%>%dplyr::filter(q.value_<=0.1)%>%pull(Nearest.gene)
dhmc_bennet_NP.genes=dhmc_bennet.NP%>%dplyr::filter(q.vlaue_<=0.1)%>%pull(Neartest.gene)
earlyAD_DCG_d5mc.overlap=lapply(earlyAD_DCGs.list,function(x)length(intersect(x,d5mc_bernstein.Dhml)))
dhmc_bennet_NFT.genes=dhmc_bennet.NFT%>%dplyr::filter(q.value_<=0.1)%>%pull(Nearest.gene)
dhmc_bennet_NP.genes=dhmc_bennet.NP%>%dplyr::filter(q.vlaue_<=0.1)%>%pull(Neartest.gene)
earlyAD_DCG_dhmc_NFT.overlap=lapply(earlyAD_DCGs.list,function(x)intersect(x,dhmc_bennet_NFT.genes))
earlyAD_DCG_dhmc_NP.overlap=lapply(earlyAD_DCGs.list,function(x)intersect(x,dhmc_bennet_NP.genes))
earlyAD_exp_dhmc_NP_DCGs=lapply(lapply(earlyAD_DCGs.list,length),function(x)round(length(dhmc_bennet_NP.genes)/19530*x,digits = 1))
earlyAD_exp_dhmc_NFT_DCGs=lapply(lapply(earlyAD_DCGs.list,length),function(x)round(length(dhmc_bennet_NFT.genes)/19530*x,digits = 1))

earlyAD_dhmc_NP.fisher_results=earlyAD_dhmc_NFT.fisher_results=earlyAD_d5mc.fisher_results=vector(mode = "list",length = length(earlyAD_DCGs.list))
names(earlyAD_dhmc_NP.fisher_results)=names(earlyAD_dhmc_NFT.fisher_results)=names(earlyAD_d5mc.fisher_results)=names(earlyAD_DCGs.list)
for(m in 1:length(earlyAD_DCGs.list)){
  earlyAD_dhmc_NP.fisher_results[[m]]=fisher.test(matrix(c(length(earlyAD_DCG_dhmc_NP.overlap[[m]]),
                                                           length(dhmc_bennet_NP.genes)-length(earlyAD_DCG_dhmc_NP.overlap[[m]]),
                                                           lapply(earlyAD_DCGs.list,length)[[m]]-length(earlyAD_DCG_dhmc_NP.overlap[[m]]),
                                                           19530-(length(dhmc_bennet_NP.genes)-length(earlyAD_DCG_dhmc_NP.overlap[[m]]))-
                                                             lapply(earlyAD_DCGs.list,length)[[m]]-length(earlyAD_DCG_dhmc_NP.overlap[[m]])),nrow = 2))
}

res_early.dhmc_NP_DCG=foreach(i=1:length(earlyAD_DCG_dhmc_NP.overlap))%do%{
  
  
  data.frame(NPgenes_in_DCGs=unlist(lapply(earlyAD_DCG_dhmc_NP.overlap[i],length)),
             DCGs=unlist(lapply(earlyAD_DCGs.list,length)[[i]]),
             Expected_dhmc_NPgenes_in_DCGs=unname(unlist(earlyAD_exp_dhmc_NP_DCGs[i])),
             Fisher_statistic=unname(unlist(lapply(earlyAD_dhmc_NP.fisher_results[i],function(s)s$estimate))),
             pval=unname(unlist(lapply(earlyAD_dhmc_NP.fisher_results[i],function(p)p$p.value))),
             stringsAsFactors = F)
}


names(res_early.dhmc_NP_DCG)=names(earlyAD_DCG_dhmc_NP.overlap)
res_early.dhmc_NP_DCG_df=data.frame(rbindlist(res_early.dhmc_NP_DCG),stringsAsFactors = F)
rownames(res_early.dhmc_NP_DCG_df)=names(earlyAD_DCG_dhmc_NP.overlap)
res_early.dhmc_NP_DCG_df$adj.p=p.adjust(p = res_early.dhmc_NP_DCG_df$pval,method = "fdr")
res_early.dhmc_NP_DCG_df$NPgenes_in_DCGs=as.numeric(res_early.dhmc_NP_DCG_df$NPgenes_in_DCGs)
res_early.dhmc_NP_DCG_df=res_early.dhmc_NP_DCG_df%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(NPgenes_in_DCGs>Expected_dhmc_NPgenes_in_DCGs,true = "Over",false = "Under"))

earlyAD_DCG_dhmc_NFT.overlap=lapply(earlyAD_DCGs.list,function(x)intersect(x,dhmc_bennet_NFT.genes))
earlyAD_DCG_dhmc_NFT.overlap=Filter(f = function(x)length(x)>1,x = earlyAD_DCG_dhmc_NFT.overlap)
for(m in names(earlyAD_DCG_dhmc_NFT.overlap)){
  earlyAD_dhmc_NFT.fisher_results[[m]]=fisher.test(matrix(c(length(earlyAD_DCG_dhmc_NFT.overlap[[m]]),
                                                            length(dhmc_bennet_NFT.genes)-length(earlyAD_DCG_dhmc_NFT.overlap[[m]]),
                                                            lapply(earlyAD_DCGs.list,length)[[m]]-length(earlyAD_DCG_dhmc_NFT.overlap[[m]]),
                                                            19530-(length(dhmc_bennet_NFT.genes)-length(earlyAD_DCG_dhmc_NFT.overlap[[m]]))-
                                                              lapply(earlyAD_DCGs.list,length)[[m]]-length(earlyAD_DCG_dhmc_NFT.overlap[[m]])),nrow = 2))
}

res_early.dhmc_NFT_DCG=foreach(i=names(earlyAD_DCG_dhmc_NFT.overlap))%do%{
  
  
  data.frame(NFTgenes_in_DCGs=unlist(lapply(earlyAD_DCG_dhmc_NFT.overlap[i],length)),
             DCGs=unlist(lapply(earlyAD_DCGs.list,length)[[i]]),
             Expected_dhmc_NFTgenes_in_DCGs=unname(unlist(earlyAD_exp_dhmc_NFT_DCGs[i])),
             Fisher_statistic=unname(unlist(lapply(earlyAD_dhmc_NFT.fisher_results[i],function(s)s$estimate))),
             pval=unname(unlist(lapply(earlyAD_dhmc_NFT.fisher_results[i],function(p)p$p.value))),
             stringsAsFactors = F)
}


names(res_early.dhmc_NFT_DCG)=names(earlyAD_DCG_dhmc_NFT.overlap)
res_early.dhmc_NFT_DCG_df=data.frame(rbindlist(res_early.dhmc_NFT_DCG),stringsAsFactors = F)
rownames(res_early.dhmc_NFT_DCG_df)=names(earlyAD_DCG_dhmc_NFT.overlap)
res_early.dhmc_NFT_DCG_df$adj.p=p.adjust(p = res_early.dhmc_NFT_DCG_df$pval,method = "fdr")
res_early.dhmc_NFT_DCG_df$NFTgenes_in_DCGs=as.numeric(res_early.dhmc_NFT_DCG_df$NFTgenes_in_DCGs)
res_early.dhmc_NFT_DCG_df=res_early.dhmc_NFT_DCG_df%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(NFTgenes_in_DCGs>Expected_dhmc_NFTgenes_in_DCGs,true = "Over",false = "Under"))

fwrite(res.dhmc_NFT_DCG_df,"EarlyAD_DhMR_NFT_earlyAD_DCG_EnrichmentSummary.txt",sep = "\t",col.names = T,row.names = F,quote = F)


#STG GW-Methylation dataset
stg_gw_methylation=fread("/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/STG-GW-Methylation-Watson/13073_2015_258_MOESM6_ESM.txt",sep = "\t",header = T,data.table = F)
stg_gw_methylation_dmr.genes=stg_gw_methylation%>%dplyr::filter(`FDR_DMR.q.value`<=0.2)%>%pull(`Closest.Gene`)%>%unique
earlyAD_DCG_stg_methylation_dmr.overlap=lapply(earlyAD_DCGs.uniq,function(x)intersect(x,stg_gw_methylation_dmr.genes))
earlyAD_exp_stg_methylation_dmr_DCGs=lapply(lapply(earlyAD_DCGs.uniq,length),function(x)round(length(stg_gw_methylation_dmr.genes)/19530*x,digits = 1))

earlAD_stg_methylation.fisher_results=vector(mode = "list",length = length(earlyAD_DCGs.uniq))

for(m in 1:length(earlyAD_DCGs.uniq)){
  earlAD_stg_methylation.fisher_results[[m]]=fisher.test(matrix(c(length(earlyAD_DCG_stg_methylation_dmr.overlap[[m]]),
                                                                  length(stg_gw_methylation_dmr.genes)-length(earlyAD_DCG_stg_methylation_dmr.overlap[[m]]),
                                                                  lapply(earlyAD_DCGs.uniq,length)[[m]]-length(earlyAD_DCG_stg_methylation_dmr.overlap[[m]]),
                                                                  19530-(length(stg_gw_methylation_dmr.genes)-length(earlyAD_DCG_stg_methylation_dmr.overlap[[m]]))-
                                                                    lapply(earlyAD_DCGs.uniq,length)[[m]]-length(earlyAD_DCG_stg_methylation_dmr.overlap[[m]])),nrow = 2))
}

res.earlyAD_stg_methylation=foreach(i=1:length(earlyAD_DCG_stg_methylation_dmr.overlap))%do%{
  
  
  data.frame(Methylated_genes_in_DCGs=unlist(lapply(earlyAD_DCG_stg_methylation_dmr.overlap[i],length)),
             DCGs=unlist(lapply(earlyAD_DCGs.uniq,length)[[i]]),
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

stg_gw_methylation=fread("/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/STG-GW-Methylation-Watson/13073_2015_258_MOESM6_ESM.txt",sep = "\t",header = T,data.table = F)
stg_gw_methylation_dmr.genes=stg_gw_methylation%>%dplyr::filter(`FDR_DMR.q.value`<=0.2)%>%pull(`Closest.Gene`)%>%unique
lateAD_DCG_stg_methylation_dmr.overlap=lapply(lateAD_DCGs.uniq,function(x)intersect(x,stg_gw_methylation_dmr.genes))
lateAD_exp_stg_methylation_dmr_DCGs=lapply(lapply(lateAD_DCGs.uniq,length),function(x)round(length(stg_gw_methylation_dmr.genes)/19530*x,digits = 1))

lateAD_stg_methylation.fisher_results=vector(mode = "list",length = length(lateAD_DCGs.uniq))

for(m in 1:length(lateAD_DCGs.uniq)){
  lateAD_stg_methylation.fisher_results[[m]]=fisher.test(matrix(c(length(lateAD_DCG_stg_methylation_dmr.overlap[[m]]),
                                                                  length(stg_gw_methylation_dmr.genes)-length(lateAD_DCG_stg_methylation_dmr.overlap[[m]]),
                                                                  lapply(lateAD_DCGs.uniq,length)[[m]]-length(lateAD_DCG_stg_methylation_dmr.overlap[[m]]),
                                                                  19530-(length(stg_gw_methylation_dmr.genes)-length(lateAD_DCG_stg_methylation_dmr.overlap[[m]]))-
                                                                    lapply(lateAD_DCGs.uniq,length)[[m]]-length(lateAD_DCG_stg_methylation_dmr.overlap[[m]])),nrow = 2))
}

res.lateAD_stg_methylation=foreach(i=1:length(lateAD_DCG_stg_methylation_dmr.overlap))%do%{
  
  
  data.frame(Methylated_genes_in_DCGs=unlist(lapply(lateAD_DCG_stg_methylation_dmr.overlap[i],length)),
             DCGs=unlist(lapply(lateAD_DCGs.uniq,length)[[i]]),
             Expected_STG_methylation_in_DCGs=unname(unlist(lateAD_exp_stg_methylation_dmr_DCGs[i])),
             Fisher_statistic=unname(unlist(lapply(lateAD_stg_methylation.fisher_results[i],function(s)s$estimate))),
             pval=unname(unlist(lapply(lateAD_stg_methylation.fisher_results[i],function(p)p$p.value))),
             stringsAsFactors = F)
}


names(res.lateAD_stg_methylation)=names(lateAD_DCG_stg_methylation_dmr.overlap)
res_lateAD_stg_methylation.df=data.frame(rbindlist(res.lateAD_stg_methylation),stringsAsFactors = F)
rownames(res_lateAD_stg_methylation.df)=names(lateAD_DCG_stg_methylation_dmr.overlap)
res_lateAD_stg_methylation.df$adj.p=p.adjust(p = res_lateAD_stg_methylation.df$pval,method = "fdr")
res_lateAD_stg_methylation.df$Methylated_genes_in_DCGs=as.numeric(res_lateAD_stg_methylation.df$Methylated_genes_in_DCGs)
res_lateAD_stg_methylation.df=res_lateAD_stg_methylation.df%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(Methylated_genes_in_DCGs>Expected_STG_methylation_in_DCGs,true = "Over",false = "Under"))






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


earlyAD.exp_Astrocytes=lapply(lapply(earlyAD_DCGs.uniq,length),function(x)round(length(zhang_celltype_ADgenes.list$ast)/19530*x,digits = 1))
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
res.earlyAD_DCG_Astrocytes_df=res.earlyAD_DCG_Astrocytes_df%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(Nr_AstrocyteMarkers_in_DCGs>Expected_AstrocyteMarkers_in_DCGs,true = "Over",false = "Under"))%>%filter(Representation=="Under")


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
res.earlyAD_DCG_Endothelial_df=res.earlyAD_DCG_Endothelial_df%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(Nr_EndothelialMarkers_in_DCGs>Expected_EndothelialMarkers_in_DCGs,true = "Over",false = "Under"))%>%filter(Representation=="Under")

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
res.earlyAD_DCG_Microglia_df=res.earlyAD_DCG_Microglia_df%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(Nr_MicrogliaMarkers_in_DCGs>Expected_MicrogliaMarkers_in_DCGs,true = "Over",false = "Under"))%>%filter(Representation=="Under")


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
res.earlyAD_DCG_Neurons_df=res.earlyAD_DCG_Neurons_df%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(Nr_NeuronsMarkers_in_DCGs>Expected_NeuronsMarkers_in_DCGs,true = "Over",false = "Under"))%>%filter(Representation=="Under")

earlyAD_DCG_Neurons.enrichedPathways=lapply(res.earlyAD_DCG_Neurons_df$NeuronsMarkers_in_DCGs,function(a)data.frame(enricher(gene = strsplit(x = a,split = ",")[[1]],pvalueCutoff = 0.05,pAdjustMethod = "BH",TERM2GENE = msigdb_c2),stringsAsFactors = F))
names(earlyAD_DCG_Neurons.enrichedPathways)=res.earlyAD_DCG_Neurons_df$`Brain-region`
earlyAD_DCG_Neurons.enrichedPathways=Filter(f = function(x)dim(x)[1]>1,earlyAD_DCG_Neurons.enrichedPathways)
earlyAD_DCG_Neurons.enrichedPathways=Filter(f = ,function(x)dim(x)[1]>1,x = lapply(earlyAD_DCG_Neurons.enrichedPathways,function(x)x[grep(pattern = "KEGG|REACTOME",x = x$ID),]))


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
res.earlyAD_DCG_Oligodendrocytes_df=res.earlyAD_DCG_Oligodendrocytes_df%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(Nr_OligodendrocytesMarkers_in_DCGs>Expected_OligodendrocytesMarkers_in_DCGs,true = "Over",false = "Under"))%>%filter(Representation=="Under")





lateAD_DCG_Astrocytes.overlap=lapply(lateAD_DCGs.uniq,intersect,zhang_celltype_ADgenes.list$ast)
lateAD_DCG_Endothelial.overlap=lapply(lateAD_DCGs.uniq,intersect,zhang_celltype_ADgenes.list$end)
lateAD_DCG_Microglia.overlap=lapply(lateAD_DCGs.uniq,intersect,zhang_celltype_ADgenes.list$mic)
lateAD_DCG_Neurons.overlap=lapply(lateAD_DCGs.uniq,intersect,zhang_celltype_ADgenes.list$neu)
lateAD_DCG_Oligodendrocytes.overlap=lapply(lateAD_DCGs.uniq,intersect,zhang_celltype_ADgenes.list$oli)


lateAD.exp_Astrocytes=lapply(lapply(lateAD_DCGs.uniq,length),function(x)round(length(zhang_celltype_ADgenes.list$ast)/19530*x,digits = 1))
lateAD.exp_Endothelial=lapply(lapply(lateAD_DCGs.uniq,length),function(x)round(length(zhang_celltype_ADgenes.list$end)/19530*x,digits = 1))
lateAD.exp_Microglia=lapply(lapply(lateAD_DCGs.uniq,length),function(x)round(length(zhang_celltype_ADgenes.list$mic)/19530*x,digits = 1))
lateAD.exp_Neurons=lapply(lapply(lateAD_DCGs.uniq,length),function(x)round(length(zhang_celltype_ADgenes.list$neu)/19530*x,digits = 1))
lateAD.exp_Oligodendrocytes=lapply(lapply(lateAD_DCGs.uniq,length),function(x)round(length(zhang_celltype_ADgenes.list$oli)/19530*x,digits = 1))


lateAD_Astrocytes.fisher_results=lateAD_Endothelial.fisher_results=lateAD_Microglia.fisher_results=lateAD_Neurons.fisher_results=lateAD_Oligodendrocytes.fisher_results=vector(mode = "list",length = length(lateAD_DCGs.uniq))
names(lateAD_Astrocytes.fisher_results)=names(lateAD_Endothelial.fisher_results)=names(lateAD_Microglia.fisher_results)=names(lateAD_Neurons.fisher_results)=names(lateAD_Oligodendrocytes.fisher_results)=names(lateAD_DCGs.uniq)

for(m in 1:length(lateAD_Astrocytes.fisher_results)){
  lateAD_Astrocytes.fisher_results[[m]]=fisher.test(matrix(c(length(lateAD_DCG_Astrocytes.overlap[[m]]),
                                                             length(zhang_celltype_ADgenes.list$ast)-length(lateAD_DCG_Astrocytes.overlap[[m]]),
                                                             lapply(lateAD_DCGs.uniq,length)[[m]]-length(lateAD_DCG_Astrocytes.overlap[[m]]),
                                                             19530-(length(zhang_celltype_ADgenes.list$ast)-length(lateAD_DCG_Astrocytes.overlap[[m]]))-
                                                               lapply(lateAD_DCGs.uniq,length)[[m]]-length(lateAD_DCG_Astrocytes.overlap[[m]])),nrow = 2))
}
for(m in 1:length(lateAD_Endothelial.fisher_results)){
  lateAD_Endothelial.fisher_results[[m]]=fisher.test(matrix(c(length(lateAD_DCG_Endothelial.overlap[[m]]),
                                                              length(zhang_celltype_ADgenes.list$end)-length(lateAD_DCG_Endothelial.overlap[[m]]),
                                                              lapply(lateAD_DCGs.uniq,length)[[m]]-length(lateAD_DCG_Endothelial.overlap[[m]]),
                                                              19530-(length(zhang_celltype_ADgenes.list$end)-length(lateAD_DCG_Endothelial.overlap[[m]]))-
                                                              lapply(lateAD_DCGs.uniq,length)[[m]]-length(lateAD_DCG_Endothelial.overlap[[m]])),nrow = 2))
}
for(m in 1:length(lateAD_Microglia.fisher_results)){
  lateAD_Microglia.fisher_results[[m]]=fisher.test(matrix(c(length(lateAD_DCG_Microglia.overlap[[m]]),
                                                            length(zhang_celltype_ADgenes.list$mic)-length(lateAD_DCG_Microglia.overlap[[m]]),
                                                            lapply(lateAD_DCGs.uniq,length)[[m]]-length(lateAD_DCG_Microglia.overlap[[m]]),
                                                            19530-(length(zhang_celltype_ADgenes.list$mic)-length(lateAD_DCG_Microglia.overlap[[m]]))-
                                                            lapply(lateAD_DCGs.uniq,length)[[m]]-length(lateAD_DCG_Microglia.overlap[[m]])),nrow = 2))
}
for(m in 1:length(lateAD_Neurons.fisher_results)){
  lateAD_Neurons.fisher_results[[m]]=fisher.test(matrix(c(length(lateAD_DCG_Neurons.overlap[[m]]),
                                                          length(zhang_celltype_ADgenes.list$neu)-length(lateAD_DCG_Neurons.overlap[[m]]),
                                                          lapply(lateAD_DCGs.uniq,length)[[m]]-length(lateAD_DCG_Neurons.overlap[[m]]),
                                                          19530-(length(zhang_celltype_ADgenes.list$neu)-length(lateAD_DCG_Neurons.overlap[[m]]))-
                                                          lapply(lateAD_DCGs.uniq,length)[[m]]-length(lateAD_DCG_Neurons.overlap[[m]])),nrow = 2))
}
for(m in 1:length(lateAD_Oligodendrocytes.fisher_results)){
  lateAD_Oligodendrocytes.fisher_results[[m]]=fisher.test(matrix(c(length(lateAD_DCG_Oligodendrocytes.overlap[[m]]),
                                                                   length(zhang_celltype_ADgenes.list$oli)-length(lateAD_DCG_Oligodendrocytes.overlap[[m]]),
                                                                   lapply(lateAD_DCGs.uniq,length)[[m]]-length(lateAD_DCG_Oligodendrocytes.overlap[[m]]),
                                                                   19530-(length(zhang_celltype_ADgenes.list$oli)-length(lateAD_DCG_Oligodendrocytes.overlap[[m]]))-
                                                                   lapply(lateAD_DCGs.uniq,length)[[m]]-length(lateAD_DCG_Oligodendrocytes.overlap[[m]])),nrow = 2))
}




res.lateAD_DCG_Astrocytes=foreach(i=1:length(lateAD_DCG_Astrocytes.overlap))%do%{
  
  
  data.frame(AstrocyteMarkers_in_DCGs=unlist(lapply(lateAD_DCG_Astrocytes.overlap[i],paste,collapse=",")),
             Nr_AstrocyteMarkers_in_DCGs=unlist(lapply(lapply(lateAD_DCG_Astrocytes.overlap,length)[i],paste,collapse=",")),
             DCGs=unlist(lapply(lateAD_DCGs.uniq,length)[[i]]),
             Expected_AstrocyteMarkers_in_DCGs=unname(unlist(lateAD.exp_Astrocytes[i])),
             Fisher_statistic=unname(unlist(lapply(lateAD_Astrocytes.fisher_results[i],function(s)s$estimate))),
             pval=unname(unlist(lapply(lateAD_Astrocytes.fisher_results[i],function(p)p$p.value))),
             stringsAsFactors = F)
}
names(res.lateAD_DCG_Astrocytes)=names(lateAD_DCG_Astrocytes.overlap)
res.lateAD_DCG_Astrocytes_df=data.frame(rbindlist(res.lateAD_DCG_Astrocytes),stringsAsFactors = F)
rownames(res.lateAD_DCG_Astrocytes_df)=names(lateAD_DCG_Astrocytes.overlap)
res.lateAD_DCG_Astrocytes_df$adj.p=p.adjust(p = res.lateAD_DCG_Astrocytes_df$pval,method = "fdr")
res.lateAD_DCG_Astrocytes_df$Expected_AstrocyteMarkers_in_DCGs=as.numeric(res.lateAD_DCG_Astrocytes_df$Expected_AstrocyteMarkers_in_DCGs)
res.lateAD_DCG_Astrocytes_df$Nr_AstrocyteMarkers_in_DCGs=as.numeric(res.lateAD_DCG_Astrocytes_df$Nr_AstrocyteMarkers_in_DCGs)
res.lateAD_DCG_Astrocytes_df=res.lateAD_DCG_Astrocytes_df%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(Nr_AstrocyteMarkers_in_DCGs>Expected_AstrocyteMarkers_in_DCGs,true = "Over",false = "Under"))%>%filter(Representation=="Under")


res.lateAD_DCG_Endothelial=foreach(i=1:length(lateAD_DCG_Endothelial.overlap))%do%{
  
  
  data.frame(EndothelialMarkers_in_DCGs=unlist(lapply(lateAD_DCG_Endothelial.overlap[i],paste,collapse=",")),
             Nr_EndothelialMarkers_in_DCGs=unlist(lapply(lapply(lateAD_DCG_Endothelial.overlap,length)[i],paste,collapse=",")),
             DCGs=unlist(lapply(lateAD_DCGs.uniq,length)[[i]]),
             Expected_EndothelialMarkers_in_DCGs=unname(unlist(lateAD.exp_Endothelial[i])),
             Fisher_statistic=unname(unlist(lapply(lateAD_Endothelial.fisher_results[i],function(s)s$estimate))),
             pval=unname(unlist(lapply(lateAD_Endothelial.fisher_results[i],function(p)p$p.value))),
             stringsAsFactors = F)
}
names(res.lateAD_DCG_Endothelial)=names(lateAD_DCG_Endothelial.overlap)
res.lateAD_DCG_Endothelial_df=data.frame(rbindlist(res.lateAD_DCG_Endothelial),stringsAsFactors = F)
rownames(res.lateAD_DCG_Endothelial_df)=names(lateAD_DCG_Endothelial.overlap)
res.lateAD_DCG_Endothelial_df$adj.p=p.adjust(p = res.lateAD_DCG_Endothelial_df$pval,method = "fdr")
res.lateAD_DCG_Endothelial_df$Expected_EndothelialMarkers_in_DCGs=as.numeric(res.lateAD_DCG_Endothelial_df$Expected_EndothelialMarkers_in_DCGs)
res.lateAD_DCG_Endothelial_df$Nr_EndothelialMarkers_in_DCGs=as.numeric(res.lateAD_DCG_Endothelial_df$Nr_EndothelialMarkers_in_DCGs)
res.lateAD_DCG_Endothelial_df=res.lateAD_DCG_Endothelial_df%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(Nr_EndothelialMarkers_in_DCGs>Expected_EndothelialMarkers_in_DCGs,true = "Over",false = "Under"))%>%filter(Representation=="Under")

res.lateAD_DCG_Microglia=foreach(i=1:length(lateAD_DCG_Microglia.overlap))%do%{
  
  
  data.frame(MicrogliaMarkers_in_DCGs=unlist(lapply(lateAD_DCG_Microglia.overlap[i],paste,collapse=",")),
             Nr_MicrogliaMarkers_in_DCGs=unlist(lapply(lapply(lateAD_DCG_Microglia.overlap,length)[i],paste,collapse=",")),
             DCGs=unlist(lapply(lateAD_DCGs.uniq,length)[[i]]),
             Expected_MicrogliaMarkers_in_DCGs=unname(unlist(lateAD.exp_Microglia[i])),
             Fisher_statistic=unname(unlist(lapply(lateAD_Microglia.fisher_results[i],function(s)s$estimate))),
             pval=unname(unlist(lapply(lateAD_Microglia.fisher_results[i],function(p)p$p.value))),
             stringsAsFactors = F)
}
names(res.lateAD_DCG_Microglia)=names(lateAD_DCG_Microglia.overlap)
res.lateAD_DCG_Microglia_df=data.frame(rbindlist(res.lateAD_DCG_Microglia),stringsAsFactors = F)
rownames(res.lateAD_DCG_Microglia_df)=names(lateAD_DCG_Microglia.overlap)
res.lateAD_DCG_Microglia_df$adj.p=p.adjust(p = res.lateAD_DCG_Microglia_df$pval,method = "fdr")
res.lateAD_DCG_Microglia_df$Expected_MicrogliaMarkers_in_DCGs=as.numeric(res.lateAD_DCG_Microglia_df$Expected_MicrogliaMarkers_in_DCGs)
res.lateAD_DCG_Microglia_df$Nr_MicrogliaMarkers_in_DCGs=as.numeric(res.lateAD_DCG_Microglia_df$Nr_MicrogliaMarkers_in_DCGs)
res.lateAD_DCG_Microglia_df=res.lateAD_DCG_Microglia_df%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(Nr_MicrogliaMarkers_in_DCGs>Expected_MicrogliaMarkers_in_DCGs,true = "Over",false = "Under"))%>%filter(Representation=="Under")

res.lateAD_DCG_Neurons=foreach(i=1:length(lateAD_DCG_Neurons.overlap))%do%{
  
  
  data.frame(NeuronsMarkers_in_DCGs=unlist(lapply(lateAD_DCG_Neurons.overlap[i],paste,collapse=",")),
             Nr_NeuronsMarkers_in_DCGs=unlist(lapply(lapply(lateAD_DCG_Neurons.overlap,length)[i],paste,collapse=",")),
             DCGs=unlist(lapply(lateAD_DCGs.uniq,length)[[i]]),
             Expected_NeuronsMarkers_in_DCGs=unname(unlist(lateAD.exp_Neurons[i])),
             Fisher_statistic=unname(unlist(lapply(lateAD_Neurons.fisher_results[i],function(s)s$estimate))),
             pval=unname(unlist(lapply(lateAD_Neurons.fisher_results[i],function(p)p$p.value))),
             stringsAsFactors = F)
}
names(res.lateAD_DCG_Neurons)=names(lateAD_DCG_Neurons.overlap)
res.lateAD_DCG_Neurons_df=data.frame(rbindlist(res.lateAD_DCG_Neurons),stringsAsFactors = F)
rownames(res.lateAD_DCG_Neurons_df)=names(lateAD_DCG_Neurons.overlap)
res.lateAD_DCG_Neurons_df$adj.p=p.adjust(p = res.lateAD_DCG_Neurons_df$pval,method = "fdr")
res.lateAD_DCG_Neurons_df$Expected_NeuronsMarkers_in_DCGs=as.numeric(res.lateAD_DCG_Neurons_df$Expected_NeuronsMarkers_in_DCGs)
res.lateAD_DCG_Neurons_df$Nr_NeuronsMarkers_in_DCGs=as.numeric(res.lateAD_DCG_Neurons_df$Nr_NeuronsMarkers_in_DCGs)
res.lateAD_DCG_Neurons_df=res.lateAD_DCG_Neurons_df%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(Nr_NeuronsMarkers_in_DCGs>Expected_NeuronsMarkers_in_DCGs,true = "Over",false = "Under"))%>%filter(Representation=="Under")

res.lateAD_DCG_Oligodendrocytes=foreach(i=1:length(lateAD_DCG_Oligodendrocytes.overlap))%do%{
  
  
  data.frame(OligodendrocytesMarkers_in_DCGs=unlist(lapply(lateAD_DCG_Oligodendrocytes.overlap[i],paste,collapse=",")),
             Nr_OligodendrocytesMarkers_in_DCGs=unlist(lapply(lapply(lateAD_DCG_Oligodendrocytes.overlap,length)[i],paste,collapse=",")),
             DCGs=unlist(lapply(lateAD_DCGs.uniq,length)[[i]]),
             Expected_OligodendrocytesMarkers_in_DCGs=unname(unlist(lateAD.exp_Oligodendrocytes[i])),
             Fisher_statistic=unname(unlist(lapply(lateAD_Oligodendrocytes.fisher_results[i],function(s)s$estimate))),
             pval=unname(unlist(lapply(lateAD_Oligodendrocytes.fisher_results[i],function(p)p$p.value))),
             stringsAsFactors = F)
}
names(res.lateAD_DCG_Oligodendrocytes)=names(lateAD_DCG_Oligodendrocytes.overlap)
res.lateAD_DCG_Oligodendrocytes_df=data.frame(rbindlist(res.lateAD_DCG_Oligodendrocytes),stringsAsFactors = F)
rownames(res.lateAD_DCG_Oligodendrocytes_df)=names(lateAD_DCG_Oligodendrocytes.overlap)
res.lateAD_DCG_Oligodendrocytes_df$adj.p=p.adjust(p = res.lateAD_DCG_Oligodendrocytes_df$pval,method = "fdr")
res.lateAD_DCG_Oligodendrocytes_df$Expected_OligodendrocytesMarkers_in_DCGs=as.numeric(res.lateAD_DCG_Oligodendrocytes_df$Expected_OligodendrocytesMarkers_in_DCGs)
res.lateAD_DCG_Oligodendrocytes_df$Nr_OligodendrocytesMarkers_in_DCGs=as.numeric(res.lateAD_DCG_Oligodendrocytes_df$Nr_OligodendrocytesMarkers_in_DCGs)
res.lateAD_DCG_Oligodendrocytes_df=res.lateAD_DCG_Oligodendrocytes_df%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(Nr_OligodendrocytesMarkers_in_DCGs>Expected_OligodendrocytesMarkers_in_DCGs,true = "Over",false = "Under"))%>%filter(Representation=="Over")




#Redhead et al. viral network drivers, loss-gain genes
redhead_nnet_driver_genes.df=read.xls("/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/Redhead_AD_HerpesVirus_Dataset/mmc2.xlsx",sheet = 3,header=T)
redhead_nnet_driver_genes.list=earlyAD_DCG_nnet_driver_genes.overlap=earlyAD_DCG_nnet_driver_genes.expected=earlyAD_DCG_nnet_driver_genes.fisher=vector(mode = "list",length = 3)
names(redhead_nnet_driver_genes.list)=names(earlyAD_DCG_nnet_driver_genes.overlap)=names(earlyAD_DCG_nnet_driver_genes.expected)=names(earlyAD_DCG_nnet_driver_genes.fisher)=unique(redhead_nnet_driver_genes.df$driver_type)

for(i in 1:3){
  redhead_nnet_driver_genes.list[[i]]=redhead_nnet_driver_genes.df%>%filter(driver_type==names(redhead_nnet_driver_genes.list)[i])%>%pull(driver_symbol)%>%droplevels%>%levels
}  
  earlyAD_DCG_nnet_driver_genes.overlap[[i]]=lapply(earlyAD_DCGs.list,intersect,redhead_nnet_driver_genes.list[[i]])  
  earlyAD_DCG_nnet_driver_genes.expected[[i]]=lapply(lapply(earlyAD_DCGs.list,length),function(x)round(length(redhead_nnet_driver_genes.list[[i]])/19530*x,digits = 1))
  earlyAD_DCG_nnet_driver_genes.fisher[[i]]=vector(mode = "list",length = 12)
  names(earlyAD_DCG_nnet_driver_genes.fisher[[i]])=names(earlyAD_DCGs.list)
#   for(m in 1:length(earlyAD_DCGs.list)){
#     earlyAD_DCG_nnet_driver_genes.fisher[[i]][[m]]=fisher.test(matrix(c(length(earlyAD_DCG_nnet_driver_genes.overlap[[i]][[m]]),
#                                                                         length(redhead_nnet_driver_genes.list[[i]])-length(earlyAD_DCG_nnet_driver_genes.overlap[[i]][[m]]),
#                                                                         lapply(earlyAD_DCGs.list,length)[[m]]-length(earlyAD_DCG_nnet_driver_genes.overlap[[i]][[m]]),
#                                                                         19530-(length(redhead_nnet_driver_genes.list[[i]])-length(earlyAD_DCG_nnet_driver_genes.overlap[[i]][[m]]))-
#                                                                         lapply(earlyAD_DCGs.list,length)[[m]]-length(earlyAD_DCG_nnet_driver_genes.overlap[[i]][[m]])),nrow = 2))}
# }
#   
# 
#   earlyAD_DCG_nnet_driver_genes.fisher[[i]]=vector(mode = "list",length = length(earlyAD_DCGs.list))
#   earlyAD_DCG_nnet_driver_genes.fisher[[i]]=

  
#eQTL analysis
AD_all_eQTLs=scan("All_eQTLs.txt",what = "char",sep = "\n")
AD_HighConf_eQTLs=scan("High_confidence_eQTLs.txt",what = "char",sep = "\n")
AD_HighConf_Disease_eQTLs=scan("High_confidence_eQTLs_high_disease_effect.txt",what = "char",sep = "\n")

names(AD_all_eQTLs)=unname(mapIds(x = EnsDb.Hsapiens.v79,keys = AD_all_eQTLs,keytype = "GENEID",column = "SYMBOL"))
names(AD_HighConf_eQTLs)=unname(mapIds(x = EnsDb.Hsapiens.v79,keys = AD_HighConf_eQTLs,keytype = "GENEID",column = "SYMBOL"))
names(AD_HighConf_Disease_eQTLs)=unname(mapIds(x = EnsDb.Hsapiens.v79,keys = AD_HighConf_Disease_eQTLs,keytype = "GENEID",column = "SYMBOL"))

earlyAD_UniqDCG_HiConf_Disease_eQTLs.KEGG=lapply(lapply(earlyAD_DCGs.uniq,intersect,unname(mapIds(x = EnsDb.Hsapiens.v79,keys = AD_HighConf_Disease_eQTLs,keytype = "GENEID",column = "SYMBOL"))),function(x)enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c2)%>%data.frame())

earlyAD_DCG_eQTLs.list=lapply(earlyAD_DCGs.uniq,intersect,names(AD_all_eQTLs))
lateAD_DCG_eQTLs.list=lapply(lateAD_DCGs.uniq,intersect,names(AD_all_eQTLs))
earlyAD_DCG_eQTLs.exp=lapply(lapply(earlyAD_DCGs.uniq,length),function(x)round(length(AD_all_eQTLs)/19529*x,digits = 1))
lateAD_DCG_eQTLs.exp=lapply(lapply(lateAD_DCGs.uniq,length),function(x)round(length(AD_all_eQTLs)/19529*x,digits = 1))

earlyAD_TF_bridged_DCL.TF_downstream_subgraph=readRDS("EarlyAD_diffcoexp/earlyAD_TF_bridged_DCL_TF_downstream_subgraph.RDS")
lateAD_TF_bridged_DCL.TF_downstream_subgraph=readRDS("LateAD_diffcoexp/lateAD_TF_bridged_DCL_TF_downstream_subgraph.RDS")
earlyAD_commonTF_downstream_DCL_eQTL.list=earlyAD_commonTF_downstream_DCL_eQTL.exp=vector(mode = "list",length=length(earlyAD_TF_bridged_DCL.TF_downstream_subgraph))
names(earlyAD_commonTF_downstream_DCL_eQTL.list)=names(earlyAD_commonTF_downstream_DCL_eQTL.exp)=names(earlyAD_TF_bridged_DCL.TF_downstream_subgraph)
lateAD_commonTF_downstream_DCL_eQTL.list=vector(mode = "list",length = length(lateAD_TF_bridged_DCL.TF_downstream_subgraph))
names(lateAD_commonTF_downstream_DCL_eQTL.list)=names(lateAD_TF_bridged_DCL.TF_downstream_subgraph)

for(i in 1:length(earlyAD_TF_bridged_DCL.TF_downstream_subgraph)){
  earlyAD_commonTF_downstream_DCL_eQTL.list[[i]]=lapply(lapply(earlyAD_TF_bridged_DCL.TF_downstream_subgraph[[i]],function(x)V(x)$name),intersect,names(AD_HighConf_eQTLs))
}
for(i in 1:length(lateAD_TF_bridged_DCL.TF_downstream_subgraph)){
  lateAD_commonTF_downstream_DCL_eQTL.list[[i]]=lapply(lapply(lateAD_TF_bridged_DCL.TF_downstream_subgraph[[i]],function(x)V(x)$name),intersect,names(AD_HighConf_eQTLs))
}

earlyAD_commonTF_downstream_DCL_eQTL.fisher=vector(mode = "list",length = length(earlyAD_TF_bridged_DCL.TF_downstream_subgraph))
names(earlyAD_commonTF_downstream_DCL_eQTL.fisher)=names(earlyAD_TF_bridged_DCL.TF_downstream_subgraph)
lateAD_commonTF_downstream_DCL_eQTL.fisher=vector(mode = "list",length = length(lateAD_TF_bridged_DCL.TF_downstream_subgraph))
names(lateAD_commonTF_downstream_DCL_eQTL.fisher)=names(lateAD_TF_bridged_DCL.TF_downstream_subgraph)


for(t in 1:length(earlyAD_TF_bridged_DCL.TF_downstream_subgraph)){
  for(i in 1:length(earlyAD_DCGs.uniq)){
    earlyAD_commonTF_downstream_DCL_eQTL.fisher[[t]][[i]]=fisher.test(matrix(c(length(earlyAD_commonTF_downstream_DCL_eQTL.list[[t]][[i]]),
                                                                               length(earlyAD_DCGs.uniq[[i]])-length(earlyAD_commonTF_downstream_DCL_eQTL.list[[t]][[i]]),
                                                                               length(names(AD_all_eQTLs))-length(earlyAD_commonTF_downstream_DCL_eQTL.list[[t]][[i]]),
                                                                               19529-(length(earlyAD_DCGs.uniq[[i]])-length(earlyAD_commonTF_downstream_DCL_eQTL.list[[t]][[i]]))-(length(names(AD_all_eQTLs))-length(earlyAD_commonTF_downstream_DCL_eQTL.list[[t]][[i]]))-length(earlyAD_commonTF_downstream_DCL_eQTL.list[[t]][[i]])),nrow = 2))
    
  }
  names(earlyAD_commonTF_downstream_DCL_eQTL.fisher[[t]])=names(earlyAD_DCGs.uniq)
}









for(t in 1:length(lateAD_TF_bridged_DCL.TF_downstream_subgraph)){
  for(i in 1:length(lateAD_DCGs.uniq)){
    lateAD_commonTF_downstream_DCL_eQTL.fisher[[t]][[i]]=fisher.test(matrix(c(length(lateAD_commonTF_downstream_DCL_eQTL.list[[t]][[i]]),
                                                                              length(lateAD_DCGs.uniq[[i]])-length(lateAD_commonTF_downstream_DCL_eQTL.list[[t]][[i]]),
                                                                              length(janssen_meta_hits2)-length(lateAD_commonTF_downstream_DCL_eQTL.list[[t]][[i]]),
                                                                              19529-(length(lateAD_DCGs.uniq[[i]])-length(lateAD_commonTF_downstream_DCL_eQTL.list[[t]][[i]]))-(length(janssen_meta_hits2)-length(lateAD_commonTF_downstream_DCL_eQTL.list[[t]][[i]]))-length(lateAD_commonTF_downstream_DCL_eQTL.list[[t]][[i]])),nrow = 2))
    
  }
  names(lateAD_commonTF_downstream_DCL_eQTL.fisher[[t]])=names(lateAD_DCGs.uniq)
}

res_earlyAD_DCG_eQTLs=vector(mode = "list",length=length(earlyAD_commonTF_downstream_DCL_eQTL.fisher))
names(res_earlyAD_DCG_eQTLs)=names(earlyAD_commonTF_downstream_DCL_eQTL.fisher)
for(t in 1:length(earlyAD_commonTF_downstream_DCL_eQTL.list)){
  

  res_earlyAD_DCG_eQTLs[[t]]=foreach(i=1:length(earlyAD_commonTF_downstream_DCL_eQTL.list[[t]]))%do%{
    
    
    data.frame(eQTLs_in_DCGs=unlist(lapply(earlyAD_commonTF_downstream_DCL_eQTL.list[[t]][i],length)),
               DCGs=unlist(lapply(lateAD_DCGs.uniq,length)[[i]]),
               Expected_eQTLs_in_DCGs=unname(unlist(earlyAD_DCG_eQTLs.exp[i])),
               Fisher_statistic=unname(unlist(lapply(earlyAD_commonTF_downstream_DCL_eQTL.fisher[[t]][i],function(s)s$estimate))),
               pval=unname(unlist(lapply(earlyAD_commonTF_downstream_DCL_eQTL.fisher[[t]][i],function(p)p$p.value))),
               stringsAsFactors = F)
  }
  res_earlyAD_DCG_eQTLs[[t]]=data.frame(rbindlist(res_earlyAD_DCG_eQTLs[[t]]),stringsAsFactors = F)
  rownames(res_earlyAD_DCG_eQTLs[[t]])=names(earlyAD_commonTF_downstream_DCL_eQTL.fisher[[t]])
  
  res_earlyAD_DCG_eQTLs[[t]]$adj.p=p.adjust(p = res_earlyAD_DCG_eQTLs[[t]]$pval,method = "fdr")
  res_earlyAD_DCG_eQTLs[[t]]$Expected_eQTLs_in_DCGs=as.numeric(res_earlyAD_DCG_eQTLs[[t]]$Expected_eQTLs_in_DCGs)
  res_earlyAD_DCG_eQTLs[[t]]$eQTLs_in_DCGs=as.numeric(res_earlyAD_DCG_eQTLs[[t]]$eQTLs_in_DCGs)
  res_earlyAD_DCG_eQTLs[[t]]=res_earlyAD_DCG_eQTLs[[t]]%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(eQTLs_in_DCGs>Expected_eQTLs_in_DCGs,true = "Over",false = "Under"))%>%dplyr::filter(Representation=="Over")
  

}

res_lateAD_DCG_eQTLs=vector(mode = "list",length=length(lateAD_commonTF_downstream_DCL_eQTL.fisher))
names(res_lateAD_DCG_eQTLs)=names(lateAD_commonTF_downstream_DCL_eQTL.fisher)
for(t in 1:length(lateAD_commonTF_downstream_DCL_eQTL.list)){
  
  
  res_lateAD_DCG_eQTLs[[t]]=foreach(i=1:length(lateAD_commonTF_downstream_DCL_eQTL.list[[t]]))%do%{
    
    
    data.frame(eQTLs_in_DCGs=unlist(lapply(lateAD_commonTF_downstream_DCL_eQTL.list[[t]][i],length)),
               DCGs=unlist(lapply(lateAD_DCGs.uniq,length)[[i]]),
               Expected_eQTLs_in_DCGs=unname(unlist(lateAD_DCG_eQTLs.exp[i])),
               Fisher_statistic=unname(unlist(lapply(lateAD_commonTF_downstream_DCL_eQTL.fisher[[t]][i],function(s)s$estimate))),
               pval=unname(unlist(lapply(lateAD_commonTF_downstream_DCL_eQTL.fisher[[t]][i],function(p)p$p.value))),
               stringsAsFactors = F)
  }
  res_lateAD_DCG_eQTLs[[t]]=data.frame(rbindlist(res_lateAD_DCG_eQTLs[[t]]),stringsAsFactors = F)
  rownames(res_lateAD_DCG_eQTLs[[t]])=names(lateAD_commonTF_downstream_DCL_eQTL.fisher[[t]])
  
  res_lateAD_DCG_eQTLs[[t]]$adj.p=p.adjust(p = res_lateAD_DCG_eQTLs[[t]]$pval,method = "fdr")
  res_lateAD_DCG_eQTLs[[t]]$Expected_eQTLs_in_DCGs=as.numeric(res_lateAD_DCG_eQTLs[[t]]$Expected_eQTLs_in_DCGs)
  res_lateAD_DCG_eQTLs[[t]]$eQTLs_in_DCGs=as.numeric(res_lateAD_DCG_eQTLs[[t]]$eQTLs_in_DCGs)
  res_lateAD_DCG_eQTLs[[t]]=res_lateAD_DCG_eQTLs[[t]]%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(eQTLs_in_DCGs>Expected_eQTLs_in_DCGs,true = "Over",false = "Under"))%>%dplyr::filter(Representation=="Over")
  
  
}

# eQTL enrichment in TF-targets
eQTLs_target.overlap=Filter(f = function(x)length(x)>0,lapply(humanRegnetwork,intersect,names(AD_HighConf_Disease_eQTLs)))

res=foreach(i = 1:length(eQTLs_target.overlap))%do%{
  fisher.test(matrix(c(length(eQTLs_target.overlap[[i]]),
                       length(humanRegnetwork[names(eQTLs_target.overlap)][[i]])-length(eQTLs_target.overlap[i]),
                       length(AD_all_eQTLs)-length(eQTLs_target.overlap[i]),
                       19529-(length(humanRegnetwork[names(eQTLs_target.overlap)][i])-length(eQTLs_target.overlap[i]))-(length(AD_all_eQTLs)-length(eQTLs_target.overlap[[i]]))-length(eQTLs_target.overlap[[i]])),nrow = 2))
}

res.df=data.frame(TFs=names(eQTLs_target.overlap),
                  TF_targets=unlist(lapply(humanRegnetwork[names(eQTLs_target.overlap)],length)),
                  eQTLs_in_TFtargets=unlist(lapply(eQTLs_target.overlap,length)),
                  exp_eQTLs_in_TFtargets=unlist(lapply(humanRegnetwork[names(eQTLs_target.overlap)],function(x)length(AD_all_eQTLs)/19529*length(x))),
                  pval=unlist(lapply(res,function(x)x$`p.value`)),
                  adj.p=p.adjust(p = unlist(lapply(res,function(x)x$`p.value`)),method = "fdr"),stringsAsFactors = F)
res.df=res.df%>%mutate(Representation=if_else(eQTLs_in_TFtargets>exp_eQTLs_in_TFtargets,true = "Over",false = "Under"))
res.df=res.df%>%dplyr::filter(`adj.p`<=0.1)%>%dplyr::filter(Representation=="Over")


#LRNs ROSMAP, Mostafavi et al.
lrn_CogDec=read.xls("/Users/sandeepamberkar/Work/Data/LRNs_ROSMAP_Mostafavi/table 3.xlsx",sheet = 2)

earlyAD_Uniq_lrn_CogDec.overlap=lapply(earlyAD_DCGs.uniq,intersect,lrn_CogDec$symbol%>%levels)
earlyAD_Uniq_lrn_CogDec.exp=lapply(earlyAD_DCGs.uniq,function(x)length(lrn_CogDec$symbol%>%levels)/(19529)*length(x))
lateAD_Uniq_lrn_CogDec.overlap=lapply(lateAD_DCGs.uniq,intersect,lrn_CogDec$symbol%>%levels)
lateAD_Uniq_lrn_CogDec.exp=lapply(lateAD_DCGs.uniq,function(x)length(lrn_CogDec$symbol%>%levels)/(19529)*length(x))

earlyAD_Uniq_lrn_CogDec.fisher=lateAD_Uniq_lrn_CogDec.fisher=vector(mode = "list",length = length(earlyAD_DCGs.uniq))


earlyAD_Uniq_lrn_CogDec.fisher=foreach(i = 1:length(earlyAD_Uniq_lrn_CogDec.overlap))%do%{
  fisher.test(matrix(c(length(earlyAD_Uniq_lrn_CogDec.overlap[[i]]),
                       length(lrn_CogDec$symbol%>%levels)-length(earlyAD_Uniq_lrn_CogDec.overlap[[i]]),
                       length(earlyAD_DCGs.uniq[[i]])-length(earlyAD_Uniq_lrn_CogDec.overlap[[i]]),
                       18517-(length(lrn_CogDec$symbol%>%levels)-length(earlyAD_Uniq_lrn_CogDec.overlap[[i]]))-(length(earlyAD_DCGs.uniq[[i]])-length(earlyAD_Uniq_lrn_CogDec.overlap[[i]]))-length(earlyAD_Uniq_lrn_CogDec.overlap[[i]])),nrow = 2))
}

lateAD_Uniq_lrn_CogDec.fisher=foreach(i = 1:length(lateAD_Uniq_lrn_CogDec.overlap))%do%{
  fisher.test(matrix(c(length(lateAD_Uniq_lrn_CogDec.overlap[[i]]),
                       length(lrn_CogDec$symbol%>%levels)-length(lateAD_Uniq_lrn_CogDec.overlap[[i]]),
                       length(lateAD_DCGs.uniq[[i]])-length(lateAD_Uniq_lrn_CogDec.overlap[[i]]),
                       19529-(length(lrn_CogDec$symbol%>%levels)-length(lateAD_Uniq_lrn_CogDec.overlap[[i]]))-(length(lateAD_DCGs.uniq[[i]])-length(lateAD_Uniq_lrn_CogDec.overlap[[i]]))-length(lateAD_Uniq_lrn_CogDec.overlap[[i]])),nrow = 2))
}

res_lrn_CogDec_DEGs_earlyAD.df=data.frame(LRN_CogDec_earlyAD_DCGs=unlist(lapply(earlyAD_Uniq_lrn_CogDec.overlap,length)),
                                          EarlyAD_DCGs_restricted=unlist(lapply(earlyAD_DCGs.uniq,length)),
                                          exp_LRN_CogDec_DEGs=unlist(earlyAD_Uniq_lrn_CogDec.exp),
                                          pval=unlist(lapply(earlyAD_Uniq_lrn_CogDec.fisher,function(x)x$`p.value`)),
                                          adj.p=p.adjust(p = unlist(lapply(earlyAD_Uniq_lrn_CogDec.fisher,function(x)x$`p.value`)),method = "fdr"),stringsAsFactors = F)
res_lrn_CogDec_DEGs_earlyAD.df=res_lrn_CogDec_DEGs_earlyAD.df%>%mutate(Representation=if_else(LRN_CogDec_earlyAD_DCGs>exp_LRN_CogDec_DEGs,true = "Over",false = "Under"))
res_lrn_CogDec_DEGs_earlyAD.df=res_lrn_CogDec_DEGs_earlyAD.df%>%dplyr::filter(pval<=0.05)%>%dplyr::filter(Representation=="Over")

res_lrn_CogDec_DEGs_lateAD.df=data.frame(LRN_CogDec_lateAD_DCGs=unlist(lapply(lateAD_Uniq_lrn_CogDec.overlap,length)),
                                         lateAD_DCGs_restricted=unlist(lapply(lateAD_DCGs.uniq,length)),
                                         exp_LRN_CogDec_DEGs=unlist(lateAD_Uniq_lrn_CogDec.exp),
                                         pval=unlist(lapply(lateAD_Uniq_lrn_CogDec.fisher,function(x)x$`p.value`)),
                                         adj.p=p.adjust(p = unlist(lapply(lateAD_Uniq_lrn_CogDec.fisher,function(x)x$`p.value`)),method = "fdr"),stringsAsFactors = F)
res_lrn_CogDec_DEGs_lateAD.df=res_lrn_CogDec_DEGs_lateAD.df%>%mutate(Representation=if_else(LRN_CogDec_lateAD_DCGs>exp_LRN_CogDec_DEGs,true = "Over",false = "Under"))
res_lrn_CogDec_DEGs_lateAD.df=res_lrn_CogDec_DEGs_lateAD.df%>%dplyr::filter(pval<=0.05)%>%dplyr::filter(Representation=="Over")

names(earlyAD_Uniq_lrn_CogDec.fisher)=names(lateAD_Uniq_lrn_CogDec.fisher)=names(earlyAD_DCGs.uniq)



#Comparison with Berchtold et al.

berchtold_hippocampus.diffcoexp=readRDS("/Users/sandeepamberkar/Work/Data/AD_GSE48350/hippocampusdiffcoexp.RDS")

berchtold_hippocampus.diffcoexp$DCGs=berchtold_hippocampus.diffcoexp$DCGs%>%dplyr::filter(!grepl(pattern = "///|\\.|\\-",x = Gene))
rand_msbb_earlyAD_DCGs.hipp=replicate(n = 1000,expr = sample(x = rownames(earlyAD_samples.exprs$Hippocampus),size = length(earlyAD_DCGs.uniq$Hippocampus),replace = F),simplify = T)
rand_btold_AD_DCGs.hipp=replicate(n = 1000,expr = sample(x = rownames(earlyAD_samples.exprs$Hippocampus),size = length(berchtold_hippocampus.diffcoexp$DCGs$Gene),replace = F),simplify = T)
rand_msbb_lateAD_DCGs.hipp=replicate(n = 1000,expr = sample(x = rownames(lateAD_samples.exprs$Hippocampus),size = length(lateAD_DCGs.uniq$Hippocampus),replace = F),simplify = T)
rand_msbb_btold_earlyAD.commonDCGs=rand_msbb_btold_lateAD.commonDCGs=list()
for(i in 1:1000){
  rand_msbb_btold_earlyAD.commonDCGs[[i]]=intersect(rand_msbb_earlyAD_DCGs.hipp[,i],rand_btold_AD_DCGs.hipp[,i])
}
for(i in 1:1000){
  rand_msbb_btold_lateAD.commonDCGs[[i]]=intersect(rand_msbb_lateAD_DCGs.hipp[,i],rand_btold_AD_DCGs.hipp[,i])
}

# fisher.test(matrix(c(427,length(earlyAD_DCGs.uniq$Hippocampus)-427,length(berchtold_hippocampus.diffcoexp$DCGs$Gene)-427,length(rownames(earlyAD_samples.exprs$Hippocampus))-427-length(earlyAD_DCGs.uniq$Hippocampus)-427-length(berchtold_hippocampus.diffcoexp$DCGs$Gene)-427),nrow = 2),alternative = "g")
# fisher.test(matrix(c(1292,length(earlyAD_DCGs.list$Hippocampus)-1292,length(berchtold_hippocampus.diffcoexp$DCGs$Gene)-1292,length(rownames(earlyAD_samples.exprs$Hippocampus))-1292-length(earlyAD_DCGs.list$Hippocampus)-1292-length(berchtold_hippocampus.diffcoexp$DCGs$Gene)-1292),nrow = 2))
# 
# fisher.test(matrix(c(1098,length(lateAD_DCGs.uniq$Hippocampus)-1098,length(berchtold_hippocampus.diffcoexp$DCGs$Gene)-1098,length(rownames(lateAD_samples.exprs$Hippocampus))-1098-length(lateAD_DCGs.uniq$Hippocampus)-1098-length(berchtold_hippocampus.diffcoexp$DCGs$Gene)-1098),nrow = 2),alternative = "g")
# fisher.test(matrix(c(1967,7894-1967,4856-1967,18320-7894-1967-4856-1967-1967),nrow = 2),alternative = "g")
 
#Test significance of overlap between early AD/late AD DCGs hippocampus & Berchtold hippocampus DCGs
library(GeneOverlap)
earlyAD.obj=newGeneOverlap(earlyAD_DCGs.list$Hippocampus,berchtold_hippocampus.diffcoexp$DCGs$Gene,genome.size=19947)
lateAD.obj=newGeneOverlap(lateAD_DCGs.list$Hippocampus,berchtold_hippocampus.diffcoexp$DCGs$Gene,genome.size=19947)

#Sample DCL
earlyAD_gene1=earlyAD_samples.exprs$Frontal_Pole[which(rownames(earlyAD_samples.exprs$Frontal_Pole)=="S100A3"),]
earlyAD_gene2=earlyAD_samples.exprs$Frontal_Pole[which(rownames(earlyAD_samples.exprs$Frontal_Pole)=="ADGRG2"),]
earlyAD_genes=data.frame(rbind(t(rbind(S100A3=earlyAD_gene1[earlyAD_samples$Frontal_Pole$SampleType=="CONTROL"],
                                       ADGRG2=earlyAD_gene2[earlyAD_samples$Frontal_Pole$SampleType=="CONTROL"])),
                               t(rbind(S100A3=earlyAD_gene1[earlyAD_samples$Frontal_Pole$SampleType=="AD"],
                                       ADGRG2=earlyAD_gene2[earlyAD_samples$Frontal_Pole$SampleType=="AD"]))),stringsAsFactors = F)
earlyAD_genes$SampleType=c(rep("No-AD",10),rep("early-AD",6))


#Compare with Tanzi rare variants

tanzi_variants=scan("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Gene_Lists/Tanzi_Ranked_SNP_Genes.txt",sep = "\n",what = "char")
tanzi_Region_genes.df=read.table("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Gene_Lists/Region_genes_top_1000.txt",header = T,sep = "\t",as.is = T)
tanzi_Region_genes.df=tanzi_Region_genes.df[which(as.numeric(unlist(lapply(strsplit(tanzi_Region_genes.df$Region_Rank,split = "_"),`[[`,3)))<=1000),]
tanzi_SNV_Region.list=union(unique(tanzi_Region_genes.df$Gene_GREAT),tanzi_variants)

tanzi_variants_earlyAD_DCGs_uniq.list=lapply(earlyAD_DCGs.uniq,intersect,tanzi_SNV_Region.list)
tanzi_variants_earlyAD_DCGs_uniq.C2=lapply(tanzi_variants_earlyAD_DCGs_uniq.list,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c2),stringsAsFactors = F))
tanzi_variants_earlyAD_DCGs_uniq.BP=lapply(tanzi_variants_earlyAD_DCGs_uniq.list,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c5.BP),stringsAsFactors = F))
tanzi_variants_earlyAD_DCGs_uniq.CC=lapply(tanzi_variants_earlyAD_DCGs_uniq.list,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c5.CC),stringsAsFactors = F))
tanzi_variants_earlyAD_DCGs_uniq.MF=lapply(tanzi_variants_earlyAD_DCGs_uniq.list,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c5.MF),stringsAsFactors = F))

tanzi_variants_earlyAD_DCGs_uniq.C2=lapply(tanzi_variants_earlyAD_DCGs_uniq.C2,function(x)x%>%filter(grepl(pattern = "KEGG|REACTOME",x = ID)&(Count>1)))

tanzi_variants_earlyAD_DCGs_uniq.FuncEnrichment=vector(mode = "list",length = length(earlyAD_DCGs.uniq))
names(tanzi_variants_earlyAD_DCGs_uniq.FuncEnrichment)=names(earlyAD_DCGs.uniq)
for(i in 1:12){
  tanzi_variants_earlyAD_DCGs_uniq.FuncEnrichment[[i]]=rbind(tanzi_variants_earlyAD_DCGs_uniq.C2[[i]]%>%filter(grepl(pattern = "KEGG|REACTOME",x = ID)&(Count>2))%>%mutate(Enriched_in_Janssen=if_else(ID%in%janssen_AD.sigGO$Gene.set,"Y","N"))%>%arrange(ID)%>%dplyr::select(1,3:9),
                                                             tanzi_variants_earlyAD_DCGs_uniq.BP[[i]]%>%filter(Count>2)%>%arrange(ID)%>%dplyr::select(1,3:9),
                                                             tanzi_variants_earlyAD_DCGs_uniq.CC[[i]]%>%filter(Count>2)%>%arrange(ID)%>%dplyr::select(1,3:9),
                                                             tanzi_variants_earlyAD_DCGs_uniq.MF[[i]]%>%filter(Count>2)%>%arrange(ID)%>%dplyr::select(1,3:9))
  tanzi_variants_earlyAD_DCGs_uniq.FuncEnrichment[[i]]$SNV_Genes=unlist(lapply(tanzi_variants_earlyAD_DCGs_uniq.FuncEnrichment[[i]]$geneID,function(x)paste(intersect(strsplit(x = x,split = "/")[[1]],tanzi_variants),collapse = ",")))
  tanzi_variants_earlyAD_DCGs_uniq.FuncEnrichment[[i]]$Region_Genes=unlist(lapply(tanzi_variants_earlyAD_DCGs_uniq.FuncEnrichment[[i]]$geneID,function(x)paste(intersect(strsplit(x = x,split = "/")[[1]],unique(tanzi_Region_genes.df$Gene_GREAT)),collapse = ",")))
  fwrite(tanzi_variants_earlyAD_DCGs_uniq.FuncEnrichment[[i]],paste(names(tanzi_variants_earlyAD_DCGs_uniq.FuncEnrichment)[i],"EarlyAD_WGS_SNV_restrictedDCGs_FunctionalEnrichment.txt",sep = "_"),sep="\t",col.names=T)
}

tanzi_variants_lateAD_DCGs_uniq.list=lapply(lateAD_DCGs.uniq,intersect,tanzi_SNV_Region.list)
tanzi_variants_lateAD_DCGs_uniq.C2=lapply(tanzi_variants_lateAD_DCGs_uniq.list,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c2),stringsAsFactors = F))
tanzi_variants_lateAD_DCGs_uniq.BP=lapply(tanzi_variants_lateAD_DCGs_uniq.list,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c5.BP),stringsAsFactors = F))
tanzi_variants_lateAD_DCGs_uniq.CC=lapply(tanzi_variants_lateAD_DCGs_uniq.list,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c5.CC),stringsAsFactors = F))
tanzi_variants_lateAD_DCGs_uniq.MF=lapply(tanzi_variants_lateAD_DCGs_uniq.list,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c5.MF),stringsAsFactors = F))

tanzi_variants_lateAD_DCGs_uniq.C2=lapply(tanzi_variants_lateAD_DCGs_uniq.C2,function(x)x%>%filter(grepl(pattern = "KEGG|REACTOME",x = ID)&(Count>1)))

tanzi_variants_lateAD_DCGs_uniq.FuncEnrichment=vector(mode = "list",length = length(lateAD_DCGs.uniq))
names(tanzi_variants_lateAD_DCGs_uniq.FuncEnrichment)=names(lateAD_DCGs.uniq)
for(i in 1:12){
  tanzi_variants_lateAD_DCGs_uniq.FuncEnrichment[[i]]=rbind(tanzi_variants_lateAD_DCGs_uniq.C2[[i]]%>%filter(grepl(pattern = "KEGG|REACTOME",x = ID)&(Count>2))%>%mutate(Enriched_in_Janssen=if_else(ID%in%janssen_AD.sigGO$Gene.set,"Y","N"))%>%arrange(ID)%>%dplyr::select(1,3:9),
                                                            tanzi_variants_lateAD_DCGs_uniq.BP[[i]]%>%filter(Count>2)%>%arrange(ID)%>%dplyr::select(1,3:9),
                                                            tanzi_variants_lateAD_DCGs_uniq.CC[[i]]%>%filter(Count>2)%>%arrange(ID)%>%dplyr::select(1,3:9),
                                                            tanzi_variants_lateAD_DCGs_uniq.MF[[i]]%>%filter(Count>2)%>%arrange(ID)%>%dplyr::select(1,3:9))
  tanzi_variants_lateAD_DCGs_uniq.FuncEnrichment[[i]]$SNV_Genes=unlist(lapply(tanzi_variants_lateAD_DCGs_uniq.FuncEnrichment[[i]]$geneID,function(x)paste(intersect(strsplit(x = x,split = "/")[[1]],tanzi_variants),collapse = ",")))
  tanzi_variants_lateAD_DCGs_uniq.FuncEnrichment[[i]]$Region_Genes=unlist(lapply(tanzi_variants_lateAD_DCGs_uniq.FuncEnrichment[[i]]$geneID,function(x)paste(intersect(strsplit(x = x,split = "/")[[1]],unique(tanzi_Region_genes.df$Gene_GREAT)),collapse = ",")))
  fwrite(tanzi_variants_lateAD_DCGs_uniq.FuncEnrichment[[i]],paste(names(tanzi_variants_lateAD_DCGs_uniq.FuncEnrichment)[i],"lateAD_WGS_SNV_restrictedDCGs_FunctionalEnrichment.txt",sep = "_"),sep="\t",col.names=T)
}

#Functional convergence at the pathway level
tanzi_SNP_genes=scan("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Gene_Lists/LarsAD_SNPgenes.txt",what = "char",sep = "\n")
tanzi_Region_genes=scan("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Gene_Lists/LarsAD_Regiongenes.txt",what = "char",sep = "\n")
tanzi_Region_genes.df=read.table("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Gene_Lists/Region_genes_top_1000.txt",header = T,sep = "\t",as.is = T)
tanzi_Region_genes.df=tanzi_Region_genes.df[which(as.numeric(unlist(lapply(strsplit(tanzi_Region_genes.df$Region_Rank,split = "_"),`[[`,3)))<=1000),]

tanzi_variants.C2=data.frame(enricher(gene = tanzi_variants,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c2),stringsAsFactors = F)
tanzi_variants.BP=data.frame(enricher(gene = tanzi_variants,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c5.BP),stringsAsFactors = F)
tanzi_variants.CC=data.frame(enricher(gene = tanzi_variants,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c5.CC),stringsAsFactors = F)
tanzi_variants.MF=data.frame(enricher(gene = tanzi_variants,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c5.MF),stringsAsFactors = F)

tanzi_variants.FuncEnrichment=rbind(tanzi_variants.C2%>%filter(grepl(pattern = "KEGG|REACTOME",x = ID))%>%mutate(Enriched_in_Janssen=if_else(ID%in%janssen_AD.sigGO$Gene.set,"Y","N"))%>%arrange(ID)%>%dplyr::select(1,3:9),
                                    tanzi_variants.BP%>%arrange(ID)%>%dplyr::select(1,3:9),
                                    tanzi_variants.CC%>%arrange(ID)%>%dplyr::select(1,3:9))
tanzi_variants.FuncEnrichment$SNP_Genes=unlist(lapply(tanzi_variants.FuncEnrichment$geneID,function(x)paste(intersect(strsplit(x = x,split = "/")[[1]],tanzi_SNP_genes),collapse = ",")))
tanzi_variants.FuncEnrichment$Region_Genes=unlist(lapply(tanzi_variants.FuncEnrichment$geneID,function(x)paste(intersect(strsplit(x = x,split = "/")[[1]],tanzi_Region_genes),collapse = ",")))

earlyAD_DCG_Uniq_convergent_functions_tanzi.list=lapply(earlyAD_DCGs.Uniq_c2,function(x)paste(intersect(x$ID,tanzi_variants.FuncEnrichment$ID),collapse = ";"))
lateAD_DCG_Uniq_convergent_functions_tanzi.list=lapply(lateAD_DCGs.Uniq_c2,function(x)paste(intersect(x$ID,tanzi_variants.FuncEnrichment$ID),collapse = ";"))
common_DCG_convergent_functions_tanzi.list=lapply(common_DCGs_c2,function(x)paste(intersect(x$ID,tanzi_variants.FuncEnrichment$ID),collapse = ";"))
earlyAD_DCG_Uniq_convergent_functions_janssen.list=lapply(earlyAD_DCGs.Uniq_c2,function(x)paste(intersect(x$ID,janssen_AD.sigGO$Gene.set),collapse = ";"))
lateAD_DCG_Uniq_convergent_functions_janssen.list=lapply(lateAD_DCGs.Uniq_c2,function(x)paste(intersect(x$ID,janssen_AD.sigGO$Gene.set),collapse = ";"))
common_DCG_convergent_functions_janssen.list=lapply(common_DCGs_c2,function(x)paste(intersect(x$ID,janssen_AD.sigGO$Gene.set),collapse = ";"))

DCG_Genetic_convergent_functions.df=data.frame(EarlyAD_DCG_uniq_WGS_AD=unlist(earlyAD_DCG_Uniq_convergent_functions_tanzi.list),
                                               LateAD_DCG_uniq_WGS_AD=unlist(lateAD_DCG_Uniq_convergent_functions_tanzi.list),
                                               EarlyAD_DCG_uniq_GWAS_AD=unlist(earlyAD_DCG_Uniq_convergent_functions_janssen.list),
                                               LateAD_DCG_uniq_GWAS_AD=unlist(lateAD_DCG_Uniq_convergent_functions_janssen.list),stringsAsFactors = F)

#Tsai snRNAseq DEGs
tsai_snRNAseq_DEGs_data.list=vector(mode = 'list',length = 6)

tsai_snRNAseq_DEGs_data.list[[1]]=read.xls('/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/Tsai-snRNAseq-ROSMAP/NoAD-earlyAD-Path.xls',sheet = 1,skip=1)
tsai_snRNAseq_DEGs_data.list[[2]]=read.xls('/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/Tsai-snRNAseq-ROSMAP/NoAD-earlyAD-Path.xls',sheet = 2,skip=1)
tsai_snRNAseq_DEGs_data.list[[3]]=read.xls('/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/Tsai-snRNAseq-ROSMAP/NoAD-earlyAD-Path.xls',sheet = 3,skip=1)
tsai_snRNAseq_DEGs_data.list[[4]]=read.xls('/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/Tsai-snRNAseq-ROSMAP/NoAD-earlyAD-Path.xls',sheet = 4,skip=1)
tsai_snRNAseq_DEGs_data.list[[5]]=read.xls('/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/Tsai-snRNAseq-ROSMAP/NoAD-earlyAD-Path.xls',sheet = 5,skip=1)
tsai_snRNAseq_DEGs_data.list[[6]]=read.xls('/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/Tsai-snRNAseq-ROSMAP/NoAD-earlyAD-Path.xls',sheet = 6,skip=1)


tsai_snRNAseq_DEGs.list=lapply(lapply(tsai_snRNAseq_DEGs_data.list,function(x)x%>%filter(.[[8]]=='TRUE'&.[[9]]=='TRUE')),function(x)droplevels(x[,1]))
earlyAD_DCG_tsai_snRNAseq_DEG.overlap=lapply(tsai_snRNAseq_DEGs.list,intersect,earlyAD_DCGs.uniq$Dorsolateral_Prefrontal_Cortex)
tsai_snRNAseq_DEGs_earlyAD_DLPFC_DCGs.list=lapply(tsai_snRNAseq_DEGs.list,intersect,earlyAD_DCGs.uniq$Dorsolateral_Prefrontal_Cortex)
tsai_snRNAseq_DEGs_earlyAD_DLPFC_DCGs.c2=lapply(tsai_snRNAseq_DEGs_earlyAD_DLPFC_DCGs.list,function(x)data.frame(enricher(gene = x,pvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = msigdb_c2),stringsAsFactors = F))
tsai_snRNAseq_DEGs_earlyAD_DLPFC_DCGs.KEGG=lapply(tsai_snRNAseq_DEGs_earlyAD_DLPFC_DCGs.c2,function(x)x%>%filter(grepl(pattern = 'KEGG',ID)))
earlyAD_DCG_tsai_snRNAseq_DEG_enrichment.df=rbindlist(res.earlyAD_DCG_tsai_snRNAseq_DEG)
rownames(earlyAD_DCG_tsai_snRNAseq_DEG_enrichment.df)