library(org.Hs.eg.db)
library(dplyr)
library(magrittr)
library(data.table)
library(DCGL)
library(clusterProfiler)

dbs <- listEnrichrDbs()
kegg_dbs=grep(pattern = "KEGG",x = dbs$libraryName,value = T)[3]
biocarta_dbs=grep(pattern = "BioCarta|Panther",x = dbs$libraryName,value = T)[4]
panther_dbs=grep(pattern = "BioCarta|Panther",x = dbs$libraryName,value = T)[5]
setwd("/Users/sandeepamberkar/Work/Data/MSMM-MSBB-HBTRC-Data/MSBB_Array19/GSE84422/Lobe_Analysis/Dyscorrelation_Landscape_Results/")
msbb_gse84422_GPL96_97_samplesToAnalyse.exprs=readRDS("msbb_gse84422_GPL96_97_samplesToAnalyse_exprs.RDS")
names(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs)=gsub(pattern = " ",replacement = "_",x = names(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs))
msbb_gse84422_GPL570_samplesToAnalyse.exprs=readRDS("msbb_gse84422_GPL570_samplesToAnalyse_exprs.RDS")
names(msbb_gse84422_GPL570_samplesToAnalyse.exprs)=gsub(pattern = " ",replacement = "_",x = names(msbb_gse84422_GPL570_samplesToAnalyse.exprs))


setwd("/Users/sandeepamberkar/Work/Data/MSMM-MSBB-HBTRC-Data/MSBB_Array19/GSE84422/Lobe_Analysis/Dyscorrelation_Landscape_Results/Early_Diffcoexp_CDR/")
early_diffcoexp_files=system("find . |grep 'CDR'|sort",intern = T)
msbb_gse84422_early_diffcoexp_results=msbb_gse84422_earlyCDR.DCGs=msbb_gse84422_earlyCDR.DCLs=msbb_gse84422_earlyCDR.DRGs=msbb_gse84422_earlyCDR.DRLs=vector(mode = "list",length = length(early_diffcoexp_files))
regnet_tf2target.HGNC=fread("/Users/sandeepamberkar/Work/Data/TF_Databases/RegNetwork_human_regulators2.txt",header = T,sep = "\t",showProgress = T,data.table = F)%>%dplyr::filter(evidence=="Experimental")%>%dplyr::select(c(regulator_symbol,target_symbol))

for(i in 1:length(msbb_gse84422_early_diffcoexp_results)){
  msbb_gse84422_early_diffcoexp_results[[i]]=readRDS(early_diffcoexp_files[i])
}
names(msbb_gse84422_early_diffcoexp_results)=names(msbb_gse84422_earlyCDR.DCGs)=names(msbb_gse84422_earlyCDR.DCLs)=names(msbb_gse84422_earlyCDR.DRGs)=names(msbb_gse84422_earlyCDR.DRLs)=gsub(pattern = "./",replacement = "",x = gsub(pattern = "_diffcoexp_CDR0_CDR1.RDS",replacement = "",x = early_diffcoexp_files))
msbb_gse84422_earlyCDR.DCGs=lapply(msbb_gse84422_early_diffcoexp_results,function(x)x$DCGs)
msbb_gse84422_earlyCDR.DCLs=lapply(msbb_gse84422_early_diffcoexp_results,function(x)x$DCLs)
msbb_gse84422_earlyCDR.DCGs_list=lapply(msbb_gse84422_earlyCDR.DCGs,function(x)x%>%dplyr::filter(q<=0.05)%>%pull(Gene))
msbb_gse84422_earlyCDR_total_DCGs=lapply(msbb_gse84422_earlyCDR.DCGs_list,length)

msbb_gse84422_earlyCDR_DCGs.KEGG=lapply(msbb_gse84422_earlyCDR.DCGs_list,function(x)enrichr(genes = x,databases = kegg_dbs)[[1]]%>%dplyr::filter(Adjusted.P.value<=0.1))

msbb_gse84422_earlyCDR.DRsort=vector(mode = "list",length = 19)
names(msbb_gse84422_earlyCDR.DRsort)=names(msbb_gse84422_earlyCDR.DCGs)
for(i in c(2:9,11:19)){
  msbb_gse84422_earlyCDR.DRsort[[i]]=DRsort(DCGs = msbb_gse84422_earlyCDR.DCGs[[i]],DCLs = msbb_gse84422_earlyCDR.DCLs[[i]],tf2target = regnet_tf2target.HGNC,expGenes =rownames(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs$Frontal_Pole))  
}
for(i in c(1,10)){
  msbb_gse84422_earlyCDR.DRsort[[i]]=DRsort(DCGs = msbb_gse84422_earlyCDR.DCGs[[i]],DCLs = msbb_gse84422_earlyCDR.DCLs[[i]],tf2target = regnet_tf2target.HGNC,expGenes =rownames(msbb_gse84422_GPL570_samplesToAnalyse.exprs$Amygdala))  
}

msbb_gse84422_earlyCDR.DRGs=lapply(msbb_gse84422_earlyCDR.DRsort,function(x)x$DRGs)
msbb_gse84422_earlyCDR.DRGs_filtered=lapply(msbb_gse84422_earlyCDR.DRsort,function(x)x$DRGs%>%dplyr::filter(DCGisTF=="TRUE"&q<=0.1))
msbb_gse84422_earlyCDR.DRGs_list=lapply(msbb_gse84422_earlyCDR.DRGs,function(x)x%>%pull(DCG))
msbb_gse84422_earlyCDR.DRLs=lapply(msbb_gse84422_earlyCDR.DRsort,function(x)x$DRLs%>%dplyr::filter(q.diffcor<=0.1))
msbb_gse84422_earlyCDR.TF_bridged_DCLs=lapply(msbb_gse84422_earlyCDR.DRsort,function(x)x$TF_bridged_DCL%>%dplyr::filter(q.diffcor<=0.1))

msbb_gse84422_earlyCDR.RS_DCGs_list=msbb_gse84422_earlyCDR.RS_DRGs_list=vector(mode = "list",length = 19)
names(msbb_gse84422_earlyCDR.RS_DCGs_list)=names(msbb_gse84422_earlyCDR.RS_DRGs_list)=names(msbb_gse84422_earlyCDR.DRGs_list)

for(i in 1:length(msbb_gse84422_earlyCDR.RS_DCGs_list)){
  msbb_gse84422_earlyCDR.RS_DCGs_list[[i]]=setdiff(msbb_gse84422_earlyCDR.DCGs_list[[i]],(Reduce(union,lapply(msbb_gse84422_earlyCDR.DCGs_list,function(x)intersect(x,msbb_gse84422_earlyCDR.DCGs_list[[i]]))[-i])))
}
msbb_gse84422_earlyCDR_RS_DCGs.KEGG=lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,function(x)enrichr(genes = x,databases = kegg_dbs)[[1]]%>%dplyr::filter(Adjusted.P.value<=0.1))
msbb_gse84422_earlyCDR_RS_DCGs.KEGG=Filter(f = function(x)dim(x)[1]>0,x = msbb_gse84422_earlyCDR_RS_DCGs.KEGG)


msbb_gse84422_earlyCDR_DCG.CompMatrix=matrix(NA,ncol=length(msbb_gse84422_earlyCDR.DCGs_list),nrow=length(msbb_gse84422_earlyCDR.DCGs_list))
for(i in 1:length(msbb_gse84422_earlyCDR.DCGs_list)){
  msbb_gse84422_earlyCDR_DCG.CompMatrix[i,]=unlist(lapply(msbb_gse84422_earlyCDR.DCGs_list,function(x)length(intersect(x,msbb_gse84422_earlyCDR.DCGs_list[[i]]))))
}
diag(msbb_gse84422_earlyCDR_DCG.CompMatrix)=0
rownames(msbb_gse84422_earlyCDR_DCG.CompMatrix)=colnames(msbb_gse84422_earlyCDR_DCG.CompMatrix)=c("AMYG","AC","CN","DLPFC","FP","HIPP","IFG","ITG","MTG","NAc","OVC","PHG","PCC","PCG","PFC","PTMN","SPL","STG","TP")

#Compare with Case-Control DCGs



msbb_gse84422_diffcoexp_results_files=list.files(path = "Lobe_Analysis/Dyscorrelation_Landscape_Results/",pattern = "diffcoexp",full.names = T)%>%grep(pattern = ".RDS",value = T)%>%sort
msbb_gse84422_diffcoexp_results=vector(mode = "list",length = length(msbb_gse84422_diffcoexp_results_files))
for(f in 1:19){
  msbb_gse84422_diffcoexp_results[[f]]=readRDS(msbb_gse84422_diffcoexp_results_files[f])  
}

msbb_gse84422.DCGs=lapply(msbb_gse84422_diffcoexp_results,function(x)x$DCGs)
msbb_gse84422.DCGs_list=lapply(msbb_gse84422.DCGs,function(x)x%>%dplyr::filter(q<=0.05)%>%pull(Gene))
msbb_gse84422.RS_DCGs_list=vector(mode = "list",length = 19)
names(msbb_gse84422.RS_DCGs_list)=names(msbb_gse84422.DCGs_list)
for(i in 1:length(msbb_gse84422.DCGs_list)){
  msbb_gse84422.RS_DCGs_list[[i]]=setdiff(msbb_gse84422.DCGs_list[[i]],(Reduce(union,lapply(msbb_gse84422.DCGs_list,function(x)intersect(x,msbb_gse84422.DCGs_list[[i]]))[-i])))
}
msbb_gse84422.DCGs_list_Entrez=lapply(msbb_gse84422.DCGs_list,function(y)unname(mapIds(x = org.Hs.eg.db,keys = y,keytype = "SYMBOL",column = "ENTREZID")))
msbb_gse84422.RS_DCGs_list_Entrez=lapply(msbb_gse84422.RS_DCGs_list,function(y)unname(mapIds(x = org.Hs.eg.db,keys = y,keytype = "SYMBOL",column = "ENTREZID")))
msbb_gse84422.DCLs=lapply(msbb_gse84422_diffcoexp_results,function(x)x$DCLs%>%dplyr::filter(q.diffcor<=0.05)%>%rownames_to_column("GenePair"))
msbb_gse84422_total_DCGs=unlist(lapply(msbb_gse84422.DCGs_list,length))
msbb_gse84422_total_DCLs=unlist(lapply(msbb_gse84422.DCLs,function(x)dim(x)[1]))

msbb_gse84422_DCG.CompMatrix=matrix(NA,ncol=length(msbb_gse84422.DCGs_list),nrow=length(msbb_gse84422.DCGs_list))
rownames(msbb_gse84422_DCG.CompMatrix)=colnames(msbb_gse84422_DCG.CompMatrix)=c("AMYG","AC","CN","DLPFC","FP","HIPP","IFG","ITG","MTG","NAc","OVC","PHG","PCC","PCG","FC","PTMN","SPL","STG","TP")

msbb_gse84422.DRsort=vector(mode = "list",length = 19)
for(i in c(2:9,11:19)){
  msbb_gse84422.DRsort[[i]]=DRsort(DCGs = msbb_gse84422.DCGs[[i]],DCLs = msbb_gse84422.DCLs[[i]],tf2target = regnet_tf2target.HGNC,expGenes =rownames(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs$Frontal_Pole))  
}


for(i in c(1,10)){
  msbb_gse84422.DRsort[[i]]=DRsort(DCGs = msbb_gse84422.DCGs[[i]],DCLs = msbb_gse84422.DCLs[[i]],tf2target = regnet_tf2target.HGNC,expGenes =rownames(msbb_gse84422_GPL570_samplesToAnalyse.exprs$Amygdala))  
}

msbb_gse84422_early_generic_DCG_jaccard.list=mapply(FUN = function(a,b)jaccard(a,b),msbb_gse84422_earlyCDR.DCGs_list,msbb_gse84422.DCGs_list)
msbb_gse84422_early_generic_RS_DCG_jaccard.list=mapply(FUN = function(a,b)jaccard(a,b),msbb_gse84422_earlyCDR.RS_DCGs_list,msbb_gse84422.RS_DCGs_list)

msbb_gse84422_earlyCDR_DCGs.random_list=msbb_gse84422_earlyCDR_RS_DCGs.random_list=msbb_gse84422_DCGs.random_list=msbb_gse84422_RS_DCGs.random_list=vector(mode = "list",length = 19)
names(msbb_gse84422_earlyCDR_DCGs.random_list)=names(msbb_gse84422_DCGs.random_list)=names(msbb_gse84422_earlyCDR_RS_DCGs.random_list)=names(msbb_gse84422_RS_DCGs.random_list)=names(msbb_gse84422.DCGs_list)

for(t in c(2:9,11:19)){
  msbb_gse84422_DCGs.random_list[[t]]=replicate(n = 10000,expr = sample(x = rownames(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs$Frontal_Pole),size = length(msbb_gse84422.DCGs_list[[t]]),replace = F),simplify = T)
  msbb_gse84422_RS_DCGs.random_list[[t]]=replicate(n = 10000,expr = sample(x = rownames(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs$Frontal_Pole),size = length(msbb_gse84422.RS_DCGs_list[[t]]),replace = F),simplify = T)
}
for(t in c(1,10)){
  msbb_gse84422_DCGs.random_list[[t]]=replicate(n = 10000,expr = sample(x = rownames(msbb_gse84422_GPL570_samplesToAnalyse.exprs$Amygdala),size = length(msbb_gse84422.DCGs_list[[t]]),replace = F),simplify = T)
  msbb_gse84422_RS_DCGs.random_list[[t]]=replicate(n = 10000,expr = sample(x = rownames(msbb_gse84422_GPL570_samplesToAnalyse.exprs$Amygdala),size = length(msbb_gse84422.RS_DCGs_list[[t]]),replace = F),simplify = T)
}

for(t in c(2:9,11:19)){
  msbb_gse84422_earlyCDR_DCGs.random_list[[t]]=replicate(n = 10000,expr = sample(x = rownames(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs$Frontal_Pole),size = length(msbb_gse84422_earlyCDR.DCGs_list[[t]]),replace = F),simplify = T)
  msbb_gse84422_earlyCDR_RS_DCGs.random_list[[t]]=replicate(n = 10000,expr = sample(x = rownames(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs$Frontal_Pole),size = length(msbb_gse84422_earlyCDR.RS_DCGs_list[[t]]),replace = F),simplify = T)
}
for(t in c(1,10)){
  msbb_gse84422_earlyCDR_DCGs.random_list[[t]]=replicate(n = 10000,expr = sample(x = rownames(msbb_gse84422_GPL570_samplesToAnalyse.exprs$Amygdala),size = length(msbb_gse84422_earlyCDR.DCGs_list[[t]]),replace = F),simplify = T)
  msbb_gse84422_earlyCDR_RS_DCGs.random_list[[t]]=replicate(n = 10000,expr = sample(x = rownames(msbb_gse84422_GPL570_samplesToAnalyse.exprs$Amygdala),size = length(msbb_gse84422_earlyCDR.RS_DCGs_list[[t]]),replace = F),simplify = T)
}

msbb_gse84422_early_generic_DCG_jaccard.random_list=msbb_gse84422_early_generic_RS_DCG_jaccard.random_list=list()

for(i in 1:10000){
  
  msbb_gse84422_early_generic_DCG_jaccard.random_list[[i]]=data.frame(random_jcd=round(mapply(FUN = function(a,b)jaccard(a,b),lapply(msbb_gse84422_earlyCDR_DCGs.random_list,function(x)x[,i]),lapply(msbb_gse84422_DCGs.random_list,function(x)x[,i])),digits = 3),stringsAsFactors = F)
  msbb_gse84422_early_generic_RS_DCG_jaccard.random_list[[i]]=data.frame(random_jcd=round(mapply(FUN = function(a,b)jaccard(a,b),lapply(msbb_gse84422_earlyCDR_RS_DCGs.random_list,function(x)x[,i]),lapply(msbb_gse84422_RS_DCGs.random_list,function(x)x[,i])),digits = 3),stringsAsFactors = F)
  
}
msbb_gse84422_early_generic_DCG_jaccard.random_df=cbind.data.frame(msbb_gse84422_early_generic_DCG_jaccard.random_list)
msbb_gse84422_early_generic_RS_DCG_jaccard.random_df=cbind.data.frame(msbb_gse84422_early_generic_RS_DCG_jaccard.random_list)

msbb_gse84422_early_generic_DCG_jaccard.pval=c(1,4,1,323,1,2,1,1,1,1,4,2,157,1,4,5,1,772,1)
msbb_gse84422_early_generic_RS_DCG_jaccard.pval=c(2,486,2195,8819,171,1798,125,3,5693,1,4538,8098,9839,319,4,275,9254,9792,8775)
msbb_gse84422_early_generic_DCG_jaccard.pval=msbb_gse84422_early_generic_DCG_jaccard.pval/10000
msbb_gse84422_early_generic_RS_DCG_jaccard.pval=msbb_gse84422_early_generic_RS_DCG_jaccard.pval/10000

msbb_gse84422_early_generic_DCG_jaccard.adjust_pval=p.adjust(p = msbb_gse84422_early_generic_DCG_jaccard.pval,method = 'fdr')
msbb_gse84422_early_generic_RS_DCG_jaccard.adjust_pval=p.adjust(p = msbb_gse84422_early_generic_RS_DCG_jaccard.pval,method = 'fdr')
names(msbb_gse84422_early_generic_DCG_jaccard.pval)=names(msbb_gse84422_early_generic_RS_DCG_jaccard.pval)=names(msbb_gse84422_early_generic_DCG_jaccard.adjust_pval)=names(msbb_gse84422_early_generic_RS_DCG_jaccard.adjust_pval)=rownames(msbb_gse84422_early_generic_DCG_jaccard.random_df)

msbb_gse84422_All_DCG_Summary.df=cbind.data.frame(list(Brain_Region=names(msbb_gse84422_earlyCDR_total_DCGs)[-c(14,16)],
                                                       EarlyCDR_DCGs=unlist(msbb_gse84422_earlyCDR_total_DCGs[-c(14,16)]),Gen_DCGs=unlist(msbb_gse84422_total_DCGs[-c(14,16)]),
                                                       EarlyCDR_DCGs_Gen_DCGs_Jaccard=unlist(msbb_gse84422_early_generic_DCG_jaccard.list[-c(14,16)]),
                                                       Jaccard_pval=msbb_gse84422_early_generic_DCG_jaccard.pval[-c(14,16)],
                                                       Jaccard_adjP=msbb_gse84422_early_generic_DCG_jaccard.adjust_pval[-c(14,16)]),
                                                       EarlyCDR_DCGs_Gen_DCGs_Jaccard_Rank=rank(-unlist(msbb_gse84422_early_generic_DCG_jaccard.list[-c(14,16)])))

msbb_gse84422_All_RS_DCG_Summary.df=cbind.data.frame(list(Brain_Region=names(msbb_gse84422_earlyCDR_total_DCGs)[-c(14,16)],
                                                       Region_specific_EarlyCDR_DCGs=unlist(lapply(msbb_gse84422_earlyCDR.RS_DCGs_list[-c(14,16)],length)),
                                                       Region_specific_Gen_DCGs=unlist(lapply(msbb_gse84422.RS_DCGs_list[-c(14,16)],length)),
                                                       EarlyCDR_RS_DCGs_RS_Gen_DCGs_Jaccard=unlist(msbb_gse84422_early_generic_RS_DCG_jaccard.list[-c(14,16)]),
                                                       Jaccard_pval=msbb_gse84422_early_generic_RS_DCG_jaccard.pval[-c(14,16)],
                                                       Jaccard_adjP=msbb_gse84422_early_generic_RS_DCG_jaccard.adjust_pval[-c(14,16)]),
                                                       EarlyCDR_RS_DCGs_RS_Gen_DCGs_Jaccard_Rank=rank(-unlist(msbb_gse84422_early_generic_RS_DCG_jaccard.list[-c(14,16)])))

fwrite(arrange(msbb_gse84422_All_DCG_Summary.df,EarlyCDR_DCGs_Gen_DCGs_Jaccard_Rank),"MSBB_All_DCG_Summary.txt",sep="\t",col.names = T,row.names = F)
fwrite(arrange(msbb_gse84422_All_RS_DCG_Summary.df,EarlyCDR_RS_DCGs_RS_Gen_DCGs_Jaccard_Rank),"MSBB_All_RS_DCG_Summary.txt",sep="\t",col.names = T,row.names = F)


#Compare with generic DCGs
msbb_early_gen_DCG_overlap.list=list()
msbb_gse84422.DCGs_list2=msbb_gse84422.DCGs_list[c(2:9,11:19)]
for(i in 1:17){
  msbb_early_gen_DCG_overlap.list[[i]]=intersect(msbb_gse84422.DCGs_list2[[i]],msbb_gse84422_earlyCDR.DCGs_list[[i]])
}

msbb_early_gen_DCG_overlap.KEGG=lapply(msbb_early_gen_DCG_overlap.list,function(x)enrichr(genes = x,databases = kegg_dbs)[[1]]%>%dplyr::filter(Adjusted.P.value<=0.1))
msbb_gse84422.random_DCGs_list=msbb_gse84422_earlyCDR.random_DCGs_list=vector(mode = "list",length = length(msbb_gse84422_earlyCDR.DCGs))
names(msbb_early_gen_DCG_overlap.list)=names(msbb_gse84422.random_DCGs_list)=names(msbb_gse84422_earlyCDR.random_DCGs_list)=names(msbb_gse84422.DCGs)[c(2:9,11:19)]




for(i in 1:length(msbb_gse84422_earlyCDR.DCGs)){
  msbb_gse84422.random_DCGs_list[[i]]=replicate(n = 10000,expr = sample(x = rownames(msbb_gse84422_GPL96_97_byRegion.exprs$`Frontal Pole`),size = length(msbb_gse84422.DCGs_list2[[i]]),replace = F))
  msbb_gse84422_earlyCDR.random_DCGs_list[[i]]=replicate(n = 10000,expr = sample(x = rownames(msbb_gse84422_GPL96_97_byRegion.exprs$`Frontal Pole`),size = length(msbb_gse84422_earlyCDR.DCGs_list[[i]]),replace = F))
}

msbb_early_gen_DCG_random_overlap.list=vector(mode = "list",length = length(msbb_gse84422.DCGs_list2))
names(msbb_early_gen_DCG_random_overlap.list)=names(msbb_gse84422.DCGs_list2)
for(i in 1:length(msbb_early_gen_DCG_random_overlap.list)){
  msbb_early_gen_DCG_random_overlap.list[[i]]=list()  
}

for(i in 1:length(msbb_early_gen_DCG_overlap.list)){
  msbb_early_gen_DCG_random_overlap.list[[i]]=foreach(n=1:10000)%do%{
    intersect(msbb_gse84422.random_DCGs_list[[i]][,n],msbb_gse84422_earlyCDR.random_DCGs_list[[i]][,n])
  }
}

names(msbb_early_gen_DCG_overlap.list)=names(msbb_gse84422_earlyCDR.DCGs_list)=names(msbb_gse84422.DCGs_list2)
msbb_early_gen_DCG_overlap.matrix=matrix(NA,nrow=10000,ncol=17)
for(c in 1:17){
  msbb_early_gen_DCG_overlap.matrix[,c]=unlist(lapply(msbb_early_gen_DCG_random_overlap.list[[c]],length))
}




rownames(msbb_early_gen_DCG_overlap.matrix)=1:10000
colnames(msbb_early_gen_DCG_overlap.matrix)=names(msbb_early_gen_DCG_overlap.list)
msbb_early_gen_DCG_overlap.df=as.data.frame(msbb_early_gen_DCG_overlap.matrix)

#In silico validation

msbb_gse84422_tanzi_SAGs.overlap=lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,intersect,tanzi_ranked_genes.AD)

msbb_exp_DCGs=lapply(lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,length),function(x)round(length(tanzi_ranked_genes.AD)/21982*x,digits = 1))
msbb_gse84422_tanzi_SAGs.overlap=lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,function(x)length(intersect(x,tanzi_ranked_genes.AD)))
msbb_gse84422_tanzi_SAGs.fisher_results=vector(mode = "list",length = 19)
for(m in 1:length(msbb_gse84422_earlyCDR.RS_DCGs_list)){
  msbb_gse84422_tanzi_SAGs.fisher_results[[m]]=fisher.test(matrix(c(length(msbb_gse84422_tanzi_SAGs.overlap[[m]]),
                                                                    length(tanzi_ranked_genes.AD)-length(msbb_gse84422_tanzi_SAGs.overlap[[m]]),
                                                                    lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,length)[[m]]-length(msbb_gse84422_tanzi_SAGs.overlap[[m]]),
                                                                    21982-(length(tanzi_ranked_genes.AD)-length(msbb_gse84422_tanzi_SAGs.overlap[[m]]))-
                                                                      lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,length)[[m]]-length(msbb_gse84422_tanzi_SAGs.overlap[[m]])),nrow = 2))
  
}
res.SAG_earlyCDR_DCG=foreach(i=1:length(msbb_gse84422_tanzi_SAGs.overlap))%do%{
  
  data.frame(SAGs_in_earlyCDR_DCGs=unlist(lapply(msbb_gse84422_tanzi_SAGs.overlap[i],paste,collapse=",")),
             earlyCDR_DCGs=unlist(lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,length)[[i]]),
             Expected_SAGs_in_earlyCDR_DCGs=unname(unlist(msbb_exp_DCGs[i])),
             Fisher_statistic=unname(unlist(lapply(msbb_gse84422_tanzi_SAGs.fisher_results[i],function(s)s$estimate))),
             pval=unname(unlist(lapply(msbb_gse84422_tanzi_SAGs.fisher_results[i],function(p)p$p.value))),
             stringsAsFactors = F)
}
names(res.SAG_earlyCDR_DCG)=names(msbb_gse84422_tanzi_SAGs.overlap)
res.SAG_earlyCDR_DCG_df=data.frame(rbindlist(res.SAG_earlyCDR_DCG),stringsAsFactors = F)
rownames(res.SAG_earlyCDR_DCG_df)=names(msbb_gse84422_tanzi_SAGs.overlap)
res.SAG_earlyCDR_DCG_df$adj.p=p.adjust(p = res.SAG_earlyCDR_DCG_df$pval,method = "fdr")
res.SAG_earlyCDR_DCG_df$SAGs_in_earlyCDR_DCGs=as.numeric(res.SAG_earlyCDR_DCG_df$SAGs_in_earlyCDR_DCGs)
res.SAG_earlyCDR_DCG_df=res.SAG_earlyCDR_DCG_df%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(SAGs_in_earlyCDR_DCGs>Expected_SAGs_in_earlyCDR_DCGs,true = "Over",false = "Under"))  
fwrite(res.SAG_earlyCDR_DCG_df,"MSBB_SAG_earlyCDR_DCG_EnrichmentSummary.txt",sep = "\t",col.names = T,row.names = F,quote = F)

dhmc_bennet.NP=read.table("/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/GW-DNA-hydroxymethylation-Bennet/1-s2.0-S155252601633059X-mmc4_NP.txt",sep = "\t",header = T,as.is = T)
dhmc_bennet.NFT=read.table("/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/GW-DNA-hydroxymethylation-Bennet/1-s2.0-S155252601633059X-mmc5_NFT.txt",sep = "\t",header = T,as.is = T)
d5mc_bernstein.Dhml=scan("/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/5hMC_Tau_AD/ddw109_Supp/Supplemental Table 6_Bernstein et al_R1.txt",sep = "\n",what = "char")

msbb_gse84422_earlyCDR_DCG_d5mc.overlap=lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,function(x)length(intersect(x,d5mc_bernstein.Dhml)))
dhmc_bennet_NFT.genes=dhmc_bennet.NFT%>%dplyr::filter(q.value_<=0.1)%>%pull(Nearest.gene)
dhmc_bennet_NP.genes=dhmc_bennet.NP%>%dplyr::filter(q.vlaue_<=0.1)%>%pull(Neartest.gene)
msbb_gse84422_earlyCDR_DCG_dhmc_NFT.overlap=lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,function(x)length(intersect(x,dhmc_bennet_NFT.genes)))
msbb_gse84422_earlyCDR_DCG_dhmc_NP.overlap=lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,function(x)length(intersect(x,dhmc_bennet_NP.genes)))
msbb_exp_dhmc_NP_DCGs=lapply(lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,length),function(x)round(length(dhmc_bennet_NP.genes)/21982*x,digits = 1))

msbb_gse84422_dhmc_NP.fisher_results=msbb_gse84422_dhmc_NFT.fisher_results=msbb_gse84422_d5mc.fisher_results=vector(mode = "list",length = length(msbb_gse84422_earlyCDR.RS_DCGs_list))
names(msbb_gse84422_dhmc_NP.fisher_results)=names(msbb_gse84422_dhmc_NFT.fisher_results)=names(msbb_gse84422_d5mc.fisher_results)=names(msbb_gse84422_earlyCDR.RS_DCGs_list)
for(m in 1:length(msbb_gse84422_earlyCDR.RS_DCGs_list)){
  msbb_gse84422_dhmc_NP.fisher_results[[m]]=fisher.test(matrix(c(length(msbb_gse84422_earlyCDR_DCG_dhmc_NP.overlap[[m]]),
                                                                 length(dhmc_bennet_NP.genes)-length(msbb_gse84422_earlyCDR_DCG_dhmc_NP.overlap[[m]]),
                                                                 lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,length)[[m]]-length(msbb_gse84422_earlyCDR_DCG_dhmc_NP.overlap[[m]]),
                                                                 21982-(length(dhmc_bennet_NP.genes)-length(msbb_gse84422_earlyCDR_DCG_dhmc_NP.overlap[[m]]))-
                                                                   lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,length)[[m]]-length(msbb_gse84422_earlyCDR_DCG_dhmc_NP.overlap[[m]])),nrow = 2))
}

res.dhmc_NP_DCG=foreach(i=1:length(msbb_gse84422_tanzi_SAGs.overlap))%do%{
  
  
  data.frame(NPgenes_in_DCGs=unlist(lapply(msbb_gse84422_earlyCDR_DCG_dhmc_NP.overlap[i],paste,collapse=",")),
             DCGs=unlist(lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,length)[[i]]),
             Expected_dhmc_NPgenes_in_DCGs=unname(unlist(msbb_exp_DCGs[i])),
             Fisher_statistic=unname(unlist(lapply(msbb_gse84422_dhmc_NP.fisher_results[i],function(s)s$estimate))),
             pval=unname(unlist(lapply(msbb_gse84422_dhmc_NP.fisher_results[i],function(p)p$p.value))),
             stringsAsFactors = F)
}


names(res.dhmc_NP_DCG)=names(msbb_gse84422_earlyCDR_DCG_dhmc_NP.overlap)
res.dhmc_NP_DCG_df=data.frame(rbindlist(res.dhmc_NP_DCG),stringsAsFactors = F)
rownames(res.dhmc_NP_DCG_df)=names(msbb_gse84422_earlyCDR_DCG_dhmc_NP.overlap)
res.dhmc_NP_DCG_df$adj.p=p.adjust(p = res.dhmc_NP_DCG_df$pval,method = "fdr")
res.dhmc_NP_DCG_df$NPgenes_in_DCGs=as.numeric(res.dhmc_NP_DCG_df$NPgenes_in_DCGs)
res.dhmc_NP_DCG_df=res.dhmc_NP_DCG_df%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(NPgenes_in_DCGs>Expected_dhmc_NPgenes_in_DCGs,true = "Over",false = "Under"))
fwrite(res.dhmc_NP_DCG_df,"MSBB_DhMR_NP_earlyCDR_DCG_EnrichmentSummary.txt",sep = "\t",col.names = T,row.names = F,quote = F)


msbb_gse84422_earlyCDR_DCG_dhmc_NFT.overlap=lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,function(x)length(intersect(x,dhmc_bennet_NFT.genes)))
for(m in 1:length(msbb_gse84422_earlyCDR.RS_DCGs_list)){
  msbb_gse84422_dhmc_NFT.fisher_results[[m]]=fisher.test(matrix(c(length(msbb_gse84422_earlyCDR_DCG_dhmc_NFT.overlap[[m]]),
                                                                  length(dhmc_bennet_NFT.genes)-length(msbb_gse84422_earlyCDR_DCG_dhmc_NFT.overlap[[m]]),
                                                                  lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,length)[[m]]-length(msbb_gse84422_earlyCDR_DCG_dhmc_NFT.overlap[[m]]),
                                                                  21982-(length(dhmc_bennet_NFT.genes)-length(msbb_gse84422_earlyCDR_DCG_dhmc_NFT.overlap[[m]]))-
                                                                    lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,length)[[m]]-length(msbb_gse84422_earlyCDR_DCG_dhmc_NFT.overlap[[m]])),nrow = 2))
}

res.dhmc_NFT_DCG=foreach(i=1:length(msbb_gse84422_earlyCDR_DCG_dhmc_NFT.overlap))%do%{
  
  
  data.frame(NFTgenes_in_DCGs=unlist(lapply(msbb_gse84422_earlyCDR_DCG_dhmc_NFT.overlap[i],paste,collapse=",")),
             DCGs=unlist(lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,length)[[i]]),
             Expected_dhmc_NFTgenes_in_DCGs=unname(unlist(msbb_exp_DCGs[i])),
             Fisher_statistic=unname(unlist(lapply(msbb_gse84422_dhmc_NFT.fisher_results[i],function(s)s$estimate))),
             pval=unname(unlist(lapply(msbb_gse84422_dhmc_NFT.fisher_results[i],function(p)p$p.value))),
             stringsAsFactors = F)
}


names(res.dhmc_NFT_DCG)=names(msbb_gse84422_earlyCDR_DCG_dhmc_NFT.overlap)
res.dhmc_NFT_DCG_df=data.frame(rbindlist(res.dhmc_NFT_DCG),stringsAsFactors = F)
rownames(res.dhmc_NFT_DCG_df)=names(msbb_gse84422_earlyCDR_DCG_dhmc_NFT.overlap)
res.dhmc_NFT_DCG_df$adj.p=p.adjust(p = res.dhmc_NFT_DCG_df$pval,method = "fdr")
res.dhmc_NFT_DCG_df$NFTgenes_in_DCGs=as.numeric(res.dhmc_NFT_DCG_df$NFTgenes_in_DCGs)
res.dhmc_NFT_DCG_df=res.dhmc_NFT_DCG_df%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(NFTgenes_in_DCGs>Expected_dhmc_NFTgenes_in_DCGs,true = "Over",false = "Under"))
fwrite(res.dhmc_NFT_DCG_df,"MSBB_DhMR_NFT_earlyCDR_DCG_EnrichmentSummary.txt",sep = "\t",col.names = T,row.names = F,quote = F)



zhang_celltype_ADgenes=read.xls('/Users/sandeepamberkar/Work/Data/BrainExpression_Datasets/Zhang_19BrainRegions_Paper/Zhang_BrainCelltype_Markers.xlsx',skip=1,sheet=3,header=T,as.is=T)
zhang_celltype_ADgenes.list=vector(mode = 'list',length = 5)
names(zhang_celltype_ADgenes.list)=sort(unique(zhang_celltype_ADgenes$Cell.type))
zhang_celltype_ADgenes.list$Astrocytes=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_ADgenes.list)[1],zhang_celltype_ADgenes$Cell.type)]
zhang_celltype_ADgenes.list$Endothelial=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_ADgenes.list)[2],zhang_celltype_ADgenes$Cell.type)]
zhang_celltype_ADgenes.list$Microglia=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_ADgenes.list)[3],zhang_celltype_ADgenes$Cell.type)]
zhang_celltype_ADgenes.list$Neurons=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_ADgenes.list)[4],zhang_celltype_ADgenes$Cell.type)]
zhang_celltype_ADgenes.list$Oligodendrocytes=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_ADgenes.list)[5],zhang_celltype_ADgenes$Cell.type)]

msbb_gse84422_earlyCDR_DCG_Astrocytes.overlap=lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,intersect,zhang_celltype_ADgenes.list$Astrocytes)
msbb_gse84422_earlyCDR_DCG_Endothelial.overlap=lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,intersect,zhang_celltype_ADgenes.list$Endothelial)
msbb_gse84422_earlyCDR_DCG_Microglia.overlap=lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,intersect,zhang_celltype_ADgenes.list$Microglia)
msbb_gse84422_earlyCDR_DCG_Neurons.overlap=lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,intersect,zhang_celltype_ADgenes.list$Neurons)
msbb_gse84422_earlyCDR_DCG_Oligodendrocytes.overlap=lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,intersect,zhang_celltype_ADgenes.list$Oligodendrocytes)

msbb_gse84422.exp_Astrocytes=lapply(lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,length),function(x)round(length(zhang_celltype_ADgenes.list$Astrocytes)/21982*x,digits = 1))
msbb_gse84422.exp_Endothelial=lapply(lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,length),function(x)round(length(zhang_celltype_ADgenes.list$Endothelial)/21982*x,digits = 1))
msbb_gse84422.exp_Microglia=lapply(lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,length),function(x)round(length(zhang_celltype_ADgenes.list$Microglia)/21982*x,digits = 1))
msbb_gse84422.exp_Neurons=lapply(lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,length),function(x)round(length(zhang_celltype_ADgenes.list$Neurons)/21982*x,digits = 1))
msbb_gse84422.exp_Oligodendrocytes=lapply(lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,length),function(x)round(length(zhang_celltype_ADgenes.list$Oligodendrocytes)/21982*x,digits = 1))

msbb_gse84422_Astrocytes.fisher_results=msbb_gse84422_Endothelial.fisher_results=msbb_gse84422_Microglia.fisher_results=msbb_gse84422_Neurons.fisher_results=msbb_gse84422_Oligodendrocytes.fisher_results=vector(mode = "list",length = length(msbb_gse84422_earlyCDR.RS_DCGs_list))
names(msbb_gse84422_Astrocytes.fisher_results)=names(msbb_gse84422_Endothelial.fisher_results)=names(msbb_gse84422_Microglia.fisher_results)=names(msbb_gse84422_Neurons.fisher_results)=names(msbb_gse84422_Oligodendrocytes.fisher_results)=names(msbb_gse84422_total_DCGs)

for(m in 1:length(msbb_gse84422_Astrocytes.fisher_results)){
  msbb_gse84422_Astrocytes.fisher_results[[m]]=fisher.test(matrix(c(length(msbb_gse84422.exp_Astrocytes[[m]]),
                                                                    length(zhang_celltype_ADgenes.list$Astrocytes)-length(msbb_gse84422_earlyCDR_DCG_Astrocytes.overlap[[m]]),
                                                                    lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,length)[[m]]-length(msbb_gse84422_earlyCDR_DCG_Astrocytes.overlap[[m]]),
                                                                    21982-(length(zhang_celltype_ADgenes.list$Astrocytes)-length(msbb_gse84422_earlyCDR_DCG_Astrocytes.overlap[[m]]))-
                                                                      lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,length)[[m]]-length(msbb_gse84422_earlyCDR_DCG_Astrocytes.overlap[[m]])),nrow = 2))
}
for(m in 1:length(msbb_gse84422_Endothelial.fisher_results)){
  msbb_gse84422_Endothelial.fisher_results[[m]]=fisher.test(matrix(c(length(msbb_gse84422.exp_Endothelial[[m]]),
                                                                     length(zhang_celltype_ADgenes.list$Endothelial)-length(msbb_gse84422_earlyCDR_DCG_Endothelial.overlap[[m]]),
                                                                     lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,length)[[m]]-length(msbb_gse84422_earlyCDR_DCG_Endothelial.overlap[[m]]),
                                                                     21982-(length(zhang_celltype_ADgenes.list$Endothelial)-length(msbb_gse84422_earlyCDR_DCG_Endothelial.overlap[[m]]))-
                                                                       lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,length)[[m]]-length(msbb_gse84422_earlyCDR_DCG_Endothelial.overlap[[m]])),nrow = 2))
}
for(m in 1:length(msbb_gse84422_Microglia.fisher_results)){
  msbb_gse84422_Microglia.fisher_results[[m]]=fisher.test(matrix(c(length(msbb_gse84422.exp_Microglia[[m]]),
                                                                   length(zhang_celltype_ADgenes.list$Microglia)-length(msbb_gse84422_earlyCDR_DCG_Microglia.overlap[[m]]),
                                                                   lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,length)[[m]]-length(msbb_gse84422_earlyCDR_DCG_Microglia.overlap[[m]]),
                                                                   21982-(length(zhang_celltype_ADgenes.list$Microglia)-length(msbb_gse84422_earlyCDR_DCG_Microglia.overlap[[m]]))-
                                                                     lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,length)[[m]]-length(msbb_gse84422_earlyCDR_DCG_Microglia.overlap[[m]])),nrow = 2))
}
for(m in 1:length(msbb_gse84422_Neurons.fisher_results)){
  msbb_gse84422_Neurons.fisher_results[[m]]=fisher.test(matrix(c(length(msbb_gse84422.exp_Neurons[[m]]),
                                                                 length(zhang_celltype_ADgenes.list$Neurons)-length(msbb_gse84422_earlyCDR_DCG_Neurons.overlap[[m]]),
                                                                 lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,length)[[m]]-length(msbb_gse84422_earlyCDR_DCG_Neurons.overlap[[m]]),
                                                                 21982-(length(zhang_celltype_ADgenes.list$Neurons)-length(msbb_gse84422_earlyCDR_DCG_Neurons.overlap[[m]]))-
                                                                   lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,length)[[m]]-length(msbb_gse84422_earlyCDR_DCG_Neurons.overlap[[m]])),nrow = 2))
}
for(m in 1:length(msbb_gse84422_Oligodendrocytes.fisher_results)){
  msbb_gse84422_Oligodendrocytes.fisher_results[[m]]=fisher.test(matrix(c(length(msbb_gse84422.exp_Oligodendrocytes[[m]]),
                                                                          length(zhang_celltype_ADgenes.list$Oligodendrocytes)-length(msbb_gse84422_earlyCDR_DCG_Oligodendrocytes.overlap[[m]]),
                                                                          lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,length)[[m]]-length(msbb_gse84422_earlyCDR_DCG_Oligodendrocytes.overlap[[m]]),
                                                                          21982-(length(zhang_celltype_ADgenes.list$Oligodendrocytes)-length(msbb_gse84422_earlyCDR_DCG_Oligodendrocytes.overlap[[m]]))-
                                                                            lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,length)[[m]]-length(msbb_gse84422_earlyCDR_DCG_Oligodendrocytes.overlap[[m]])),nrow = 2))
}

res.earlyCDR_DCG_Astrocytes=foreach(i=1:length(msbb_gse84422_earlyCDR_DCG_Astrocytes.overlap))%do%{
  
  
  data.frame(AstrocyteMarkers_in_earlyCDR_DCGs=unlist(lapply(msbb_gse84422_earlyCDR_DCG_Astrocytes.overlap[i],paste,collapse=",")),
             Nr_AstrocyteMarkers_in_earlyCDR_DCGs=unlist(lapply(lapply(msbb_gse84422_earlyCDR_DCG_Astrocytes.overlap,length)[i],paste,collapse=",")),
             DCGs=unlist(lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,length)[[i]]),
             Expected_AstrocyteMarkers_in_earlyCDR_DCGs=unname(unlist(msbb_gse84422.exp_Astrocytes[i])),
             Fisher_statistic=unname(unlist(lapply(msbb_gse84422_Astrocytes.fisher_results[i],function(s)s$estimate))),
             pval=unname(unlist(lapply(msbb_gse84422_Astrocytes.fisher_results[i],function(p)p$p.value))),
             stringsAsFactors = F)
}
names(res.earlyCDR_DCG_Astrocytes)=names(msbb_gse84422_earlyCDR_DCG_Astrocytes.overlap)
res.earlyCDR_DCG_Astrocytes_df=data.frame(rbindlist(res.earlyCDR_DCG_Astrocytes),stringsAsFactors = F)
rownames(res.earlyCDR_DCG_Astrocytes_df)=names(msbb_gse84422_earlyCDR_DCG_Astrocytes.overlap)
res.earlyCDR_DCG_Astrocytes_df$adj.p=p.adjust(p = res.earlyCDR_DCG_Astrocytes_df$pval,method = "fdr")
res.earlyCDR_DCG_Astrocytes_df$Expected_AstrocyteMarkers_in_earlyCDR_DCGs=as.numeric(res.earlyCDR_DCG_Astrocytes_df$Expected_AstrocyteMarkers_in_earlyCDR_DCGs)
res.earlyCDR_DCG_Astrocytes_df$Nr_AstrocyteMarkers_in_earlyCDR_DCGs=as.numeric(res.earlyCDR_DCG_Astrocytes_df$Nr_AstrocyteMarkers_in_earlyCDR_DCGs)
res.earlyCDR_DCG_Astrocytes_df=res.earlyCDR_DCG_Astrocytes_df%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(Nr_AstrocyteMarkers_in_earlyCDR_DCGs>Expected_AstrocyteMarkers_in_earlyCDR_DCGs,true = "Over",false = "Under"))


res.earlyCDR_DCG_Endothelial=foreach(i=1:length(msbb_gse84422_earlyCDR_DCG_Endothelial.overlap))%do%{
  
  
  data.frame(EndothelialMarkers_in_earlyCDR_DCGs=unlist(lapply(msbb_gse84422_earlyCDR_DCG_Endothelial.overlap[i],paste,collapse=",")),
             Nr_EndothelialMarkers_in_earlyCDR_DCGs=unlist(lapply(lapply(msbb_gse84422_earlyCDR_DCG_Endothelial.overlap,length)[i],paste,collapse=",")),
             DCGs=unlist(lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,length)[[i]]),
             Expected_EndothelialMarkers_in_earlyCDR_DCGs=unname(unlist(msbb_gse84422.exp_Endothelial[i])),
             Fisher_statistic=unname(unlist(lapply(msbb_gse84422_Endothelial.fisher_results[i],function(s)s$estimate))),
             pval=unname(unlist(lapply(msbb_gse84422_Endothelial.fisher_results[i],function(p)p$p.value))),
             stringsAsFactors = F)
}
names(res.earlyCDR_DCG_Endothelial)=names(msbb_gse84422_earlyCDR_DCG_Endothelial.overlap)
res.earlyCDR_DCG_Endothelial_df=data.frame(rbindlist(res.earlyCDR_DCG_Endothelial),stringsAsFactors = F)
rownames(res.earlyCDR_DCG_Endothelial_df)=names(msbb_gse84422_earlyCDR_DCG_Endothelial.overlap)
res.earlyCDR_DCG_Endothelial_df$adj.p=p.adjust(p = res.earlyCDR_DCG_Endothelial_df$pval,method = "fdr")
res.earlyCDR_DCG_Endothelial_df$Expected_EndothelialMarkers_in_earlyCDR_DCGs=as.numeric(res.earlyCDR_DCG_Endothelial_df$Expected_EndothelialMarkers_in_earlyCDR_DCGs)
res.earlyCDR_DCG_Endothelial_df$Nr_EndothelialMarkers_in_earlyCDR_DCGs=as.numeric(res.earlyCDR_DCG_Endothelial_df$Nr_EndothelialMarkers_in_earlyCDR_DCGs)
res.earlyCDR_DCG_Endothelial_df=res.earlyCDR_DCG_Endothelial_df%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(Nr_EndothelialMarkers_in_earlyCDR_DCGs>Expected_EndothelialMarkers_in_earlyCDR_DCGs,true = "Over",false = "Under"))

res.earlyCDR_DCG_Microglia=foreach(i=1:length(msbb_gse84422_earlyCDR_DCG_Microglia.overlap))%do%{
  
  
  data.frame(MicrogliaMarkers_in_earlyCDR_DCGs=unlist(lapply(msbb_gse84422_earlyCDR_DCG_Microglia.overlap[i],paste,collapse=",")),
             Nr_MicrogliaMarkers_in_earlyCDR_DCGs=unlist(lapply(lapply(msbb_gse84422_earlyCDR_DCG_Microglia.overlap,length)[i],paste,collapse=",")),
             DCGs=unlist(lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,length)[[i]]),
             Expected_MicrogliaMarkers_in_earlyCDR_DCGs=unname(unlist(msbb_gse84422.exp_Microglia[i])),
             Fisher_statistic=unname(unlist(lapply(msbb_gse84422_Microglia.fisher_results[i],function(s)s$estimate))),
             pval=unname(unlist(lapply(msbb_gse84422_Microglia.fisher_results[i],function(p)p$p.value))),
             stringsAsFactors = F)
}
names(res.earlyCDR_DCG_Microglia)=names(msbb_gse84422_earlyCDR_DCG_Microglia.overlap)
res.earlyCDR_DCG_Microglia_df=data.frame(rbindlist(res.earlyCDR_DCG_Microglia),stringsAsFactors = F)
rownames(res.earlyCDR_DCG_Microglia_df)=names(msbb_gse84422_earlyCDR_DCG_Microglia.overlap)
res.earlyCDR_DCG_Microglia_df$adj.p=p.adjust(p = res.earlyCDR_DCG_Microglia_df$pval,method = "fdr")
res.earlyCDR_DCG_Microglia_df$Expected_MicrogliaMarkers_in_earlyCDR_DCGs=as.numeric(res.earlyCDR_DCG_Microglia_df$Expected_MicrogliaMarkers_in_earlyCDR_DCGs)
res.earlyCDR_DCG_Microglia_df$Nr_MicrogliaMarkers_in_earlyCDR_DCGs=as.numeric(res.earlyCDR_DCG_Microglia_df$Nr_MicrogliaMarkers_in_earlyCDR_DCGs)
res.earlyCDR_DCG_Microglia_df=res.earlyCDR_DCG_Microglia_df%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(Nr_MicrogliaMarkers_in_earlyCDR_DCGs>Expected_MicrogliaMarkers_in_earlyCDR_DCGs,true = "Over",false = "Under"))

res.earlyCDR_DCG_Neurons=foreach(i=1:length(msbb_gse84422_earlyCDR_DCG_Neurons.overlap))%do%{
  
  
  data.frame(NeuronsMarkers_in_earlyCDR_DCGs=unlist(lapply(msbb_gse84422_earlyCDR_DCG_Neurons.overlap[i],paste,collapse=",")),
             Nr_NeuronsMarkers_in_earlyCDR_DCGs=unlist(lapply(lapply(msbb_gse84422_earlyCDR_DCG_Neurons.overlap,length)[i],paste,collapse=",")),
             DCGs=unlist(lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,length)[[i]]),
             Expected_NeuronsMarkers_in_earlyCDR_DCGs=unname(unlist(msbb_gse84422.exp_Neurons[i])),
             Fisher_statistic=unname(unlist(lapply(msbb_gse84422_Neurons.fisher_results[i],function(s)s$estimate))),
             pval=unname(unlist(lapply(msbb_gse84422_Neurons.fisher_results[i],function(p)p$p.value))),
             stringsAsFactors = F)
}
names(res.earlyCDR_DCG_Neurons)=names(msbb_gse84422_earlyCDR_DCG_Neurons.overlap)
res.earlyCDR_DCG_Neurons_df=data.frame(rbindlist(res.earlyCDR_DCG_Neurons),stringsAsFactors = F)
rownames(res.earlyCDR_DCG_Neurons_df)=names(msbb_gse84422_earlyCDR_DCG_Neurons.overlap)
res.earlyCDR_DCG_Neurons_df$adj.p=p.adjust(p = res.earlyCDR_DCG_Neurons_df$pval,method = "fdr")
res.earlyCDR_DCG_Neurons_df$Expected_NeuronsMarkers_in_earlyCDR_DCGs=as.numeric(res.earlyCDR_DCG_Neurons_df$Expected_NeuronsMarkers_in_earlyCDR_DCGs)
res.earlyCDR_DCG_Neurons_df$Nr_NeuronsMarkers_in_earlyCDR_DCGs=as.numeric(res.earlyCDR_DCG_Neurons_df$Nr_NeuronsMarkers_in_earlyCDR_DCGs)
res.earlyCDR_DCG_Neurons_df=res.earlyCDR_DCG_Neurons_df%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(Nr_NeuronsMarkers_in_earlyCDR_DCGs>Expected_NeuronsMarkers_in_earlyCDR_DCGs,true = "Over",false = "Under"))

res.earlyCDR_DCG_Oligodendrocytes=foreach(i=1:length(msbb_gse84422_earlyCDR_DCG_Oligodendrocytes.overlap))%do%{
  
  
  data.frame(OligodendrocytesMarkers_in_earlyCDR_DCGs=unlist(lapply(msbb_gse84422_earlyCDR_DCG_Oligodendrocytes.overlap[i],paste,collapse=",")),
             Nr_OligodendrocytesMarkers_in_earlyCDR_DCGs=unlist(lapply(lapply(msbb_gse84422_earlyCDR_DCG_Oligodendrocytes.overlap,length)[i],paste,collapse=",")),
             DCGs=unlist(lapply(msbb_gse84422_earlyCDR.RS_DCGs_list,length)[[i]]),
             Expected_OligodendrocytesMarkers_in_earlyCDR_DCGs=unname(unlist(msbb_gse84422.exp_Oligodendrocytes[i])),
             Fisher_statistic=unname(unlist(lapply(msbb_gse84422_Oligodendrocytes.fisher_results[i],function(s)s$estimate))),
             pval=unname(unlist(lapply(msbb_gse84422_Oligodendrocytes.fisher_results[i],function(p)p$p.value))),
             stringsAsFactors = F)
}
names(res.earlyCDR_DCG_Oligodendrocytes)=names(msbb_gse84422_earlyCDR_DCG_Oligodendrocytes.overlap)
res.earlyCDR_DCG_Oligodendrocytes_df=data.frame(rbindlist(res.earlyCDR_DCG_Oligodendrocytes),stringsAsFactors = F)
rownames(res.earlyCDR_DCG_Oligodendrocytes_df)=names(msbb_gse84422_earlyCDR_DCG_Oligodendrocytes.overlap)
res.earlyCDR_DCG_Oligodendrocytes_df$adj.p=p.adjust(p = res.earlyCDR_DCG_Oligodendrocytes_df$pval,method = "fdr")
res.earlyCDR_DCG_Oligodendrocytes_df$Expected_OligodendrocytesMarkers_in_earlyCDR_DCGs=as.numeric(res.earlyCDR_DCG_Oligodendrocytes_df$Expected_OligodendrocytesMarkers_in_earlyCDR_DCGs)
res.earlyCDR_DCG_Oligodendrocytes_df$Nr_OligodendrocytesMarkers_in_earlyCDR_DCGs=as.numeric(res.earlyCDR_DCG_Oligodendrocytes_df$Nr_OligodendrocytesMarkers_in_earlyCDR_DCGs)
res.earlyCDR_DCG_Oligodendrocytes_df=res.earlyCDR_DCG_Oligodendrocytes_df%>%rownames_to_column("Brain-region")%>%mutate(Representation=if_else(Nr_OligodendrocytesMarkers_in_earlyCDR_DCGs>Expected_OligodendrocytesMarkers_in_earlyCDR_DCGs,true = "Over",false = "Under"))

#Comparison with generic DCGs
msbb_gse84422_diffcoexp_results_files=list.files(path = ".",pattern = "diffcoexp")%>%grep(pattern = ".RDS",value = T)%>%sort
msbb_gse84422_diffcoexp_results=vector(mode = "list",length = length(msbb_gse84422_diffcoexp_results_files))
names(msbb_gse84422_diffcoexp_results)=gsub(pattern = " ",unlist(lapply(lapply(msbb_gse84422_diffcoexp_results_files,function(y)strsplit(x = y,split = "_")[[1]]),`[[`,1)),replacement = "_")
for(f in 1:19){
  msbb_gse84422_diffcoexp_results[[f]]=readRDS(msbb_gse84422_diffcoexp_results_files[f])  
}

msbb_gse84422.DCGs=lapply(msbb_gse84422_diffcoexp_results,function(x)x$DCGs)
msbb_gse84422.DCGs_list=lapply(msbb_gse84422.DCGs,function(x)x%>%dplyr::filter(q<=0.05)%>%pull(Gene))
msbb_gse84422.RS_DCGs_list=vector(mode = "list",length = 19)
names(msbb_gse84422.RS_DCGs_list)=names(msbb_gse84422.DCGs_list)

for(i in 1:length(msbb_gse84422.DCGs_list)){
  msbb_gse84422.RS_DCGs_list[[i]]=setdiff(msbb_gse84422.DCGs_list[[i]],(Reduce(union,lapply(msbb_gse84422.DCGs_list,function(x)intersect(x,msbb_gse84422.DCGs_list[[i]]))[-i])))
}



msbb_array_earlyCDR_DRrank_TDD.files=msbb_array_earlyCDR_DRrank_TED.files=msbb_array_earlyCDR_DRrank.TDD=msbb_array_earlyCDR_DRrank.TED=vector(mode = "list",length = 19)
names(msbb_array_earlyCDR_DRrank_TDD.files)=names(msbb_array_earlyCDR_DRrank_TED.files)=names(msbb_array_earlyCDR_DRrank.TDD)=names(msbb_array_earlyCDR_DRrank.TED)=names(msbb_gse84422_earlyCDR.DCGs)

msbb_array_earlyCDR_DRrank_TDD.files=list.files(path = "Early_Diffcoexp_CDR/DRrank_EarlyCDR_TDD_Results/",pattern = "TDD",full.names = T)
msbb_array_earlyCDR_DRrank_TED.files=list.files(path = "Early_Diffcoexp_CDR/DRrank_EarlyCDR_TED_Results/",pattern = "TED",full.names = T)
for(t in 1:19){
  msbb_array_earlyCDR_DRrank.TDD[[t]]=readRDS(msbb_array_earlyCDR_DRrank_TDD.files[[t]])
  msbb_array_earlyCDR_DRrank.TED[[t]]=readRDS(msbb_array_earlyCDR_DRrank_TED.files[[t]])
}

msbb_array_genDCG_DRrank.TDD=msbb_array_genDCG_DRrank.TED=vector(mode = "list",length = 19)
names(msbb_array_genDCG_DRrank.TDD)=names(msbb_array_genDCG_DRrank.TED)=names(msbb_gse84422.DCGs_list)
msbb_array_genDCG_DRrank_TED.files=list.files(path = "./Gen_Diffcoexp/DRrank_Gen_DCG_TED/",full.names = T)
msbb_array_genDCG_DRrank_TDD.files=list.files(path = "./Gen_Diffcoexp/DRrank_Gen_DCG_TDD/",full.names = T)
for(i in 1:19){
  msbb_array_genDCG_DRrank.TED[[i]]=readRDS(msbb_array_genDCG_DRrank_TED.files[[i]])
  msbb_array_genDCG_DRrank.TDD[[i]]=readRDS(msbb_array_genDCG_DRrank_TDD.files[[i]])
}
msbb_array_genDCG_DRrank.TED=lapply(msbb_array_genDCG_DRrank.TED,function(x)x%>%dplyr::filter(p<=0.05)%>%mutate(Rank=with_order(order_by = score,fun = row_number,x = -score))%>%arrange(desc(-Rank))%>%dplyr::select(c(Gene,Rank,score,p,p.adj)))
