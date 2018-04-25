library(diffcoexp)
library(org.Hs.eg.db)
library(dplyr)
library(magrittr)
library(enrichR)
library(data.table)
library(DCGL)

#setwd("/Users/sandeepamberkar/Work/Data/MSMM-MSBB-HBTRC-Data/MSBB_Array19/GSE84422/")
setwd("/shared/hidelab2/user/md4zsa/Work/Data/MSBB_Array19/GSE84422")
msbb_gse84422_GPL96_97_samplesToAnalyse.exprs=readRDS("msbb_gse84422_GPL96_97_samplesToAnalyse_exprs.RDS")
names(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs)=gsub(pattern = " ",replacement = "_",x = names(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs))
msbb_gse84422_GPL570_samplesToAnalyse.exprs=readRDS("msbb_gse84422_GPL570_samplesToAnalyse_exprs.RDS")
names(msbb_gse84422_GPL570_samplesToAnalyse.exprs)=gsub(pattern = " ",replacement = "_",x = names(msbb_gse84422_GPL570_samplesToAnalyse.exprs))
regnet_tf2target.HGNC=fread("/shared/hidelab2/user/md4zsa/Work/Data/TF_Databases/RegNetwork_human_regulators2.txt",header = T,sep = "\t",showProgress = T,data.table = F)%>%dplyr::filter(evidence=="Experimental")%>%dplyr::select(c(regulator_symbol,target_symbol))
#regnet_tf2target.HGNC=fread("/Users/sandeepamberkar/Work/Data/TF_Databases/RegNetwork_human_regulators2.txt",header = T,sep = "\t",showProgress = T,data.table = F)%>%filter(evidence=="Experimental")%>%dplyr::select(c(regulator_symbol,target_symbol))

msbb_gse84422_diffcoexp_results_files=list.files(path = ".",pattern = "diffcoexp")%>%grep(pattern = ".RDS",value = T)%>%sort
msbb_gse84422_diffcoexp_results=vector(mode = "list",length = length(msbb_gse84422_diffcoexp_results_files))
names(msbb_gse84422_diffcoexp_results)=gsub(pattern = " ",unlist(lapply(lapply(msbb_gse84422_diffcoexp_results_files,function(y)strsplit(x = y,split = "_")[[1]]),`[[`,1)),replacement = "_")
for(f in 1:19){
  msbb_gse84422_diffcoexp_results[[f]]=readRDS(msbb_gse84422_diffcoexp_results_files[f])  
}
msbb_gse84422.DCGs=lapply(msbb_gse84422_diffcoexp_results,function(x)x$DCGs)
msbb_gse84422.DCGs_list=lapply(msbb_gse84422.DCGs,function(x)x%>%dplyr::filter(q<=0.05)%>%pull(Gene))
msbb_gse84422.DCLs=lapply(msbb_gse84422_diffcoexp_results,function(x)x$DCLs)
msbb_gse84422.DCLs_filtered=lapply(msbb_gse84422_diffcoexp_results,function(x)x$DCLs%>%filter(q.diffcor<=0.05))


msbb_array_DRrank.TDD=msbb_array_DRrank.TED=vector(mode = "list",length = 19)
names(msbb_array_DRrank.TDD)=names(msbb_array_DRrank.TED)=names(msbb_gse84422_diffcoexp_results)
for(i in c(1,10)){
  msbb_array_DRrank.TDD[[i]]=DRrank(DCGs=msbb_gse84422.DCGs[[i]],DCLs=msbb_gse84422.DCLs[[i]],tf2target=regnet_tf2target.HGNC,expGenes=rownames(msbb_gse84422_GPL570_samplesToAnalyse.exprs$Amygdala),rank.method="TDD",Nperm=100)
  saveRDS(msbb_array_DRrank.TDD[[i]],paste(names(msbb_array_DRrank.TDD)[i],"DRrank_TDD.RDS",sep = "_"))
}
for(i in c(2:9,11:19)){
  msbb_array_DRrank.TDD[[i]]=DRrank(DCGs=msbb_gse84422.DCGs[[i]],DCLs=msbb_gse84422.DCLs[[i]],tf2target=regnet_tf2target.HGNC,expGenes=rownames(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs[[1]]),rank.method="TDD",Nperm=100)
  saveRDS(msbb_array_DRrank.TDD[[i]],paste(names(msbb_array_DRrank.TDD)[i],"DRrank_TDD.RDS",sep = "_"))
}
