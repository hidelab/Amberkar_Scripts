library(dplyr)
library(magrittr)
library(data.table)
library(DCGL)

setwd("/shared/hidelab2/user/md4zsa/Work/Data/MSBB_Array19/GSE84422/")
earlyAD_diffcoexp_files=list.files(path = ".",pattern = "earlyAD_diffcoexp",full.names = T)
earlyAD_samples=readRDS("./msbb_gse84422_GPL96_97_earlyAD_samplesToAnalyse.RDS")
earlyAD_samples.exprs=readRDS("./msbb_gse84422_GPL96_97_earlyAD_samplesToAnalyse_exprs.RDS")

regnet_tf2target.HGNC=fread("/shared/hidelab2/user/md4zsa/Work/Data/TF_Databases/RegNetwork_human_regulators2.txt",header = T,sep = "\t",showProgress = T,data.table = F)%>%dplyr::filter(evidence=="Experimental")%>%dplyr::select(c(regulator_symbol,target_symbol))

earlyAD_diffcoexp_results=vector(mode = "list",length = length(earlyAD_diffcoexp_files))
for(i in 1:length(earlyAD_diffcoexp_files)){
  earlyAD_diffcoexp_results[[i]]=readRDS(earlyAD_diffcoexp_files[[i]])
}
names(earlyAD_diffcoexp_results)=names(earlyAD_samples.exprs)

earlyAD_DCGs=lapply(earlyAD_diffcoexp_results,function(x)x$DCGs)
earlyAD_DCLs=lapply(earlyAD_diffcoexp_results,function(x)x$DCLs)

earlyAD_DRrank.TED=earlyAD_DRrank.TED=vector(mode = "list",length = length(earlyAD_samples.exprs))
names(earlyAD_DRrank.TED)=names(earlyAD_DRrank.TED)=names(earlyAD_diffcoexp_results)
for(i in 1:length(earlyAD_samples.exprs)){
  earlyAD_DRrank.TED[[i]]=DRrank(DCGs=earlyAD_DCGs[[i]],DCLs=earlyAD_DCLs[[i]],tf2target=regnet_tf2target.HGNC,expGenes=rownames(earlyAD_samples.exprs$Frontal_Pole),rank.method="TED",Nperm=1000)
  saveRDS(earlyAD_DRrank.TED[[i]],paste(names(earlyAD_DRrank.TED)[i],"earlyAD_DRrank_TED.RDS",sep = "_"))
}
