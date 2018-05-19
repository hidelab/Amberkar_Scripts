library(diffcoexp)
library(dplyr)
library(magrittr)
library(data.table)
library(DCGL)

setwd("/shared/hidelab2/user/md4zsa/Work/Data/MSBB_Array19/GSE84422/")
msbb_gse84422_GPL96_97_samplesToAnalyse.exprs=readRDS("msbb_gse84422_GPL96_97_samplesToAnalyse_exprs.RDS")
names(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs)=gsub(pattern = " ",replacement = "_",x = names(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs))
msbb_gse84422_GPL570_samplesToAnalyse.exprs=readRDS("msbb_gse84422_GPL570_samplesToAnalyse_exprs.RDS")
names(msbb_gse84422_GPL570_samplesToAnalyse.exprs)=gsub(pattern = " ",replacement = "_",x = names(msbb_gse84422_GPL570_samplesToAnalyse.exprs))
regnet_tf2target.HGNC=fread("/shared/hidelab2/user/md4zsa/Work/Data/TF_Databases/RegNetwork_human_regulators2.txt",header = T,sep = "\t",showProgress = T,data.table = F)%>%dplyr::filter(evidence=="Experimental")%>%dplyr::select(c(regulator_symbol,target_symbol))


early_diffcoexp_files=system("find . |grep 'CDR0'|sort",intern = T)
msbb_gse84422_early_diffcoexp_results=msbb_gse84422_earlyCDR.DCGs=msbb_gse84422_earlyCDR.DCLs=msbb_gse84422_earlyCDR.DRGs=msbb_gse84422_earlyCDR.DRLs=vector(mode = "list",length = length(early_diffcoexp_files))



for(i in 1:length(msbb_gse84422_early_diffcoexp_results)){
  msbb_gse84422_early_diffcoexp_results[[i]]=readRDS(early_diffcoexp_files[i])
}
names(msbb_gse84422_early_diffcoexp_results)=names(msbb_gse84422_earlyCDR.DCGs)=names(msbb_gse84422_earlyCDR.DCLs)=names(msbb_gse84422_earlyCDR.DRGs)=names(msbb_gse84422_earlyCDR.DRLs)=gsub(pattern = "./",replacement = "",x = gsub(pattern = "_diffcoexp_CDR0_CDR1.RDS",replacement = "",x = early_diffcoexp_files))
msbb_gse84422_earlyCDR.DCGs=lapply(msbb_gse84422_early_diffcoexp_results,function(x)x$DCGs)
msbb_gse84422_earlyCDR.DCLs=lapply(msbb_gse84422_early_diffcoexp_results,function(x)x$DCLs)

msbb_array_earlyCDR_DRrank.TDD=msbb_array_earlyCDR_DRrank.TED=vector(mode = "list",length = 19)
names(msbb_array_earlyCDR_DRrank.TDD)=names(msbb_array_earlyCDR_DRrank.TED)=names(msbb_gse84422_early_diffcoexp_results)
# for(i in c(1,10)){
#   msbb_array_earlyCDR_DRrank.TDD[[i]]=DRrank(DCGs=msbb_gse84422_earlyCDR.DCGs[[i]],DCLs=msbb_gse84422_earlyCDR.DCLs[[i]],tf2target=regnet_tf2target.HGNC,expGenes=rownames(msbb_gse84422_GPL570_samplesToAnalyse.exprs$Amygdala),rank.method="TDD",Nperm=1000)
#   saveRDS(msbb_array_earlyCDR_DRrank.TDD[[i]],paste(names(msbb_array_earlyCDR_DRrank.TDD)[i],"earlyCDR_DRrank_TDD.RDS",sep = "_"))
# }
for(i in c(4:9,11:19)){
  msbb_array_earlyCDR_DRrank.TDD[[i]]=DRrank(DCGs=msbb_gse84422_earlyCDR.DCGs[[i]],DCLs=msbb_gse84422_earlyCDR.DCLs[[i]],tf2target=regnet_tf2target.HGNC,expGenes=rownames(msbb_gse84422_GPL570_samplesToAnalyse.exprs$Amygdala),rank.method="TDD",Nperm=1000)
  saveRDS(msbb_array_earlyCDR_DRrank.TDD[[i]],paste(names(msbb_array_earlyCDR_DRrank.TDD)[i],"earlyCDR_DRrank_TDD.RDS",sep = "_"))
}