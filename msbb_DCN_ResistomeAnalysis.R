library(igraph)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(gdata)

msbb_diffcorr_fdr10=msbb_resistant_DCN_Nstringent=msbb_resistant_DCN_stringent=vector(mode = 'list',length = 4)
tanzi_select_protectiveGenes=read.xls('../../../../Collaborations/Tanzi_WGS/Protective_exonic_SNPs.xlsx',verbose = T)
names(msbb_diffcorr_fdr10)=names(msbb_resistant_DCN_Nstringent)=names(msbb_resistant_DCN_stringent)=c('FP','IFG','PHG','STG')

msbb_diffcorr_fdr10$FP=msbb_diffcorr_fdr10$FP=fread('FP/FP_allResults_DiffCorr_FDR10.txt',sep = '\t',header = T,data.table = F,showProgress = T)
msbb_diffcorr_fdr10$IFG=msbb_diffcorr_fdr10$IFG=fread('IFG/IFG_allResults_DiffCorr_FDR10.txt',sep = '\t',header = T,data.table = F,showProgress = T)
msbb_diffcorr_fdr10$PHG=msbb_diffcorr_fdr10$PHG=fread('PHG/PHG_allResults_DiffCorr_FDR10.txt',sep = '\t',header = T,data.table = F,showProgress = T)
msbb_diffcorr_fdr10$STG=msbb_diffcorr_fdr10$STG=fread('STG/STG_allResults_DiffCorr_FDR10.txt',sep = '\t',header = T,data.table = F,showProgress = T)

