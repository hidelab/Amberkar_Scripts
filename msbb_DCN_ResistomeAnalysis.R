library(igraph)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(gdata)
library(data.table)

msbb_diffcorr_fdr10=msbb_Resistant_DCN_Nstringent=msbb_Resistant_DCN_stringent=vector(mode = 'list',length = 4)
tanzi_select_protectiveGenes=read.xls('/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Protective_exonic_SNPs.xlsx',verbose = T)
names(msbb_diffcorr_fdr10)=names(msbb_Resistant_DCN_Nstringent)=names(msbb_Resistant_DCN_stringent)=c('FP','IFG','PHG','STG')
msbb_Resistant_DCN_Nstringent$FP=msbb_Resistant_DCN_Nstringent$IFG=msbb_Resistant_DCN_Nstringent$PHG=msbb_Resistant_DCN_Nstringent$STG=vector(mode = 'list',length = 2)
msbb_Resistant_DCN_stringent$FP=msbb_Resistant_DCN_stringent$IFG=msbb_Resistant_DCN_stringent$PHG=msbb_Resistant_DCN_stringent$STG=vector(mode = 'list',length = 2)
names(msbb_Resistant_DCN_Nstringent$FP)=names(msbb_Resistant_DCN_Nstringent$IFG)=names(msbb_Resistant_DCN_Nstringent$PHG)=names(msbb_Resistant_DCN_Nstringent$STG)=c('Normal','AD')
names(msbb_Resistant_DCN_stringent$FP)=names(msbb_Resistant_DCN_stringent$IFG)=names(msbb_Resistant_DCN_stringent$PHG)=names(msbb_Resistant_DCN_stringent$STG)=c('Normal','AD')

msbb_diffcorr_fdr10$FP=msbb_diffcorr_fdr10$FP=fread('/Users/sandeepamberkar/Work/Data/MSMM-MSBB-HBTRC-Data/MSMM_RNAseq/MSMM_RNAseq_FinalRelease2/FP/FP_allResults_DiffCorr_FDR10.txt',sep = '\t',header = T,data.table = F,showProgress = T)
msbb_diffcorr_fdr10$IFG=msbb_diffcorr_fdr10$IFG=fread('/Users/sandeepamberkar/Work/Data/MSMM-MSBB-HBTRC-Data/MSMM_RNAseq/MSMM_RNAseq_FinalRelease2/IFG/IFG_allResults_DiffCorr_FDR10.txt',sep = '\t',header = T,data.table = F,showProgress = T)
msbb_diffcorr_fdr10$PHG=msbb_diffcorr_fdr10$PHG=fread('/Users/sandeepamberkar/Work/Data/MSMM-MSBB-HBTRC-Data/MSMM_RNAseq/MSMM_RNAseq_FinalRelease2/PHG/PHG_allResults_DiffCorr_FDR10.txt',sep = '\t',header = T,data.table = F,showProgress = T)
msbb_diffcorr_fdr10$STG=msbb_diffcorr_fdr10$STG=fread('/Users/sandeepamberkar/Work/Data/MSMM-MSBB-HBTRC-Data/MSMM_RNAseq/MSMM_RNAseq_FinalRelease2/STG/STG_allResults_DiffCorr_FDR10.txt',sep = '\t',header = T,data.table = F,showProgress = T)

msbb_Resistant_DCN_Nstringent$FP$Normal=msbb_Resistant_DCN_Nstringent$FP$AD=graph.data.frame(msbb_diffcorr_fdr10$FP[which((msbb_diffcorr_fdr10$FP$Gene.A%in%tanzi_select_protectiveGenes$Gene)|(msbb_diffcorr_fdr10$FP$Gene.B%in%tanzi_select_protectiveGenes$Gene)),c(1:2)])
msbb_Resistant_DCN_Nstringent$IFG$Normal=msbb_Resistant_DCN_Nstringent$IFG$AD=graph.data.frame(msbb_diffcorr_fdr10$IFG[which((msbb_diffcorr_fdr10$IFG$Gene.A%in%tanzi_select_protectiveGenes$Gene)|(msbb_diffcorr_fdr10$IFG$Gene.B%in%tanzi_select_protectiveGenes$Gene)),c(1:2)])
msbb_Resistant_DCN_Nstringent$PHG$Normal=msbb_Resistant_DCN_Nstringent$PHG$AD=graph.data.frame(msbb_diffcorr_fdr10$PHG[which((msbb_diffcorr_fdr10$PHG$Gene.A%in%tanzi_select_protectiveGenes$Gene)|(msbb_diffcorr_fdr10$PHG$Gene.B%in%tanzi_select_protectiveGenes$Gene)),c(1:2)])
msbb_Resistant_DCN_Nstringent$STG$Normal=msbb_Resistant_DCN_Nstringent$STG$AD=graph.data.frame(msbb_diffcorr_fdr10$STG[which((msbb_diffcorr_fdr10$STG$Gene.A%in%tanzi_select_protectiveGenes$Gene)|(msbb_diffcorr_fdr10$STG$Gene.B%in%tanzi_select_protectiveGenes$Gene)),c(1:2)])

E(msbb_Resistant_DCN_Nstringent$FP$Normal)$weight=msbb_diffcorr_fdr10$FP[which((msbb_diffcorr_fdr10$FP$Gene.A%in%tanzi_select_protectiveGenes$Gene)|(msbb_diffcorr_fdr10$FP$Gene.B%in%tanzi_select_protectiveGenes$Gene)),3]+1
E(msbb_Resistant_DCN_Nstringent$IFG$Normal)$weight=msbb_diffcorr_fdr10$IFG[which((msbb_diffcorr_fdr10$IFG$Gene.A%in%tanzi_select_protectiveGenes$Gene)|(msbb_diffcorr_fdr10$IFG$Gene.B%in%tanzi_select_protectiveGenes$Gene)),3]+1
E(msbb_Resistant_DCN_Nstringent$PHG$Normal)$weight=msbb_diffcorr_fdr10$PHG[which((msbb_diffcorr_fdr10$PHG$Gene.A%in%tanzi_select_protectiveGenes$Gene)|(msbb_diffcorr_fdr10$PHG$Gene.B%in%tanzi_select_protectiveGenes$Gene)),3]+1
E(msbb_Resistant_DCN_Nstringent$STG$Normal)$weight=msbb_diffcorr_fdr10$STG[which((msbb_diffcorr_fdr10$STG$Gene.A%in%tanzi_select_protectiveGenes$Gene)|(msbb_diffcorr_fdr10$STG$Gene.B%in%tanzi_select_protectiveGenes$Gene)),3]+1

E(msbb_Resistant_DCN_Nstringent$FP$AD)$weight=msbb_diffcorr_fdr10$FP[which((msbb_diffcorr_fdr10$FP$Gene.A%in%tanzi_select_protectiveGenes$Gene)|(msbb_diffcorr_fdr10$FP$Gene.B%in%tanzi_select_protectiveGenes$Gene)),6]+1
E(msbb_Resistant_DCN_Nstringent$IFG$AD)$weight=msbb_diffcorr_fdr10$IFG[which((msbb_diffcorr_fdr10$IFG$Gene.A%in%tanzi_select_protectiveGenes$Gene)|(msbb_diffcorr_fdr10$IFG$Gene.B%in%tanzi_select_protectiveGenes$Gene)),6]+1
E(msbb_Resistant_DCN_Nstringent$PHG$AD)$weight=msbb_diffcorr_fdr10$PHG[which((msbb_diffcorr_fdr10$PHG$Gene.A%in%tanzi_select_protectiveGenes$Gene)|(msbb_diffcorr_fdr10$PHG$Gene.B%in%tanzi_select_protectiveGenes$Gene)),6]+1
E(msbb_Resistant_DCN_Nstringent$STG$AD)$weight=msbb_diffcorr_fdr10$STG[which((msbb_diffcorr_fdr10$STG$Gene.A%in%tanzi_select_protectiveGenes$Gene)|(msbb_diffcorr_fdr10$STG$Gene.B%in%tanzi_select_protectiveGenes$Gene)),6]+1

write.graph(graph = msbb_Resistant_DCN_Nstringent$FP$Normal,file = 'msbb_Resist_DCN_NS_FP_Normal.gml',format = 'gml')
write.graph(graph = msbb_Resistant_DCN_Nstringent$IFG$Normal,file = 'msbb_Resist_DCN_NS_IFG_Normal.gml',format = 'gml')
write.graph(graph = msbb_Resistant_DCN_Nstringent$PHG$Normal,file = 'msbb_Resist_DCN_NS_PHG_Normal.gml',format = 'gml')
write.graph(graph = msbb_Resistant_DCN_Nstringent$STG$Normal,file = 'msbb_Resist_DCN_NS_STG_Normal.gml',format = 'gml')

write.graph(graph = msbb_Resistant_DCN_Nstringent$FP$AD,file = 'msbb_Resist_DCN_NS_FP_AD.gml',format = 'gml')
write.graph(graph = msbb_Resistant_DCN_Nstringent$IFG$AD,file = 'msbb_Resist_DCN_NS_IFG_AD.gml',format = 'gml')
write.graph(graph = msbb_Resistant_DCN_Nstringent$PHG$AD,file = 'msbb_Resist_DCN_NS_PHG_AD.gml',format = 'gml')
write.graph(graph = msbb_Resistant_DCN_Nstringent$STG$AD,file = 'msbb_Resist_DCN_NS_STG_AD.gml',format = 'gml')

msbb_Resistant_DCN_stringent$FP$Normal=msbb_Resistant_DCN_stringent$FP$AD=graph.data.frame(msbb_diffcorr_fdr10$FP[which((msbb_diffcorr_fdr10$FP$Gene.A%in%tanzi_select_protectiveGenes$Gene)&(msbb_diffcorr_fdr10$FP$Gene.B%in%tanzi_select_protectiveGenes$Gene)),c(1:2)])
msbb_Resistant_DCN_stringent$IFG$Normal=msbb_Resistant_DCN_stringent$IFG$AD=graph.data.frame(msbb_diffcorr_fdr10$IFG[which((msbb_diffcorr_fdr10$IFG$Gene.A%in%tanzi_select_protectiveGenes$Gene)&(msbb_diffcorr_fdr10$IFG$Gene.B%in%tanzi_select_protectiveGenes$Gene)),c(1:2)])
msbb_Resistant_DCN_stringent$PHG$Normal=msbb_Resistant_DCN_stringent$PHG$AD=graph.data.frame(msbb_diffcorr_fdr10$PHG[which((msbb_diffcorr_fdr10$PHG$Gene.A%in%tanzi_select_protectiveGenes$Gene)&(msbb_diffcorr_fdr10$PHG$Gene.B%in%tanzi_select_protectiveGenes$Gene)),c(1:2)])
msbb_Resistant_DCN_stringent$STG$Normal=msbb_Resistant_DCN_stringent$STG$AD=graph.data.frame(msbb_diffcorr_fdr10$STG[which((msbb_diffcorr_fdr10$STG$Gene.A%in%tanzi_select_protectiveGenes$Gene)&(msbb_diffcorr_fdr10$STG$Gene.B%in%tanzi_select_protectiveGenes$Gene)),c(1:2)])

E(msbb_Resistant_DCN_stringent$FP$Normal)$weight=msbb_diffcorr_fdr10$FP[which((msbb_diffcorr_fdr10$FP$Gene.A%in%tanzi_select_protectiveGenes$Gene)&(msbb_diffcorr_fdr10$FP$Gene.B%in%tanzi_select_protectiveGenes$Gene)),3]+1
E(msbb_Resistant_DCN_stringent$IFG$Normal)$weight=msbb_diffcorr_fdr10$IFG[which((msbb_diffcorr_fdr10$IFG$Gene.A%in%tanzi_select_protectiveGenes$Gene)&(msbb_diffcorr_fdr10$IFG$Gene.B%in%tanzi_select_protectiveGenes$Gene)),3]+1
E(msbb_Resistant_DCN_stringent$PHG$Normal)$weight=msbb_diffcorr_fdr10$PHG[which((msbb_diffcorr_fdr10$PHG$Gene.A%in%tanzi_select_protectiveGenes$Gene)&(msbb_diffcorr_fdr10$PHG$Gene.B%in%tanzi_select_protectiveGenes$Gene)),3]+1
E(msbb_Resistant_DCN_stringent$STG$Normal)$weight=msbb_diffcorr_fdr10$STG[which((msbb_diffcorr_fdr10$STG$Gene.A%in%tanzi_select_protectiveGenes$Gene)&(msbb_diffcorr_fdr10$STG$Gene.B%in%tanzi_select_protectiveGenes$Gene)),3]+1

E(msbb_Resistant_DCN_stringent$FP$AD)$weight=msbb_diffcorr_fdr10$FP[which((msbb_diffcorr_fdr10$FP$Gene.A%in%tanzi_select_protectiveGenes$Gene)&(msbb_diffcorr_fdr10$FP$Gene.B%in%tanzi_select_protectiveGenes$Gene)),6]+1
E(msbb_Resistant_DCN_stringent$IFG$AD)$weight=msbb_diffcorr_fdr10$IFG[which((msbb_diffcorr_fdr10$IFG$Gene.A%in%tanzi_select_protectiveGenes$Gene)&(msbb_diffcorr_fdr10$IFG$Gene.B%in%tanzi_select_protectiveGenes$Gene)),6]+1
E(msbb_Resistant_DCN_stringent$PHG$AD)$weight=msbb_diffcorr_fdr10$PHG[which((msbb_diffcorr_fdr10$PHG$Gene.A%in%tanzi_select_protectiveGenes$Gene)&(msbb_diffcorr_fdr10$PHG$Gene.B%in%tanzi_select_protectiveGenes$Gene)),6]+1
E(msbb_Resistant_DCN_stringent$STG$AD)$weight=msbb_diffcorr_fdr10$STG[which((msbb_diffcorr_fdr10$STG$Gene.A%in%tanzi_select_protectiveGenes$Gene)&(msbb_diffcorr_fdr10$STG$Gene.B%in%tanzi_select_protectiveGenes$Gene)),6]+1

write.graph(graph = msbb_Resistant_DCN_stringent$FP$Normal,file = 'msbb_Resist_DCN_S_FP_Normal.gml',format = 'gml')
write.graph(graph = msbb_Resistant_DCN_stringent$IFG$Normal,file = 'msbb_Resist_DCN_S_IFG_Normal.gml',format = 'gml')
write.graph(graph = msbb_Resistant_DCN_stringent$PHG$Normal,file = 'msbb_Resist_DCN_S_PHG_Normal.gml',format = 'gml')
write.graph(graph = msbb_Resistant_DCN_stringent$STG$Normal,file = 'msbb_Resist_DCN_S_STG_Normal.gml',format = 'gml')

write.graph(graph = msbb_Resistant_DCN_stringent$FP$AD,file = 'msbb_Resist_DCN_S_FP_AD.gml',format = 'gml')
write.graph(graph = msbb_Resistant_DCN_stringent$IFG$AD,file = 'msbb_Resist_DCN_S_IFG_AD.gml',format = 'gml')
write.graph(graph = msbb_Resistant_DCN_stringent$PHG$AD,file = 'msbb_Resist_DCN_S_PHG_AD.gml',format = 'gml')
write.graph(graph = msbb_Resistant_DCN_stringent$STG$AD,file = 'msbb_Resist_DCN_S_STG_AD.gml',format = 'gml')