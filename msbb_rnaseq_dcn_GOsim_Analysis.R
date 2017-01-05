library(org.Hs.eg.db)
library(igraph)
library(GOSemSim)
require(parallel)

setwd('/shared/hidelab2/user/md4zsa/Work/Data/MSMM_RNAseq/MSMM_RNAseq_FinalRelease2')
fp_dcn=readRDS('FP/DCN/MSBB_RNAseq_FP.DCN.RDS')
ifg_dcn=readRDS('IFG/DCN/MSBB_RNAseq_IFG_FDR01.DCN.RDS')
phg_dcn=readRDS('PHGG/DCN/MSBB_RNAseq_PHG_FDR01.DCN.RDS')
stg_dcn=readRDS('STG/DCN/MSBB_RNAseq_STG_FDR01.DCN.RDS')
hsGO.BP=godata(OrgDb='org.Hs.eg.db',keytype='SYMBOL',ont='BP')
hsGO.CC=godata(OrgDb='org.Hs.eg.db',keytype='SYMBOL',ont='CC')
hsGO.MF=godata(OrgDb='org.Hs.eg.db',keytype='SYMBOL',ont='MF')
msbb_dcn.clusters=list()
msbb_dcn.clusters[[1]]=V(fp_dcn)$name
msbb_dcn.clusters[[2]]=V(ifg_dcn)$name
msbb_dcn.clusters[[3]]=V(phg_dcn)$name
msbb_dcn.clusters[[4]]=V(stg_dcn)$name
names(msbb_dcn.clusters)=c('FP','IFG','PHG','STG')
cat(paste("Computing GO.BP semsim ...\n",sep = ""))
msbb_dcn.clusters_GOBP=lapply(msbb_dcn.clusters,mgeneSim,semData=hsGO.BP, measure="Wang", combine="BMA", verbose=T)
saveRDS(msbb_dcn.clusters_GOBP)
cat(paste("Computing GO.CC semsim ...\n",sep = ""))
msbb_dcn.clusters_GOCC=lapply(msbb_dcn.clusters,mgeneSim,semData=hsGO.CC, measure="Wang", combine="BMA", verbose=T)
saveRDS(msbb_dcn.clusters_GOCC)
cat(paste("Computing GO.MF semsim ...\n",sep = ""))
msbb_dcn.clusters_GOMF=lapply(msbb_dcn.clusters,mgeneSim,semData=hsGO.MF, measure="Wang", combine="BMA", verbose=T)
saveRDS(msbb_dcn.clusters_GOMF)

