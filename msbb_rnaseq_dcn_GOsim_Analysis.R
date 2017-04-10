library(org.Hs.eg.db)
library(igraph)
library(GOSemSim)
require(parallel)
library()

setwd('/shared/hidelab2/user/md4zsa/Work/Data/MSMM_RNAseq/MSMM_RNAseq_FinalRelease2')
fp_dcn=readRDS('FP/DCN/FP_DCN_FDR01.RDS') 
ifg_dcn=readRDS('IFG/DCN/IFG_DCN_FDR01.RDS') 
phg_dcn=readRDS('PHG/DCN/PHG_DCN_FDR01.RDS') 
stg_dcn=readRDS('STG/DCN/STG_DCN_FDR01.RDS') 

fp_plq_dcn=readRDS('FP/PLQ_DCN/FP_PLQGenes_FDR01.DCN.RDS') 
ifg_plq_dcn=readRDS('IFG/PLQ_DCN/IFG_PLQGenes_FDR01.DCN.RDS') 
phg_plq_dcn=readRDS('PHG/PLQ_DCN/PHG_PLQGenes_FDR01.DCN.RDS') 
stg_plq_dcn=readRDS('STG/PLQ_DCN/STG_PLQGenes_FDR01.DCN.RDS')

hsGO=vector(mode = "list",length = 3)
names(hsGO)=c('BP','CC','MF')
hsGO$BP=godata(OrgDb='org.Hs.eg.db',keytype='SYMBOL',ont='BP')
hsGO$CC=godata(OrgDb='org.Hs.eg.db',keytype='SYMBOL',ont='CC')
hsGO$MF=godata(OrgDb='org.Hs.eg.db',keytype='SYMBOL',ont='MF')
msbb_dcn.clusters=msbb_plq_dcn.clusters=list()
msbb_dcn.clusters[[1]]=V(fp_dcn$High)$name
msbb_dcn.clusters[[2]]=V(ifg_dcn$High)$name
msbb_dcn.clusters[[3]]=V(phg_dcn$High)$name
msbb_dcn.clusters[[4]]=V(stg_dcn$High)$name

msbb_plq_dcn.clusters[[1]]=V(fp_plq_dcn$High)$name
msbb_plq_dcn.clusters[[2]]=V(ifg_plq_dcn$High)$name
msbb_plq_dcn.clusters[[3]]=V(phg_plq_dcn$High)$name
msbb_plq_dcn.clusters[[4]]=V(stg_plq_dcn$High)$name

names(msbb_dcn.clusters)=names(msbb_plq_dcn.clusters)=c('FP','IFG','PHG','STG')
msbb_dcn.clusters=lapply(msbb_dcn.clusters,grep,pattern='^LOC|^LINC',invert=T,value=T)
msbb_semsim=vector(mode = "list",length = 4)
names(msbb_semsim)=c('FP','IFG','PHG','STG')
msbb_semsim$FP=msbb_semsim$IFG=msbb_semsim$PHG=msbb_semsim$STG=vector(mode = "list",length = 3)
names(msbb_semsim$FP)=names(msbb_semsim$IFG)=names(msbb_semsim$PHG)=names(msbb_semsim$STG)=c("BP","CC","MF")

nc = 8
blocksize=50

for(r in 4:1){
  #cat(paste("In brain region ...",names(msbb_dcn.clusters)[r],"\n",sep=" "))
  for (ont in 1:length(hsGO)){
    #cat(paste("Computing semantic similarity for ...",names(hsGO)[ont],"\n",sep=" "))    
    number_of_network_genes=length(msbb_dcn.clusters[[r]])
    i<-0
    start<-i*blocksize+1
    end<-min((i+1)*blocksize, number_of_network_genes)
    go_semsim=c()
    while(start < number_of_network_genes){
      cat(paste("In brain region,",names(msbb_dcn.clusters)[r],"computing semantic similarity for",names(hsGO)[ont],"processing block ...",i,"\n",sep=" "))
      input=msbb_dcn.clusters[[r]][start:end]
      pb = txtProgressBar(min=0,max=length(input),style=3,initial=0)
      res = mgeneSim(genes = input,semData = hsGO[[ont]],measure = 'Wang',drop = 'IEA',combine = 'BMA',verbose = T)
      res2=rowSums(res)
      go_semsim=append(go_semsim,res2,after = length(go_semsim))
      i<-i+1
      start<-i*blocksize+1
      end<-min((i+1)*blocksize, number_of_network_genes)
      }
      msbb_semsim[[r]][[ont]]=go_semsim      
  }
  
}






cat(paste("Computing semsim ...\n",sep = ""))
msbb_dcn.clusters_GOBP=lapply(msbb_dcn.clusters[-1],mgeneSim,semData=hsGO.BP, measure="Wang", combine="BMA", verbose=T)
saveRDS(msbb_dcn.clusters_GOBP)
cat(paste("Computing GO.CC semsim ...\n",sep = ""))
msbb_dcn.clusters_GOCC=lapply(msbb_dcn.clusters,mgeneSim,semData=hsGO.CC, measure="Wang", combine="BMA", verbose=T)
saveRDS(msbb_dcn.clusters_GOCC)
cat(paste("Computing GO.MF semsim ...\n",sep = ""))
msbb_dcn.clusters_GOMF=lapply(msbb_dcn.clusters,mgeneSim,semData=hsGO.MF, measure="Wang", combine="BMA", verbose=T)
saveRDS(msbb_dcn.clusters_GOMF)

