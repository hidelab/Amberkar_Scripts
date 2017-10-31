library(igraph)
library(foreach)
library(cocor)
library(synapseClient)
library(data.table)
library(org.Hs.eg.db)
# library(doParallel)
library(DCGL)
library(dplyr)
library(magrittr)
source("/shared/hidelab2/shared/wenbin/diffcoexp/diffcoexp.R")
synapseLogin()
# cl=makeCluster(8)
# registerDoParallel(cl)
#Define functions to be used
mapIds2<-function(IDs,IDFrom,IDTo){
  idmap=mapIds(x = org.Hs.eg.db,keys = IDs,column = IDTo,keytype = IDFrom,multiVals = "first")
  na_vec=names(idmap[is.na(idmap)==T])
  idmap=idmap[is.na(idmap)==F]
  idmap_df=data.frame("From"=names(idmap),"To"=unlist(unname(idmap)),stringsAsFactors = F)
  return(list(map=idmap_df,noMap=na_vec))
}

#Set working directory
setwd("/shared/hidelab2/user/md4zsa/Work/Data/AMP-AD_RNAseq_ReSeq/Normalised_covariate_corrected_NoResiduals/MAYO")
#Download data from Synapse
mayo_reseq_data_pointer<-synGet(id='syn8466826')
mayo_reseq_data=fread(mayo_reseq_data_pointer@filePath,sep = "\t",header = T,stringsAsFactors = F,showProgress = T,data.table = F)
#Read ROSMAP covariates
mayo_covariates_pointer <- synGet(id='syn8466814')
mayo_covariates=read.table(mayo_covariates_pointer@filePath,header = T,sep = "\t",stringsAsFactors = F)

#Collapse Ensembl IDs to gene symbols. For multiple Ensembl IDs for same gene, compute average expression
ensembl_geneSymbol_map=mapIds2(IDs = mayo_reseq_data$ensembl_gene_id,IDFrom = "ENSEMBL",IDTo = "SYMBOL")
mayo_reseq_data2=mayo_reseq_data[-which(mayo_reseq_data$ensembl_gene_id%in%mapIds2(IDs = mayo_reseq_data$ensembl_gene_id,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[2]]),]
mayo_reseq_data2$gene_symbol=mapIds2(IDs = mayo_reseq_data2$ensembl_gene_id,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[1]][,2]
mayo_reseq_data2.agg=aggregate(x=mayo_reseq_data2[,-c(1,(length(colnames(mayo_reseq_data2))))],by=list(Symbol=mayo_reseq_data2$gene_symbol),mean)
rownames(mayo_reseq_data2.agg)=mayo_reseq_data2.agg$Symbol
mayo_reseq_data2.agg=mayo_reseq_data2.agg[,-1]

#Segregate data into TCX and CER brain regions
mayo_reseq_tcx_data=mayo_reseq_data2.agg[,which(colnames(mayo_reseq_data2.agg)%in%mayo_covariates$SampleID[grep(pattern = "TCX.AD|TCX.Control",mayo_covariates$BrainRegion.Diagnosis)])]
mayo_reseq_cer_data=mayo_reseq_data2.agg[,which(colnames(mayo_reseq_data2.agg)%in%mayo_covariates$SampleID[grep(pattern = "CER.AD|CER.Control",mayo_covariates$BrainRegion.Diagnosis)])]

#Read TF-Target interactions, preprocess data
regnet_tf2target=fread("/shared/hidelab2/user/md4zsa/Work/Data/TF_Databases/RegNetwork_human_regulators2.txt",header = T,sep = "\t",showProgress = T,data.table = F)%>%filter(evidence=="Experimental")%>%select(c(regulator_symbol,target_symbol))

#Segragate Control and AD samples
cer_c_counts=mayo_reseq_cer_data[,mayo_covariates$SampleID[grep(pattern = "CER.Control",mayo_covariates$BrainRegion.Diagnosis)]]
cer_t_counts=mayo_reseq_cer_data[,mayo_covariates$SampleID[grep(pattern = "CER.AD",mayo_covariates$BrainRegion.Diagnosis)]]
#tcx_DCp=DCp(exprs.1 = tcx_c_counts,exprs.2 = tcx_t_counts,r.method = "spearman",link.method = "qth",cutoff = 0.05,N = 1000)
cer_DiffCoexp=diffcoexp(exprs.1 = cer_c_counts,exprs.2 = cer_t_counts,rth=0.6, qth=0.1, r.diffth=0.1, q.diffth=0.1)
cer_DRsort=DRsort(DCGs = cer_DiffCoexp$DCGs,DCLs = cer_DiffCoexp$DCLs,tf2target = regnet_tf2target,expGenes = rownames(mayo_reseq_tcx_data))
saveRDS(cer_DiffCoexp,"CER_DiffCoexp_noLFC_Results.RDS")
saveRDS(cer_DRsort,"CER_DiffCoexp_DRsort_noLFC.res.RDS")

# DCecutoff = 0.05
# tcx_DCe.DCG <- tcx_DCe$DCGs[tcx_DCe$DCGs[, "q"] < DCecutoff, ]
# tcx_DCe.DCG <- data.frame(DCG = rownames(tcx_DCe.DCG), tcx_DCe.DCG)
# DCG<-tcx_DCe.DCG
# x<-tcx_DCe$DCLs
# x<-subset(x, subset=(Gene.1 %in% rownames(DCG) | Gene.2 %in% rownames(DCG) ))
# expGenes<-rownames(tcx_DCe$DCGs)
# tcx_DRsort.res<- DRsort(DCG, x, regnet_tf2target, expGenes)
# saveRDS(tcx_DRsort.res,"TCX_DRsort_Results_q005.RDS")

proc.time()
cat(paste("Completed!"))
stopCluster(cl)
