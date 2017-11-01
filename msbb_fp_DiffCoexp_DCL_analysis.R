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
setwd("/shared/hidelab2/user/md4zsa/Work/Data/AMP-AD_RNAseq_ReSeq/Normalised_covariate_corrected_NoResiduals/MSMM")
#Download data from Synapse
msbb_reseq_data_pointer<-synGet(id='syn8485027')
msbb_reseq_data=fread("MSSM_FP_STG_PHG_IFG_netResidualExpression.tsv",sep = "\t",header = T,stringsAsFactors = F,data.table = F,showProgress = T)
msbb_reseq_covariates_pointer<-synGet(id='syn11024323')
msbb_covariates=fread(msbb_reseq_covariates_pointer@filePath,header = T,sep = "\t",stringsAsFactors = F)


#Collapse Ensembl IDs to gene symbols. For multiple Ensembl IDs for same gene, compute average expression
ensembl_geneSymbol_map=mapIds2(IDs = msbb_reseq_data$ensembl_gene_id,IDFrom = "ENSEMBL",IDTo = "SYMBOL")
msbb_reseq_data2=msbb_reseq_data[-which(msbb_reseq_data$ensembl_gene_id%in%mapIds2(IDs = msbb_reseq_data$ensembl_gene_id,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[2]]),]
msbb_reseq_data2$gene_symbol=mapIds2(IDs = msbb_reseq_data2$ensembl_gene_id,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[1]][,2]
msbb_reseq_data2.agg=aggregate(x=msbb_reseq_data2[,-c(1,dim(msbb_reseq_data2)[2])],by=list(Symbol=msbb_reseq_data2$gene_symbol),mean)
rownames(msbb_reseq_data2.agg)=msbb_reseq_data2.agg$Symbol
msbb_reseq_data2.agg=msbb_reseq_data2.agg[,-1]

#Read TF-Target interactions, preprocess data
regnet_tf2target=fread("/shared/hidelab2/user/md4zsa/Work/Data/TF_Databases/RegNetwork_human_regulators2.txt",header = T,sep = "\t",showProgress = T,data.table = F)%>%filter(evidence=="Experimental")%>%select(c(regulator_symbol,target_symbol))

#Segragate Control and AD samples
msbb_fp_c_counts=msbb_reseq_data2.agg[,msbb_covariates$SampleID[msbb_covariates$BrainRegion.Diagnosis=="FP.CONTROL"]]
msbb_fp_t_counts=msbb_reseq_data2.agg[,msbb_covariates$SampleID[msbb_covariates$BrainRegion.Diagnosis=="FP.AD"]]

msbb_fp_DiffCoexp=diffcoexp(exprs.1 = msbb_fp_c_counts,exprs.2 = msbb_fp_t_counts,rth=0.6, qth=0.1, r.diffth=0.1, q.diffth=0.1)
msbb_fp_DRsort=DRsort(DCGs = msbb_fp_DiffCoexp$DCGs,DCLs = msbb_fp_DiffCoexp$DCLs,tf2target = regnet_tf2target,expGenes = rownames(msbb_reseq_data2.agg))
saveRDS(msbb_fp_DiffCoexp,"MSBB_FP_DiffCoexp_noLFC_bugfix.res.RDS")
saveRDS(msbb_fp_DRsort,"MSBB_FP_DRsort_noLFC_bugfix.res.RDS")

proc.time()
cat(paste("Completed!"))
stopCluster(cl)