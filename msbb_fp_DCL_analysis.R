library(igraph)
library(foreach)
library(cocor)
library(synapseClient)
library(data.table)
library(org.Hs.eg.db)
library(doParallel)
library(DCGL)

synapseLogin()
cl=makeCluster(8)
registerDoParallel(cl)
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
msbb_reseq_data=fread(msbb_reseq_data_pointer@filePath,sep = "\t",header = T,stringsAsFactors = F,data.table = F,showProgress = T)
msbb_reseq_covariates_pointer<-synGet(id='syn11024323')
msbb_covariates=fread(msbb_reseq_covariates_pointer@filePath,header = T,sep = "\t",stringsAsFactors = F)

#Collapse Ensembl IDs to gene symbols. For multiple Ensembl IDs for same gene, compute average expression
ensembl_geneSymbol_map=mapIds2(IDs = msbb_reseq_data$ensembl_gene_id,IDFrom = "ENSEMBL",IDTo = "SYMBOL")
msbb_reseq_data2=msbb_reseq_data[-which(msbb_reseq_data$ensembl_gene_id%in%mapIds2(IDs = msbb_reseq_data$ensembl_gene_id,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[2]]),]
msbb_reseq_data2$gene_symbol=mapIds2(IDs = msbb_reseq_data2$ensembl_gene_id,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[1]][,2]
msbb_reseq_data2.agg=aggregate(x=msbb_reseq_data2[,-c(1,'gene_symbol')],by=list(Symbol=msbb_reseq_data2$gene_symbol),mean)
rownames(msbb_reseq_data2.agg)=msbb_reseq_data2.agg$Symbol
msbb_reseq_data2.agg=msbb_reseq_data2.agg[,-1]

#Read TF-Target interactions, preprocess data
regnet_tf2target=fread("/shared/hidelab2/user/md4zsa/Work/Data/TF_Databases/RegNetwork_human_regulators2.txt",header = T,sep = "\t",showProgress = T,data.table = F)%>%filter(evidence=="Experimental")%>%select(c(regulator_symbol,target_symbol))

#Segragate Control and AD samples
msbb_fp_c_counts=msbb_reseq_data2.agg[,msbb_covariates$SampleID[msbb_covariates$BrainRegion.Diagnosis=="FP.CONTROL"]]
msbb_fp_t_counts=msbb_reseq_data2.agg[,msbb_covariates$SampleID[msbb_covariates$BrainRegion.Diagnosis=="FP.AD"]]

msbb_fp_DCe=DCe(exprs.1 = msbb_fp_c_counts,exprs.2 = msbb_fp_t_counts,r.method = "spearman",p = 0.05,link.method = "qth",cutoff = 0.05)
#msbb_fp_DCsum.res=DCsum(msbb_fp_DCp,msbb_fp_DCe,DCpcutoff=0.05,DCecutoff=0.05)
DCecutoff = 0.25
msbb_fp_DCe.DCG <- msbb_fp_DCe$DCGs[msbb_fp_DCe$DCGs[, "q"] < DCecutoff, ]
msbb_fp_DCe.DCG <- data.frame(DCG = rownames(msbb_fp_DCe.DCG), msbb_fp_DCe.DCG)
DCG <-msbb_fp_DCe.DCG
x<-msbb_fp_DCe$DCLs
x<-subset(x, subset=(Gene.1 %in% rownames(DCG) | Gene.2 %in% rownames(DCG) ))
expGenes<-rownames(msbb_fp_DCe$DCGs)
msbb_fp_DCe_DRsort.res<- DRsort(DCG, x, regnet_tf2target, expGenes)
saveRDS(msbb_fp_DCe_DRsort.res,"MSBB_FP_DCe_DRsort.res.RDS")
saveRDS(msbb_fp_DCe,"MSBB_FP_DCe_res.RDS")

proc.time()
cat(paste("Completed!"))
stopCluster(cl)