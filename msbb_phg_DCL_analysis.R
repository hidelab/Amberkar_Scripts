library(igraph)
library(foreach)
library(cocor)
library(synapseClient)
library(data.table)
library(org.Hs.eg.db)
library(doParallel)
library(DCGL)

synapseLogin(username = "s.amberkar@sheffield.ac.uk",apiKey = "evb/5m+/10KmKAOP2vS1G6+a20iWAQlDosD9UfoQhvvFUdip/R/kZCzuk3jYecQ7zti5F4ZePz8djJQ8PoRC6Q==",rememberMe = T)
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
msbb_reseq_data=fread(msbb_reseq_data_pointer@filePath,sep = "\t",header = T,stringsAsFactors = F,showProgress = T)
msbb_reseq_covariates_pointer<-synGet(id='syn11024323')
msbb_covariates=read.table(msbb_reseq_covariates_pointer@filePath,header = T,sep = "\t",stringsAsFactors = F)

#Collapse Ensembl IDs to gene symbols. For multiple Ensembl IDs for same gene, compute average expression
ensembl_geneSymbol_map=mapIds2(IDs = msbb_reseq_data$ensembl_gene_id,IDFrom = "ENSEMBL",IDTo = "SYMBOL")
msbb_reseq_data2=msbb_reseq_data[-which(msbb_reseq_data$ensembl_gene_id%in%mapIds2(IDs = msbb_reseq_data$ensembl_gene_id,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[2]]),]
msbb_reseq_data2$gene_symbol=mapIds2(IDs = msbb_reseq_data2$ensembl_gene_id,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[1]][,2]
msbb_reseq_data2.agg=aggregate(x=msbb_reseq_data2[,-c(1,634)],by=list(Symbol=msbb_reseq_data2$gene_symbol),mean)
rownames(msbb_reseq_data2.agg)=msbb_reseq_data2.agg$Symbol
msbb_reseq_data2.agg=msbb_reseq_data2.agg[,-1]

#Read TF-Target interactions, preprocess data
regnet_tf2target=fread("/shared/hidelab2/user/md4zsa/Work/Data/TF_Databases/RegNetwork_human_regulators2.txt",header = T,sep = "\t",showProgress = T,data.table = F)%>%filter(evidence=="Experimental")%>%select(c(regulator_symbol,target_symbol))

#Segragate Control and AD samples
msbb_phg_c_counts=msbb_reseq_data2.agg[,msbb_covariates$SampleID[msbb_covariates$BrainRegion.Diagnosis=="PHG.CONTROL"]]
msbb_phg_t_counts=msbb_reseq_data2.agg[,msbb_covariates$SampleID[msbb_covariates$BrainRegion.Diagnosis=="PHG.AD"]]

msbb_phg_DCe=DCe(exprs.1 = msbb_phg_c_counts,exprs.2 = msbb_phg_t_counts,r.method = "spearman",p = 0.05,link.method = "qth",cutoff = 0.05)
#msbb_phg_DCsum.res=DCsum(msbb_phg_DCp,msbb_phg_DCe,DCpcutoff=0.05,DCecutoff=0.05)
DCecutoff = 0.25
msbb_phg_DCe.DCG <- msbb_phg_DCe$DCGs[msbb_phg_DCe$DCGs[, "q"] < DCecutoff, ]
msbb_phg_DCe.DCG <- data.frame(DCG = rownames(msbb_phg_DCe.DCG), msbb_phg_DCe.DCG)
DCG <-msbb_phg_DCe.DCG
x<-msbb_phg_DCe$DCLs
x<-subset(x, subset=(Gene.1 %in% rownames(DCG) | Gene.2 %in% rownames(DCG) ))
expGenes<-rownames(msbb_phg_DCe$DCGs)
msbb_phg_DCe_DRsort.res<- DRsort(DCG, x, regnet_tf2target, expGenes)
saveRDS(msbb_phg_DCe_DRsort.res,"MSBB_PHG_DCe_DRsort.res.RDS")
saveRDS(msbb_phg_DCe,"MSBB_PHG_DCe_res.RDS")

proc.time()
cat(paste("Completed!"))
stopCluster(cl)