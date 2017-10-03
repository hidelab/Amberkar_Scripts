library(igraph)
library(foreach)
library(cocor)
library(synapseClient)
library(data.table)
library(org.Hs.eg.db)
library(doParallel)
library(DCGL)
library(dplyr)
library(magrittr)

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
setwd("/shared/hidelab2/user/md4zsa/Work/Data/AMP-AD_RNAseq_ReSeq/Normalised_covariate_corrected_NoResiduals")
#Download data from Synapse
rosmap_reseq_data_pointer<-synGet(id='syn8456719')
rosmap_reseq_data=fread(rosmap_reseq_data_pointer@filePath,sep = "\t",header = T,stringsAsFactors = F,showProgress = T)

#Read ROSMAP covariates
rosmap_covariates=read.table("ROSMAP/ROSMAP_DLPFC_Covariates.tsv",header = T,sep = "\t",stringsAsFactors = F)

#Collapse Ensembl IDs to gene symbols. For multiple Ensembl IDs for same gene, compute average expression
ensembl_geneSymbol_map=mapIds2(IDs = rosmap_reseq_data$ensembl_gene_id,IDFrom = "ENSEMBL",IDTo = "SYMBOL")
rosmap_reseq_data2=rosmap_reseq_data[-which(rosmap_reseq_data$ensembl_gene_id%in%mapIds2(IDs = rosmap_reseq_data$ensembl_gene_id,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[2]]),]
rosmap_reseq_data2$gene_symbol=mapIds2(IDs = rosmap_reseq_data2$ensembl_gene_id,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[1]][,2]
rosmap_reseq_data2.agg=aggregate(x=rosmap_reseq_data2[,-c(1,634)],by=list(Symbol=rosmap_reseq_data2$gene_symbol),mean)
rownames(rosmap_reseq_data2.agg)=rosmap_reseq_data2.agg$Symbol
rosmap_reseq_data2.agg=rosmap_reseq_data2.agg[,-1]

#Read TF-Target interactions, preprocess data
regnet_tf2target=fread("/shared/hidelab2/user/md4zsa/Work/Data/TF_Databases/RegNetwork_human_regulators2.txt",header = T,sep = "\t",showProgress = T,data.table = F)%>%filter(evidence=="Experimental")%>%select(c(regulator_symbol,target_symbol))

#Segragate Control and AD samples
rosmap_c_counts=rosmap_reseq_data2.agg[,rosmap_covariates$SampleID[rosmap_covariates$Diagnosis=="CONTROL"]]
rosmap_t_counts=rosmap_reseq_data2.agg[,rosmap_covariates$SampleID[rosmap_covariates$Diagnosis=="AD"]]
rosmap_DCp=DCp(exprs.1 = rosmap_c_counts,exprs.2 = rosmap_t_counts,r.method = "spearman",link.method = "qth",cutoff = 0.05,N = 1000)
rosmap_DCe=DCe(exprs.1 = rosmap_c_counts,exprs.2 = rosmap_t_counts,r.method = "spearman",p = 0.05,link.method = "qth",cutoff = 0.05)
rosmap_DCsum.res=DCsum(rosmap_DCp,rosmap_DCe,DCpcutoff=0.05,DCecutoff=0.05)
rosmap_DRsort.res=DRsort(DCGs = rosmap_DCsum.res$DCGs,DCLs = rosmap_DCsum.res$DCLs,tf2target = regnet_tf2target,expGenes = rosmap_reseq_data2.agg)
saveRDS(rosmap_DCp,"ROSMAP_DCp_AnalysisResults.RDS")
saveRDS(rosmap_DCe,"ROSMAP_DCe_AnalysisResults.RDS")
saveRDS(rosmap_DCsum.res,"ROSMAP_DCsum_AnalysisResults.RDS")
saveRDS(rosmap_DRsort.res,"ROSMAP_DRsort_AnalysisResults.RDS")

proc.time()
cat(paste("Completed!"))
stopCluster(cl)