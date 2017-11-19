library(igraph)
library(foreach)
library(synapseClient)
library(data.table)
library(org.Hs.eg.db)
library(doParallel)
library(DCGL)
library(dplyr)
library(magrittr)
library(devtools)
library(diffcoexp)
synapseLogin()
#cl=makeCluster(8)
#registerDoParallel(cl)
#Define functions to be used
mapIds2<-function(IDs,IDFrom,IDTo){
  idmap=mapIds(x = org.Hs.eg.db,keys = IDs,column = IDTo,keytype = IDFrom,multiVals = "first")
  na_vec=names(idmap[is.na(idmap)==T])
  idmap=idmap[is.na(idmap)==F]
  idmap_df=data.frame("From"=names(idmap),"To"=unlist(unname(idmap)),stringsAsFactors = F)
  return(list(map=idmap_df,noMap=na_vec))
}
#Set working directory
setwd("/shared/hidelab2/user/md4zsa/Work/Data/AMP-AD_RNAseq_ReSeq/Normalised_covariate_corrected_NoResiduals/ROSMAP")
#Download data from Synapse
rosmap_reseq_data_pointer<-synGet(id='syn8456719')
rosmap_reseq_data=fread(rosmap_reseq_data_pointer@filePath,sep = "\t",header = T,stringsAsFactors = F,showProgress = T)

#Read ROSMAP covariates
rosmap_covariates=read.table("ROSMAP_DLPFC_Covariates.tsv",header = T,sep = "\t",stringsAsFactors = F)

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
#rosmap_DCp=DCp(exprs.1 = rosmap_c_counts,exprs.2 = rosmap_t_counts,r.method = "spearman",link.method = "qth",cutoff = 0.05,N = 1000)
rosmap_diffcoexp=diffcoexp(exprs.1 = rosmap_c_counts,exprs.2 = rosmap_t_counts,rth=0.6, qth=0.1, r.diffth=0.1, q.diffth=0.1)
rosmap_diffcoexp_DRsort=DRsort(DCGs = rosmap_diffcoexp$DCGs,DCLs = rosmap_diffcoexp$DCLs,tf2target = regnet_tf2target,expGenes = rownames(rosmap_reseq_data2.agg))
saveRDS(rosmap_diffcoexp,"ROSMAP_DiffCoexp_noLFC_bugfix2.res.RDS")
saveRDS(rosmap_diffcoexp_DRsort,"ROSMAP_DiffCoexp_DRsort_noLFC_bugfix2.res.RDS")

# DCecutoff = 0.05
# rosmap_DCe.DCG <- rosmap_DCe$DCGs[rosmap_DCe$DCGs[, "q"] < DCecutoff, ]
# rosmap_DCe.DCG <- data.frame(DCG = rownames(rosmap_DCe.DCG), rosmap_DCe.DCG)
# DCG <-rosmap_DCe.DCG
# x<-rosmap_DCe$DCLs
# x<-subset(x, subset=(Gene.1 %in% rownames(DCG) | Gene.2 %in% rownames(DCG) ))
# expGenes<-rownames(rosmap_DCe$DCGs)
# rosmap_DRsort.res<- DRsort(DCG, x, regnet_tf2target, expGenes)
# saveRDS(rosmap_DRsort.res,"ROSMAP_DCe_DRsort_q005.res.RDS")
#saveRDS(rosmap_DCe,"ROSMAP_DCe_AnalysisResults.RDS")

proc.time()
cat(paste("Completed!"))
stopCluster(cl)