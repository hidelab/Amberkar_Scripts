library(igraph)
library(synapseClient)
library(data.table)
library(org.Hs.eg.db)
library(dplyr)
library(magrittr)
library(devtools)
library(DESeq2)
library(EnsDb.Hsapiens.v79)
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
rosmap_reseq_data_pointer<-synGet(id='syn8456637')
rosmap_reseq_data=fread(rosmap_reseq_data_pointer@filePath,sep = "\t",header = T,stringsAsFactors = F,showProgress = T)
#Read ROSMAP covariates
rosmap_covariates=read.table("ROSMAP_DLPFC_Covariates.tsv",header = T,sep = "\t",stringsAsFactors = F)

#Collapse Ensembl IDs to gene symbols. For multiple Ensembl IDs for same gene, compute average expression
rosmap_reseq_data$gene_symbol=unname(mapIds(x = EnsDb.Hsapiens.v79,keys=rosmap_reseq_data$ensembl_gene_id,keytype = "GENEID",column = "SYMBOL",multiVals = "first"))
rosmap_reseq_data.agg=aggregate.data.frame(x=rosmap_reseq_data[,-c('ensembl_gene_id','gene_symbol')],by=list(Symbol=rosmap_reseq_data$gene_symbol),mean)
rownames(rosmap_reseq_data.agg)=rosmap_reseq_data.agg$Symbol
rosmap_reseq_data.agg=rosmap_reseq_data.agg[,-1]


rosmap_covariates$SampleType="OTHER"
rosmap_covariates$SampleType[rosmap_covariates$cogdx==1]="CONTROL"
rosmap_covariates$SampleType[rosmap_covariates$cogdx==2]="COGDX2"
rosmap_covariates$SampleType[rosmap_covariates$cogdx==4]="COGDX4"
rosmap_covariates$SampleType[rosmap_covariates$cogdx==5]="COGDX5"
rosmap_coldata=cbind.data.frame(condition=rosmap_covariates$SampleType,type="paired-end")
rownames(rosmap_coldata)=colnames(rosmap_reseq_data.agg)

rosmap_reseq_counts.list=vector(mode = "list",length = 3)
names(rosmap_reseq_counts.list)=c("cogdx2","cogdx4","cogdx5")
rosmap_reseq_counts.list$cogdx2=rosmap_reseq_data.agg[,rownames(rosmap_coldata[which(rosmap_covariates$SampleType%in%c("CONTROL","COGDX2")),])]
rosmap_reseq_counts.list$cogdx4=rosmap_reseq_data.agg[,rownames(rosmap_coldata[which(rosmap_covariates$SampleType%in%c("CONTROL","COGDX4")),])]
rosmap_reseq_counts.list$cogdx5=rosmap_reseq_data.agg[,rownames(rosmap_coldata[which(rosmap_covariates$SampleType%in%c("CONTROL","COGDX5")),])]

rosmap_reseq_counts.list$cogdx2=data.frame(apply(rosmap_reseq_counts.list$cogdx2,2,function(x)round(x = x,digits = 0)))
rosmap_reseq_counts.list$cogdx4=data.frame(apply(rosmap_reseq_counts.list$cogdx4,2,function(x)round(x = x,digits = 0)))
rosmap_reseq_counts.list$cogdx5=data.frame(apply(rosmap_reseq_counts.list$cogdx5,2,function(x)round(x = x,digits = 0)))

rosmap_reseq_counts.list=lapply(rosmap_reseq_counts.list,function(x){colnames(x) <- gsub(pattern = "X",replacement = "",x = colnames(x));x})
rosmap_reseq_counts.dds=vector(mode = "list",length = 3)
names(rosmap_reseq_counts.dds)=c("cogdx2","cogdx4","cogdx5")
for(i in 1:3){
  rosmap_reseq_counts.dds[[i]]=DESeqDataSetFromMatrix(countData = rosmap_reseq_counts.list[[i]],colData = rosmap_coldata[rownames(rosmap_coldata)%in%colnames(rosmap_reseq_counts.list[[i]]),],design = ~condition)
}
rosmap_reseq_counts.dds=lapply(rosmap_reseq_counts.dds,DESeq)
rosmap_reseq_counts.dds=lapply(rosmap_reseq_counts.dds,function(x){levels(x$condition) <- levels(x$condition)[2:1];x})
rosmap_reseq_counts.res=lapply(rosmap_reseq_counts.dds,results)

saveRDS(rosmap_reseq_counts.dds,"rosmap_reseq_counts_dds.RDS")
saveRDS(rosmap_reseq_counts.res,"rosmap_reseq_counts_res.RDS")
