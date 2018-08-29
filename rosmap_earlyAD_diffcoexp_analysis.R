library(igraph)
library(synapseClient)
library(data.table)
library(org.Hs.eg.db)
library(dplyr)
library(magrittr)
library(devtools)
synapseLogin()

#Set working directory
setwd("/shared/hidelab2/user/md4zsa/Work/Data/AMP-AD_RNAseq_ReSeq/Normalised_covariate_corrected_NoResiduals/ROSMAP")
#Download data from Synapse - adjusted for network analysis
rosmap_reseq_data_pointer<-synGet(id='syn8456719')
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

rosmap_control.counts=rosmap_reseq_data.agg[,rosmap_covariates%>%dplyr::filter(SampleType=="CONTROL")%>%dplyr::pull(SampleID)]
rosmap_cogdx2.counts=rosmap_reseq_data.agg[,rosmap_covariates%>%dplyr::filter(SampleType=="COGDX2")%>%dplyr::pull(SampleID)]
rosmap_diffcoexp=diffcoexp(exprs.1 = rosmap_control.counts,exprs.2 = rosmap_cogdx2.counts,r.method="spearman",rth=0.6,q.diffth=0.1,q.dcgth=0.1)
saveRDS(rosmap_diffcoexp,"ROSMAP_earlyAD_diffcoexp_results.RDS")


