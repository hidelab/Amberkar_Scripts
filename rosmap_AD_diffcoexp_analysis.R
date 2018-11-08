library(igraph)
library(synapseClient)
library(data.table)
library(org.Hs.eg.db)
library(dplyr)
library(magrittr)
library(diffcoexp)
#library(EnsDb.Hsapiens.v79)

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
#Download data from Synapse
# rosmap_reseq_data_pointer<-synGet(id='syn11807272')
# rosmap_reseq_data=fread(rosmap_reseq_data_pointer@filePath,sep = "\t",header = T,stringsAsFactors = F,showProgress = T,data.table = F)
# rosmap_reseq_data=rosmap_reseq_data%>%dplyr::mutate(GeneSymbol=unname(mapIds(x = EnsDb.Hsapiens.v79,keys = ensembl_gene_id,column = "SYMBOL",keytype = "GENEID")))
# rosmap_reseq_data.agg=aggregate.data.frame(x = rosmap_reseq_data[,-c(1,634)],by=list(symbol=rosmap_reseq_data$GeneSymbol),mean)
# rownames(rosmap_reseq_data.agg)=rosmap_reseq_data.agg$symbol
# rosmap_reseq_data2=rosmap_reseq_data.agg[,-1]

#Read ROSMAP covariates
rosmap_covariates=synGet("syn11024258")
rosmap_covariates.df=fread(rosmap_covariates@filePath,sep = "\t",header = T,stringsAsFactors = F,showProgress = T,data.table = F)
# rosmap_covariates.df=rosmap_covariates.df%>%mutate(diagnosis=if_else((cogdx=4 & braaksc>=4 & ceradsc <= 2),true = "AD",false = "OTHER"))
# rosmap_covariates.df=rosmap_covariates.df%>%mutate(diagnosis=if_else((cogdx=1 & braaksc<=3 & ceradsc >= 3),true = "CONTROL",false = "OTHER"))

rosmap_reseq_data2=readRDS("rosmap_reseq_data_agg.RDS")
c_rosmap_exprs=rosmap_reseq_data2[,rosmap_covariates.df%>%dplyr::filter(Diagnosis=="CONTROL")%>%pull(SampleID)]
t_rosmap_exprs=rosmap_reseq_data2[,rosmap_covariates.df%>%dplyr::filter(Diagnosis=="AD")%>%pull(SampleID)]
rosmap.diffcoexp=diffcoexp(exprs.1=c_rosmap_exprs,exprs.2=t_rosmap_exprs,r.method="spearman",rth=0.6,q.diffth=0.1,q.dcgth=0.1)
saveRDS(rosmap.diffcoexp,"rosmap_diffcoexp_CTRL_AD_out.RDS")
