library(foreach)
library(doParallel)

msbb_rnaseq2016_data=read.table("AMP-AD_MSBB_MSSM_IlluminaHiSeq2500_normalized_counts_September_2016.txt",sep="\t",header = T,as.is = T)
msbb_rnaseq_covariates=read.csv("MSBB_RNAseq_covariates.csv",header = T,as.is = T)
msbb_rnaseq_clinical_covariates=read.csv("MSBB_clinical.csv",header = T,as.is = T)
colnames(msbb_rnaseq2016_data)=unlist(lapply(strsplit(x = colnames(msbb_rnaseq2016_data),split = "X"),`[[`,2))
msbb_ensembl_symbol=data.frame(Ensembl=names(mapIds(x = org.Hs.eg.db,keys=rownames(msbb_rnaseq2016_data),column = "SYMBOL",keytype = "ENSEMBL",multiVals = "first")),Symbol=mapIds(x = org.Hs.eg.db,keys=rownames(msbb_rnaseq2016_data),column = "SYMBOL",keytype = "ENSEMBL",multiVals = "first"),stringsAsFactors = F)
msbb_rnaseq_covariates.merged=merge(x = msbb_rnaseq_clinical_covariates,y=msbb_rnaseq_covariates,by=c("individualIdentifier","individualIdentifier"))
msbb_rnaseq_covariates.merged2=msbb_rnaseq_covariates.merged[grep(pattern = "unmapped|resequenced",x = msbb_rnaseq_covariates.merged$fileName,invert = T),]


msbb_rnaseq2016_byRegion=msbb_rnaseq_covariates.merged_final=msbb_rnaseq2016_PLQGenes=vector(mode = "list",length = 4)
names(msbb_rnaseq2016_byRegion)=names(msbb_rnaseq_covariates.merged_final)=names(msbb_rnaseq2016_PLQGenes)=c("FP","IFG","PHG","STG") 
msbb_rnaseq2016_data2=msbb_rnaseq2016_data[-which(rownames(msbb_rnaseq2016_data)%in%msbb_ensembl_symbol$Ensembl[which(is.na(msbb_ensembl_symbol$Symbol)==T)]),]
msbb_rnaseq2016_data2=data.frame(cbind(Symbol=msbb_ensembl_symbol$Symbol[-which(is.na(msbb_ensembl_symbol$Symbol)==T)],msbb_rnaseq2016_data2),stringsAsFactors = F)
msbb_rnaseq2016_data2.agg=aggregate(x=msbb_rnaseq2016_data2[,-1],by=list(geneSymbol=msbb_rnaseq2016_data2$Symbol),mean)
rm_genes=grep(pattern = "_|\\.|^RP|-",msbb_rnaseq2016_data2.agg$geneSymbol,invert = F)
rownames(msbb_rnaseq2016_data2.agg)=msbb_rnaseq2016_data2.agg$geneSymbol
msbb_rnaseq2016_data2.agg=msbb_rnaseq2016_data2.agg[-rm_genes,-1]
colnames(msbb_rnaseq2016_data2.agg)=gsub(pattern = "X",replacement = "",x = colnames(msbb_rnaseq2016_data2.agg))

msbb_rnaseq_covariates.merged_final$FP=msbb_rnaseq_covariates.merged2[which(unlist(lapply(strsplit(x=msbb_rnaseq_covariates.merged2$sampleIdentifier,split = "_"),`[[`,3))%in%colnames(msbb_rnaseq2016_byRegion$FP)),c(1:11,13,15:16)]
msbb_rnaseq_covariates.merged_final$IFG=msbb_rnaseq_covariates.merged2[which(unlist(lapply(strsplit(x=msbb_rnaseq_covariates.merged2$sampleIdentifier,split = "_"),`[[`,3))%in%colnames(msbb_rnaseq2016_byRegion$IFG)),c(1:11,13,15:16)]
msbb_rnaseq_covariates.merged_final$PHG=msbb_rnaseq_covariates.merged2[which(unlist(lapply(strsplit(x=msbb_rnaseq_covariates.merged2$sampleIdentifier,split = "_"),`[[`,3))%in%colnames(msbb_rnaseq2016_byRegion$PHG)),c(1:11,13,15:16)]
msbb_rnaseq_covariates.merged_final$STG=msbb_rnaseq_covariates.merged2[which(unlist(lapply(strsplit(x=msbb_rnaseq_covariates.merged2$sampleIdentifier,split = "_"),`[[`,3))%in%colnames(msbb_rnaseq2016_byRegion$STG)),c(1:11,13,15:16)]

msbb_rnaseq2016_byRegion$FP=msbb_rnaseq2016_data2.agg[,which(colnames(msbb_rnaseq2016_data2.agg)%in%unlist(lapply(strsplit(x = msbb_rnaseq_covariates.merged2$sampleIdentifier[which(msbb_rnaseq_covariates.merged2$BrodmannArea=="BM10")],split = "_"),`[[`,3)))]
msbb_rnaseq2016_byRegion$IFG=msbb_rnaseq2016_data2.agg[,which(colnames(msbb_rnaseq2016_data2.agg)%in%unlist(lapply(strsplit(x = msbb_rnaseq_covariates.merged2$sampleIdentifier[which(msbb_rnaseq_covariates.merged2$BrodmannArea=="BM44")],split = "_"),`[[`,3)))]
msbb_rnaseq2016_byRegion$PHG=msbb_rnaseq2016_data2.agg[,which(colnames(msbb_rnaseq2016_data2.agg)%in%unlist(lapply(strsplit(x = msbb_rnaseq_covariates.merged2$sampleIdentifier[which(msbb_rnaseq_covariates.merged2$BrodmannArea=="BM36")],split = "_"),`[[`,3)))]
msbb_rnaseq2016_byRegion$STG=msbb_rnaseq2016_data2.agg[,which(colnames(msbb_rnaseq2016_data2.agg)%in%unlist(lapply(strsplit(x = msbb_rnaseq_covariates.merged2$sampleIdentifier[which(msbb_rnaseq_covariates.merged2$BrodmannArea=="BM22")],split = "_"),`[[`,3)))]

lowPlaque_samples=highPlaque_samples=vector(mode = "list",length = 4)
names(lowPlaque_samples)=names(highPlaque_samples)=names(msbb_rnaseq2016_byRegion)
lowPlaque_samples=lapply(msbb_rnaseq_covariates.merged_final,function(x)x$sampleIdentifier[which(x$PlaqueMean<=1)])
highPlaque_samples=lapply(msbb_rnaseq_covariates.merged_final,function(x)x$sampleIdentifier[which(x$PlaqueMean>=15)])

i=0
blocksize=5
blocks=200#length(rownames(msbb_rnaseq2016_byRegion[[1]]))/5
start<-i*blocksize+1
end<-min((i+1)*blocksize, blocks)

indices=data.frame(cbind(seq(1,16556,by = 55)[seq(1,301,by=2)],seq(1,16556,by = 55)[seq(0,302,by=2)]),stringsAsFactors = F)
colnames(indices)=c("Start","End")
indices[151,]=c(16501,16505)

for (n in 1:dim(indices)[1]){
  for (r in 1:length(msbb_rnaseq2016_byRegion)){
    cat(paste("Processing brain region ",names(msbb_rnaseq2016_byRegion)[r]," ...\n",sep=""))
    tmp2=foreach(i=indices$Start[n]:indices$End[n],.combine = rbind)%dopar%{
      cat(paste("Correlating for block ", indices$Start[n]," ...\n",sep=""))
      pb = txtProgressBar(min=0,max=length(indices$Start[n]:indices$End[n]),style=3,initial=0)
      cat("\n")
      x1=unlist(unname(msbb_rnaseq2016_byRegion[[r]][i,which(colnames(msbb_rnaseq2016_byRegion[[r]])%in%unlist(lapply(strsplit(x = msbb_rnaseq_covariates.merged_final[[r]]$sampleIdentifier,split = "_"),`[[`,3)))]))
      y1=rep(msbb_rnaseq_covariates.merged_final[[r]]$PlaqueMean,length(i))
      rho=cor.test(x = x1,y = y1,method = "spearman")
      rho.sp=unname(rho$estimate)
      rho.p=rho$p.value
      gene=rownames(msbb_rnaseq2016_byRegion[[r]])[i]
      padj=p.adjust(rho.p,method = "fdr")
      tmp1=c(gene,rho.sp,rho.p,padj)
      close(pb)
    }
    res=data.frame(tmp2[which(tmp2[,3]<=0.05),],stringsAsFactors = F)
    msbb_rnaseq2016_PLQGenes[[r]]=data.frame(Genes=res$X1,Rho.Spearman=round(as.numeric(res$X2),digits = 5),Pval=round(as.numeric(res$X3),digits = 5),FDR=round(as.numeric(res$X4),digits = 5),stringsAsFactors = F)
    
  }
}





# i<-0
# start<-i*blocksize+1
# end<-min((i+1)*blocksize, number_of_combinations)
# GenePlaqueCorr <- function(){
#   
#   x1=geneCounts
#   y1=plaque
#   rho=cor.test(x = x1,y = y1,method = "spearman")
#   rho.p=rho$p.value
#   genes=rownames(geneCounts)
#   tmp1[j,1:3]=cbind(genes,rho.sp,rho.p)
#   tmp2=data.frame(Genes=tmp1[,1],Rho.Spearman=as.numeric(tmp1[,2]),Pval=as.numeric(tmp1[,3]),stringsAsFactors = F)
#   return(tmp2)
#   
# }

