library(org.Hs.eg.db)
library(dplyr)
library(magrittr)
library(data.table)
library(diffcoexp)
library(GEOquery)

# 
# 
# 
# msbb_gse84422_series_matrix=getGEO(GEO = "GSE84422",GSEMatrix = T)
# msbb_gse84422_series_matrix.GPL96=msbb_gse84422_series_matrix$`GSE84422-GPL96_series_matrix.txt.gz`
# msbb_gse84422_series_matrix.GPL97=msbb_gse84422_series_matrix$`GSE84422-GPL97_series_matrix.txt.gz`
# msbb_gse84422_series_matrix.GPL570=msbb_gse84422_series_matrix$`GSE84422-GPL570_series_matrix.txt.gz`
# 
# msbb_gse84422.fData=vector(mode = "list",length = 3)
# names(msbb_gse84422.fData)=c("GPL96","GPL97","GPL570")
# msbb_gse84422.fData$GPL96=fData(msbb_gse84422_series_matrix.GPL96)[,c(1:2,11:13)]
# msbb_gse84422.fData$GPL97=fData(msbb_gse84422_series_matrix.GPL97)[,c(1:2,11:13)]
# msbb_gse84422.fData$GPL570=fData(msbb_gse84422_series_matrix.GPL570)[,c(1:2,11:13)]
# 
# 
# msbb_gse84422.pData=vector(mode = "list",length = 3)
# names(msbb_gse84422.pData)=c("GPL96","GPL97","GPL570")
# msbb_gse84422.pData$GPL96=pData(phenoData(msbb_gse84422_series_matrix.GPL96))
# msbb_gse84422.pData$GPL97=pData(phenoData(msbb_gse84422_series_matrix.GPL97))
# msbb_gse84422.pData$GPL570=pData(phenoData(msbb_gse84422_series_matrix.GPL570))
# msbb_gse84422.pData$GPL96$pseudoSampleID=msbb_gse84422.pData$GPL97$pseudoSampleID=paste("pSample",1:dim(msbb_gse84422.pData$GPL96)[1],sep = "")
# msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$`brain region:ch1`=gsub(pattern = "^ ",replacement = "",x = x$`brain region:ch1`);x})
# msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$`clinical dementia rating:ch1`=as.numeric(gsub(pattern = "^ ",replacement = "",x = x$`clinical dementia rating:ch1`));x})
# msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$`braak neurofibrillary tangle score:ch1`=as.numeric(gsub(pattern = "^ ",replacement = "",x = x$`braak neurofibrillary tangle score:ch1`));x})
# msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$`average neuritic plaque density:ch1`=as.numeric(gsub(pattern = "^ ",replacement = "",x = x$`average neuritic plaque density:ch1`));x})
# msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$`neuropathological category:ch1`=gsub(pattern = "^ ",replacement = "",x = x$`neuropathological category:ch1`);x})
# msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$`SampleType`="OTHER";x})
# msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$SampleType[which(x$`neuropathological category:ch1`=="Normal"&x$`clinical dementia rating:ch`<=0.5&x$`braak neurofibrillary tangle score:ch1`<=3)]="CONTROL";x})
# msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$SampleType[which((x$`neuropathological category:ch1`=="Probable AD"|x$`neuropathological category:ch1`=="Possible AD"|x$`neuropathological category:ch1`=="definite AD")&x$`clinical dementia rating:ch`>=1&x$`braak neurofibrillary tangle score:ch1`>=4)]="AD";x})
# msbb_gse84422_GPL96_97_samples.Control=msbb_gse84422.pData$GPL96$pseudoSampleID[which(msbb_gse84422.pData$GPL96$`neuropathological category:ch1`=="Normal"&msbb_gse84422.pData$GPL96$`clinical dementia rating:ch`<=0.5&msbb_gse84422.pData$GPL96$`braak neurofibrillary tangle score:ch1`<=3)]
# msbb_gse84422_GPL96_97_samples.AD=msbb_gse84422.pData$GPL96$pseudoSampleID[which(msbb_gse84422.pData$GPL96$`neuropathological category:ch1`=="Normal"&msbb_gse84422.pData$GPL96$`clinical dementia rating:ch`<=0.5&msbb_gse84422.pData$GPL96$`braak neurofibrillary tangle score:ch1`<=3)]
# 
# 
# 
# msbb_gse84422.exprs=vector(mode = "list",length = 3)
# names(msbb_gse84422.exprs)=c("GPL96","GPL97","GPL570")
# msbb_gse84422.exprs$GPL96=exprs(msbb_gse84422_series_matrix.GPL96)
# msbb_gse84422.exprs$GPL97=exprs(msbb_gse84422_series_matrix.GPL97)
# msbb_gse84422.exprs$GPL570=data.frame(exprs(msbb_gse84422_series_matrix.GPL570))
# 
# colnames(msbb_gse84422.exprs$GPL96)=colnames(msbb_gse84422.exprs$GPL97)=msbb_gse84422.pData$GPL96$pseudoSampleID
# 
# msbb_gse84422_exprs.GPL96_97=rbind.data.frame(msbb_gse84422.exprs$GPL96,msbb_gse84422.exprs$GPL97)
# msbb_gse84422_exprs.GPL96_97$GeneSymbol=c(msbb_gse84422.fData$GPL96$`Gene Symbol`,msbb_gse84422.fData$GPL97$`Gene Symbol`)
# msbb_gse84422_exprs_GPL96_97.agg=aggregate.data.frame(x=msbb_gse84422_exprs.GPL96_97[,-which(colnames(msbb_gse84422_exprs.GPL96_97)=="GeneSymbol")],by=list(symbol=msbb_gse84422_exprs.GPL96_97$GeneSymbol),mean)
# rownames(msbb_gse84422_exprs_GPL96_97.agg)=msbb_gse84422_exprs_GPL96_97.agg$symbol
# msbb_gse84422_exprs_GPL96_97.agg=msbb_gse84422_exprs_GPL96_97.agg[,-which(colnames(msbb_gse84422_exprs_GPL96_97.agg)=="symbol")]
# 
# msbb_gse84422.exprs$GPL570$GeneSymbol=msbb_gse84422.fData$GPL570$`Gene Symbol`
# msbb_gse84422.exprs$GPL570=aggregate.data.frame(x=msbb_gse84422.exprs$GPL570[,-which(colnames(msbb_gse84422.exprs$GPL570)=="GeneSymbol")],by=list(symbol=msbb_gse84422.exprs$GPL570$GeneSymbol),mean)
# rownames(msbb_gse84422.exprs$GPL570)=msbb_gse84422.exprs$GPL570$symbol
# msbb_gse84422.exprs$GPL570=msbb_gse84422.exprs$GPL570[,-which(colnames(msbb_gse84422.exprs$GPL570)=="symbol")]
# 
# msbb_gse84422_GPL96_97_byRegion.exprs=vector(mode = "length",length = 17)
# names(msbb_gse84422_GPL96_97_byRegion.exprs)=unique(msbb_gse84422.pData$GPL96$`brain region:ch1`)
# msbb_gse84422_GPL96_97_samplesToAnalyse=lapply(unique(msbb_gse84422.pData$GPL96$`brain region:ch1`),function(y)msbb_gse84422.pData$GPL96%>%filter(`brain region:ch1`==y&SampleType!="OTHER")%>%select(c(SampleType,pseudoSampleID)))
# names(msbb_gse84422_GPL96_97_samplesToAnalyse)=unique(msbb_gse84422.pData$GPL96$`brain region:ch1`)
# msbb_gse84422_GPL570_samplesToAnalyse=lapply(unique(msbb_gse84422.pData$GPL570$`brain region:ch1`),function(y)msbb_gse84422.pData$GPL570%>%filter(`brain region:ch1`==y&SampleType!="OTHER")%>%select(c(SampleType,geo_accession)))
# names(msbb_gse84422_GPL570_samplesToAnalyse)=unique(msbb_gse84422.pData$GPL570$`brain region:ch1`)
# 
# 
# msbb_gse84422_GPL96_97_samplesToAnalyse.exprs=lapply(msbb_gse84422_GPL96_97_samplesToAnalyse,function(y)msbb_gse84422_exprs_GPL96_97.agg[,colnames(msbb_gse84422_exprs_GPL96_97.agg)%in%y$pseudoSampleID])
# msbb_gse84422_GPL570_samplesToAnalyse.exprs=lapply(msbb_gse84422_GPL570_samplesToAnalyse, function(y)msbb_gse84422.exprs$GPL570[,colnames(msbb_gse84422.exprs$GPL570)%in%y$geo_accession])
# 
# saveRDS(msbb_gse84422_GPL570_samplesToAnalyse,"msbb_gse84422_GPL570_samplesToAnalyse.RDS")
# saveRDS(msbb_gse84422_GPL570_samplesToAnalyse.exprs,"msbb_gse84422_GPL570_samplesToAnalyse_exprs.RDS")
# saveRDS(msbb_gse84422_GPL96_97_samplesToAnalyse,"msbb_gse84422_GPL96_97_samplesToAnalyse.RDS")
# saveRDS(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs,"msbb_gse84422_GPL96_97_samplesToAnalyse_exprs.RDS")
setwd("/shared/hidelab2/user/md4zsa/Work/Data/MSBB_Array19/GSE84422")
msbb_gse84422_GPL96_97_samplesToAnalyse.exprs=readRDS("msbb_gse84422_GPL96_97_samplesToAnalyse_exprs.RDS")
msbb_gse84422_GPL96_97_samplesToAnalyse=readRDS("msbb_gse84422_GPL96_97_samplesToAnalyse.RDS")
msbb_gse84422_GPL570_samplesToAnalyse.exprs=readRDS("msbb_gse84422_GPL570_samplesToAnalyse_exprs.RDS")
msbb_gse84422_GPL570_samplesToAnalyse=readRDS("msbb_gse84422_GPL570_samplesToAnalyse.RDS")

#Diffcoexp analysis
for(i in 1:2){
  c_exprs=msbb_gse84422_GPL570_samplesToAnalyse.exprs[[i]][,msbb_gse84422_GPL570_samplesToAnalyse[[i]]$SampleType=="CONTROL"]
  d_exprs=msbb_gse84422_GPL570_samplesToAnalyse.exprs[[i]][,msbb_gse84422_GPL570_samplesToAnalyse[[i]]$SampleType=="AD"]
  diffcoexp_out=diffcoexp(exprs.1 = c_exprs[1:15000,],exprs.2 = d_exprs[1:15000,],r.method="spearman",rth=0.6,q.diffth=0.1,q.dcgth=0.1)
  saveRDS(diffcoexp_out,paste(names(msbb_gse84422_GPL570_samplesToAnalyse.exprs)[i],"diffcoexp_results.RDS",sep = "_"))
}

for(i in 1:17){
  c_exprs=msbb_gse84422_GPL96_97_samplesToAnalyse.exprs[[i]][,msbb_gse84422_GPL96_97_samplesToAnalyse[[i]]$SampleType=="CONTROL"]
  d_exprs=msbb_gse84422_GPL96_97_samplesToAnalyse.exprs[[i]][,msbb_gse84422_GPL96_97_samplesToAnalyse[[i]]$SampleType=="AD"]
  diffcoexp_out=diffcoexp(exprs.1 = c_exprs,exprs.2 = d_exprs,r.method="spearman",rth=0.6,q.diffth=0.1,q.dcgth=0.1)
  saveRDS(diffcoexp_out,paste(names(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs)[i],"diffcoexp_results.RDS",sep = "_"))
}
saveRDS(c_exprs,"CONTROL_exprs.RDS")
saveRDS(d_exprs,"AD_exprs.RDS")
warnings