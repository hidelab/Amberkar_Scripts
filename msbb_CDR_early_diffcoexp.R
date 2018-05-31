library(org.Hs.eg.db)
library(dplyr)
library(magrittr)
library(data.table)
library(diffcoexp)
library(GEOquery)

setwd("/shared/hidelab2/user/md4zsa/Work/Data/MSBB_Array19/GSE84422/")
load("msbb_gse84422_Lobe_DEG_DEP.RData")

msbb_gse84422_series_matrix=readRDS("msbb_gse84422_series_matrix.RDS")


msbb_gse84422_series_matrix.GPL570=msbb_gse84422_series_matrix$`GSE84422-GPL570_series_matrix.txt.gz`
msbb_gse84422.pData$GPL570=pData(phenoData(msbb_gse84422_series_matrix.GPL570))

msbb_gse84422.pData=vector(mode = "list",length = 3)
names(msbb_gse84422.pData)=c("GPL96","GPL97","GPL570")
msbb_gse84422.pData$GPL96=pData(phenoData(msbb_gse84422_series_matrix$`GSE84422-GPL96_series_matrix.txt.gz`))
msbb_gse84422.pData$GPL97=pData(phenoData(msbb_gse84422_series_matrix$`GSE84422-GPL97_series_matrix.txt.gz`))
msbb_gse84422.pData$GPL570=pData(phenoData(msbb_gse84422_series_matrix$`GSE84422-GPL570_series_matrix.txt.gz`))


msbb_gse84422.fData=vector(mode = "list",length = 3)
names(msbb_gse84422.fData)=c("GPL96","GPL97","GPL570")
msbb_gse84422.fData$GPL96=fData(msbb_gse84422_series_matrix$`GSE84422-GPL96_series_matrix.txt.gz`)[,c(1:2,11:13)]
msbb_gse84422.fData$GPL97=fData(msbb_gse84422_series_matrix$`GSE84422-GPL97_series_matrix.txt.gz`)[,c(1:2,11:13)]
msbb_gse84422.fData$GPL570=fData(msbb_gse84422_series_matrix$`GSE84422-GPL570_series_matrix.txt.gz`)[,c(1:2,11:13)]

msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$`brain region:ch1`=gsub(pattern = "^ ",replacement = "",x = x$`brain region:ch1`);x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$`clinical dementia rating:ch1`=as.numeric(gsub(pattern = "^ ",replacement = "",x = x$`clinical dementia rating:ch1`));x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$`braak neurofibrillary tangle score:ch1`=as.numeric(gsub(pattern = "^ ",replacement = "",x = x$`braak neurofibrillary tangle score:ch1`));x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$`average neuritic plaque density:ch1`=as.numeric(gsub(pattern = "^ ",replacement = "",x = x$`average neuritic plaque density:ch1`));x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$`neuropathological category:ch1`=gsub(pattern = "^ ",replacement = "",x = x$`neuropathological category:ch1`);x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$`SampleTypeCDR`="OTHER";x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$SampleTypeCDR[which(x$`clinical dementia rating:ch`==0)]="CDR0";x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$SampleTypeCDR[which(x$`clinical dementia rating:ch`==0.5)]="CDR05";x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$SampleTypeCDR[which(x$`clinical dementia rating:ch`==1)]="CDR1";x})

msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$`SampleTypeCDR`="OTHER";x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$SampleTypeCDR[which(x$`clinical dementia rating:ch`==0)]="CDR0";x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$SampleTypeCDR[which(x$`clinical dementia rating:ch`==0.5)]="CDR05";x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$SampleTypeCDR[which(x$`clinical dementia rating:ch`==1)]="CDR1";x})

msbb_gse84422_exprs.GPL570=data.frame(exprs(msbb_gse84422_series_matrix.GPL570))
msbb_gse84422_exprs.GPL570$GeneSymbol=msbb_gse84422.fData$GPL570$`Gene Symbol`
msbb_gse84422_exprs.GPL570$ENTREZ_GENE_ID=msbb_gse84422.fData$GPL570$ENTREZ_GENE_ID
msbb_gse84422_exprs_GPL570.agg=aggregate.data.frame(x=msbb_gse84422_exprs.GPL570[,-which(colnames(msbb_gse84422_exprs.GPL570)=="GeneSymbol"|colnames(msbb_gse84422_exprs.GPL570)=="ENTREZ_GENE_ID")],by=list(symbol=msbb_gse84422_exprs.GPL570$GeneSymbol),mean)
msbb_gse84422_exprs_GPL570.agg_Entrez=aggregate.data.frame(x=msbb_gse84422_exprs.GPL570[,-which(colnames(msbb_gse84422_exprs.GPL570)=="GeneSymbol"|colnames(msbb_gse84422_exprs.GPL570)=="ENTREZ_GENE_ID")],by=list(entrez=msbb_gse84422_exprs.GPL570$ENTREZ_GENE_ID),mean)
rownames(msbb_gse84422_exprs_GPL570.agg)=msbb_gse84422_exprs_GPL570.agg$symbol
rownames(msbb_gse84422_exprs_GPL570.agg_Entrez)=msbb_gse84422_exprs_GPL570.agg_Entrez$entrez
msbb_gse84422_exprs_GPL570.agg=msbb_gse84422_exprs_GPL570.agg[,-which(colnames(msbb_gse84422_exprs_GPL570.agg)=="symbol")]
msbb_gse84422_GPL570_byRegion.exprs=msbb_gse84422_GPL570_byRegion.exprs_Entrez=vector(mode = "list",length = 2)
names(msbb_gse84422_GPL570_byRegion.exprs)=names(msbb_gse84422_GPL570_byRegion.exprs_Entrez)=unique(msbb_gse84422.pData$GPL570$`brain region:ch1`)
msbb_gse84422_GPL570_samplesByRegion=lapply(unique(msbb_gse84422.pData$GPL570$`brain region:ch1`),function(y)msbb_gse84422.pData$GPL570%>%filter(`brain region:ch1`==y)%>%pull(geo_accession))
msbb_gse84422_GPL570_byRegion.exprs=lapply(msbb_gse84422_GPL570_samplesByRegion, function(x)msbb_gse84422_exprs_GPL570.agg[,x])
names(msbb_gse84422_GPL570_byRegion.exprs)=unique(msbb_gse84422.pData$GPL570$`brain region:ch1`)

for(i in 1:length(msbb_gse84422_GPL96_97_byRegion.exprs)){
  c_exprs=msbb_gse84422_GPL96_97_byRegion.exprs[[i]][1:1000,msbb_gse84422.pData$GPL96$pseudoSampleID[which(msbb_gse84422.pData$GPL96$`brain region:ch1`==names(msbb_gse84422_GPL96_97_byRegion.exprs)[i]&msbb_gse84422.pData$GPL96$SampleTypeCDR=="CDR0")]]
  d_exprs=msbb_gse84422_GPL96_97_byRegion.exprs[[i]][1:1000,msbb_gse84422.pData$GPL96$pseudoSampleID[which(msbb_gse84422.pData$GPL96$`brain region:ch1`==names(msbb_gse84422_GPL96_97_byRegion.exprs)[i]&msbb_gse84422.pData$GPL96$SampleTypeCDR=="CDR1")]]
  diffcoexp_out=diffcoexp(exprs.1 = c_exprs[,],exprs.2 = d_exprs[,],r.method="spearman",rth=0.6,q.diffth=0.1,q.dcgth=0.1)
  fn_prefix=gsub(pattern = " ",replacement = "_",names(msbb_gse84422_GPL96_97_byRegion.exprs)[i])
  saveRDS(diffcoexp_out,paste(fn_prefix,"diffcoexp_CDR0_CDR1.RDS",sep="_"))
  proc.time()
}

for(i in 1:length(msbb_gse84422_GPL570_byRegion.exprs)){
  c_exprs=msbb_gse84422_GPL570_byRegion.exprs[[i]][1:15000,msbb_gse84422.pData$GPL570$geo_accession[which(msbb_gse84422.pData$GPL570$`brain region:ch1`==names(msbb_gse84422_GPL570_byRegion.exprs)[i]&msbb_gse84422.pData$GPL570$SampleTypeCDR=="CDR0")]]
  d_exprs=msbb_gse84422_GPL570_byRegion.exprs[[i]][1:15000,msbb_gse84422.pData$GPL570$geo_accession[which(msbb_gse84422.pData$GPL570$`brain region:ch1`==names(msbb_gse84422_GPL570_byRegion.exprs)[i]&msbb_gse84422.pData$GPL570$SampleTypeCDR=="CDR1")]]
  diffcoexp_out=diffcoexp(exprs.1 = c_exprs[,],exprs.2 = d_exprs[,],r.method="spearman",rth=0.6,q.diffth=0.1,q.dcgth=0.1)
  fn_prefix=gsub(pattern = " ",replacement = "_",names(msbb_gse84422_GPL570_byRegion.exprs)[i])
  saveRDS(diffcoexp_out,paste(fn_prefix,"diffcoexp_CDR0_CDR1.RDS",sep="_"))
  proc.time()
}
