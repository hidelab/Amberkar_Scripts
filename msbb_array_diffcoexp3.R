library(org.Hs.eg.db)
library(dplyr)
library(magrittr)
library(data.table)
library(diffcoexp)
library(GEOquery)

# 
# 
# 
msbb_gse84422_series_matrix=getGEO(GEO = "GSE84422",GSEMatrix = T)
msbb_gse84422_series_matrix.GPL96=msbb_gse84422_series_matrix$`GSE84422-GPL96_series_matrix.txt.gz`
msbb_gse84422_series_matrix.GPL97=msbb_gse84422_series_matrix$`GSE84422-GPL97_series_matrix.txt.gz`
msbb_gse84422_series_matrix.GPL570=msbb_gse84422_series_matrix$`GSE84422-GPL570_series_matrix.txt.gz`
# 
msbb_gse84422.fData=vector(mode = "list",length = 3)
names(msbb_gse84422.fData)=c("GPL96","GPL97","GPL570")
msbb_gse84422.fData$GPL96=fData(msbb_gse84422_series_matrix.GPL96)[,c(1:2,11:13)]
msbb_gse84422.fData$GPL97=fData(msbb_gse84422_series_matrix.GPL97)[,c(1:2,11:13)]
msbb_gse84422.fData$GPL570=fData(msbb_gse84422_series_matrix.GPL570)[,c(1:2,11:13)]
# 
# 
msbb_gse84422.pData=vector(mode = "list",length = 3)
names(msbb_gse84422.pData)=c("GPL96","GPL97","GPL570")
msbb_gse84422.pData$GPL96=pData(phenoData(msbb_gse84422_series_matrix.GPL96))
msbb_gse84422.pData$GPL97=pData(phenoData(msbb_gse84422_series_matrix.GPL97))
msbb_gse84422.pData$GPL570=pData(phenoData(msbb_gse84422_series_matrix.GPL570))
msbb_gse84422.pData$GPL96$pseudoSampleID=msbb_gse84422.pData$GPL97$pseudoSampleID=paste("pSample",1:dim(msbb_gse84422.pData$GPL96)[1],sep = "")
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$`brain region:ch1`=gsub(pattern = "^ ",replacement = "",x = x$`brain region:ch1`);x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$`clinical dementia rating:ch1`=as.numeric(gsub(pattern = "^ ",replacement = "",x = x$`clinical dementia rating:ch1`));x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$`braak neurofibrillary tangle score:ch1`=as.numeric(gsub(pattern = "^ ",replacement = "",x = x$`braak neurofibrillary tangle score:ch1`));x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$`average neuritic plaque density:ch1`=as.numeric(gsub(pattern = "^ ",replacement = "",x = x$`average neuritic plaque density:ch1`));x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$`neuropathological category:ch1`=gsub(pattern = "^ ",replacement = "",x = x$`neuropathological category:ch1`);x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$`SampleType`="OTHER";x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$`SampleTypeCDR`="OTHER";x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$`SampleTypeBraak`="OTHER";x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$`SampleTypeCERAD`="OTHER";x})

msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$SampleTypeCERAD[which(x$`neuropathological category:ch1`=="Normal")]="Normal";x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$SampleTypeCERAD[which(x$`neuropathological category:ch1`=="Probable AD")]="Probable_AD";x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$SampleTypeCERAD[which(x$`neuropathological category:ch1`=="Possible AD")]="Possible_AD";x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$SampleTypeCERAD[which(x$`neuropathological category:ch1`=="definite AD")]="Definite_AD";x})

msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$SampleTypeBraak[which(x$`braak neurofibrillary tangle score:ch1`==0)]="Braak0";x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$SampleTypeBraak[which(x$`braak neurofibrillary tangle score:ch1`>=1&x$`braak neurofibrillary tangle score:ch1`<=2)]="Braak12";x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$SampleTypeBraak[which(x$`braak neurofibrillary tangle score:ch1`>=3&x$`braak neurofibrillary tangle score:ch1`<=4)]="Braak34";x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$SampleTypeBraak[which(x$`braak neurofibrillary tangle score:ch1`>=5)]="Braak56";x})

msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$SampleTypeCDR[which(x$`clinical dementia rating:ch`==0)]="CDR0";x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$SampleTypeCDR[which(x$`clinical dementia rating:ch`==0.5)]="CDR05";x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$SampleTypeCDR[which(x$`clinical dementia rating:ch`==1)]="CDR1";x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$SampleType[which(x$`neuropathological category:ch1`=="Normal"&x$`clinical dementia rating:ch`<=0.5&x$`braak neurofibrillary tangle score:ch1`<=3)]="CONTROL";x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$SampleType[which((x$`neuropathological category:ch1`=="Probable AD"|x$`neuropathological category:ch1`=="Possible AD"|x$`neuropathological category:ch1`=="definite AD")&x$`clinical dementia rating:ch`>=1&x$`braak neurofibrillary tangle score:ch1`>=4)]="AD";x})
msbb_gse84422_GPL96_97_samples.Control=msbb_gse84422.pData$GPL96$pseudoSampleID[which(msbb_gse84422.pData$GPL96$`neuropathological category:ch1`=="Normal"&msbb_gse84422.pData$GPL96$`clinical dementia rating:ch`<=0.5&msbb_gse84422.pData$GPL96$`braak neurofibrillary tangle score:ch1`<=3)]
msbb_gse84422_GPL96_97_samples.AD=msbb_gse84422.pData$GPL96$pseudoSampleID[which(msbb_gse84422.pData$GPL96$`neuropathological category:ch1`=="Normal"&msbb_gse84422.pData$GPL96$`clinical dementia rating:ch`<=0.5&msbb_gse84422.pData$GPL96$`braak neurofibrillary tangle score:ch1`<=3)]
msbb_gse84422_GPL96_97_samples.CDR0=msbb_gse84422.pData$GPL96$pseudoSampleID[which(msbb_gse84422.pData$GPL96$`clinical dementia rating:ch`==0)]
msbb_gse84422_GPL96_97_samples.CDR05=msbb_gse84422.pData$GPL96$pseudoSampleID[which(msbb_gse84422.pData$GPL96$`clinical dementia rating:ch`==0.5)]

# 
# 
# 
msbb_gse84422.exprs=vector(mode = "list",length = 3)
names(msbb_gse84422.exprs)=c("GPL96","GPL97","GPL570")
msbb_gse84422.exprs$GPL96=exprs(msbb_gse84422_series_matrix.GPL96)
msbb_gse84422.exprs$GPL97=exprs(msbb_gse84422_series_matrix.GPL97)
msbb_gse84422.exprs$GPL570=data.frame(exprs(msbb_gse84422_series_matrix.GPL570))
# 
colnames(msbb_gse84422.exprs$GPL96)=colnames(msbb_gse84422.exprs$GPL97)=msbb_gse84422.pData$GPL96$pseudoSampleID
# 
msbb_gse84422_exprs.GPL96_97=rbind.data.frame(msbb_gse84422.exprs$GPL96,msbb_gse84422.exprs$GPL97)
msbb_gse84422_exprs.GPL96_97$GeneSymbol=c(msbb_gse84422.fData$GPL96$`Gene Symbol`,msbb_gse84422.fData$GPL97$`Gene Symbol`)
msbb_gse84422_exprs_GPL96_97.agg=aggregate.data.frame(x=msbb_gse84422_exprs.GPL96_97[,-which(colnames(msbb_gse84422_exprs.GPL96_97)=="GeneSymbol")],by=list(symbol=msbb_gse84422_exprs.GPL96_97$GeneSymbol),mean)
rownames(msbb_gse84422_exprs_GPL96_97.agg)=msbb_gse84422_exprs_GPL96_97.agg$symbol
msbb_gse84422_exprs_GPL96_97.agg=msbb_gse84422_exprs_GPL96_97.agg[,-which(colnames(msbb_gse84422_exprs_GPL96_97.agg)=="symbol")]
# 
msbb_gse84422.exprs$GPL570$GeneSymbol=msbb_gse84422.fData$GPL570$`Gene Symbol`
msbb_gse84422.exprs$GPL570=aggregate.data.frame(x=msbb_gse84422.exprs$GPL570[,-which(colnames(msbb_gse84422.exprs$GPL570)=="GeneSymbol")],by=list(symbol=msbb_gse84422.exprs$GPL570$GeneSymbol),mean)
rownames(msbb_gse84422.exprs$GPL570)=msbb_gse84422.exprs$GPL570$symbol
msbb_gse84422.exprs$GPL570=msbb_gse84422.exprs$GPL570[,-which(colnames(msbb_gse84422.exprs$GPL570)=="symbol")]
# 
msbb_gse84422_GPL96_97_byRegion.exprs=vector(mode = "list",length = 17)
names(msbb_gse84422_GPL96_97_byRegion.exprs)=unique(msbb_gse84422.pData$GPL96$`brain region:ch1`)
msbb_gse84422_GPL96_97_samplesToAnalyse=lapply(unique(msbb_gse84422.pData$GPL96$`brain region:ch1`),function(y)msbb_gse84422.pData$GPL96%>%filter(`brain region:ch1`==y&SampleType!="OTHER")%>%select(c(SampleType,pseudoSampleID)))
msbb_gse84422_GPL96_97_samplesToAnalyse.CDR=lapply(unique(msbb_gse84422.pData$GPL96$`brain region:ch1`),function(y)msbb_gse84422.pData$GPL96%>%filter(`brain region:ch1`==y&SampleTypeCDR!="OTHER")%>%select(c(SampleTypeCDR,pseudoSampleID)))
msbb_gse84422_GPL96_97_samplesToAnalyse.Braak=lapply(unique(msbb_gse84422.pData$GPL96$`brain region:ch1`),function(y)msbb_gse84422.pData$GPL96%>%filter(`brain region:ch1`==y&SampleTypeBraak!="OTHER")%>%select(c(SampleTypeBraak,pseudoSampleID)))
msbb_gse84422_GPL96_97_samplesToAnalyse.CERAD=lapply(unique(msbb_gse84422.pData$GPL96$`brain region:ch1`),function(y)msbb_gse84422.pData$GPL96%>%filter(`brain region:ch1`==y&SampleTypeCERAD!="OTHER")%>%select(c(SampleTypeCERAD,pseudoSampleID)))
names(msbb_gse84422_GPL96_97_samplesToAnalyse)=names(msbb_gse84422_GPL96_97_samplesToAnalyse.CERAD)=names(msbb_gse84422_GPL96_97_samplesToAnalyse.CDR)=names(msbb_gse84422_GPL96_97_samplesToAnalyse.Braak)=names(msbb_gse84422_GPL96_97_samplesToAnalyse.CDR)=unique(msbb_gse84422.pData$GPL96$`brain region:ch1`)

msbb_gse84422_GPL570_samplesToAnalyse=lapply(unique(msbb_gse84422.pData$GPL570$`brain region:ch1`),function(y)msbb_gse84422.pData$GPL570%>%filter(`brain region:ch1`==y&SampleType!="OTHER")%>%select(c(SampleType,geo_accession)))
msbb_gse84422_GPL570_samplesToAnalyse.CDR=lapply(unique(msbb_gse84422.pData$GPL570$`brain region:ch1`),function(y)msbb_gse84422.pData$GPL570%>%filter(`brain region:ch1`==y&SampleTypeCDR!="OTHER")%>%select(c(SampleTypeCDR,geo_accession)))

names(msbb_gse84422_GPL570_samplesToAnalyse)=names(msbb_gse84422_GPL570_samplesToAnalyse.CDR)=unique(msbb_gse84422.pData$GPL570$`brain region:ch1`)
# 
# 


#Regroup brain regions by lobes
lobe_bm_area.map=matrix(nrow=17,ncol=2)
lobe_bm_area.map[1,]=c("Frontal_Lobe","Frontal Pole")
lobe_bm_area.map[2,]=c("Frontal_Lobe","Anterior Cingulate")
lobe_bm_area.map[3,]=c("Frontal_Lobe","Prefrontal Cortex")
lobe_bm_area.map[4,]=c("Occipetal_Lobe","Occipital Visual Cortex")
lobe_bm_area.map[5,]=c("Temporal_Lobe","Inferior Temporal Gyrus")
lobe_bm_area.map[6,]=c("Temporal_Lobe","Middle Temporal Gyrus")
lobe_bm_area.map[7,]=c("Temporal_Lobe","Superior Temporal Gyrus")
lobe_bm_area.map[8,]=c("Parietal_Lobe","Posterior Cingulate Cortex")
lobe_bm_area.map[9,]=c("Temporal_Lobe","Parahippocampal Gyrus")
lobe_bm_area.map[10,]=c("Temporal_Lobe","Temporal Pole")
lobe_bm_area.map[11,]=c("Frontal_Lobe","Precentral Gyrus")
lobe_bm_area.map[12,]=c("Frontal_Lobe","Inferior Frontal Gyrus")
lobe_bm_area.map[13,]=c("Frontal_Lobe","Dorsolateral Prefrontal Cortex")
lobe_bm_area.map[14,]=c("Parietal_Lobe","Superior Parietal Lobule")
lobe_bm_area.map[16,]=c("Temporal_Lobe","Hippocampus")
lobe_bm_area.map[15,]=c("Dorsal_striatum","Putamen")
lobe_bm_area.map[17,]=c("Dorsal_striatum","Caudate Nucleus")
lobe_bm_area.map=data.frame(Lobe=lobe_bm_area.map[,1],BrainRegion=lobe_bm_area.map[,2],stringsAsFactors = F)

msbb_gse84422.byLobe=vector(mode = "list",length = length(unique(lobe_bm_area.map$Lobe)))
names(msbb_gse84422.byLobe)=unique(lobe_bm_area.map$Lobe)
msbb_gse84422.byLobe=foreach(i=1:length(names(msbb_gse84422.byLobe)))%dopar%{
  array_byLobe=do.call("cbind",msbb_gse84422_exprs_GPL96_97.agg[which(names(msbb_gse84422_exprs_GPL96_97.agg)%in%lobe_bm_area.map$BrainRegion[lobe_bm_area.map$Lobe==names(msbb_gse84422.byLobe)[i]])])  
}

msbb_gse84422_GPL96_97_byRegion.exprs=vector(mode = "list",length = 17)
names(msbb_gse84422_GPL96_97_byRegion.exprs)=unique(msbb_gse84422.pData$GPL96$`brain region:ch1`)
msbb_gse84422_GPL96_97_samplesToAnalyse=lapply(unique(msbb_gse84422.pData$GPL96$`brain region:ch1`),function(y)msbb_gse84422.pData$GPL96%>%filter(`brain region:ch1`==y&SampleType!="OTHER")%>%select(c(SampleType,pseudoSampleID)))
names(msbb_gse84422_GPL96_97_samplesToAnalyse)=unique(msbb_gse84422.pData$GPL96$`brain region:ch1`)




msbb_gse84422_GPL96_97_samplesToAnalyse.exprs=lapply(msbb_gse84422_GPL96_97_samplesToAnalyse,function(y)msbb_gse84422_exprs_GPL96_97.agg[,colnames(msbb_gse84422_exprs_GPL96_97.agg)%in%y$pseudoSampleID])
msbb_gse84422_GPL570_samplesToAnalyse.exprs=lapply(msbb_gse84422_GPL570_samplesToAnalyse, function(y)msbb_gse84422.exprs$GPL570[,colnames(msbb_gse84422.exprs$GPL570)%in%y$geo_accession])
msbb_gse84422_GPL96_97_samplesToAnalyse_CDR.exprs=lapply(msbb_gse84422_GPL96_97_samplesToAnalyse.CDR,function(y)msbb_gse84422_exprs_GPL96_97.agg[,colnames(msbb_gse84422_exprs_GPL96_97.agg)%in%y$pseudoSampleID])
msbb_gse84422_GPL96_97_samplesToAnalyse_CERAD.exprs=lapply(msbb_gse84422_GPL96_97_samplesToAnalyse.CERAD,function(y)msbb_gse84422_exprs_GPL96_97.agg[,colnames(msbb_gse84422_exprs_GPL96_97.agg)%in%y$pseudoSampleID])
msbb_gse84422_GPL96_97_samplesToAnalyse_Braak.exprs=lapply(msbb_gse84422_GPL96_97_samplesToAnalyse.Braak,function(y)msbb_gse84422_exprs_GPL96_97.agg[,colnames(msbb_gse84422_exprs_GPL96_97.agg)%in%y$pseudoSampleID])
msbb_gse84422_GPL570_samplesToAnalyse_CDR.exprs=lapply(msbb_gse84422_GPL570_samplesToAnalyse.CDR, function(y)msbb_gse84422.exprs$GPL570[,colnames(msbb_gse84422.exprs$GPL570)%in%y$geo_accession])

names(msbb_gse84422_GPL96_97_samplesToAnalyse_CDR.exprs)=names(msbb_gse84422_GPL96_97_samplesToAnalyse_CERAD.exprs)=names(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs)
names(msbb_gse84422_GPL570_samplesToAnalyse_CDR.exprs)=names(msbb_gse84422_GPL570_samplesToAnalyse.exprs)
names(msbb_gse84422_GPL96_97_samplesToAnalyse_Braak.exprs)=names(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs)
saveRDS(msbb_gse84422_GPL570_samplesToAnalyse.CDR,"msbb_gse84422_GPL570_samplesToAnalyse_CDR.RDS")
saveRDS(msbb_gse84422_GPL96_97_samplesToAnalyse.CDR,"msbb_gse84422_GPL96_97_samplesToAnalyse_CDR.RDS")
#saveRDS(msbb_gse84422_GPL570_samplesToAnalyse,"msbb_gse84422_GPL570_samplesToAnalyse.RDS")
# saveRDS(msbb_gse84422_GPL570_samplesToAnalyse.exprs,"msbb_gse84422_GPL570_samplesToAnalyse_exprs.RDS")
saveRDS(msbb_gse84422_GPL570_samplesToAnalyse_CDR.exprs,"msbb_gse84422_GPL570_samplesToAnalyse_CDR_exprs.RDS")
# saveRDS(msbb_gse84422_GPL96_97_samplesToAnalyse,"msbb_gse84422_GPL96_97_samplesToAnalyse.RDS")
# saveRDS(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs,"msbb_gse84422_GPL96_97_samplesToAnalyse_exprs.RDS")
saveRDS(msbb_gse84422_GPL96_97_samplesToAnalyse_CDR.exprs,"msbb_gse84422_GPL96_97_samplesToAnalyse_CDR_exprs.RDS")
setwd("/shared/hidelab2/user/md4zsa/Work/Data/msbb_gse8442219/GSE84422")
msbb_gse84422_GPL96_97_samplesToAnalyse.exprs=readRDS("msbb_gse84422_GPL96_97_samplesToAnalyse_exprs.RDS")
msbb_gse84422_GPL96_97_samplesToAnalyse=readRDS("msbb_gse84422_GPL96_97_samplesToAnalyse.RDS")
msbb_gse84422_GPL570_samplesToAnalyse.exprs=readRDS("msbb_gse84422_GPL570_samplesToAnalyse_exprs.RDS")
msbb_gse84422_GPL570_samplesToAnalyse=readRDS("msbb_gse84422_GPL570_samplesToAnalyse.RDS")

#Diffcoexp analysis
for(i in 1:17){
  c_exprs=msbb_gse84422_GPL96_97_samplesToAnalyse.exprs[[i]][,msbb_gse84422_GPL96_97_samplesToAnalyse[[i]]$SampleType=="CONTROL"]
  d_exprs=msbb_gse84422_GPL96_97_samplesToAnalyse.exprs[[i]][,msbb_gse84422_GPL96_97_samplesToAnalyse[[i]]$SampleType=="AD"]
  diffcoexp_out=diffcoexp(exprs.1 = c_exprs,exprs.2 = d_exprs,r.method="spearman",rth=0.6,q.diffth=0.1,q.dcgth=0.1)
  saveRDS(diffcoexp_out,paste(names(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs)[i],"diffcoexp_results.RDS",sep = "_"))
  proc.time()
}

for(i in 1:2){
  c_exprs=msbb_gse84422_GPL570_samplesToAnalyse.exprs[[i]][,msbb_gse84422_GPL570_samplesToAnalyse[[i]]$SampleType=="CONTROL"]
  d_exprs=msbb_gse84422_GPL570_samplesToAnalyse.exprs[[i]][,msbb_gse84422_GPL570_samplesToAnalyse[[i]]$SampleType=="AD"]
  diffcoexp_out=diffcoexp(exprs.1 = c_exprs[,],exprs.2 = d_exprs[,],r.method="spearman",rth=0.6,q.diffth=0.1,q.dcgth=0.1)
  saveRDS(diffcoexp_out,paste(names(msbb_gse84422_GPL570_samplesToAnalyse.exprs)[i],"diffcoexp_results.RDS",sep = "_"))
}


#Differential Expression analysis
msbb_gse84422_GPL96_97.DEGs=vector(mode = "list",length = 17)
msbb_gse84422_GPL570.DEGs=vector(mode = "list",length = 2)
names(msbb_gse84422_GPL96_97.DEGs)=names(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs)
names(msbb_gse84422_GPL570.DEGs)=names(msbb_gse84422_GPL570_samplesToAnalyse.exprs)
for(i in 1:17){
  c_exprs=msbb_gse84422_GPL96_97_samplesToAnalyse.exprs[[i]][,msbb_gse84422_GPL96_97_samplesToAnalyse[[i]]$SampleType=="CONTROL"]
  d_exprs=msbb_gse84422_GPL96_97_samplesToAnalyse.exprs[[i]][,msbb_gse84422_GPL96_97_samplesToAnalyse[[i]]$SampleType=="AD"]
  group=factor(c(rep("CONTROL",dim(c_exprs)[2]),rep("AD",dim(d_exprs)[2])))
  design_df=model.matrix(~0+group)
  colnames(design_df)=levels(group)
  contrasts_matrix=makeContrasts(CONTROL-AD,levels = design_df)
  fit=lmFit(object = cbind.data.frame(c_exprs,d_exprs),design = design_df)
  fit2=contrasts.fit(fit,contrasts_matrix)
  fit2=eBayes(fit2,trend = T)
  msbb_gse84422_GPL96_97.DEGs[[i]]=topTable(fit2,coef = 1,number = 20442,p.value = 0.1,adjust.method = "BH")%>%rownames_to_column("Gene")
}

for(i in 1:2){
  c_exprs=msbb_gse84422_GPL570_samplesToAnalyse.exprs[[i]][,msbb_gse84422_GPL570_samplesToAnalyse[[i]]$SampleType=="CONTROL"]
  d_exprs=msbb_gse84422_GPL570_samplesToAnalyse.exprs[[i]][,msbb_gse84422_GPL570_samplesToAnalyse[[i]]$SampleType=="AD"]
  group=factor(c(rep("CONTROL",dim(c_exprs)[2]),rep("AD",dim(d_exprs)[2])))
  design_df=model.matrix(~0+group)
  colnames(design_df)=levels(group)
  contrasts_matrix=makeContrasts(CONTROL-AD,levels = design_df)
  fit=lmFit(object = cbind.data.frame(c_exprs,d_exprs),design = design_df)
  fit2=contrasts.fit(fit,contrasts_matrix)
  fit2=eBayes(fit2,trend = T)
  msbb_gse84422_GPL570.DEGs[[i]]=topTable(fit2,coef = 1,number = 23521,p.value = 0.1,adjust.method = "BH")%>%rownames_to_column("Gene")
}



msbb_gse84422_sample_breakdown_GPL96_97.NP1=lapply(unique(msbb_gse84422.pData$GPL96$`neuropathological category:ch1`),function(x)msbb_gse84422.pData$GPL96%>%filter(`neuropathological category:ch1`==x)%>%select(c(`brain region:ch1`,`average neuritic plaque density:ch1`))%>%table)
msbb_gse84422_sample_breakdown_GPL570.NP1=lapply(unique(msbb_gse84422.pData$GPL570$`neuropathological category:ch1`),function(x)msbb_gse84422.pData$GPL570%>%filter(`neuropathological category:ch1`==x)%>%pull(`brain region:ch1`)%>%table)
msbb_gse84422_sample_breakdown_GPL96_97.Braak=lapply(c(2,4,6),function(x)msbb_gse84422.pData$GPL96%>%filter(`braak neurofibrillary tangle score:ch1`==x)%>%select(c('geo_accession','brain region:ch1','pseudoSampleID')))
msbb_gse84422_sample_breakdown_GPL570.Braak=lapply(c(2,4,6),function(x)msbb_gse84422.pData$GPL570%>%filter(`braak neurofibrillary tangle score:ch1`==x)%>%select(c('geo_accession','brain region:ch1')))
names(msbb_gse84422_sample_breakdown_GPL96_97.Braak)=names(msbb_gse84422_sample_breakdown_GPL570.Braak)=c(2,4,6)
names(msbb_gse84422_sample_breakdown_GPL96_97.NP1)=names(msbb_gse84422_sample_breakdown_GPL570.NP1)=unique(msbb_gse84422.pData$GPL570$`neuropathological category:ch1`)


msbb_gse84422_GPL96_97_samplesToAnalyse_Braak=vector(mode = "list",length = 3)
names(msbb_gse84422_GPL96_97_samplesToAnalyse_Braak)=c(2,4,6)

msbb_gse84422_GPL96_97_samplesToAnalyse_Braak$`2`=lapply(unique(msbb_gse84422_sample_breakdown_GPL96_97.Braak$`2`$`brain region:ch1`),function(x)msbb_gse84422_sample_breakdown_GPL96_97.Braak$`2`$pseudoSampleID[msbb_gse84422_sample_breakdown_GPL96_97.Braak$`2`$`brain region:ch1`==x])
msbb_gse84422_GPL96_97_samplesToAnalyse_Braak$`4`=lapply(unique(msbb_gse84422_sample_breakdown_GPL96_97.Braak$`4`$`brain region:ch1`),function(x)msbb_gse84422_sample_breakdown_GPL96_97.Braak$`4`$pseudoSampleID[msbb_gse84422_sample_breakdown_GPL96_97.Braak$`4`$`brain region:ch1`==x])
msbb_gse84422_GPL96_97_samplesToAnalyse_Braak$`6`=lapply(unique(msbb_gse84422_sample_breakdown_GPL96_97.Braak$`6`$`brain region:ch1`),function(x)msbb_gse84422_sample_breakdown_GPL96_97.Braak$`6`$pseudoSampleID[msbb_gse84422_sample_breakdown_GPL96_97.Braak$`6`$`brain region:ch1`==x])
names(msbb_gse84422_GPL96_97_samplesToAnalyse_Braak$`2`)=names(msbb_gse84422_GPL96_97_samplesToAnalyse_Braak$`4`)=names(msbb_gse84422_GPL96_97_samplesToAnalyse_Braak$`6`)=names(msbb_gse84422_GPL96_97_samplesToAnalyse)

braak_c_exprs=msbb_gse84422_exprs_GPL96_97.agg[,msbb_gse84422_GPL96_97_samplesToAnalyse_Braak$`2`$`Frontal Pole`]
braak_d1_exprs=msbb_gse84422_exprs_GPL96_97.agg[,msbb_gse84422_GPL96_97_samplesToAnalyse_Braak$`4`$`Frontal Pole`]
braak_d2_exprs=msbb_gse84422_exprs_GPL96_97.agg[,msbb_gse84422_GPL96_97_samplesToAnalyse_Braak$`6`$`Frontal Pole`]

braak_c_d1.diffcoexp=diffcoexp(exprs.1 = braak_c_exprs[1:1000,],exprs.2 = braak_d1_exprs[1:1000,],rth = 0.6,qth = 0.1,q.diffth = 0.1,r.method = "spearman",q.dcgth = 0.1)
braak_c_d2.diffcoexp=diffcoexp(exprs.1 = braak_c_exprs[1:1000,],exprs.2 = braak_d2_exprs[1:1000,],rth = 0.6,qth = 0.1,q.diffth = 0.1,r.method = "spearman",q.dcgth = 0.1)
braak_d1_d2.diffcoexp=diffcoexp(exprs.1 = braak_d1_exprs[1:1000,],exprs.2 = braak_d2_exprs[1:1000,],rth = 0.6,qth = 0.1,q.diffth = 0.1,r.method = "spearman",q.dcgth = 0.1)



#POC for Early Diffcoexp, Anterior Cingulate
c_exprs=msbb_gse84422_GPL96_97_samplesToAnalyse.exprs$`Anterior Cingulate`[which(rownames(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs$`Anterior Cingulate`)%in%zhang_celltype_ADgenes.list$Astrocytes),msbb_gse84422_GPL96_97_samplesToAnalyse$`Anterior Cingulate`$SampleType=="CONTROL"]
d_exprs=msbb_gse84422_GPL96_97_samplesToAnalyse.exprs$`Anterior Cingulate`[which(rownames(msbb_gse84422_GPL96_97_samplesToAnalyse.exprs$`Anterior Cingulate`)%in%zhang_celltype_ADgenes.list$Astrocytes),msbb_gse84422_GPL96_97_samplesToAnalyse$`Anterior Cingulate`$SampleType=="AD"]
diffcoexp_out=diffcoexp(exprs.1 = c_exprs[,],exprs.2 = d_exprs[,],r.method="spearman",rth=0.6,q.diffth=0.2,q.dcgth=0.1)


c_exprs.CERAD=msbb_gse84422_GPL96_97_samplesToAnalyse_CERAD.exprs$`Anterior Cingulate`[which(rownames(msbb_gse84422_GPL96_97_samplesToAnalyse_CERAD.exprs$`Anterior Cingulate`)%in%zhang_celltype_ADgenes.list$Astrocytes),msbb_gse84422_GPL96_97_samplesToAnalyse.CERAD$`Anterior Cingulate`$SampleTypeCERAD=="Normal"]
d_exprs.CERAD=msbb_gse84422_GPL96_97_samplesToAnalyse_CERAD.exprs$`Anterior Cingulate`[which(rownames(msbb_gse84422_GPL96_97_samplesToAnalyse_CERAD.exprs$`Anterior Cingulate`)%in%zhang_celltype_ADgenes.list$Astrocytes),msbb_gse84422_GPL96_97_samplesToAnalyse.CERAD$`Anterior Cingulate`$SampleTypeCERAD=="Probable_AD"]
diffcoexp_out.CERAD=diffcoexp(exprs.1 = c_exprs.CERAD[,],exprs.2 = d_exprs.CERAD[,],r.method="spearman",rth=0.6,q.diffth=0.2,q.dcgth=0.1)

c_exprs.CDR=msbb_gse84422_GPL96_97_samplesToAnalyse_CDR.exprs$`Anterior Cingulate`[which(rownames(msbb_gse84422_GPL96_97_samplesToAnalyse_CDR.exprs$`Anterior Cingulate`)%in%zhang_celltype_ADgenes.list$Astrocytes),msbb_gse84422_GPL96_97_samplesToAnalyse.CDR$`Anterior Cingulate`$SampleTypeCDR=="CDR0"]
d_exprs.CDR=msbb_gse84422_GPL96_97_samplesToAnalyse_CDR.exprs$`Anterior Cingulate`[which(rownames(msbb_gse84422_GPL96_97_samplesToAnalyse_CDR.exprs$`Anterior Cingulate`)%in%zhang_celltype_ADgenes.list$Astrocytes),msbb_gse84422_GPL96_97_samplesToAnalyse.CDR$`Anterior Cingulate`$SampleTypeCDR=="CDR1"]
diffcoexp_out.CDR=diffcoexp(exprs.1 = c_exprs.CDR[,],exprs.2 = d_exprs.CDR[,],r.method="spearman",rth=0.6,q.diffth=0.2,q.dcgth=0.1)

c_exprs.Braak=msbb_gse84422_GPL96_97_samplesToAnalyse_Braak.exprs$`Anterior Cingulate`[which(rownames(msbb_gse84422_GPL96_97_samplesToAnalyse_Braak.exprs$`Anterior Cingulate`)%in%zhang_celltype_ADgenes.list$Astrocytes),msbb_gse84422_GPL96_97_samplesToAnalyse.Braak$`Anterior Cingulate`$SampleTypeBraak=="Braak12"]
d_exprs.Braak=msbb_gse84422_GPL96_97_samplesToAnalyse_Braak.exprs$`Anterior Cingulate`[which(rownames(msbb_gse84422_GPL96_97_samplesToAnalyse_Braak.exprs$`Anterior Cingulate`)%in%zhang_celltype_ADgenes.list$Astrocytes),msbb_gse84422_GPL96_97_samplesToAnalyse.Braak$`Anterior Cingulate`$SampleTypeBraak=="Braak34"]
diffcoexp_out.Braak=diffcoexp(exprs.1 = c_exprs.Braak[,],exprs.2 = d_exprs.Braak[,],r.method="spearman",rth=0.6,q.diffth=0.2,q.dcgth=0.1)
