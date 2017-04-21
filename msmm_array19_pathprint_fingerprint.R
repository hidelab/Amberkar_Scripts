#library(limma)
library(Biobase)
#library(jetset)
library(pathprint)
library(metaArray)
library(doMC)
library(pheatmap)

setwd("/shared/hidelab2/user/md4zsa/Work/Data/MSMM_Array/MSBB_Array19/Normalised_Data")
remove_ZeroSumPathways=function(x,y)list(x[-y,])
getEntropy <- function(mat, index){
  if (index > 2 | index < 1)
    stop("Indicate 1 for rows or 2 for columns")
  d <- apply(as.matrix(mat), index, function(x){discretize(x, numBins = 3, r=c(-1,1))})
  entropy.vec <- apply(d, 2, entropy)
  return(entropy.vec)
}

msbb_array19.files=list.files(pattern = "*.tsv",full.names = T)
msbb_array19.files=msbb_array19.files[-c(18:19)]
msbb_array19=mapply(FUN = read.table,msbb_array19.files,MoreArgs = list(header=T,sep="\t",as.is=T))
names(msbb_array19)=c("AC","CN","DLPFC","FP","HP","IFG","ITG","MTG","OVC","PHG","PCC","PCG","PFC","PTMN","SPL","STG","TP")
msbb_array19.covariates=read.delim2("AMP-AD_MSBB_MSSM_covariates_mRNA_AffymetrixU133AB.tsv",header = T,as.is = T)
msbb_array19.covariates$Age[which(msbb_array19.covariates$Age=="89+")]=90 #For simplicity, all samples with AOD as 89+ were labeled 90
msbb_array19.covariates$Age=as.integer(msbb_array19.covariates$Age)
msbb_array19.covariates$pH=as.numeric(msbb_array19.covariates$pH)
msbb_array19.covariates$CDR=as.numeric(msbb_array19.covariates$CDR)
msbb_array19.covariates$PLQ_Mn=as.numeric(msbb_array19.covariates$PLQ_Mn)
msbb_array19.covariates$BrainBank=paste("X",msbb_array19.covariates$BrainBank,sep="")

msbb_array19.2=lapply(msbb_array19,function(x){rownames(x)<-x$ID;x})
msbb_array19.2_CDR=msbb_array19.fingerprint=msbb_array19_fingerprint.HighEntropy05_Fingerprint=msbb_array19_fingerprint_HighEntropy05_Fingerprint.Core=vector(mode = "list",length = 7)
names(msbb_array19.2_CDR)=names(msbb_array19.fingerprint)=sort(unique(msbb_array19.covariates$CDR))
for(c in 1:7){
  msbb_array19.2_CDR[[c]]=lapply(msbb_array19,function(x)x[,which(colnames(x)%in%msbb_array19.covariates$BrainBank[which(msbb_array19.covariates$CDR==names(msbb_array19.2_CDR)[c])])])
  msbb_array19.2_CDR[[c]]=lapply(msbb_array19,function(x){rownames(x)<-x$ID;x})
  msbb_array19.fingerprint[[c]]=lapply(lapply(msbb_array19.2_CDR[[c]],function(x)x[,-c(1:4)]),exprs2fingerprint,platform = "GPL96",species = "human",progressBar = T)
  msbb_array19_fingerprint.ZeroSumPathways=lapply(msbb_array19.fingerprint[[c]],function(x)which(rowSums(x)==0))
  msbb_array19_fingerprint.noZeromSumPathways=mapply(remove_ZeroSumPathways,msbb_array19.fingerprint[[c]],msbb_array19_fingerprint.ZeroSumPathways)
  msbb_array19_fingerprint.HighEntropy05_Indices=lapply(lapply(msbb_array19_fingerprint.noZeromSumPathways,getEntropy,1),function(x)unname(which(x>0.5)))
  msbb_array19_fingerprint.HighEntropy05_Fingerprint[[c]]=mapply(FUN = function(x,y)x[y,],msbb_array19_fingerprint.noZeromSumPathways,msbb_array19_fingerprint.HighEntropy05_Indices)
  msbb_array19_fingerprint_HighEntropy05_Fingerprint.Core[[c]]=Reduce(intersect,lapply(msbb_array19_fingerprint.HighEntropy05_Fingerprint[[c]],rownames))
  msbb_array19_fingerprint_HighEntropy05.CoreFingerprint=lapply(msbb_array19_fingerprint.HighEntropy05_Fingerprint[[c]],function(x)x[which(rownames(x)%in%msbb_array19_fingerprint_HighEntropy05_Fingerprint.Core[[c]]),])
  
}
saveRDS(msbb_array19.fingerprint,"msbb_CDR_Fingerprint.RDS")
saveRDS(msbb_array19_fingerprint.HighEntropy05_Fingerprint,"msbb_CDR_HE05_Fingerprint.RDS")



# CoreFingerprint.df=data.frame(msbb_array19_fingerprint_HighEntropy05.CoreFingerprint,stringsAsFactors = F)
# msbb_array19.corr_PLQ_PP_SigPathways=mapply(FUN = function(x,y)x[y,1],msbb_array19.corr_PLQ_PP,msbb_array19.corr_PLQ_PP.IndicesPval005)
# 
# msbb_array19.corr_PLQ_PP_SigPathways.CompMatrix=matrix(NA,nrow = 17,ncol = 17)
# colnames(msbb_array19.corr_PLQ_PP_SigPathways.CompMatrix)=rownames(msbb_array19.corr_PLQ_PP_SigPathways.CompMatrix)=names(msbb_array19)
# for(i in 1:17){
#   msbb_array19.corr_PLQ_PP_SigPathways.CompMatrix[i,]=unname(unlist(lapply(lapply(msbb_array19.corr_PLQ_PP_SigPathways,intersect,msbb_array19.corr_PLQ_PP_SigPathways[[i]]),length))/sum(unlist(lapply(msbb_array19.corr_PLQ_PP_SigPathways,length))))
# }
# 
# msbb_PLQ_PP_CoreFingerprint_Overlap=data.frame(BrainRegion=names(lapply(lapply(msbb_array19.corr_PLQ_PP_SigPathways,intersect,rownames(CoreFingerprint.df)),length)),OverlapPathways=as.numeric(unlist(unname(lapply(lapply(msbb_array19.corr_PLQ_PP_SigPathways,intersect,rownames(CoreFingerprint.df)),length)))),stringsAsFactors = F)
# 
# msbb_array19_PLQ_PP_FingerprintUnion=Reduce(union,lapply(msbb_array19_fingerprint.HighEntropy05_Fingerprint,rownames))

# msbb_array19.HighPLQ_Samples=lapply(msbb_array19_fingerprint_HighEntropy05.CoreFingerprint,function(x)colnames(x)[which(colnames(x)%in%msbb_array19.covariates$BrainBank[which(msbb_array19.covariates$PLQ_Mn>15)])])
# msbb_array19_fingerprint_HighEntropy05.CoreFingerprint_HighPLQ=mapply(FUN = function(x,y)x[,y],msbb_array19_fingerprint_HighEntropy05.CoreFingerprint,msbb_array19.HighPLQ_Samples)
# HighPLQ_CoreFingerprint.df=data.frame(msbb_array19_fingerprint_HighEntropy05.CoreFingerprint_HighPLQ,stringsAsFactors = F)
# 
# msbb_array19.LowPLQ_Samples=lapply(msbb_array19_fingerprint_HighEntropy05.CoreFingerprint,function(x)colnames(x)[which(colnames(x)%in%msbb_array19.covariates$BrainBank[which(msbb_array19.covariates$PLQ_Mn<1)])])
# msbb_array19_fingerprint_HighEntropy05.CoreFingerprint_LowPLQ=mapply(FUN = function(x,y)x[,y],msbb_array19_fingerprint_HighEntropy05.CoreFingerprint,msbb_array19.LowPLQ_Samples)
# LowPLQ_CoreFingerprint.df=data.frame(msbb_array19_fingerprint_HighEntropy05.CoreFingerprint_LowPLQ,stringsAsFactors = F)

# anno_df=data.frame(BrainRegion=unlist(lapply(strsplit(x = colnames(CoreFingerprint.df),split = "\\."),`[[`,1)),stringsAsFactors = F)
# anno_df2=data.frame(BrainRegion=unlist(lapply(strsplit(x = colnames(HighPLQ_CoreFingerprint.df),split = "\\."),`[[`,1)),stringsAsFactors = F)
# anno_df3=data.frame(BrainRegion=unlist(lapply(strsplit(x = colnames(LowPLQ_CoreFingerprint.df),split = "\\."),`[[`,1)),stringsAsFactors = F)
# rownames(anno_df)=colnames(CoreFingerprint.df)
# rownames(anno_df2)=colnames(HighPLQ_CoreFingerprint.df)
# rownames(anno_df3)=colnames(LowPLQ_CoreFingerprint.df)
# 
# msbb_array19.NTr_PLQ_phenoData=msbb_array19.covariates[,c(1,6,10,12)]
# msbb_array19.phenoDataIndices=lapply(lapply(msbb_array19.2,colnames),function(x)which(msbb_array19.NTr_PLQ_phenoData$BrainBank%in%x))
# msbb_array19.NTr_PLQ_phenoData2=lapply(msbb_array19.phenoDataIndices,function(x)msbb_array19.NTr_PLQ_phenoData[x,])
# msbb_array19.NTr_PLQ_phenoData2=lapply(msbb_array19.NTr_PLQ_phenoData2,function(x){rownames(x)<-x$BrainBank;x})
# msbb_array19.NTr_PLQ_phenoData2=lapply(msbb_array19.NTr_PLQ_phenoData2,`[`,c(-1))
# 
# msbb_array19.eset=msbb_array19.DEG_PLQ=msbb_array19.DEG_NTR=list()
# control.plq=lapply(msbb_array19.NTr_PLQ_phenoData2,function(x)which(x$PLQ_Mn<1))
# disease.plq=lapply(msbb_array19.NTr_PLQ_phenoData2,function(x)which(x$PLQ_Mn>10))
# for (i in 1:length(names(msbb_array19.NTr_PLQ_phenoData2))){
#   design.df=matrix(NA,nrow=sum(length(control.plq[[i]]),length(disease.plq[[i]])),ncol=2)
#   rownames(design.df)=c(rownames(msbb_array19.NTr_PLQ_phenoData2[[i]][control.plq[[i]],]),rownames(msbb_array19.NTr_PLQ_phenoData2[[i]][disease.plq[[i]],]))
#   colnames(design.df)=c("control.plq","disease.plq")
#   design.df[which(rownames(design.df)%in%rownames(msbb_array19.NTr_PLQ_phenoData2[[i]][control.plq[[i]],])),1]=1
#   design.df[-which(rownames(design.df)%in%rownames(msbb_array19.NTr_PLQ_phenoData2[[i]][control.plq[[i]],])),1]=0
#   design.df[which(rownames(design.df)%in%rownames(msbb_array19.NTr_PLQ_phenoData2[[i]][disease.plq[[i]],])),2]=1
#   design.df[-which(rownames(design.df)%in%rownames(msbb_array19.NTr_PLQ_phenoData2[[i]][disease.plq[[i]],])),2]=0
#   phenoData=AnnotatedDataFrame(data=msbb_array19.NTr_PLQ_phenoData2[[i]][which(rownames(msbb_array19.NTr_PLQ_phenoData2[[i]])%in%rownames(design.df)),])
#   msbb_array19.eset[[i]]=ExpressionSet(assayData=as.matrix(msbb_array19.2[[i]][,which(colnames(msbb_array19.2[[i]])%in%rownames(design.df))]),phenoData=phenoData)  
#   fit=lmFit(msbb_array19.eset[[i]],design=design.df)
#   fit=eBayes(fit)
#   contMatrix=makeContrasts(CtrlvsDisease=disease.plq-control.plq,levels=design.df)
#   fit2=contrasts.fit(fit,contMatrix)
#   fit2=eBayes(fit2)
#   msbb_array19.DEG_PLQ[[i]]=topTableF(fit2,adjust.method="BH",number=Inf)
# }
# 
# 
# names(msbb_array19.eset)=names(msbb_array19.DEG_PLQ)=names(msbb_array19)

#PathPrint on PLQ bins

