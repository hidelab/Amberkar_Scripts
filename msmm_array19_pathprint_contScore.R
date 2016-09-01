library(limma)
library(Biobase)
library(jetset)
library(pathprint)
library(metaArray)
library(doMC)
library(pheatmap)
library(entropy)
getEntropy <- function(mat, index){
  if (index > 2 | index < 1)
    stop("Indicate 1 for rows or 2 for columns")
  d <- apply(as.matrix(mat), index, function(x){discretize(x, numBins = 3, r=c(-1,1))})
  entropy.vec <- apply(d, 2, entropy)
  return(entropy.vec)
}
setwd("/Users/sandeepamberkar/Work/Data/MSMM-MSBB-HBTRC-Data/MSBB_Array19/Normalised_Data/")
msbb_array19.files=list.files(pattern = "*.tsv",full.names = T)
msbb_array19.files=msbb_array19.files[-c(18:19)]
msbb_array19=mapply(FUN = read.table,msbb_array19.files,MoreArgs = list(header=T,sep="\t",as.is=T))
names(msbb_array19)=c("AC","CN","DLPFC","FP","HP","IFG","ITG","MTG","OVC","PHG","PCC","PCG","PFC","PTMN","SPL","STG","TP")

# hgu133a.jscores=jscores(chip = 'hgu133a')
# hgu133a.jscores_filtered=hgu133a.jscores[-which(is.na(hgu133a.jscores$specificity)|is.na(hgu133a.jscores$coverage)|is.na(hgu133a.jscores$robust)==T),]
# hgu133a.jscores_os_filtered=hgu133a.jscores_filtered[which(hgu133a.jscores_filtered$overall>=0.25),]
# 
# uniq_entrez=unique(hgu133a.jscores_filtered$EntrezID)
# msbb_array19.2=list()
# for (m in 1:length(msbb_array19)){
#   tmp=matrix(NA,nrow=length(uniq_entrez),ncol=length(colnames(msbb_array19[[m]]))-4)
#   cat(paste("Working on region",names(msbb_array19)[m],"...\n",sep=" "))
#   cat(paste("****************************************************************************\n"))
#   for (e in 1:length(uniq_entrez)){
#     cat(paste("Computing median for gene",uniq_entrez[e],"...\n",sep=" "))
#     match=which(msbb_array19[[m]]$ENTREZ_GENE_ID%in%uniq_entrez[e])
#     cat(paste("Found",length(match),"probes for gene",uniq_entrez[e],"...\n",sep=" "))
#     tmp[e,]=apply(msbb_array19[[m]][match,-c(1:4)],2,median)
# 
#   }
#   msbb_array19.2[[m]]=data.frame(tmp)
#   rownames(msbb_array19.2[[m]])=uniq_entrez
#   colnames(msbb_array19.2[[m]])=colnames(msbb_array19[[m]][-c(1:4)])
# }
# msbb_array19=msbb_array19.2
msbb_array19.covariates=read.delim2("../AMP-AD_MSBB_MSSM_covariates_mRNA_AffymetrixU133AB.tsv",header = T,as.is = T)
msbb_array19.covariates$Age[which(msbb_array19.covariates$Age=="89+")]=90 #For simplicity, all samples with AOD as 89+ were labeled 90
msbb_array19.covariates$Age=as.integer(msbb_array19.covariates$Age)
msbb_array19.covariates$pH=as.numeric(msbb_array19.covariates$pH)
msbb_array19.covariates$CDR=as.numeric(msbb_array19.covariates$CDR)
msbb_array19.covariates$PLQ_Mn=as.numeric(msbb_array19.covariates$PLQ_Mn)
msbb_array19.covariates$BrainBank=paste("X",msbb_array19.covariates$BrainBank,sep="")
msbb_array19.covariates_noZeroNTrSum=msbb_array19.covariates[-which(as.numeric(msbb_array19.covariates$NTrSum)==0),]
msbb_array19.covariates_noZeroPLQMn=msbb_array19.covariates[-which(as.numeric(msbb_array19.covariates$PLQ_Mn)==0),]

braak1=msbb_array19.covariates$BrainBank[which(msbb_array19.covariates_noZeroNTrSum$Braak==1)] # PCG and PTM had few (3) samples
braak2=msbb_array19.covariates$BrainBank[which(msbb_array19.covariates_noZeroNTrSum$Braak==2)]
braak3=msbb_array19.covariates$BrainBank[which(msbb_array19.covariates_noZeroNTrSum$Braak==3)]
braak4=msbb_array19.covariates$BrainBank[which(msbb_array19.covariates_noZeroNTrSum$Braak==4)]
#braak5=paste("X",msbb_array19.covariates$BrainBank[which(msbb_array19.covariates$Braak==5)],sep=""), Dropped as HP,PFC,SPL had no samples while other regions had single sample
#There was only a single sample of Braak=9 hence, it was dropped from the analysis!
braak6=msbb_array19.covariates$BrainBank[which(msbb_array19.covariates_noZeroNTrSum$Braak==6)]

braak1_indices=lapply(lapply(msbb_array19,colnames),function(x)which(x%in%braak1))
braak2_indices=lapply(lapply(msbb_array19,colnames),function(x)which(x%in%braak2))
braak3_indices=lapply(lapply(msbb_array19,colnames),function(x)which(x%in%braak3))
braak4_indices=lapply(lapply(msbb_array19,colnames),function(x)which(x%in%braak4))
braak6_indices=lapply(lapply(msbb_array19,colnames),function(x)which(x%in%braak6))

braak1_samples=braak2_samples=braak3_samples=braak4_samples=braak6_samples=list()
for (i in 1:length(names(braak1_indices))){
  braak1_samples[[i]]=msbb_array19[[i]][,c(4,braak1_indices[[i]])]
}
for (i in 1:length(names(braak2_indices))){
  braak2_samples[[i]]=msbb_array19[[i]][,c(4,braak2_indices[[i]])]
}
for (i in 1:length(names(braak3_indices))){
  braak3_samples[[i]]=msbb_array19[[i]][,c(4,braak3_indices[[i]])]
}
for (i in 1:length(names(braak4_indices))){
  braak4_samples[[i]]=msbb_array19[[i]][,c(4,braak4_indices[[i]])]
}
for (i in 1:length(names(braak6_indices))){
  braak6_samples[[i]]=msbb_array19[[i]][,c(4,braak6_indices[[i]])]
}
names(braak1_samples)=names(braak1_indices)
names(braak2_samples)=names(braak2_indices)
names(braak3_samples)=names(braak3_indices)
names(braak4_samples)=names(braak4_indices)
names(braak6_samples)=names(braak6_indices)

braak1_samples.agg=lapply(lapply(braak1_samples,`[`,-c(1)),aggregate,by=list(EntrezID=msbb_array19[[1]]$ENTREZ_GENE_ID),mean)
for(i in 1:length(names(braak1_samples.agg))){
  braak1_samples.agg[[i]]=braak1_samples.agg[[i]][-1,]
  rownames(braak1_samples.agg[[i]])=braak1_samples.agg[[i]][,1]
  braak1_samples.agg[[i]]=braak1_samples.agg[[i]][,-1]
}
braak2_samples.agg=lapply(lapply(braak2_samples,`[`,-c(1)),aggregate,by=list(EntrezID=msbb_array19[[1]]$ENTREZ_GENE_ID),mean)
for(i in 1:length(names(braak2_samples.agg))){
  braak2_samples.agg[[i]]=braak2_samples.agg[[i]][-1,]
  rownames(braak2_samples.agg[[i]])=braak2_samples.agg[[i]][,1]
  braak2_samples.agg[[i]]=braak2_samples.agg[[i]][,-1]
}
braak3_samples.agg=lapply(lapply(braak3_samples,`[`,-c(1)),aggregate,by=list(EntrezID=msbb_array19[[1]]$ENTREZ_GENE_ID),mean)
for(i in 1:length(names(braak3_samples.agg))){
  braak3_samples.agg[[i]]=braak3_samples.agg[[i]][-1,]
  rownames(braak3_samples.agg[[i]])=braak3_samples.agg[[i]][,1]
  braak3_samples.agg[[i]]=braak3_samples.agg[[i]][,-1]
}
braak4_samples.agg=lapply(lapply(braak4_samples,`[`,-c(1)),aggregate,by=list(EntrezID=msbb_array19[[1]]$ENTREZ_GENE_ID),mean)
for(i in 1:length(names(braak4_samples.agg))){
  braak4_samples.agg[[i]]=braak4_samples.agg[[i]][-1,]
  rownames(braak4_samples.agg[[i]])=braak4_samples.agg[[i]][,1]
  braak4_samples.agg[[i]]=braak4_samples.agg[[i]][,-1]
}
braak6_samples.agg=lapply(lapply(braak6_samples,`[`,-c(1)),aggregate,by=list(EntrezID=msbb_array19[[1]]$ENTREZ_GENE_ID),mean)
for(i in 1:length(names(braak6_samples.agg))){
  braak6_samples.agg[[i]]=braak6_samples.agg[[i]][-1,]
  rownames(braak6_samples.agg[[i]])=braak6_samples.agg[[i]][,1]
  braak6_samples.agg[[i]]=braak6_samples.agg[[i]][,-1]
}


# braak1_samples=lapply(lapply(msbb_array19.agg,colnames),intersect,braak1)[-c(9,15)]
# braak2_samples=lapply(lapply(msbb_array19.agg,colnames),intersect,braak2)[-c(9,15)]
# braak3_samples=lapply(lapply(msbb_array19.agg,colnames),intersect,braak3)[-c(9,15)]
# braak4_samples=lapply(lapply(msbb_array19.agg,colnames),intersect,braak4)[-c(9,15)]
# braak6_samples=lapply(lapply(msbb_array19.agg,colnames),intersect,braak6)[-c(9,15)]


msbb_array19_Braak1_sampleData=msbb_array19_Braak2_sampleData=msbb_array19_Braak3_sampleData=msbb_array19_Braak4_sampleData=msbb_array19_Braak6_sampleData=list()
#Remove genes with NA
# for (r in 1:length(names(msbb_array19))){
#   msbb_array19_Braak1_sampleData[[r]]=msbb_array19.agg.braak1[[r]][-which(unlist(lapply(apply(msbb_array19.agg.braak1[[r]],1,function(x)x[which(is.na(x)==T)]),length))>0),]
#   msbb_array19_Braak2_sampleData[[r]]=msbb_array19.agg.braak2[[r]][-which(unlist(lapply(apply(msbb_array19.agg.braak2[[r]],1,function(x)x[which(is.na(x)==T)]),length))>0),]
#   msbb_array19_Braak3_sampleData[[r]]=msbb_array19.agg.braak3[[r]][-which(unlist(lapply(apply(msbb_array19.agg.braak3[[r]],1,function(x)x[which(is.na(x)==T)]),length))>0),]
#   msbb_array19_Braak4_sampleData[[r]]=msbb_array19.agg.braak4[[r]][-which(unlist(lapply(apply(msbb_array19.agg.braak4[[r]],1,function(x)x[which(is.na(x)==T)]),length))>0),]
#   msbb_array19_Braak6_sampleData[[r]]=msbb_array19.agg.braak6[[r]][-which(unlist(lapply(apply(msbb_array19.agg.braak6[[r]],1,function(x)x[which(is.na(x)==T)]),length))>0),]
#   
# }


# for (m in 1:length(braak1_samples)){
#   msbb_array19_Braak1_sampleData[[m]]=braak1_samples.agg[[m]][,c(1,which(colnames(braak1_samples.agg[[m]])%in%unlist(braak1_samples[[m]])))]
# }
# for (i in 1:length(braak2_samples)){
#   msbb_array19_Braak2_sampleData[[i]]=msbb_array19.agg.braak2[[i]][,c(1,which(colnames(msbb_array19.agg.braak2[[i]])%in%unlist(braak2_samples[[i]])))]
# }
# for (j in 1:length(braak3_samples)){
#   msbb_array19_Braak3_sampleData[[j]]=msbb_array19.agg.braak3[[j]][,c(1,which(colnames(msbb_array19.agg.braak3[[j]])%in%unlist(braak3_samples[[j]])))]
# }
# for (k in 1:length(braak4_samples)){
#   msbb_array19_Braak4_sampleData[[k]]=msbb_array19.agg.braak4[[k]][,c(1,which(colnames(msbb_array19.agg.braak4[[k]])%in%unlist(braak4_samples[[k]])))]
# }
# for (l in 1:length(braak6_samples)){
#   msbb_array19_Braak6_sampleData[[l]]=msbb_array19.agg.braak6[[l]][,c(1,which(colnames(msbb_array19.agg.braak6[[l]])%in%unlist(braak6_samples[[l]])))]
# }

# names(msbb_array19_Braak1_sampleData)=names(msbb_array19.agg.braak1)
# names(msbb_array19_Braak2_sampleData)=names(msbb_array19.agg.braak2)
# names(msbb_array19_Braak3_sampleData)=names(msbb_array19.agg.braak3)
# names(msbb_array19_Braak4_sampleData)=names(msbb_array19.agg.braak4)
# names(msbb_array19_Braak6_sampleData)=names(msbb_array19.agg.braak6)

msbb_array19_braak1.SCE=mapply(single.chip.enrichment,braak1_samples.agg,MoreArgs = list(geneset = pathprint.Hs.gs,transformation = "squared.rank","mean",T,T))
msbb_array19_braak2.SCE=mapply(single.chip.enrichment,braak2_samples.agg,MoreArgs = list(geneset = pathprint.Hs.gs,transformation = "squared.rank","mean",T,T))  
msbb_array19_braak3.SCE=mapply(single.chip.enrichment,braak3_samples.agg,MoreArgs = list(geneset = pathprint.Hs.gs,transformation = "squared.rank","mean",T,T))  
msbb_array19_braak4.SCE=mapply(single.chip.enrichment,braak4_samples.agg,MoreArgs = list(geneset = pathprint.Hs.gs,transformation = "squared.rank","mean",T,T))  
msbb_array19_braak6.SCE=mapply(single.chip.enrichment,braak6_samples.agg,MoreArgs = list(geneset = pathprint.Hs.gs,transformation = "squared.rank","mean",T,T))  

# pb <- txtProgressBar(min = 0,max = length(pathprint.Hs.gs),style = 1)
# for (mb in 1:length(names(msbb_array19_Braak6_sampleData.SCE))){
#     msbb_array19_Braak6_sampleData.POE=foreach(i=names(pathprint.Hs.gs)) %do%{
#     sample=msbb_array19_Braak6_sampleData.SCE[[mb]][i,]
#     sample.POE=sample
#     sample.POE[]=NA
#     sample.valid=sample[!(is.na(sample))]
#     fit=fit.em(sample.valid,cl=rep(0,length(sample.valid)))
#     setTxtProgressBar(pb,match(i,names(pathprint.Hs.gs)))
#     sample.POE[!(is.na(sample))]=fit$expr
#   }
# }
#   
msbb_array19_braak1.indices=lapply(lapply(msbb_array19_braak1.SCE,colnames),function(x)which(x%in%msbb_array19.covariates$BrainBank))
msbb_array19_braak2.indices=lapply(lapply(msbb_array19_braak2.SCE,colnames),function(x)which(x%in%msbb_array19.covariates$BrainBank))
msbb_array19_braak3.indices=lapply(lapply(msbb_array19_braak3.SCE,colnames),function(x)which(x%in%msbb_array19.covariates$BrainBank))
msbb_array19_braak4.indices=lapply(lapply(msbb_array19_braak4.SCE,colnames),function(x)which(x%in%msbb_array19.covariates$BrainBank))
msbb_array19_braak6.indices=lapply(lapply(msbb_array19_braak6.SCE,colnames),function(x)which(x%in%msbb_array19.covariates$BrainBank))

msbb_array19_braak1.NTrSum=lapply(msbb_array19_braak1.indices,function(x)msbb_array19.covariates_noZeroNTrSum[x,12])
msbb_array19_braak2.NTrSum=lapply(msbb_array19_braak2.indices,function(x)msbb_array19.covariates_noZeroNTrSum[x,12])
msbb_array19_braak3.NTrSum=lapply(msbb_array19_braak3.indices,function(x)msbb_array19.covariates_noZeroNTrSum[x,12])
msbb_array19_braak4.NTrSum=lapply(msbb_array19_braak4.indices,function(x)msbb_array19.covariates_noZeroNTrSum[x,12])
msbb_array19_braak6.NTrSum=lapply(msbb_array19_braak6.indices,function(x)msbb_array19.covariates_noZeroNTrSum[x,12])

msbb_array19_braak1.PLQ_Mn=lapply(msbb_array19_braak1.indices,function(x)msbb_array19.covariates_noZeroNTrSum[x,10])
msbb_array19_braak2.PLQ_Mn=lapply(msbb_array19_braak2.indices,function(x)msbb_array19.covariates_noZeroNTrSum[x,10])
msbb_array19_braak3.PLQ_Mn=lapply(msbb_array19_braak3.indices,function(x)msbb_array19.covariates_noZeroNTrSum[x,10])
msbb_array19_braak4.PLQ_Mn=lapply(msbb_array19_braak4.indices,function(x)msbb_array19.covariates_noZeroNTrSum[x,10])
msbb_array19_braak6.PLQ_Mn=lapply(msbb_array19_braak6.indices,function(x)msbb_array19.covariates_noZeroNTrSum[x,10])

msbb_array19_braak1.SCE=msbb_array19_braak1.SCE[-c(9,15)]
msbb_array19_braak2.SCE=msbb_array19_braak2.SCE[-c(9,15)]
msbb_array19_braak3.SCE=msbb_array19_braak3.SCE[-c(9,15)]
msbb_array19_braak4.SCE=msbb_array19_braak4.SCE[-c(9,15)]
msbb_array19_braak6.SCE=msbb_array19_braak6.SCE[-c(9,15)]

msbb_array19_braak1.corr_NTr_PP=msbb_array19_braak2.corr_NTr_PP=msbb_array19_braak3.corr_NTr_PP=msbb_array19_braak4.corr_NTr_PP=msbb_array19_braak6.corr_NTr_PP=list()
for (i in 1:length(names(msbb_array19_braak1.SCE))){
    msbb_array19_braak1.corr_NTr_PP[[i]]=apply(msbb_array19_braak1.SCE[[i]],1,cor,msbb_array19_braak1.NTrSum[[i]])    
}
for (i in 1:length(names(msbb_array19_braak2.SCE))){
  msbb_array19_braak2.corr_NTr_PP[[i]]=apply(msbb_array19_braak2.SCE[[i]],1,cor,msbb_array19_braak2.NTrSum[[i]])    
}
for (l in 1:length(names(msbb_array19_braak3.SCE))){
  msbb_array19_braak3.corr_NTr_PP[[l]]=apply(msbb_array19_braak3.SCE[[l]],1,cor,msbb_array19_braak3.NTrSum[[l]])
}
for (i in 1:length(names(msbb_array19_braak4.NTrSum))){
  msbb_array19_braak4.corr_NTr_PP[[i]]=apply(msbb_array19_braak4.SCE[[i]],1,cor,msbb_array19_braak4.NTrSum[[i]])    
}
for (i in 1:length(names(msbb_array19_braak6.NTrSum))){
  msbb_array19_braak6.corr_NTr_PP[[i]]=apply(msbb_array19_braak6.SCE[[i]],1,cor,msbb_array19_braak6.NTrSum[[i]])    
}
names(msbb_array19_braak1.corr_NTr_PP)=names(msbb_array19_braak1.SCE)
names(msbb_array19_braak2.corr_NTr_PP)=names(msbb_array19_braak2.SCE)
names(msbb_array19_braak3.corr_NTr_PP)=names(msbb_array19_braak3.SCE)
names(msbb_array19_braak4.corr_NTr_PP)=names(msbb_array19_braak4.SCE)
names(msbb_array19_braak6.corr_NTr_PP)=names(msbb_array19_braak6.SCE)

#Check correlation with amyloid plaques using the PLQ_Mn

msbb_array19.BraakAggCorr=list()
for (j in 1:length(names(msbb_array19_braak1.SCE))){
  msbb_array19.BraakAggCorr[[j]]=data.frame(msbb_array19_braak1.corr_NTr_PP[[j]],
                                            msbb_array19_braak2.corr_NTr_PP[[j]],
                                            msbb_array19_braak3.corr_NTr_PP[[j]],
                                            msbb_array19_braak4.corr_NTr_PP[[j]],
                                            msbb_array19_braak6.corr_NTr_PP[[j]],
                                            stringsAsFactors = F)
  colnames(msbb_array19.BraakAggCorr[[j]])=c("Braak1","Braak2","Braak3","Braak4","Braak6")
}
names(msbb_array19.BraakAggCorr)=names(msbb_array19_braak1.SCE)

for (j in 1:length(names(msbb_array19.BraakAggCorr))){
  msbb_array19.BraakAggCorr[[j]]=lapply(apply(sapply(msbb_array19.BraakAggCorr[[j]],is.na),2,function(x)which(x==T)),function(x)msbb_array19.BraakAggCorr[[j]][-x,])
  
}

#Correlate with amyloid plaues using the PLQ_Mn covariate
msbb_array19_braak1.corr_PLQ_PP=msbb_array19_braak2.corr_PLQ_PP=msbb_array19_braak3.corr_PLQ_PP=msbb_array19_braak4.corr_PLQ_PP=msbb_array19_braak6.corr_PLQ_PP=list()
for (i in 1:length(names(msbb_array19_braak1.SCE))){
  msbb_array19_braak1.corr_PLQ_PP[[i]]=apply(msbb_array19_braak1.SCE[[i]],1,cor,msbb_array19_braak1.NTrSum[[i]])    
}
for (i in 1:length(names(msbb_array19_braak2.SCE))){
  msbb_array19_braak2.corr_PLQ_PP[[i]]=apply(msbb_array19_braak2.SCE[[i]],1,cor,msbb_array19_braak2.NTrSum[[i]])    
}
for (l in 1:length(names(msbb_array19_braak3.SCE))){
  msbb_array19_braak3.corr_PLQ_PP[[l]]=apply(msbb_array19_braak3.SCE[[l]],1,cor,msbb_array19_braak3.NTrSum[[l]])
}
for (i in 1:length(names(msbb_array19_braak4.SCE))){
  msbb_array19_braak4.corr_PLQ_PP[[i]]=apply(msbb_array19_braak4.SCE[[i]],1,cor,msbb_array19_braak4.NTrSum[[i]])    
}
for (i in 1:length(names(msbb_array19_braak6.SCE))){
  msbb_array19_braak6.corr_PLQ_PP[[i]]=apply(msbb_array19_braak6.SCE[[i]],1,cor,msbb_array19_braak6.NTrSum[[i]])    
}
names(msbb_array19_braak1.corr_PLQ_PP)=names(msbb_array19_braak1.SCE)
names(msbb_array19_braak2.corr_PLQ_PP)=names(msbb_array19_braak2.SCE)
names(msbb_array19_braak3.corr_PLQ_PP)=names(msbb_array19_braak3.SCE)
names(msbb_array19_braak4.corr_PLQ_PP)=names(msbb_array19_braak4.SCE)
names(msbb_array19_braak6.corr_PLQ_PP)=names(msbb_array19_braak6.SCE)

#Check correlation with amyloid plaques using the PLQ_Mn

msbb_array19.BraakAgg_PLQCorr=list()
for (j in 1:length(names(msbb_array19_braak1.SCE))){
  msbb_array19.BraakAggCorr[[j]]=data.frame(msbb_array19_braak1.corr_PLQ_PP[[j]],
                                            msbb_array19_braak2.corr_PLQ_PP[[j]],
                                            msbb_array19_braak3.corr_PLQ_PP[[j]],
                                            msbb_array19_braak4.corr_PLQ_PP[[j]],
                                            msbb_array19_braak6.corr_PLQ_PP[[j]],
                                            stringsAsFactors = F)
  colnames(msbb_array19.BraakAggCorr[[j]])=c("Braak1","Braak2","Braak3","Braak4","Braak6")
}
names(msbb_array19.BraakAggCorr)=names(msbb_array19_braak1.SCE)

for (j in 1:length(names(msbb_array19.BraakAggCorr))){
  msbb_array19.BraakAggCorr[[j]]=lapply(apply(sapply(msbb_array19.BraakAggCorr[[j]],is.na),2,function(x)which(x==T)),function(x)msbb_array19.BraakAggCorr[[j]][-x,])
  
}
save(msbb_array19.BraakAggCorr,file="msbb_array19_PathPrint_ContCorrAnalysis.RData")
