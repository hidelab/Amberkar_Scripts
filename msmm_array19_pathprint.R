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
hgu133a.jscores=jscores(chip = 'hgu133a')
hgu133a.jscores_filtered=hgu133a.jscores[-which(is.na(hgu133a.jscores$specificity)|is.na(hgu133a.jscores$coverage)|is.na(hgu133a.jscores$robust)==T),]
hgu133a.jscores_os_filtered=hgu133a.jscores_filtered[which(hgu133a.jscores_filtered$overall>=0.25),]
#uniq_entrez=unique(hgu133a.jscores_os_0.2$EntrezID)
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
#   
#Read covariate info and do minor preprocessing
msbb_array19.covariates=read.delim2("../AMP-AD_MSBB_MSSM_covariates_mRNA_AffymetrixU133AB.tsv",header = T,as.is = T)
msbb_array19.covariates$Age[which(msbb_array19.covariates$Age=="89+")]=90 #For simplicity, all samples with AOD as 89+ were labeled 90
msbb_array19.covariates$Age=as.integer(msbb_array19.covariates$Age)
msbb_array19.covariates$pH=as.numeric(msbb_array19.covariates$pH)
msbb_array19.covariates$CDR=as.numeric(msbb_array19.covariates$CDR)
msbb_array19.covariates$PLQ_Mn=as.numeric(msbb_array19.covariates$PLQ_Mn)
msbb_array19.covariates$BrainBank=paste("X",msbb_array19.covariates$BrainBank,sep="")
#braak0=paste("X",msbb_array19.covariates$BrainBank[which(msbb_array19.covariates$Braak==0)],sep="")
braak1=msbb_array19.covariates$BrainBank[which(msbb_array19.covariates$Braak==1)] # PCG and PTM had few (3) samples
braak2=msbb_array19.covariates$BrainBank[which(msbb_array19.covariates$Braak==2)]
braak3=msbb_array19.covariates$BrainBank[which(msbb_array19.covariates$Braak==3)]
braak4=msbb_array19.covariates$BrainBank[which(msbb_array19.covariates$Braak==4)]
#braak5=paste("X",msbb_array19.covariates$BrainBank[which(msbb_array19.covariates$Braak==5)],sep=""), Dropped as HP,PFC,SPL had no samples while other regions had single sample
#There was only a single sample of Braak=9 hence, it was dropped from the analysis!
braak6=msbb_array19.covariates$BrainBank[which(msbb_array19.covariates$Braak==6)]


#braak_controls_length=lapply(lapply(lapply(msbb_array19,colnames),intersect,braak0_samples.region),length)
#Select Braak-score specific samples
#braak0_samples=lapply(lapply(msbb_array19,colnames),intersect,braak0)
braak1_samples=lapply(lapply(msbb_array19,colnames),intersect,braak1)
braak2_samples=lapply(lapply(msbb_array19,colnames),intersect,braak2)
braak3_samples=lapply(lapply(msbb_array19,colnames),intersect,braak3)
braak4_samples=lapply(lapply(msbb_array19,colnames),intersect,braak4)
braak6_samples=lapply(lapply(msbb_array19,colnames),intersect,braak6)


#Extract data for Braak-score specific samples from median probe intensity computation
# msbb_array19.braak2raw=msbb_array19.2[c(names(braak2_samples))]
# msbb_array19.braak3raw=msbb_array19.2[c(names(braak3_samples))]
# msbb_array19.braak4raw=msbb_array19.2[c(names(braak4_samples))]
# msbb_array19.braak6raw=msbb_array19.2[c(names(braak6_samples))]
#Extract data for Braak-score specific samples
#msbb_array19.braak0raw=msbb_array19[c(names(msbb_array19))]
msbb_array19.braak1raw=msbb_array19[c(names(braak1_samples))]
msbb_array19.braak2raw=msbb_array19[c(names(braak2_samples))]
msbb_array19.braak3raw=msbb_array19[c(names(braak3_samples))]
msbb_array19.braak4raw=msbb_array19[c(names(braak4_samples))]
msbb_array19.braak6raw=msbb_array19[c(names(braak6_samples))]

#msbb_array19_Braak0_sampleData=list()
msbb_array19_Braak1_sampleData=msbb_array19_Braak2_sampleData=msbb_array19_Braak3_sampleData=msbb_array19_Braak4_sampleData=msbb_array19_Braak6_sampleData=list()
for (m in 1:length(braak1_samples)){
  msbb_array19_Braak1_sampleData[[m]]=msbb_array19.braak1raw[[m]][,c(1,which(colnames(msbb_array19.braak1raw[[m]])%in%unlist(braak1_samples[[m]])))]
}
for (i in 1:length(braak2_samples)){
  msbb_array19_Braak2_sampleData[[i]]=msbb_array19.braak2raw[[i]][,c(1,which(colnames(msbb_array19.braak2raw[[i]])%in%unlist(braak2_samples[[i]])))]
}
for (j in 1:length(braak3_samples)){
  msbb_array19_Braak3_sampleData[[j]]=msbb_array19.braak3raw[[j]][,c(1,which(colnames(msbb_array19.braak3raw[[j]])%in%unlist(braak3_samples[[j]])))]
}
for (k in 1:length(braak4_samples)){
  msbb_array19_Braak4_sampleData[[k]]=msbb_array19.braak4raw[[k]][,c(1,which(colnames(msbb_array19.braak4raw[[k]])%in%unlist(braak4_samples[[k]])))]
}
for (l in 1:length(braak6_samples)){
  msbb_array19_Braak6_sampleData[[l]]=msbb_array19.braak6raw[[l]][,c(1,which(colnames(msbb_array19.braak6raw[[l]])%in%unlist(braak6_samples[[l]])))]
}
# for (m in 1:length(braak0_samples)){
#    msbb_array19_Braak0_sampleData[[m]]=msbb_array19.braak0raw[[m]][,c(1,which(colnames(msbb_array19.braak0raw[[m]])%in%unlist(braak0_samples[[m]])))]
# }

#names(msbb_array19_Braak0_sampleData)=names(braak0_samples.region)
names(msbb_array19_Braak1_sampleData)=names(braak1_samples)
names(msbb_array19_Braak2_sampleData)=names(braak2_samples)
names(msbb_array19_Braak3_sampleData)=names(braak3_samples)
names(msbb_array19_Braak4_sampleData)=names(braak4_samples)
names(msbb_array19_Braak6_sampleData)=names(braak6_samples)
#Remove the ID column from the datatables
#msbb_array19_Braak0_sampleData=lapply(msbb_array19_Braak0_sampleData,function(x){rownames(x)<-x$ID;x})
msbb_array19_Braak1_sampleData=lapply(msbb_array19_Braak1_sampleData,function(x){rownames(x)<-x$ID;x})
msbb_array19_Braak2_sampleData=lapply(msbb_array19_Braak2_sampleData,function(x){rownames(x)<-x$ID;x})
msbb_array19_Braak3_sampleData=lapply(msbb_array19_Braak3_sampleData,function(x){rownames(x)<-x$ID;x})
msbb_array19_Braak4_sampleData=lapply(msbb_array19_Braak4_sampleData,function(x){rownames(x)<-x$ID;x})
msbb_array19_Braak6_sampleData=lapply(msbb_array19_Braak6_sampleData,function(x){rownames(x)<-x$ID;x})

#Filter probes from Jetset dataset
#Braak0.filteredProbesList=lapply(lapply(msbb_array19_Braak0_sampleData,rownames),intersect,rownames(hgu133a.jscores_filtered))
Braak1.filteredProbesList=lapply(lapply(msbb_array19_Braak1_sampleData,rownames),intersect,rownames(hgu133a.jscores_filtered))
Braak2.filteredProbesList=lapply(lapply(msbb_array19_Braak2_sampleData,rownames),intersect,rownames(hgu133a.jscores_filtered))
Braak3.filteredProbesList=lapply(lapply(msbb_array19_Braak3_sampleData,rownames),intersect,rownames(hgu133a.jscores_filtered))
Braak4.filteredProbesList=lapply(lapply(msbb_array19_Braak4_sampleData,rownames),intersect,rownames(hgu133a.jscores_filtered))
Braak6.filteredProbesList=lapply(lapply(msbb_array19_Braak6_sampleData,rownames),intersect,rownames(hgu133a.jscores_filtered))

#msbb_array19_Braak0_sampleData.filteredProbes=list()
msbb_array19_Braak1_sampleData.filteredProbes=msbb_array19_Braak2_sampleData.filteredProbes=msbb_array19_Braak3_sampleData.filteredProbes=msbb_array19_Braak4_sampleData.filteredProbes=msbb_array19_Braak6_sampleData.filteredProbes=list()
for (m in 1:length(msbb_array19_Braak1_sampleData)){
  msbb_array19_Braak1_sampleData.filteredProbes[[m]]=msbb_array19_Braak1_sampleData[[m]][c(which(rownames(msbb_array19_Braak2_sampleData[[m]])%in%Braak1.filteredProbesList[[m]])),]
}
for (i in 1:length(msbb_array19_Braak2_sampleData)){
  msbb_array19_Braak2_sampleData.filteredProbes[[i]]=msbb_array19_Braak2_sampleData[[i]][c(which(rownames(msbb_array19_Braak2_sampleData[[i]])%in%Braak2.filteredProbesList[[i]])),]
}
for (j in 1:length(msbb_array19_Braak3_sampleData)){
  msbb_array19_Braak3_sampleData.filteredProbes[[j]]=msbb_array19_Braak3_sampleData[[j]][c(which(rownames(msbb_array19_Braak3_sampleData[[j]])%in%Braak3.filteredProbesList[[j]])),]
}
for (k in 1:length(msbb_array19_Braak4_sampleData)){
  msbb_array19_Braak4_sampleData.filteredProbes[[k]]=msbb_array19_Braak4_sampleData[[k]][c(which(rownames(msbb_array19_Braak4_sampleData[[k]])%in%Braak4.filteredProbesList[[k]])),]
}
for (l in 1:length(msbb_array19_Braak6_sampleData)){
  msbb_array19_Braak6_sampleData.filteredProbes[[l]]=msbb_array19_Braak6_sampleData[[l]][c(which(rownames(msbb_array19_Braak6_sampleData[[l]])%in%Braak6.filteredProbesList[[l]])),]
}
# for (m in 1:length(msbb_array19_Braak0_sampleData)){
#    msbb_array19_Braak6_sampleData.filteredProbes[[m]]=msbb_array19_Braak0_sampleData[[m]][c(which(rownames(msbb_array19_Braak0_sampleData[[m]])%in%Braak0.filteredProbesList[[m]])),]
#  }
#names(msbb_array19_Braak0_sampleData.filteredProbes)=names(msbb_array19_Braak0_sampleData)
names(msbb_array19_Braak1_sampleData.filteredProbes)=names(msbb_array19_Braak1_sampleData)
names(msbb_array19_Braak2_sampleData.filteredProbes)=names(msbb_array19_Braak2_sampleData)
names(msbb_array19_Braak3_sampleData.filteredProbes)=names(msbb_array19_Braak3_sampleData)
names(msbb_array19_Braak4_sampleData.filteredProbes)=names(msbb_array19_Braak4_sampleData)
names(msbb_array19_Braak6_sampleData.filteredProbes)=names(msbb_array19_Braak6_sampleData)

#msbb_array19_Braak0_sampleData.filteredProbes=lapply(msbb_array19_Braak0_sampleData.filteredProbes,`[`,-c(1:2))
msbb_array19_Braak1_sampleData.filteredProbes=lapply(msbb_array19_Braak1_sampleData.filteredProbes,`[`,-c(1:2))
msbb_array19_Braak2_sampleData.filteredProbes=lapply(msbb_array19_Braak2_sampleData.filteredProbes,`[`,-c(1:2))
msbb_array19_Braak3_sampleData.filteredProbes=lapply(msbb_array19_Braak3_sampleData.filteredProbes,`[`,-c(1:2))
msbb_array19_Braak4_sampleData.filteredProbes=lapply(msbb_array19_Braak4_sampleData.filteredProbes,`[`,-c(1:2))
msbb_array19_Braak6_sampleData.filteredProbes=lapply(msbb_array19_Braak6_sampleData.filteredProbes,`[`,-c(1:2))

#msbb_array19_Braak0_sampleData.filteredProbes=lapply(msbb_array19_Braak0_sampleData.filteredProbes,function(x){x[c("ID")]<-NULL;x})
msbb_array19_Braak1_sampleData.filteredProbes=lapply(msbb_array19_Braak1_sampleData.filteredProbes,function(x){x[c("ID")]<-NULL;x})
msbb_array19_Braak2_sampleData.filteredProbes=lapply(msbb_array19_Braak2_sampleData.filteredProbes,function(x){x[c("ID")]<-NULL;x})
# msbb_array19_Braak2_sampleData=lapply(msbb_array19_Braak2_sampleData,function(x){x[c("GB_ACC")]<-NULL;x})
msbb_array19_Braak3_sampleData.filteredProbes=lapply(msbb_array19_Braak3_sampleData.filteredProbes,function(x){x[c("ID")]<-NULL;x})
# msbb_array19_Braak3_sampleData=lapply(msbb_array19_Braak3_sampleData,function(x){x[c("GB_ACC")]<-NULL;x})
msbb_array19_Braak4_sampleData.filteredProbes=lapply(msbb_array19_Braak4_sampleData.filteredProbes,function(x){x[c("ID")]<-NULL;x})
# msbb_array19_Braak4_sampleData=lapply(msbb_array19_Braak4_sampleData,function(x){x[c("GB_ACC")]<-NULL;x})
msbb_array19_Braak6_sampleData.filteredProbes=lapply(msbb_array19_Braak6_sampleData.filteredProbes,function(x){x[c("ID")]<-NULL;x})
# msbb_array19_Braak6_sampleData=lapply(msbb_array19_Braak6_sampleData,function(x){x[c("GB_ACC")]<-NULL;x})
for (d in 1:length(msbb_array19_Braak1_sampleData.filteredProbes)){
  write.table(msbb_array19_Braak1_sampleData.filteredProbes[[d]],paste("msbb_B1_jsFiltered025_",names(msbb_array19_Braak1_sampleData.filteredProbes)[d],".txt",sep=""),sep = "\t",col.names = F,row.names = F,quote = F)
}
for (d in 1:length(msbb_array19_Braak2_sampleData.filteredProbes)){
  write.table(msbb_array19_Braak2_sampleData.filteredProbes[[d]],paste("msbb_B2_jsFiltered025_",names(msbb_array19_Braak2_sampleData.filteredProbes)[d],".txt",sep=""),sep = "\t",col.names = F,row.names = F,quote = F)
}
for (d in 1:length(msbb_array19_Braak3_sampleData.filteredProbes)){
  write.table(msbb_array19_Braak3_sampleData.filteredProbes[[d]],paste("msbb_B3_jsFiltered025_",names(msbb_array19_Braak3_sampleData.filteredProbes)[d],".txt",sep=""),sep = "\t",col.names = F,row.names = F,quote = F)
}
for (d in 1:length(msbb_array19_Braak4_sampleData.filteredProbes)){
  write.table(msbb_array19_Braak4_sampleData.filteredProbes[[d]],paste("msbb_B4_jsFiltered025_",names(msbb_array19_Braak4_sampleData.filteredProbes)[d],".txt",sep=""),sep = "\t",col.names = F,row.names = F,quote = F)
}
for (d in 1:length(msbb_array19_Braak6_sampleData.filteredProbes)){
  write.table(msbb_array19_Braak6_sampleData.filteredProbes[[d]],paste("msbb_B6_jsFiltered025_",names(msbb_array19_Braak6_sampleData.filteredProbes)[d],".txt",sep=""),sep = "\t",col.names = F,row.names = F,quote = F)
}
#####################################################################################################################################
msbb_pathprint.Braak=list()
msbb_pathprint.Braak[[1]]=msbb_pathprint.Braak[[2]]=msbb_pathprint.Braak[[3]]=msbb_pathprint.Braak[[4]]=msbb_pathprint.Braak[[5]]=list()
msbb_array19_Braak1_sampleData.filteredProbes.SCE=mapply(single.chip.enrichment,msbb_array19_Braak1_sampleData.filteredProbes,MoreArgs = list(geneset = "pathprint.Hs.gs",progressBar = T))
msbb_array19_Braak2_sampleData.filteredProbes.SCE=mapply(single.chip.enrichment,msbb_array19_Braak2_sampleData.filteredProbes,MoreArgs = list(geneset = "pathprint.Hs.gs",progressBar = T))
msbb_array19_Braak3_sampleData.filteredProbes.SCE=mapply(single.chip.enrichment,msbb_array19_Braak3_sampleData.filteredProbes,MoreArgs = list(geneset = "pathprint.Hs.gs",progressBar = T))
msbb_array19_Braak4_sampleData.filteredProbes.SCE=mapply(single.chip.enrichment,msbb_array19_Braak4_sampleData.filteredProbes,MoreArgs = list(geneset = "pathprint.Hs.gs",progressBar = T))
msbb_array19_Braak6_sampleData.filteredProbes.SCE=mapply(single.chip.enrichment,msbb_array19_Braak6_sampleData.filteredProbes,MoreArgs = list(geneset = "pathprint.Hs.gs",progressBar = T))

#msbb_array19_Braak0_sampleData.filteredProbes.fingerprint=mapply(exprs2fingerprint,msbb_array19_Braak0_sampleData.filteredProbes,MoreArgs = list(platform = "GPL96",species = "human",progressBar = T))
msbb_array19_Braak1_sampleData.filteredProbes.fingerprint=mapply(exprs2fingerprint,msbb_array19_Braak1_sampleData.filteredProbes,MoreArgs = list(platform = "GPL96",species = "human",progressBar = T))
msbb_array19_Braak2_sampleData.filteredProbes.fingerprint=mapply(exprs2fingerprint,msbb_array19_Braak2_sampleData.filteredProbes,MoreArgs = list(platform = "GPL96",species = "human",progressBar = T))
msbb_array19_Braak3_sampleData.filteredProbes.fingerprint=mapply(exprs2fingerprint,msbb_array19_Braak3_sampleData.filteredProbes,MoreArgs = list(platform = "GPL96",species = "human",progressBar = T))
msbb_array19_Braak4_sampleData.filteredProbes.fingerprint=mapply(exprs2fingerprint,msbb_array19_Braak4_sampleData.filteredProbes,MoreArgs = list(platform = "GPL96",species = "human",progressBar = T))
msbb_array19_Braak6_sampleData.filteredProbes.fingerprint=mapply(exprs2fingerprint,msbb_array19_Braak6_sampleData.filteredProbes,MoreArgs = list(platform = "GPL96",species = "human",progressBar = T))
#msbb_pathprint.Braak[[1]]=msbb_array19_Braak0_sampleData.filteredProbes.fingerprint
msbb_pathprint.Braak[[1]]=msbb_array19_Braak1_sampleData.filteredProbes.fingerprint
msbb_pathprint.Braak[[2]]=msbb_array19_Braak2_sampleData.filteredProbes.fingerprint
msbb_pathprint.Braak[[3]]=msbb_array19_Braak3_sampleData.filteredProbes.fingerprint
msbb_pathprint.Braak[[4]]=msbb_array19_Braak4_sampleData.filteredProbes.fingerprint
msbb_pathprint.Braak[[5]]=msbb_array19_Braak6_sampleData.filteredProbes.fingerprint
names(msbb_pathprint.Braak)=c("Braak1","Braak2","Braak3","Braak4","Braak6")
msbb_pathprint.Braak.AggRegion=msbb_pathprint.Braak.AnnVector=msbb_pathprint.Braak.rowEntropy=msbb_pathprint.Braak.colEntropy=msbb_pathprint.Braak.AggRegionSelect=list()
for (t in 1:length(names(msbb_array19))){
  msbb_pathprint.Braak.AggRegion[[t]]=data.frame(lapply(lapply(lapply(msbb_pathprint.Braak,`[`,t),data.frame),function(x)x[c(1:633),])[[1]],
                                                 lapply(lapply(lapply(msbb_pathprint.Braak,`[`,t),data.frame),function(x)x[c(1:633),])[[2]],
                                                 lapply(lapply(lapply(msbb_pathprint.Braak,`[`,t),data.frame),function(x)x[c(1:633),])[[3]],
                                                 lapply(lapply(lapply(msbb_pathprint.Braak,`[`,t),data.frame),function(x)x[c(1:633),])[[4]],
                                                 lapply(lapply(lapply(msbb_pathprint.Braak,`[`,t),data.frame),function(x)x[c(1:633),])[[5]])
  msbb_pathprint.Braak.AnnVector[[t]]=data.frame(BraakSample=c(rep("Braak1",length(colnames(lapply(lapply(lapply(msbb_pathprint.Braak,`[`,t),data.frame),function(x)x[c(1:633),])[[1]]))),
                                                               rep("Braak2",length(colnames(lapply(lapply(lapply(msbb_pathprint.Braak,`[`,t),data.frame),function(x)x[c(1:633),])[[2]]))),
                                                               rep("Braak3",length(colnames(lapply(lapply(lapply(msbb_pathprint.Braak,`[`,t),data.frame),function(x)x[c(1:633),])[[3]]))),
                                                               rep("Braak4",length(colnames(lapply(lapply(lapply(msbb_pathprint.Braak,`[`,t),data.frame),function(x)x[c(1:633),])[[4]]))),
                                                               rep("Braak6",length(colnames(lapply(lapply(lapply(msbb_pathprint.Braak,`[`,t),data.frame),function(x)x[c(1:633),])[[5]])))),
                                                 BrainRegion=c(rep(names(msbb_array19)[t],dim(msbb_pathprint.Braak.AggRegion[[t]])[2])),stringsAsFactors = F)
  rownames(msbb_pathprint.Braak.AnnVector[[t]])=colnames(msbb_pathprint.Braak.AggRegion[[t]])
}

names(msbb_pathprint.Braak.AggRegion)=names(msbb_array19)
msbb_pathprint.Braak.rowEntropy=lapply(msbb_pathprint.Braak.AggRegion,getEntropy,1)
msbb_pathprint.Braak.rowEntropy=lapply(lapply(msbb_pathprint.Braak.rowEntropy,function(x)x[which(x>=0.5)]),names)
msbb_pathprint.Braak.colEntropy=lapply(msbb_pathprint.Braak.AggRegion,getEntropy,2)
msbb_pathprint.Braak.colEntropy=lapply(lapply(msbb_pathprint.Braak.colEntropy,function(x)x[which(x>=0.5)]),names)
names(msbb_pathprint.Braak.colEntropy)=names(msbb_pathprint.Braak.rowEntropy)=names(msbb_array19)

for (p in 1:length(msbb_array19)){
  msbb_pathprint.Braak.AggRegionSelect[[p]]=msbb_pathprint.Braak.AggRegion[[p]][msbb_pathprint.Braak.rowEntropy[[p]],]
  pheatmap(msbb_pathprint.Braak.AggRegionSelect[[p]],
           annotation=msbb_pathprint.Braak.AnnVector[[p]],
           main=paste("Pathprint analysis",names(msbb_pathprint.Braak.AggRegion)[p],sep=" "),
           color=c("Blue","White","Red"),
           cluster_rows=T,cluster_cols=T,
           clustering_distance_rows="manhattan",
           clustering_distance_cols="manhattan",
           legend_breaks=c(-1, 0, 1),
           fontsize=8,fontsize_row=8,fontsize_col=8,
           cellwidth=8,cellheight=8)
}
names(msbb_pathprint.Braak.AggRegionSelect)=names(msbb_array19)
#lapply(lapply(msbb_pathprint.Braak.AggRegionSelect,colnames),intersect,lapply(lapply(msbb_pathprint.Braak.AnnVector,function(x)x[which(x$BraakSample=="Braak6"),]),rownames)[[1]])[[1]]
#rowSums(msbb_pathprint.Braak.AggRegionSelect[[1]][which(colnames(msbb_pathprint.Braak.AggRegionSelect[[1]])%in%lapply(lapply(msbb_pathprint.Braak.AggRegionSelect,colnames),intersect,lapply(lapply(msbb_pathprint.Braak.AnnVector,function(x)x[which(x$BraakSample=="Braak6"),]),rownames)[[1]])[[1]]),])

msbb_pathprint.Braak_AllActivePathways=list()
msbb_pathprint_Braak1_activePathways=msbb_pathprint_Braak2_activePathways=msbb_pathprint_Braak3_activePathways=msbb_pathprint_Braak4_activePathways=msbb_pathprint_Braak6_activePathways=list()
msbb_pathprint.Braak_AllActivePathways[[1]]=msbb_pathprint_Braak1_activePathways
msbb_pathprint.Braak_AllActivePathways[[2]]=msbb_pathprint_Braak2_activePathways
msbb_pathprint.Braak_AllActivePathways[[3]]=msbb_pathprint_Braak3_activePathways
msbb_pathprint.Braak_AllActivePathways[[4]]=msbb_pathprint_Braak4_activePathways
msbb_pathprint.Braak_AllActivePathways[[5]]=msbb_pathprint_Braak6_activePathways
names(msbb_pathprint.Braak_AllActivePathways)=c("Braak1","Braak2","Braak3","Braak4","Braak6")

for (bb in 1:length(names(msbb_pathprint.Braak_AllActivePathways))){
  for (ap in 1:length(msbb_pathprint.Braak.AggRegionSelect)){
    active_pathways=sort(rowSums(msbb_pathprint.Braak.AggRegionSelect[[ap]][which(colnames(msbb_pathprint.Braak.AggRegionSelect[[ap]])%in%lapply(lapply(msbb_pathprint.Braak.AggRegionSelect,colnames),intersect,lapply(lapply(msbb_pathprint.Braak.AnnVector,function(x)x[which(x$BraakSample==names(msbb_pathprint.Braak_AllActivePathways)[bb]),]),rownames)[[ap]])[[ap]])]))
    active_patients=sort(colSums(msbb_pathprint.Braak.AggRegionSelect[[ap]][which(colnames(msbb_pathprint.Braak.AggRegionSelect[[ap]])%in%lapply(lapply(msbb_pathprint.Braak.AggRegionSelect,colnames),intersect,lapply(lapply(msbb_pathprint.Braak.AnnVector,function(x)x[which(x$BraakSample==names(msbb_pathprint.Braak_AllActivePathways)[bb]),]),rownames)[[ap]])[[ap]])]))
    active_pathwaysInpatients=intersect(names(active_pathways),rownames(msbb_pathprint.Braak.AggRegionSelect[[ap]][,names(active_patients)]))
    df.active=data.frame(active_pathwaysInpatients,abs(unname(unlist(active_pathways))),stringsAsFactors = F)
    df.active2=df.active[which(df.active[,2]>3),]
    colnames(df.active2)=c(paste(names(msbb_pathprint.Braak.AggRegionSelect)[ap],"Pathways",sep="_"),"Nr.of.Samples")
    
    # up=sort(rowSums(msbb_pathprint.Braak.AggRegionSelect[[i]][which(colnames(msbb_pathprint.Braak.AggRegionSelect[[i]])%in%lapply(lapply(msbb_pathprint.Braak.AggRegionSelect,colnames),intersect,lapply(lapply(msbb_pathprint.Braak.AnnVector,function(x)x[which(x$BraakSample=="Braak6"),]),rownames)[[i]])[[i]])]),decreasing = F)
    # df.up=data.frame(names(up),abs(unname(unlist(up))),stringsAsFactors = F)
    # colnames(df.up)=c(paste(names(msbb_pathprint.Braak.AggRegionSelect)[i],"Pathways",sep="_"),"Nr.of.Samples")
    #
    if (names(msbb_pathprint.Braak_AllActivePathways)[bb]=="Braak1"){
      msbb_pathprint_Braak1_activePathways[[ap]]=df.active2  
    }
    if (names(msbb_pathprint.Braak_AllActivePathways)[bb]=="Braak2"){
      msbb_pathprint_Braak2_activePathways[[ap]]=df.active2  
    }
    if (names(msbb_pathprint.Braak_AllActivePathways)[bb]=="Braak3"){
      msbb_pathprint_Braak3_activePathways[[ap]]=df.active2  
    }
    if (names(msbb_pathprint.Braak_AllActivePathways)[bb]=="Braak4"){
      msbb_pathprint_Braak4_activePathways[[ap]]=df.active2  
    }
    if (names(msbb_pathprint.Braak_AllActivePathways)[bb]=="Braak6"){
      msbb_pathprint_Braak6_activePathways[[ap]]=df.active2  
    }
    #msbb_pathprint_Braak6_selectPathwaysDown[[i]]=df.down
  }
}
names(msbb_pathprint_Braak1_activePathways)=names(msbb_pathprint.Braak.AggRegionSelect)
names(msbb_pathprint_Braak2_activePathways)=names(msbb_pathprint.Braak.AggRegionSelect)
names(msbb_pathprint_Braak3_activePathways)=names(msbb_pathprint.Braak.AggRegionSelect)
names(msbb_pathprint_Braak4_activePathways)=names(msbb_pathprint.Braak.AggRegionSelect)
names(msbb_pathprint_Braak6_activePathways)=names(msbb_pathprint.Braak.AggRegionSelect)

msbb_pathprint.Braak_AllActivePathways[[1]]=msbb_pathprint_Braak1_activePathways
msbb_pathprint.Braak_AllActivePathways[[2]]=msbb_pathprint_Braak2_activePathways
msbb_pathprint.Braak_AllActivePathways[[3]]=msbb_pathprint_Braak3_activePathways
msbb_pathprint.Braak_AllActivePathways[[4]]=msbb_pathprint_Braak4_activePathways
msbb_pathprint.Braak_AllActivePathways[[5]]=msbb_pathprint_Braak6_activePathways

#PCG,PTMN had no significant active pathways, hence were removed
msbb_pathprint.Braak_AllActivePathways[[1]]=msbb_pathprint.Braak_AllActivePathways[[1]][-c(12,14)]
msbb_pathprint.Braak_AllActivePathways[[2]]=msbb_pathprint.Braak_AllActivePathways[[2]][-c(12,14)]
msbb_pathprint.Braak_AllActivePathways[[3]]=msbb_pathprint.Braak_AllActivePathways[[3]][-c(12,14)]
msbb_pathprint.Braak_AllActivePathways[[4]]=msbb_pathprint.Braak_AllActivePathways[[4]][-c(12,14)]
msbb_pathprint.Braak_AllActivePathways[[5]]=msbb_pathprint.Braak_AllActivePathways[[5]][-c(12,14)]

# msbb_pathprint.Braak_AllActiveMedianPathways=list()
# for (mm in 1:length(names(msbb_pathprint.Braak_AllActivePathways[[1]]))){
#   msbb_pathprint.Braak_AllActiveMedianPathways[[mm]]=sapply(sapply(lapply(lapply(msbb_pathprint.Braak_AllActivePathways,`[`,.id=names(msbb_pathprint.Braak_AllActivePathways[[1]])[mm]),as.data.frame),`[`,2),median)
#   
# }
# names(msbb_pathprint.Braak_AllActiveMedianPathways)=names(msbb_array19)[-c(12,14)]
# msbb_pathprint.Braak_AllActiveMedianPathways=data.frame(msbb_pathprint.Braak_AllActiveMedianPathways,stringsAsFactors = F)

#Compute common and differential pathways each region
msbb_pathprint.Braak1.AggRegionSelect_AllActivePathwayNames=lapply(lapply(msbb_pathprint.Braak_AllActivePathways[[1]],`[[`,1),sort,decreasing=F)
msbb_pathprint.Braak2.AggRegionSelect_AllActivePathwayNames=lapply(lapply(msbb_pathprint.Braak_AllActivePathways[[2]],`[[`,1),sort,decreasing=F)
msbb_pathprint.Braak3.AggRegionSelect_AllActivePathwayNames=lapply(lapply(msbb_pathprint.Braak_AllActivePathways[[3]],`[[`,1),sort,decreasing=F)
msbb_pathprint.Braak4.AggRegionSelect_AllActivePathwayNames=lapply(lapply(msbb_pathprint.Braak_AllActivePathways[[4]],`[[`,1),sort,decreasing=F)
msbb_pathprint.Braak6.AggRegionSelect_AllActivePathwayNames=lapply(lapply(msbb_pathprint.Braak_AllActivePathways[[5]],`[[`,1),sort,decreasing=F)

#msbb_pathprint.Braak_AllActiveMedianPathways_CDMatrix=matrix(NA,nrow=length(names(msbb_pathprint.Braak_AllActiveMedianPathways)),ncol=length(names(msbb_pathprint.Braak_AllActiveMedianPathways)))
#colnames(msbb_pathprint.Braak_AllActiveMedianPathways_CDMatrix)=rownames(msbb_pathprint.Braak_AllActiveMedianPathways_CDMatrix)=names(msbb_pathprint.Braak_AllActiveMedianPathways)

msbb_pathprint.AllActivePathways_Braak1.CDMatrix=matrix(NA,nrow=length(names(msbb_pathprint.Braak_AllActivePathways[[1]])),ncol=length(names(msbb_pathprint.Braak_AllActivePathways[[1]])))
msbb_pathprint.AllActivePathways_Braak2.CDMatrix=matrix(NA,nrow=length(names(msbb_pathprint.Braak_AllActivePathways[[2]])),ncol=length(names(msbb_pathprint.Braak_AllActivePathways[[2]])))
msbb_pathprint.AllActivePathways_Braak3.CDMatrix=matrix(NA,nrow=length(names(msbb_pathprint.Braak_AllActivePathways[[3]])),ncol=length(names(msbb_pathprint.Braak_AllActivePathways[[3]])))
msbb_pathprint.AllActivePathways_Braak4.CDMatrix=matrix(NA,nrow=length(names(msbb_pathprint.Braak_AllActivePathways[[4]])),ncol=length(names(msbb_pathprint.Braak_AllActivePathways[[4]])))
msbb_pathprint.AllActivePathways_Braak6.CDMatrix=matrix(NA,nrow=length(names(msbb_pathprint.Braak_AllActivePathways[[5]])),ncol=length(names(msbb_pathprint.Braak_AllActivePathways[[5]])))

colnames(msbb_pathprint.AllActivePathways_Braak1.CDMatrix)=rownames(msbb_pathprint.AllActivePathways_Braak1.CDMatrix)=names(msbb_pathprint.Braak_AllActivePathways[[1]])
colnames(msbb_pathprint.AllActivePathways_Braak2.CDMatrix)=rownames(msbb_pathprint.AllActivePathways_Braak2.CDMatrix)=names(msbb_pathprint.Braak_AllActivePathways[[2]])
colnames(msbb_pathprint.AllActivePathways_Braak3.CDMatrix)=rownames(msbb_pathprint.AllActivePathways_Braak3.CDMatrix)=names(msbb_pathprint.Braak_AllActivePathways[[3]])
colnames(msbb_pathprint.AllActivePathways_Braak4.CDMatrix)=rownames(msbb_pathprint.AllActivePathways_Braak4.CDMatrix)=names(msbb_pathprint.Braak_AllActivePathways[[4]])
colnames(msbb_pathprint.AllActivePathways_Braak6.CDMatrix)=rownames(msbb_pathprint.AllActivePathways_Braak6.CDMatrix)=names(msbb_pathprint.Braak_AllActivePathways[[5]])

# for (i in 1:length(names(msbb_pathprint.Braak_AllActiveMedianPathways))){
#     msbb_pathprint.Braak_AllActiveMedianPathways_CDMatrix[i,]=unlist(unname(lapply(lapply(msbb_pathprint.Braak.AggRegionSelect_PathwayNames,intersect,msbb_pathprint.Braak.AggRegionSelect_PathwayNames[[i]]),length)))
#     msbb_pathprint.Braak_AllActiveMedianPathways_CDMatrix[,i]=unlist(unname(lapply(lapply(msbb_pathprint.Braak.AggRegionSelect_PathwayNames,setdiff,msbb_pathprint.Braak.AggRegionSelect_PathwayNames[[i]]),length)))
# }
for (i in 1:length(colnames(msbb_pathprint.AllActivePathways_Braak1.CDMatrix))){
  msbb_pathprint.AllActivePathways_Braak1.CDMatrix[i,]=unlist(unname(lapply(lapply(msbb_pathprint.Braak1.AggRegionSelect_AllActivePathwayNames,intersect,msbb_pathprint.Braak1.AggRegionSelect_AllActivePathwayNames[[i]]),length)))
  msbb_pathprint.AllActivePathways_Braak1.CDMatrix[,i]=unlist(unname(lapply(lapply(msbb_pathprint.Braak1.AggRegionSelect_AllActivePathwayNames,setdiff,msbb_pathprint.Braak1.AggRegionSelect_AllActivePathwayNames[[i]]),length)))
}
for (i in 1:length(colnames(msbb_pathprint.AllActivePathways_Braak2.CDMatrix))){
  msbb_pathprint.AllActivePathways_Braak2.CDMatrix[i,]=unlist(unname(lapply(lapply(msbb_pathprint.Braak2.AggRegionSelect_AllActivePathwayNames,intersect,msbb_pathprint.Braak2.AggRegionSelect_AllActivePathwayNames[[i]]),length)))
  msbb_pathprint.AllActivePathways_Braak2.CDMatrix[,i]=unlist(unname(lapply(lapply(msbb_pathprint.Braak2.AggRegionSelect_AllActivePathwayNames,setdiff,msbb_pathprint.Braak2.AggRegionSelect_AllActivePathwayNames[[i]]),length)))
}
for (i in 1:length(colnames(msbb_pathprint.AllActivePathways_Braak3.CDMatrix))){
  msbb_pathprint.AllActivePathways_Braak3.CDMatrix[i,]=unlist(unname(lapply(lapply(msbb_pathprint.Braak3.AggRegionSelect_AllActivePathwayNames,intersect,msbb_pathprint.Braak3.AggRegionSelect_AllActivePathwayNames[[i]]),length)))
  msbb_pathprint.AllActivePathways_Braak3.CDMatrix[,i]=unlist(unname(lapply(lapply(msbb_pathprint.Braak3.AggRegionSelect_AllActivePathwayNames,setdiff,msbb_pathprint.Braak3.AggRegionSelect_AllActivePathwayNames[[i]]),length)))
}
for (i in 1:length(colnames(msbb_pathprint.AllActivePathways_Braak4.CDMatrix))){
  msbb_pathprint.AllActivePathways_Braak4.CDMatrix[i,]=unlist(unname(lapply(lapply(msbb_pathprint.Braak4.AggRegionSelect_AllActivePathwayNames,intersect,msbb_pathprint.Braak4.AggRegionSelect_AllActivePathwayNames[[i]]),length)))
  msbb_pathprint.AllActivePathways_Braak4.CDMatrix[,i]=unlist(unname(lapply(lapply(msbb_pathprint.Braak4.AggRegionSelect_AllActivePathwayNames,setdiff,msbb_pathprint.Braak4.AggRegionSelect_AllActivePathwayNames[[i]]),length)))
}
for (i in 1:length(colnames(msbb_pathprint.AllActivePathways_Braak6.CDMatrix))){
  msbb_pathprint.AllActivePathways_Braak6.CDMatrix[i,]=unlist(unname(lapply(lapply(msbb_pathprint.Braak6.AggRegionSelect_AllActivePathwayNames,intersect,msbb_pathprint.Braak6.AggRegionSelect_AllActivePathwayNames[[i]]),length)))
  msbb_pathprint.AllActivePathways_Braak6.CDMatrix[,i]=unlist(unname(lapply(lapply(msbb_pathprint.Braak6.AggRegionSelect_AllActivePathwayNames,setdiff,msbb_pathprint.Braak6.AggRegionSelect_AllActivePathwayNames[[i]]),length)))
}