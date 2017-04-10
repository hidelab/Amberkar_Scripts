library(limma)
library(Biobase)
library(jetset)
library(pathprint)
library(metaArray)
library(doMC)
library(pheatmap)
library(entropy)
library(Rarity)
getEntropy <- function(mat, index){
  if (index > 2 | index < 1)
    stop("Indicate 1 for rows or 2 for columns")
  d <- apply(as.matrix(mat), index, function(x){discretize(x, numBins = 3, r=c(-1,1))})
  entropy.vec <- apply(d, 2, entropy)
  return(entropy.vec)
}
setwd("/Users/sandeepamberkar/Work/Data/MSMM-MSBB-HBTRC-Data/msbb_array19/Normalised_Data/")
msbb_array19.files=list.files(pattern = "*.tsv",full.names = T)
msbb_array19.files=msbb_array19.files[-c(18:19)]
msbb_array19=mapply(FUN = read.table,msbb_array19.files,MoreArgs = list(header=T,sep="\t",as.is=T))
names(msbb_array19)=c("AC","CN","DLPFC","FP","HP","IFG","ITG","MTG","OVC","PHG","PCC","PCG","PFC","PTMN","SPL","STG","TP")
msbb_array19.2=lapply(msbb_array19,function(x){rownames(x)<-x$ID;x})
msbb_array19.fingerprint=lapply(msbb_array19.2,exprs2fingerprint,platform = "GPL96",species = "human",progressBar = T)
#Select regions that show known regions to be affected!!!
#msbb_array19=msbb_array19[c("PHG","HP","AC","STG","SPL","DLPFC","OVC")]
#msbb_array19.SampleperBrainRegion=list()
msbb_array19.SampleperBrainRegion=list()
# for (a in 1:length(msbb_array19.covariates$BrainBank)){
#   msbb_array19.SampleperBrainRegion[[a]]=paste(names(unlist(lapply(lapply(lapply(lapply(msbb_array19,`[`,-c(1:4)),colnames),sort),function(x)which(msbb_array19.covariates$BrainBank[a]%in%x)))),collapse = ",")
# }
msbb_array19.covariates=read.delim2("../AMP-AD_MSBB_MSSM_covariates_mRNA_AffymetrixU133AB.tsv",header = T,as.is = T)
msbb_array19.covariates$Age[which(msbb_array19.covariates$Age=="89+")]=90 #For simplicity, all samples with AOD as 89+ were labeled 90
msbb_array19.covariates$Age=as.integer(msbb_array19.covariates$Age)
msbb_array19.covariates$pH=as.numeric(msbb_array19.covariates$pH)
msbb_array19.covariates$CDR=as.numeric(msbb_array19.covariates$CDR)
msbb_array19.covariates$PLQ_Mn=as.numeric(msbb_array19.covariates$PLQ_Mn)
msbb_array19.covariates$BrainBank=paste("X",msbb_array19.covariates$BrainBank,sep="")
for (a in 1:length(msbb_array19.covariates$BrainBank)){
  msbb_array19.SampleperBrainRegion[[a]]=paste(names(unlist(lapply(lapply(lapply(lapply(msbb_array19,`[`,-c(1:4)),colnames),sort),function(x)which(msbb_array19.covariates$BrainBank[a]%in%x)))),collapse = ",")
  names(msbb_array19.SampleperBrainRegion)[a]=msbb_array19.covariates$BrainBank[a]
}

noRegion=names(unlist(lapply(msbb_array19.SampleperBrainRegion,function(x)x[which(x=="")])))
msbb_array19.SampleperBrainRegion=msbb_array19.SampleperBrainRegion[-which(names(msbb_array19.SampleperBrainRegion)%in%noRegion)]
msbbFiltered_indices=lapply(lapply(msbb_array19,colnames),function(x)which(x%in%names(msbb_array19.SampleperBrainRegion)))

msbb_array19.covariates_noZeroNTrSum=msbb_array19.covariates[-union(which(as.numeric(msbb_array19.covariates$NTrSum)==0),which(msbb_array19.covariates$BrainBank%in%noRegion)),12]
names(msbb_array19.covariates_noZeroNTrSum)=msbb_array19.covariates[-union(which(as.numeric(msbb_array19.covariates$NTrSum)==0),which(msbb_array19.covariates$BrainBank%in%noRegion)),1]
msbb_array19.covariates_noZeroPLQ=msbb_array19.covariates[-union(which(as.numeric(msbb_array19.covariates$NTrSum)==0),which(msbb_array19.covariates$BrainBank%in%noRegion)),10]
names(msbb_array19.covariates_noZeroPLQ)=msbb_array19.covariates[-union(which(as.numeric(msbb_array19.covariates$NTrSum)==0),which(msbb_array19.covariates$BrainBank%in%noRegion)),1]

msbb_array19.agg=list()
multiID.u133a=msbb_array19[[1]]$ENTREZ_GENE_ID[grep(pattern="///",msbb_array19[[1]]$ENTREZ_GENE_ID)]
u133a_universe=c(msbb_array19[[1]]$ENTREZ_GENE_ID[-grep(pattern="///",msbb_array19[[1]]$ENTREZ_GENE_ID)],gsub(pattern=" ",replacement="",lapply(strsplit(x=multiID.u133a,split="///"),`[`,1)))
msbb_array19=lapply(msbb_array19,function(x){x$ENTREZ_GENE_ID<-u133a_universe;x})
msbb_array19.agg=lapply(lapply(msbb_array19,`[`,-c(1:4)),aggregate,by=list(EntrezID=msbb_array19[[1]]$ENTREZ_GENE_ID),mean)

for (i in 1:length(names(msbb_array19.agg))){
  msbb_array19.agg[[i]]=msbb_array19.agg[[i]][-1,]
  rownames(msbb_array19.agg[[i]])=msbb_array19.agg[[i]][,1]
  msbb_array19.agg[[i]]=msbb_array19.agg[[i]][,-1]
}
msbb_array19.SCE=mapply(single.chip.enrichment,msbb_array19.agg,MoreArgs = list(geneset = pathprint.Hs.gs,transformation = "log.rank","median",F,T))
msbb_array19.SCE2=mapply(single.chip.enrichment,msbb_array19.agg,MoreArgs = list(geneset = pathprint.Hs.gs,transformation = "squared.rank","mean",T,T))
# msbb_array19.SampleperBrainRegion=list()
# for (a in 1:length(msbb_array19.covariates$BrainBank)){
#   msbb_array19.SampleperBrainRegion[[a]]=paste(names(unlist(lapply(lapply(lapply(lapply(msbb_array19,`[`,-c(1:4)),colnames),sort),function(x)which(msbb_array19.covariates$BrainBank[a]%in%x)))),collapse = ",")
# }


msbb_array19.corr_NTr_PP=msbb_array19.corr_PLQ_PP=msbb_array19.corr_TrackGenes_NTr=msbb_array19.corr_TrackGenes_PLQ=list()
#Correlate Pathprint SCE with nr.of.Tangles
tmp=matrix(NA,nrow=633,ncol=3)
for (i in 1:length(names(msbb_array19.SCE))){
  for (j in 1:633){
    x=msbb_array19.SCE[[i]][j,which(colnames(msbb_array19.SCE[[i]])%in%names(msbb_array19.covariates_noZeroNTrSum))]
    y=msbb_array19.covariates_noZeroNTrSum[which(names(msbb_array19.covariates_noZeroNTrSum)%in%colnames(msbb_array19.SCE[[i]]))]
    rho.sp=unname(cor.test(x,y,method = "spearman")$estimate)
    rho.p=cor.test(x,y)$p.value
    pathway=rownames(msbb_array19.SCE[[i]][,which(colnames(msbb_array19.SCE[[i]])%in%names(msbb_array19.covariates_noZeroNTrSum))])[j]
    tmp[j,1:3]=cbind(pathway,rho.sp,rho.p)
    msbb_array19.corr_NTr_PP[[i]]=data.frame(Pathway=tmp[,1],Rho.Spearman=as.numeric(tmp[,2]),Pval=as.numeric(tmp[,3]),stringsAsFactors = F)
  }
}
names(msbb_array19.corr_NTr_PP)=names(msbb_array19.SCE)

#Correlate Pathprint SCE with nr.of.Plaques
tmp2=matrix(NA,nrow=633,ncol=3)
for (i in 1:length(names(msbb_array19.SCE2))){
  for (j in 1:633){
    x=msbb_array19.SCE2[[i]][j,which(colnames(msbb_array19.SCE2[[i]])%in%names(msbb_array19.covariates_noZeroPLQ))]
    y=msbb_array19.covariates_noZeroPLQ[which(names(msbb_array19.covariates_noZeroPLQ)%in%colnames(msbb_array19.SCE2[[i]]))]
    rho.sp=unname(cor.test(x,y,method = "spearman")$estimate)
    rho.p=cor.test(x,y)$p.value
    pathway=rownames(msbb_array19.SCE2[[i]][,which(colnames(msbb_array19.SCE2[[i]])%in%names(msbb_array19.covariates_noZeroPLQ))])[j]
    tmp2[j,1:3]=cbind(pathway,rho.sp,rho.p)
    msbb_array19.corr_PLQ_PP[[i]]=data.frame(Pathway=tmp2[,1],Rho.Spearman=as.numeric(tmp2[,2]),Pval=as.numeric(tmp2[,3]),stringsAsFactors = F)
  }
}
names(msbb_array19.corr_PLQ_PP)=names(msbb_array19.SCE2)

#Correlate genes with nr.of.Tangles
msbb_array19.corr_TrackGenes_PLQ=msbb_array19.corr_TrackGenes_NTr=list()
tmp3=matrix(NA,nrow=dim(msbb_array19.agg[[1]])[1],ncol=3)
for (i in 1:length(names(msbb_array19.agg))){
  
  for (j in 1:dim(msbb_array19.agg[[1]])[1]){
      cat(paste("Processing gene - ",rownames(msbb_array19.agg[[i]])[j],"for brain region & NTr",names(msbb_array19.agg)[i],"\n",sep=" "))
      x1=unname(unlist(msbb_array19.agg[[i]][j,which(colnames(msbb_array19.agg[[i]])%in%names(msbb_array19.covariates_noZeroNTrSum))]))
      y1=msbb_array19.covariates_noZeroNTrSum[which(names(msbb_array19.covariates_noZeroNTrSum)%in%colnames(msbb_array19.agg[[i]]))]
      rho.sp=unname(cor.test(x1,y1,method = "spearman")$estimate)
      rho.p=cor.test(x1,y1)$p.value
      genes1=rownames(msbb_array19.agg[[i]][,which(colnames(msbb_array19.agg[[i]])%in%names(msbb_array19.covariates_noZeroNTrSum))])[j]
      tmp3[j,1:3]=cbind(genes1,rho.sp,rho.p)
      msbb_array19.corr_TrackGenes_NTr[[i]]=data.frame(Genes=tmp3[,1],Rho.Spearman=as.numeric(tmp3[,2]),Pval=as.numeric(tmp3[,3]),stringsAsFactors = F)
    }
}
names(msbb_array19.corr_TrackGenes_NTr)=names(msbb_array19.agg)   

#Correlate genes with nr.of.Plaques
tmp4=matrix(NA,nrow=dim(msbb_array19.agg[[1]])[1],ncol=3)
for (i in 1:length(names(msbb_array19.agg))){
  
  for (j in 1:dim(msbb_array19.agg[[1]])[1]){
    cat(paste("Processing gene - ",rownames(msbb_array19.agg[[i]])[j],"for brain region & PLQ",names(msbb_array19.agg)[i],"\n",sep=" "))
    x2=unname(unlist(msbb_array19.agg[[i]][j,which(colnames(msbb_array19.agg[[i]])%in%names(msbb_array19.covariates_noZeroPLQ))]))
    y2=msbb_array19.covariates_noZeroPLQ[which(names(msbb_array19.covariates_noZeroPLQ)%in%colnames(msbb_array19.agg[[i]]))]
    rho.sp=unname(cor.test(x2,y2,method = "spearman")$estimate)
    rho.p=cor.test(x2,y2)$p.value
    genes2=rownames(msbb_array19.agg[[i]][,which(colnames(msbb_array19.agg[[i]])%in%names(msbb_array19.covariates_noZeroPLQ))])[j]
    tmp4[j,1:3]=cbind(genes2,rho.sp,rho.p)
    #tmp4=tmp4[order(as.numeric(tmp4[which(tmp4[,3]<0.05),3]),decreasing = F),]
    msbb_array19.corr_TrackGenes_PLQ[[i]]=data.frame(Genes=tmp4[,1],Rho.Spearman=as.numeric(tmp4[,2]),Pval=as.numeric(tmp4[,3]),stringsAsFactors = F)
  }
}
names(msbb_array19.corr_TrackGenes_NTr)=names(msbb_array19.agg) 
names(msbb_array19.corr_TrackGenes_PLQ)=names(msbb_array19.agg)
names(msbb_array19.corr_NTr_PP)=names(msbb_array19.SCE)
names(msbb_array19.corr_PLQ_PP)=names(msbb_array19.SCE2)

msbb_array19.corr_NTr_PP.pval005=lapply(lapply(msbb_array19.corr_NTr_PP,function(x)x[which(x$Pval<0.05),3]),sort,decreasing=F)
msbb_array19.corr_PLQ_PP.pval005=lapply(lapply(msbb_array19.corr_PLQ_PP,function(x)x[which(x$Pval<0.05),3]),sort,decreasing=F)
msbb_array19.corr_NTr_PP.IndicesPval005=msbb_array19.corr_PLQ_PP.IndicesPval005=list()
for (i in 1:length(names(msbb_array19.corr_NTr_PP.pval005))){
  msbb_array19.corr_NTr_PP.IndicesPval005[[i]]=which(msbb_array19.corr_NTr_PP[[i]]$Pval%in%msbb_array19.corr_NTr_PP.pval005[[i]])
  msbb_array19.corr_PLQ_PP.IndicesPval005[[i]]=which(msbb_array19.corr_PLQ_PP[[i]]$Pval%in%msbb_array19.corr_PLQ_PP.pval005[[i]])
}
names(msbb_array19.corr_NTr_PP.IndicesPval005)=names(msbb_array19.corr_NTr_PP.pval005)
names(msbb_array19.corr_PLQ_PP.IndicesPval005)=names(msbb_array19.corr_PLQ_PP.pval005)

#Filter genes for 5% significance
msbb_array19.corr_TrackGenes_NTr_Indices005=lapply(lapply(msbb_array19.corr_TrackGenes_NTr,`[[`,3),function(x)which(x<=0.05))
msbb_array19.corr_TrackGenes_PLQ_Indices005=lapply(lapply(msbb_array19.corr_TrackGenes_PLQ,`[[`,3),function(x)which(x<0.05))
msbb_array19.corr_TrackGenes_NTr_Genes005=msbb_array19.corr_TrackGenes_PLQ_Genes005=list()
for (i in length(names(msbb_array19.corr_TrackGenes_NTr))){
  msbb_array19.corr_TrackGenes_NTr_Genes005[[i]]=msbb_array19.corr_TrackGenes_NTr[[i]][msbb_array19.corr_TrackGenes_NTr_Indices005[[i]],1]
}
for (i in 1:length(names(msbb_array19.corr_TrackGenes_PLQ))){
  msbb_array19.corr_TrackGenes_PLQ_Genes005[[i]]=msbb_array19.corr_TrackGenes_PLQ[[i]][msbb_array19.corr_TrackGenes_PLQ_Indices005[[i]],1]
}
names(msbb_array19.corr_TrackGenes_PLQ_Genes005)=names(msbb_array19.corr_TrackGenes_PLQ)  
names(msbb_array19.corr_TrackGenes_NTr_Genes005)=names(msbb_array19.corr_TrackGenes_NTr)  
save(msbb_array19.corr_TrackGenes_PLQ_Genes005,file="msbb_corr_PLQ_TrackGenes005.RData") 
save(msbb_array19.corr_TrackGenes_NTr_Genes005,file="msbb_corr_NTr_TrackGenes005.RData") 
# msbb_array19.corr_TrackGenes_PLQ_GenesMultipleIDsIndices=lapply(msbb_array19.corr_TrackGenes_PLQ_Genes005,grep,pattern="///")
# msbb_array19.corr_TrackGenes_NTr_GenesMultipleIDsIndices=lapply(msbb_array19.corr_TrackGenes_NTr_Genes005,grep,pattern="///")
# 
# for (a in 1:length(names(msbb_array19.corr_TrackGenes_PLQ_GenesMultipleIDsIndices))){
#   msbb_array19.corr_TrackGenes_PLQ_Genes005[[a]][msbb_array19.corr_TrackGenes_PLQ_GenesMultipleIDsIndices[[a]]]=as.numeric(unlist(lapply(strsplit(x = msbb_array19.corr_TrackGenes_PLQ_Genes005[[a]][msbb_array19.corr_TrackGenes_PLQ_GenesMultipleIDsIndices[[a]]],split = "///"),`[`,1)))
# }
# for (a in 1:length(names(msbb_array19.corr_TrackGenes_NTr_GenesMultipleIDsIndices))){
#   msbb_array19.corr_TrackGenes_NTr_Genes005[[a]][msbb_array19.corr_TrackGenes_NTr_GenesMultipleIDsIndices[[a]]]=as.numeric(unlist(lapply(strsplit(x = msbb_array19.corr_TrackGenes_NTr_Genes005[[a]][msbb_array19.corr_TrackGenes_NTr_GenesMultipleIDsIndices[[a]]],split = "///"),`[`,1)))
# }

#Enlist top pathways based on Spearman Rank for each region
msbb_array19.corr_NTr_PP.IndicesCorr025=lapply(lapply(msbb_array19.corr_NTr_PP,`[`,c(2:3)),function(x)which((x$Rho.Spearman > 0.25|x$Rho.Spearman < -0.25)&x$Pval<0.05))
msbb_array19.corr_PLQ_PP.IndicesCorr025=lapply(lapply(msbb_array19.corr_PLQ_PP,`[`,c(2:3)),function(x)which((x$Rho.Spearman > 0.25|x$Rho.Spearman < -0.25)&x$Pval<0.05))
msbb_array19.corr_NTr_PP_Corr025=msbb_array19.corr_PLQ_PP_Corr025=list()
for (r in 1:length(names(msbb_array19.corr_NTr_PP.IndicesCorr025))){
  msbb_array19.corr_NTr_PP_Corr025[[r]]=msbb_array19.corr_NTr_PP[[r]][msbb_array19.corr_NTr_PP.IndicesCorr025[[r]],][order(msbb_array19.corr_NTr_PP[[r]][msbb_array19.corr_NTr_PP.IndicesCorr025[[r]],2]),]
}
for (r in 1:length(names(msbb_array19.corr_PLQ_PP.IndicesCorr025))){
  msbb_array19.corr_PLQ_PP_Corr025[[r]]=msbb_array19.corr_PLQ_PP[[r]][msbb_array19.corr_PLQ_PP.IndicesCorr025[[r]],][order(msbb_array19.corr_PLQ_PP[[r]][msbb_array19.corr_PLQ_PP.IndicesCorr025[[r]],2]),]
}
names(msbb_array19.corr_NTr_PP_Corr025)=names(msbb_array19.corr_PLQ_PP_Corr025)=names(msbb_array19.corr_PLQ_PP.IndicesCorr025)
msbb_array19.corr_NTr_PP_Corr025_AllPathways=unique(unlist(lapply(lapply(msbb_array19.corr_NTr_PP_Corr025,`[[`,1),sort)))
msbb_array19.corr_PLQ_PP_Corr025_AllPathways=unique(unlist(lapply(lapply(msbb_array19.corr_PLQ_PP_Corr025,`[[`,1),sort)))

msbb_array19.corr_NTr_PP_Corr025_AllPathways_CountMatrix=matrix(NA,nrow=length(msbb_array19.corr_NTr_PP_Corr025_AllPathways),ncol=length(names(msbb_array19.corr_NTr_PP_Corr025)))
msbb_array19.corr_PLQ_PP_Corr025_AllPathways_CountMatrix=matrix(0,nrow=length(msbb_array19.corr_PLQ_PP_Corr025_AllPathways),ncol=length(names(msbb_array19.corr_PLQ_PP_Corr025)))
rownames(msbb_array19.corr_PLQ_PP_Corr025_AllPathways_CountMatrix)=msbb_array19.corr_PLQ_PP_Corr025_AllPathways
colnames(msbb_array19.corr_PLQ_PP_Corr025_AllPathways_CountMatrix)=names(msbb_array19.corr_PLQ_PP_Corr025)
for (r in 1:length(rownames(msbb_array19.corr_PLQ_PP_Corr025_AllPathways_CountMatrix))){
  nd=unlist(lapply(lapply(msbb_array19.corr_PLQ_PP_Corr025,`[[`,1),function(x)which(x%in%msbb_array19.corr_PLQ_PP_Corr025_AllPathways[r])))
  msbb_array19.corr_PLQ_PP_Corr025_AllPathways_CountMatrix[r,which(colnames(msbb_array19.corr_PLQ_PP_Corr025_AllPathways_CountMatrix)%in%names(nd))]=unname(unlist(lapply(nd,length)))
  
}

msbb_array19.corr_NTr_PP_Corr025_AllPathways_CountMatrix=matrix(0,nrow=length(msbb_array19.corr_NTr_PP_Corr025_AllPathways),ncol=length(names(msbb_array19.corr_NTr_PP_Corr025)))
rownames(msbb_array19.corr_NTr_PP_Corr025_AllPathways_CountMatrix)=msbb_array19.corr_NTr_PP_Corr025_AllPathways
colnames(msbb_array19.corr_NTr_PP_Corr025_AllPathways_CountMatrix)=names(msbb_array19.corr_NTr_PP_Corr025)
for (r in 1:length(rownames(msbb_array19.corr_NTr_PP_Corr025_AllPathways_CountMatrix))){
  nd=unlist(lapply(lapply(msbb_array19.corr_NTr_PP_Corr025,`[[`,1),function(x)which(x%in%msbb_array19.corr_NTr_PP_Corr025_AllPathways[r])))
  msbb_array19.corr_NTr_PP_Corr025_AllPathways_CountMatrix[r,which(colnames(msbb_array19.corr_NTr_PP_Corr025_AllPathways_CountMatrix)%in%names(nd))]=unname(unlist(lapply(nd,length)))
  
}
#Plot correlation plots for NTr
for (i in 1:length(names(msbb_array19.SCE))){
  for (j in 1:5){
    x1=msbb_array19.SCE[[i]][msbb_array19.corr_NTr_PP.IndicesPval005[[i]][j],names(msbb_array19.SCE[[i]][msbb_array19.corr_NTr_PP.IndicesPval005[[i]][j],which(colnames(msbb_array19.SCE[[i]])%in%names(msbb_array19.covariates_noZeroNTrSum))])]
    y1=msbb_array19.covariates_noZeroNTrSum[which(names(msbb_array19.covariates_noZeroNTrSum)%in%colnames(msbb_array19.SCE[[i]][msbb_array19.corr_NTr_PP.IndicesPval005[[i]],]))]  
    tiff(filename = paste("correlation_plots/corplot_NTr_PP",names(msbb_array19.SCE)[i],"TopPathway",j,".tiff",sep="_"),height = 800,width=800,units = "px",res = 200)
    corPlot(data.frame(x1,y1),method = "spearman",xlab = paste(names(msbb_array19.SCE)[i],"Rank",j,rownames(msbb_array19.SCE[[i]][msbb_array19.corr_NTr_PP.IndicesPval005[[i]],which(colnames(msbb_array19.SCE[[i]])%in%names(msbb_array19.covariates_noZeroNTrSum))])[j],sep="_"),
            cex.axis=0.5,cex.lab=1,cex.main=0.5,cex.sub=0.5,ylab="NTr",col="blue4",ties.method = "average")
    dev.off()
    }
  
}

#Plot correlation plots for PLQ
for (i in 1:length(names(msbb_array19.SCE2))){
for (j in 1:5){
  x1=msbb_array19.SCE2[[i]][msbb_array19.corr_PLQ_PP.IndicesPval005[[i]][j],which(colnames(msbb_array19.SCE2[[i]])%in%names(msbb_array19.covariates_noZeroNTrSum))]
  y1=msbb_array19.covariates_noZeroNTrSum[which(names(msbb_array19.covariates_noZeroNTrSum)%in%colnames(msbb_array19.SCE2[[i]]))]  
  tiff(filename = paste("correlation_plots/corplot_PLQ_PP",names(msbb_array19.SCE2)[i],"TopPathway",j,".tiff",sep="_"),height = 1000,width=1000,units = "px",res = 200)
  corPlot(data.frame(x1,y1),method = "spearman",xlab = paste(names(msbb_array19.SCE2)[i],"Rank",j,rownames(msbb_array19.SCE2[[i]][msbb_array19.corr_PLQ_PP.IndicesPval005[[i]],which(colnames(msbb_array19.SCE2[[i]])%in%names(msbb_array19.covariates_noZeroPLQ))])[j],sep="_"),cex.axis=0.5,cex.lab=1,cex.main=0.5,cex.sub=0.5,ylab="PLQ",col="blue4",ties.method = "average")
  dev.off()
}

}
#Map filtered tracking genes for plaques and tangles to gene symbol and merge tables
msbb_array19.corr_TrackGenes_NTr2_Genes005=msbb_array19.corr_TrackGenes_NTr_Genes005
msbb_array19.corr_TrackGenes_NTr2_Genes005=lapply(msbb_array19.corr_TrackGenes_NTr2_Genes005,as.character)
msbb_array19.corr_TrackGenes_NTr_Genes005.geneSYMBOL=lapply(msbb_array19.corr_TrackGenes_NTr2_Genes005,select,x=org.Hs.eg.db,keytype = "ENTREZID",columns = "SYMBOL")
msbb_array19.corr_TrackGenes_NTr_Genes005.gpRKEGG=lapply(lapply(msbb_array19.corr_TrackGenes_NTr_Genes005.geneSYMBOL,`[[`,2),gprofiler,organism = "hsapiens",significant = T,max_p_value=1,correction_method = "bonferroni",src_filter = "KEGG")
msbb_array19.corr_TrackGenes_NTr_Genes005.cpRKEGG=lapply(lapply(msbb_array19.corr_TrackGenes_NTr_Genes005,enrichKEGG,organism = "human",pvalueCutoff = 0.05,pAdjustMethod = "BH"),summary)

msbb_array19.corr_TrackGenes_PLQ2_Genes005=msbb_array19.corr_TrackGenes_PLQ_Genes005
msbb_array19.corr_TrackGenes_PLQ2_Genes005=lapply(msbb_array19.corr_TrackGenes_PLQ2_Genes005,as.character)
msbb_array19.corr_TrackGenes_PLQ_Genes005.geneSYMBOL=lapply(msbb_array19.corr_TrackGenes_PLQ2_Genes005,select,x=org.Hs.eg.db,keytype = "ENTREZID",columns = "SYMBOL")
msbb_array19.corr_TrackGenes_PLQ_Genes005.gpRKEGG=lapply(lapply(msbb_array19.corr_TrackGenes_PLQ_Genes005.geneSYMBOL,`[[`,2),gprofiler,organism = "hsapiens",significant = T,max_p_value=1,correction_method = "bonferroni",src_filter = "KEGG")
msbb_array19.corr_TrackGenes_PLQ_Genes005.cpRKEGG=lapply(lapply(msbb_array19.corr_TrackGenes_PLQ_Genes005,enrichKEGG,organism = "human",pvalueCutoff = 0.05,pAdjustMethod = "bonferroni"),summary)

msbb_pathprint.NTr_PP.CDMatrix=matrix(NA,nrow=length(names(msbb_array19.corr_NTr_PP)),ncol=length(names(msbb_array19.corr_NTr_PP)))
rownames(msbb_pathprint.NTr_PP.CDMatrix)=colnames(msbb_pathprint.NTr_PP.CDMatrix)=names(msbb_array19.corr_NTr_PP.IndicesPval005)
for (m in 1:length(colnames(msbb_pathprint.NTr_PP.CDMatrix))){
  msbb_pathprint.NTr_PP.CDMatrix[m,]=unname(unlist(lapply(lapply(msbb_array19.corr_NTr_PP.IndicesPval005,intersect,msbb_array19.corr_NTr_PP.IndicesPval005[[m]]),length)))
  #msbb_pathprint.NTr_PP.CDMatrix[,m]=unname(unlist(lapply(lapply(msbb_array19.corr_NTr_PP.IndicesPval005,setdiff,msbb_array19.corr_NTr_PP.IndicesPval005[[m]]),length)))
}


msbb_pathprint.PLQ_PP.CDMatrix=matrix(NA,nrow=length(names(msbb_array19.corr_PLQ_PP)),ncol=length(names(msbb_array19.corr_PLQ_PP)))
rownames(msbb_pathprint.PLQ_PP.CDMatrix)=colnames(msbb_pathprint.PLQ_PP.CDMatrix)=names(msbb_array19.corr_PLQ_PP.IndicesPval005)
for (m in 1:length(colnames(msbb_pathprint.PLQ_PP.CDMatrix))){
  msbb_pathprint.PLQ_PP.CDMatrix[m,]=unname(unlist(lapply(lapply(msbb_array19.corr_PLQ_PP.IndicesPval005,intersect,msbb_array19.corr_PLQ_PP.IndicesPval005[[m]]),length)))
  #msbb_pathprint.PLQ_PP.CDMatrix[,m]=unname(unlist(lapply(lapply(msbb_array19.corr_PLQ_PP.IndicesPval005,setdiff,msbb_array19.corr_PLQ_PP.IndicesPval005[[m]]),length)))
}

msbb_pathprint.NTr_TrackGenes.CDMatrix=matrix(NA,nrow=length(names(msbb_array19.corr_TrackGenes_NTr_Genes005.gpRKEGG)),ncol=length(names(msbb_array19.corr_TrackGenes_NTr_Genes005.gpRKEGG)))
rownames(msbb_pathprint.NTr_TrackGenes.CDMatrix)=colnames(msbb_pathprint.NTr_TrackGenes.CDMatrix)=names(msbb_array19.corr_TrackGenes_NTr_Genes005.gpRKEGG)
for (n in 1:length(colnames(msbb_pathprint.NTr_TrackGenes.CDMatrix))){
  msbb_pathprint.NTr_TrackGenes.CDMatrix[n,]=unname(unlist(lapply(lapply(lapply(msbb_array19.corr_TrackGenes_NTr_Genes005.gpRKEGG,`[[`,12),intersect,lapply(msbb_array19.corr_TrackGenes_NTr_Genes005.gpRKEGG,`[[`,12)[[n]]),length)))
}

msbb_pathprint.PLQ_TrackGenes.CDMatrix=matrix(NA,nrow=length(names(msbb_array19.corr_TrackGenes_PLQ_Genes005)),ncol=length(names(msbb_array19.corr_TrackGenes_PLQ_Genes005)))
rownames(msbb_pathprint.PLQ_TrackGenes.CDMatrix)=colnames(msbb_pathprint.PLQ_TrackGenes.CDMatrix)=names(msbb_array19.corr_TrackGenes_PLQ_Genes005)
for (n in 1:length(colnames(msbb_pathprint.PLQ_TrackGenes.CDMatrix))){
  msbb_pathprint.PLQ_TrackGenes.CDMatrix[n,]=unlist(lapply(lapply(msbb_array19.corr_TrackGenes_PLQ_Genes005,intersect,msbb_array19.corr_TrackGenes_PLQ_Genes005[[n]]),length))
}
pheatmap(msbb_pathprint.PLQ_TrackGenes.CDMatrix,cluster_rows = F,cluster_cols = F,main = "PLQ_PP,KEGG enrichment CDMatrix")

msbb_pathprint.NTr_TrackGenes.CDMatrix=matrix(NA,nrow=length(names(msbb_array19.corr_TrackGenes_NTr_Genes005)),ncol=length(names(msbb_array19.corr_TrackGenes_NTr_Genes005)))
rownames(msbb_pathprint.NTr_TrackGenes.CDMatrix)=colnames(msbb_pathprint.NTr_TrackGenes.CDMatrix)=names(msbb_array19.corr_TrackGenes_NTr_Genes005)
for (n in 1:length(colnames(msbb_pathprint.NTr_TrackGenes.CDMatrix))){
  msbb_pathprint.NTr_TrackGenes.CDMatrix[n,]=unlist(lapply(lapply(msbb_array19.corr_TrackGenes_NTr_Genes005,intersect,msbb_array19.corr_TrackGenes_NTr_Genes005[[n]]),length))
}
pheatmap(msbb_pathprint.NTr_TrackGenes.CDMatrix,cluster_rows = F,cluster_cols = F,main = "NTr_Genes CDMatrix")

msbb_pathprint.NTr_PLQ_PP.CDMatrix=matrix(NA,nrow=length(names(msbb_array19.corr_PLQ_PP)),ncol=length(names(msbb_array19.corr_PLQ_PP)))
rownames(msbb_pathprint.NTr_PLQ_PP.CDMatrix)=colnames(msbb_pathprint.NTr_PLQ_PP.CDMatrix)=names(msbb_array19.corr_PLQ_PP.IndicesPval005)
for (m in 1:length(colnames(msbb_pathprint.NTr_PLQ_PP.CDMatrix))){
  msbb_pathprint.NTr_PLQ_PP.CDMatrix[m,]=unname(unlist(lapply(lapply(msbb_array19.corr_PLQ_PP.IndicesPval005,intersect,msbb_array19.corr_NTr_PP.IndicesPval005[[m]]),length)))
  msbb_pathprint.NTr_PLQ_PP.CDMatrix[,m]=unname(unlist(lapply(lapply(msbb_array19.corr_PLQ_PP.IndicesPval005,intersect,msbb_array19.corr_NTr_PP.IndicesPval005[[m]]),length)))
}


write.table(msbb_array19.corr_NTr_PP[[16]][msbb_array19.corr_NTr_PP.IndicesPval005[[16]][which(msbb_array19.corr_NTr_PP.IndicesPval005[[16]]%in%unlist(msbb_array19.corr_NTr_PP.IndicesPval005[c(6:8,10,12)]))],],"NTr_PP_CommPathways_Gyri.txt",sep = "\t",col.names = T,row.names = F,quote = F)
comm_pathways_gyri_indices=intersect(intersect(intersect(intersect(intersect(unlist(msbb_array19.corr_NTr_PP.IndicesPval005[6]),unlist(msbb_array19.corr_NTr_PP.IndicesPval005[7])),unlist(msbb_array19.corr_NTr_PP.IndicesPval005[8])),unlist(msbb_array19.corr_NTr_PP.IndicesPval005[10])),unlist(msbb_array19.corr_NTr_PP.IndicesPval005[12])),unlist(msbb_array19.corr_NTr_PP.IndicesPval005[16]))

#Select samples based on specific brain regions
select.samples=paste(names(unlist(lapply(lapply(lapply(lapply(msbb_array19.SampleperBrainRegion,strsplit,split = ","),unlist),function(x)which(x%in%c("PHG","HP","AC","STG","SPL","DLPFC","OVC"))),function(z)which(length(z)==7)))),collapse = ",")
select.samplesIndices=lapply(lapply(msbb_array19.SCE[c('PHG','HP','AC','STG','SPL','DLPFC','OVC')],colnames),function(x)which(x%in%select.samples))
msbb_array19.selectSamplesData=msbb_array19.SCE[c('PHG','HP','AC','STG','SPL','DLPFC','OVC')]
for (s in 1:length(names(select.samplesIndices))){
  df=lapply(msbb_array19.selectSamplesData,data.frame)
  select.samplesData[[s]]=df[[s]][,select.samplesIndices[[s]]]
}
rm(df)

for(s in 1:length(names(select.samplesData))){
  write.table(select.samplesData[[s]],file = paste("msbb_C3DmatR1/msbb_array19_selectSamples",names(select.samplesData)[s],".txt",sep = "_"),sep = "\t",col.names = F,row.names = F,quote = F)
}
write(rownames(select.samplesData[[1]]),"msbb_C3DmatR1/geneid.txt",sep = "\n")


msbb_array19.corr_TrackGenes_NTr_Genes005.geneSymbolMerged=msbb_array19.corr_TrackGenes_PLQ_Genes005.geneSymbolMerged=list()
for (m1 in 1:length(names(msbb_array19.corr_TrackGenes_NTr))){
  msbb_array19.corr_TrackGenes_NTr_Genes005.geneSymbolMerged[[m1]]=merge.data.frame(msbb_array19.corr_TrackGenes_NTr_Genes005.geneSYMBOL[[m1]],msbb_array19.corr_TrackGenes_NTr[[m1]][msbb_array19.corr_TrackGenes_NTr_Indices005[[m1]],],by.x = "ENTREZID",by.y="Genes")
}
for (m1 in 1:length(names(msbb_array19.corr_TrackGenes_PLQ))){
  msbb_array19.corr_TrackGenes_PLQ_Genes005.geneSymbolMerged[[m1]]=merge.data.frame(msbb_array19.corr_TrackGenes_PLQ_Genes005.geneSYMBOL[[m1]],msbb_array19.corr_TrackGenes_PLQ[[m1]][msbb_array19.corr_TrackGenes_PLQ_Indices005[[m1]],],by.x = "ENTREZID",by.y="Genes")
}
save(msbb_array19.corr_TrackGenes_NTr_Genes005.geneSymbolMerged,file = "Mt.Sinai_array_TangleCorrGenes005.RData")
save(msbb_array19.corr_TrackGenes_PLQ_Genes005.geneSymbolMerged,file = "Mt.Sinai_array_PLQCorrGenes005.RData")
#Collate correlations for all brain regions, remove NAs and plot heatmaps each for tangles and plaques
# msbb_array19.corr_TrackGenes_PLQ.df=as.data.frame(msbb_array19.corr_TrackGenes_PLQ)
# msbb_array19.corr_TrackGenes_PLQ.df=msbb_array19.corr_TrackGenes_PLQ.df[-unique(unname(unlist(apply(apply(msbb_array19.corr_TrackGenes_PLQ.df,1,is.na),1,function(x)which(x==T))))),]
# msbb_array19.corr_TrackGenes_NTr.df=as.data.frame(msbb_array19.corr_TrackGenes_NTr)
# msbb_array19.corr_TrackGenes_NTr.df=msbb_array19.corr_TrackGenes_NTr.df[-unique(unname(unlist(apply(apply(msbb_array19.corr_TrackGenes_NTr.df,1,is.na),1,function(x)which(x==T))))),]
# 
# #Apply filter of correlation coefficient 0.3 >|< -0.3 for tangles and plaques
# tiff(filename = "Heatmap_msbb_array19_AggCorr0.3_tangles.tiff")
# pheatmap(msbb_array19.AggCorr.NTrdf[unique(unname(unlist(apply(msbb_array19.AggCorr.NTrdf,1,function(x)which(x > 0.3|x < -0.3))))),],cluster_rows = T,cluster_cols = T,clustering_method = "average",main="PathPrint- Correlated pathways with tangles/brain region")
# dev.off()
# tiff(filename = "Heatmap_msbb_array19_AggCorr0.3_plaques.tiff")
# pheatmap(msbb_array19.AggCorr.PLQdf[unique(unname(unlist(apply(msbb_array19.AggCorr.PLQdf,1,function(x)which(x > 0.3|x < -0.3))))),],cluster_rows = T,cluster_cols = T,clustering_method = "average",main="PathPrint- Correlated pathways with plaques/brain region")
# dev.off()
# 
# #Check correlation with amyloid plaques using the PLQ_Mn
# 
# msbb_array19.AggCorr=list()
# for (j in 1:length(names(msbb_array19.SCE))){
#   msbb_array19.AggCorr[[j]]=data.frame(msbb_array19.corr_NTr_PP[[j]],
#                                        msbb_array19.corr_PLQ_PP[[j]],     
#                                        stringsAsFactors = F)
#   colnames(msbb_array19.AggCorr[[j]])=c("NTr","PLQ")
# }
# names(msbb_array19.AggCorr)=names(msbb_array19.SCE)
# 
# for (j in 1:length(names(msbb_array19.AggCorr))){
#   msbb_array19.AggCorr[[j]]=msbb_array19.AggCorr[[j]][-apply(sapply(msbb_array19.AggCorr[[j]],is.na),2,function(x)which(x==T)),]
#   
# }
msbb_array19.AggCorr_Pathways=lapply(msbb_array19.AggCorr,rownames)
core_pathways=union(union(union(union(union(union(union(union(union(union(union(union(union(union(union(union(msbb_array19.AggCorr_Pathways[[1]],
                                                                                                              msbb_array19.AggCorr_Pathways[[2]]),
                                                                                                              msbb_array19.AggCorr_Pathways[[3]]),
                                                                                                              msbb_array19.AggCorr_Pathways[[4]]),
                                                                                                              msbb_array19.AggCorr_Pathways[[5]]),
                                                                                                              msbb_array19.AggCorr_Pathways[[6]]),
                                                                                                              msbb_array19.AggCorr_Pathways[[7]]),
                                                                                                              msbb_array19.AggCorr_Pathways[[8]]),
                                                                                                              msbb_array19.AggCorr_Pathways[[9]]),
                                                                                                              msbb_array19.AggCorr_Pathways[[10]]),
                                                                                                              msbb_array19.AggCorr_Pathways[[11]]),
                                                                                                              msbb_array19.AggCorr_Pathways[[12]]),
                                                                                                              msbb_array19.AggCorr_Pathways[[13]]),
                                                                                                              msbb_array19.AggCorr_Pathways[[14]]),
                                                                                                              msbb_array19.AggCorr_Pathways[[15]]),
                                                                                                              msbb_array19.AggCorr_Pathways[[16]]),
                                                                                                              msbb_array19.AggCorr_Pathways[[17]])



msbb_array19.AggCorr_PLQ_filtIndices=lapply(lapply(msbb_array19.AggCorr,`[`,2),function(x)which(x > 0.25|x < -0.25))
msbb_array19.AggCorr_NTr_filtIndices=lapply(lapply(msbb_array19.AggCorr,`[`,1),function(x)which(x > 0.25|x < -0.25))

for (i in 1:length(names(msbb_array19.AggCorr))){
  
}
for (j in 1:length(names(msbb_array19.AggCorr))){
  tiff(filename = paste("msbb_array19_NTr_PLQ",names(msbb_array19.AggCorr)[[j]],"PathPrint.tiff",sep="_"))
  scatterplot(NTr~PLQ,data = msbb_array19.AggCorr[[j]],
              xlab = "PLQ",ylab = "NTr",legend.coords = "topleft",ellipse=T,
              main= paste(names(msbb_array19.AggCorr)[[j]],"Tangles vs Plaques correlation, PathPrint results",sep=" "))
  dev.off()
}
# save(msbb_array19.AggCorr,file="msbb_array19_PathPrint_ContCorrAnalysis.RData")
df_anno=vector(mode = 'list',length=17)
names(df_anno)=names(msbb_array19)
for(i in 1:17){
  df_anno[[i]]=msbb_array19.covariates[msbb_array19.covariates$BrainBank%in%colnames(msbb_array19.SCE[[i]][msbb_array19.corr_PLQ_PP.IndicesPval005[[i]],]),c(7:12)]
  rownames(df_anno[[i]])=colnames(msbb_array19.SCE[[i]][msbb_array19.corr_PLQ_PP.IndicesPval005[[i]],])
}
for (i in 1:17){
  jpeg(filename = paste('Heatmap',names(msbb_array19.SCE)[i],'PLQ_PathPrint.jpg',sep = '_'),quality = 300,units = 'px',width = 1500,height = 1500)
  pheatmap(msbb_array19.SCE[[i]],cluster_rows = T,clustering_distance_rows = 'euclidean',cluster_cols = T,clustering_distance_cols= 'euclidean',annotation = df_anno[[i]],main = paste(names(msbb_array19.SCE)[i],'PLQ_PathPrint',sep = '_'),height = 1500,width = 1500,cellheight = 10,cellwidth = 10)
  dev.off()
}
