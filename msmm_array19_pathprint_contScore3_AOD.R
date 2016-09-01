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
setwd("/Users/sandeepamberkar/Work/Data/MSMM-MSBB-HBTRC-Data/msbb_array19/Normalised_Data/")
msbb_array19.files=list.files(pattern = "*.tsv",full.names = T)
msbb_array19.files=msbb_array19.files[-c(18:19)]
msbb_array19=mapply(FUN = read.table,msbb_array19.files,MoreArgs = list(header=T,sep="\t",as.is=T))
names(msbb_array19)=c("AC","CN","DLPFC","FP","HP","IFG","ITG","MTG","OVC","PHG","PCC","PCG","PFC","PTMN","SPL","STG","TP")
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

#Segregate dataset as per AOD
msbb_array19.AOD=list()
msbb_array19.AOD[[1]]=msbb_array19.covariates[which(msbb_array19.covariates$Age>=60 & msbb_array19.covariates$Age<=80),1]
#msbb_array19.AOD[[2]]=msbb_array19.covariates[which(msbb_array19.covariates$Age>=71 & msbb_array19.covariates$Age<=80),1]
msbb_array19.AOD[[2]]=msbb_array19.covariates[which(msbb_array19.covariates$Age>=81 & msbb_array19.covariates$Age<=89),1]
msbb_array19.AOD[[3]]=msbb_array19.covariates[which(msbb_array19.covariates$Age>=90),1]
names(msbb_array19.AOD)=c("AOD1","AOD2","AOD3")
#msbb_array19.covariates_noZeroNTrSum=msbb_array19.covariates[-union(which(as.numeric(msbb_array19.covariates$NTrSum)==0),which(msbb_array19.covariates$BrainBank%in%noRegion)),]
msbb_array19.NTrSum=list()
msbb_array19.NTrSum[[1]]=msbb_array19.covariates[which(msbb_array19.covariates$BrainBank%in%msbb_array19.AOD[[1]]),12]
names(msbb_array19.NTrSum[[1]])=msbb_array19.covariates[which(msbb_array19.covariates$BrainBank%in%msbb_array19.AOD[[1]]),1]
msbb_array19.NTrSum[[2]]=msbb_array19.covariates[which(msbb_array19.covariates$BrainBank%in%msbb_array19.AOD[[2]]),12]
names(msbb_array19.NTrSum[[2]])=msbb_array19.covariates[which(msbb_array19.covariates$BrainBank%in%msbb_array19.AOD[[2]]),1]
msbb_array19.NTrSum[[3]]=msbb_array19.covariates[which(msbb_array19.covariates$BrainBank%in%msbb_array19.AOD[[3]]),12]
names(msbb_array19.NTrSum[[3]])=msbb_array19.covariates[which(msbb_array19.covariates$BrainBank%in%msbb_array19.AOD[[3]]),1]
# msbb_array19.NTrSum[[4]]=msbb_array19.covariates[which(msbb_array19.covariates$BrainBank%in%msbb_array19.AOD[[4]]),12]
# names(msbb_array19.NTrSum[[4]])=msbb_array19.covariates[which(msbb_array19.covariates$BrainBank%in%msbb_array19.AOD[[4]]),1]
names(msbb_array19.NTrSum)=c("AOD1.NTrSum","AOD2.NTrSum","AOD3.NTrSum")
msbb_array19.PLQMn=list()
msbb_array19.PLQMn[[1]]=msbb_array19.covariates[which(msbb_array19.covariates$BrainBank%in%msbb_array19.AOD[[1]]),12]
names(msbb_array19.PLQMn[[1]])=msbb_array19.covariates[which(msbb_array19.covariates$BrainBank%in%msbb_array19.AOD[[1]]),1]
msbb_array19.PLQMn[[2]]=msbb_array19.covariates[which(msbb_array19.covariates$BrainBank%in%msbb_array19.AOD[[2]]),12]
names(msbb_array19.PLQMn[[2]])=msbb_array19.covariates[which(msbb_array19.covariates$BrainBank%in%msbb_array19.AOD[[2]]),1]
msbb_array19.PLQMn[[3]]=msbb_array19.covariates[which(msbb_array19.covariates$BrainBank%in%msbb_array19.AOD[[3]]),12]
names(msbb_array19.PLQMn[[3]])=msbb_array19.covariates[which(msbb_array19.covariates$BrainBank%in%msbb_array19.AOD[[3]]),1]
# msbb_array19.PLQMn[[4]]=msbb_array19.covariates[which(msbb_array19.covariates$BrainBank%in%msbb_array19.AOD[[4]]),12]
# names(msbb_array19.PLQMn[[4]])=msbb_array19.covariates[which(msbb_array19.covariates$BrainBank%in%msbb_array19.AOD[[4]]),1]
names(msbb_array19.PLQMn)=c("AOD1.PLQMn","AOD2.PLQMn","AOD3.PLQMn")

msbb_array19.agg=list()
msbb_array19.agg=lapply(lapply(msbb_array19,`[`,-c(1:4)),aggregate,by=list(EntrezID=msbb_array19[[1]]$ENTREZ_GENE_ID),median)
for (i in 1:length(names(msbb_array19.agg))){
  msbb_array19.agg[[i]]=msbb_array19.agg[[i]][-1,]
  rownames(msbb_array19.agg[[i]])=msbb_array19.agg[[i]][,1]
  msbb_array19.agg[[i]]=msbb_array19.agg[[i]][,-1]
}
msbb_array19.SCE=mapply(single.chip.enrichment,msbb_array19.agg,MoreArgs = list(geneset = pathprint.Hs.gs,transformation = "log.rank","median",F,T))
# msbb_array19.SampleperBrainRegion=list()
# for (a in 1:length(msbb_array19.covariates$BrainBank)){
#   msbb_array19.SampleperBrainRegion[[a]]=paste(names(unlist(lapply(lapply(lapply(lapply(msbb_array19,`[`,-c(1:4)),colnames),sort),function(x)which(msbb_array19.covariates$BrainBank[a]%in%x)))),collapse = ",")
# }

msbb_array19.AODindices=list()
msbb_array19.AODindices[[1]]=lapply(lapply(msbb_array19.SCE,colnames),function(x)which(x%in%msbb_array19.AOD[[1]]))
msbb_array19.AODindices[[2]]=lapply(lapply(msbb_array19.SCE,colnames),function(x)which(x%in%msbb_array19.AOD[[2]]))
msbb_array19.AODindices[[3]]=lapply(lapply(msbb_array19.SCE,colnames),function(x)which(x%in%msbb_array19.AOD[[3]]))
names(msbb_array19.AODindices)=names(msbb_array19.AOD)

msbb_array19.corr_NTr_PP=msbb_array19.corr_PLQ_PP=msbb_array19.corr_NTrGenes=msbb_array19.corr_PLQGenes=list()
msbb_array19.corr_NTr_PP[[1]]=msbb_array19.corr_NTr_PP[[2]]=msbb_array19.corr_NTr_PP[[3]]=list()
msbb_array19.corr_PLQ_PP[[1]]=msbb_array19.corr_PLQ_PP[[2]]=msbb_array19.corr_PLQ_PP[[3]]=list()
msbb_array19.corr_NTrGenes[[1]]=msbb_array19.corr_NTrGenes[[2]]=msbb_array19.corr_NTrGenes[[3]]=list()
msbb_array19.corr_PLQGenes[[1]]=msbb_array19.corr_PLQGenes[[2]]=msbb_array19.corr_PLQGenes[[3]]=list()
names(msbb_array19.corr_NTr_PP)=names(msbb_array19.corr_PLQ_PP)=names(msbb_array19.corr_NTrGenes)=names(msbb_array19.corr_PLQGenes)=names(msbb_array19.AODindices)
tmp=matrix(NA,nrow=633,ncol=3)
tmp[,1]=rownames(msbb_array19.SCE[[1]])
colnames(tmp)=c("Pathway","PCCMC","PValue")
for (i in 1:length(names(msbb_array19.corr_NTr_PP))){
  for (j in 1:length(names(msbb_array19.SCE))){
    for (k in 1:633){
      
      tmp[k,2:3]=cbind(PCCMC=unname(cor.test(x=msbb_array19.SCE[[j]][k,msbb_array19.AODindices[[i]][[j]]],y = msbb_array19.NTrSum[[i]][which(names(msbb_array19.NTrSum[[i]])%in%colnames(msbb_array19.SCE[[j]])[msbb_array19.AODindices[[i]][[j]]])])$estimate),
                      PValue=cor.test(x=msbb_array19.SCE[[j]][k,msbb_array19.AODindices[[i]][[j]]],y = msbb_array19.NTrSum[[i]][which(names(msbb_array19.NTrSum[[i]])%in%colnames(msbb_array19.SCE[[j]])[msbb_array19.AODindices[[i]][[j]]])])$p.value)
      #msbb_array19.corr_NTr_PP[[i]][[j]]=data.frame(tmp[which(tmp[,3]<=0.05),],stringsAsFactors = F)
      msbb_array19.corr_NTr_PP[[i]][[j]]=data.frame(tmp,stringsAsFactors = F)
    }
    
  }
}
names(msbb_array19.corr_NTr_PP[[1]])=names(msbb_array19.corr_NTr_PP[[2]])=names(msbb_array19.corr_NTr_PP[[3]])=names(msbb_array19.SCE)
# # for (i in 1:length(names(msbb_array19.SCE))){
# #   msbb_array19.corr_NTr_PP[[i]]=apply(msbb_array19.SCE[[i]][which(rownames(msbb_array19.SCE[[i]])%in%core_pathways),],1,cor,msbb_array19.NTrSum[[i]])    
# # }
# #[which(rownames(msbb_array19.SCE[[i]])%in%core_pathways),]
# for (i in 1:length(names(msbb_array19.SCE))){
#   msbb_array19.corr_PLQ_PP[[i]]=apply(msbb_array19.SCE[[i]],1,cor,msbb_array19.PLQ_Mn[[i]])    
# }
# for (t1 in 1:length(names(msbb_array19.agg))){
#   msbb_array19.corr_TrackGenes_NTr[[t1]]=apply(msbb_array19.agg[[t1]],1,cor,msbb_array19.NTrSum[[t1]])    
# }
# for (t2 in 1:length(names(msbb_array19.agg))){
#   msbb_array19.corr_TrackGenes_PLQ[[t2]]=apply(msbb_array19.agg[[t2]],1,cor,msbb_array19.PLQ_Mn[[t2]])    
# }
# names(msbb_array19.corr_TrackGenes_PLQ)=names(msbb_array19.agg)
# names(msbb_array19.corr_TrackGenes_NTr)=names(msbb_array19.agg)
# names(msbb_array19.corr_NTr_PP)=names(msbb_array19.SCE)
# names(msbb_array19.corr_PLQ_PP)=names(msbb_array19.SCE)
# 
# #Collate correlations for all brain regions, remove NAs and plot heatmaps each for tangles and plaques
# msbb_array19.corr_TrackGenes_PLQ.df=as.data.frame(msbb_array19.corr_TrackGenes_PLQ)
# msbb_array19.corr_TrackGenes_PLQ.df=msbb_array19.corr_TrackGenes_PLQ.df[-unique(unname(unlist(apply(apply(msbb_array19.corr_TrackGenes_PLQ.df,1,is.na),1,function(x)which(x==T))))),]
# msbb_array19.corr_TrackGenes_NTr.df=as.data.frame(msbb_array19.corr_TrackGenes_NTr)
# msbb_array19.corr_TrackGenes_NTr.df=msbb_array19.corr_TrackGenes_NTr.df[-unique(unname(unlist(apply(apply(msbb_array19.corr_TrackGenes_NTr.df,1,is.na),1,function(x)which(x==T))))),]
# 
# msbb_array19.corr_TrackGenes_NTr.genes=lapply(msbb_array19.corr_TrackGenes_NTr.df_indices0.3,names)
# msbb_array19.corr_TrackGenes_PLQ.genes=lapply(msbb_array19.corr_TrackGenes_PLQ.df_indices0.3,names)
# msbb_array19.corr_TrackGenes_NTr.genesDF=data.frame(Entrez=unlist(msbb_array19.corr_TrackGenes_NTr.genes),
#                                                     group=c(rep(names(msbb_array19.corr_TrackGenes_NTr.genes)[1],length(msbb_array19.corr_TrackGenes_NTr.genes[[1]])),
#                                                             rep(names(msbb_array19.corr_TrackGenes_NTr.genes)[2],length(msbb_array19.corr_TrackGenes_NTr.genes[[2]])),
#                                                             rep(names(msbb_array19.corr_TrackGenes_NTr.genes)[3],length(msbb_array19.corr_TrackGenes_NTr.genes[[3]])),
#                                                             rep(names(msbb_array19.corr_TrackGenes_NTr.genes)[4],length(msbb_array19.corr_TrackGenes_NTr.genes[[4]])),
#                                                             rep(names(msbb_array19.corr_TrackGenes_NTr.genes)[5],length(msbb_array19.corr_TrackGenes_NTr.genes[[5]])),
#                                                             rep(names(msbb_array19.corr_TrackGenes_NTr.genes)[6],length(msbb_array19.corr_TrackGenes_NTr.genes[[6]])),
#                                                             rep(names(msbb_array19.corr_TrackGenes_NTr.genes)[7],length(msbb_array19.corr_TrackGenes_NTr.genes[[7]])),
#                                                             rep(names(msbb_array19.corr_TrackGenes_NTr.genes)[8],length(msbb_array19.corr_TrackGenes_NTr.genes[[8]])),
#                                                             rep(names(msbb_array19.corr_TrackGenes_NTr.genes)[9],length(msbb_array19.corr_TrackGenes_NTr.genes[[9]])),
#                                                             rep(names(msbb_array19.corr_TrackGenes_NTr.genes)[10],length(msbb_array19.corr_TrackGenes_NTr.genes[[10]])),
#                                                             rep(names(msbb_array19.corr_TrackGenes_NTr.genes)[11],length(msbb_array19.corr_TrackGenes_NTr.genes[[11]])),
#                                                             rep(names(msbb_array19.corr_TrackGenes_NTr.genes)[12],length(msbb_array19.corr_TrackGenes_NTr.genes[[12]])),
#                                                             rep(names(msbb_array19.corr_TrackGenes_NTr.genes)[13],length(msbb_array19.corr_TrackGenes_NTr.genes[[13]])),
#                                                             rep(names(msbb_array19.corr_TrackGenes_NTr.genes)[14],length(msbb_array19.corr_TrackGenes_NTr.genes[[14]])),
#                                                             rep(names(msbb_array19.corr_TrackGenes_NTr.genes)[15],length(msbb_array19.corr_TrackGenes_NTr.genes[[15]])),
#                                                             rep(names(msbb_array19.corr_TrackGenes_NTr.genes)[16],length(msbb_array19.corr_TrackGenes_NTr.genes[[16]])),
#                                                             rep(names(msbb_array19.corr_TrackGenes_NTr.genes)[17],length(msbb_array19.corr_TrackGenes_NTr.genes[[17]]))),
#                                                     stringsAsFactors = F)
# msbb_array19.AggCorr.NTrdf=as.data.frame(msbb_array19.corr_NTr_PP)
# msbb_array19.AggCorr.NTrdf=msbb_array19.AggCorr.NTrdf[-unique(unname(unlist(apply(apply(msbb_array19.AggCorr.NTrdf,1,is.na),1,function(x)which(x==T))))),]
# msbb_array19.AggCorr.PLQdf=as.data.frame(msbb_array19.corr_PLQ_PP)
# msbb_array19.AggCorr.PLQdf=msbb_array19.AggCorr.PLQdf[-unique(unname(unlist(apply(apply(msbb_array19.AggCorr.PLQdf,1,is.na),1,function(x)which(x==T))))),]
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
# msbb_array19.AggCorr_Pathways=lapply(msbb_array19.AggCorr,rownames)
# core_pathways=union(union(union(union(union(union(union(union(union(union(union(union(union(union(union(union(msbb_array19.AggCorr_Pathways[[1]],
#                                                                                                               msbb_array19.AggCorr_Pathways[[2]]),
#                                                                                                         msbb_array19.AggCorr_Pathways[[3]]),
#                                                                                                   msbb_array19.AggCorr_Pathways[[4]]),
#                                                                                             msbb_array19.AggCorr_Pathways[[5]]),
#                                                                                       msbb_array19.AggCorr_Pathways[[6]]),
#                                                                                 msbb_array19.AggCorr_Pathways[[7]]),
#                                                                           msbb_array19.AggCorr_Pathways[[8]]),
#                                                                     msbb_array19.AggCorr_Pathways[[9]]),
#                                                               msbb_array19.AggCorr_Pathways[[10]]),
#                                                         msbb_array19.AggCorr_Pathways[[11]]),
#                                                   msbb_array19.AggCorr_Pathways[[12]]),
#                                             msbb_array19.AggCorr_Pathways[[13]]),
#                                       msbb_array19.AggCorr_Pathways[[14]]),
#                                 msbb_array19.AggCorr_Pathways[[15]]),
#                           msbb_array19.AggCorr_Pathways[[16]]),
#                     msbb_array19.AggCorr_Pathways[[17]])
# 
# 
# 
# msbb_array19.AggCorr_PLQ_filtIndices=lapply(lapply(msbb_array19.AggCorr,`[`,2),function(x)which(x > 0.25|x < -0.25))
# msbb_array19.AggCorr_NTr_filtIndices=lapply(lapply(msbb_array19.AggCorr,`[`,1),function(x)which(x > 0.25|x < -0.25))
# 
# for (i in 1:length(names(msbb_array19.AggCorr))){
#   
# }
# for (j in 1:length(names(msbb_array19.AggCorr))){
#   tiff(filename = paste("msbb_array19_NTr_PLQ",names(msbb_array19.AggCorr)[[j]],"PathPrint.tiff",sep="_"))
#   scatterplot(NTr~PLQ,data = msbb_array19.AggCorr[[j]],
#               xlab = "PLQ",ylab = "NTr",legend.coords = "topleft",ellipse=T,
#               main= paste(names(msbb_array19.AggCorr)[[j]],"Tangles vs Plaques correlation, PathPrint results",sep=" "))
#   dev.off()
# }
 save(msbb_array19.AggCorr,file="msbb_array19_PathPrint_ContCorrAnalysis_AOD.RData")
