#library(limma)
library(Biobase)
#library(jetset)
library(pathprint)
library(metaArray)
library(doMC)
library(pheatmap)
library(org.Hs.eg.db)

setwd("/shared/hidelab2/user/md4zsa/Work/Data/MSMM_Array/MSBB_Array19/Normalised_Data")
single.chip.enrichment2=function(exprs, geneset, statistic = "mean",normalizedScore = FALSE, progressBar = TRUE){
  # if (!(transformation %in% c("rank", "squared.rank", "log.rank"))) 
  #   stop("transformation should be rank, squared.rank or log.rank")
  # if (!(statistic %in% c("mean", "median"))) 
  #   stop("transformation should be mean or median")
  # if ((normalizedScore == TRUE & !(statistic == "mean"))) 
  #   stop("Parameteric normalization can only be used for statistic = mean")
  Ns <- ncol(exprs)
  Ng <- nrow(exprs)
  gene.names <- rownames(exprs)
  geneset.names <- names(geneset)
  # exprs <- apply(exprs, 2, rank, ties.method = "average")
  # if (transformation == "log.rank") {
  #   exprs <- log(exprs)
  # }
  # else if (transformation == "squared.rank") {
  #   exprs <- exprs^2
  # }
   if (progressBar == TRUE) {
     pb <- txtProgressBar(min = 0, max = length(geneset), 
                          style = 3)
   }
  score.matrix <- matrix(0, nrow = length(geneset), ncol = Ns)
  for (i in 1:length(geneset)) {
    overlap <- intersect(geneset[[i]], gene.names)
    if (length(overlap) == 0) {
      score.matrix[i, ] <- NA
    }
    else {
      if (statistic == "mean") {
        score.matrix[i, ] <- apply(exprs, 2, function(x) {
          mean(x[overlap])
        })
        # if (normalizedScore == TRUE) {  
        # 
        #   n <- length(overlap)
        #   if (transformation == "rank") {
        #     E.mean <- mean(1:Ng)
        #     E.sd <- ((sd(1:Ng)/(n^0.5))) * (((Ng - n)/(Ng - 
        #                                                  1))^0.5)
        #   }
        #   else if (transformation == "log.rank") {
        #     E.mean <- mean(log(1:Ng))
        #     E.sd <- ((sd(log(1:Ng))/(n^0.5))) * (((Ng - 
        #                                              n)/(Ng - 1))^0.5)
        #   }
        #   else if (transformation == "squared.rank") {
        #     E.mean <- mean((1:Ng)^2)
        #     E.sd <- ((sd((1:Ng)^2)/(n^0.5))) * (((Ng - 
        #                                             n)/(Ng - 1))^0.5)
        #   }
        #   score.matrix[i, ] <- sapply(score.matrix[i, 
        #                                            ], pnorm, mean = E.mean, sd = E.sd) - 0.5
        # }
      }
      else if (statistic == "median") {
        score.matrix[i, ] <- apply(exprs, 2, function(x) {
          median(x[overlap])
        })
      }
    }
    if (progressBar == TRUE) {
      setTxtProgressBar(pb, i)
    }
  }
  colnames(score.matrix) <- colnames(exprs)
  rownames(score.matrix) <- geneset.names
  return(score.matrix)
}
remove_ZeroSumPathways=function(x,y)list(x[-y,])
getEntropy <- function(mat, index){
  if (index > 2 | index < 1)
    stop("Indicate 1 for rows or 2 for columns")
  d <- apply(as.matrix(mat), index, function(x){discretize(x, numBins = 3, r=c(-1,1))})
  entropy.vec <- apply(d, 2, entropy)
  return(entropy.vec)
}
ds = c("chipframe", "genesets","pathprint.Hs.gs" ,"platform.thresholds")
data(list = ds)

msbb_array19.files=list.files(pattern = "*.tsv",full.names = T)
msbb_array19.files=msbb_array19.files[-c(18:19)]
msbb_array19=mapply(FUN = read.table,msbb_array19.files,MoreArgs = list(header=T,sep="\t",as.is=T))
names(msbb_array19)=c("AC","CN","DLPFC","FP","HP","IFG","ITG","MTG","OVC","PHG","PCC","PCG","PFC","PTMN","SPL","STG","TP")
msbb_array19.covariates=read.delim2("AMP-AD_MSBB_MSSM_covariates_mRNA_AffymetrixU133AB.txt",header = T,as.is = T)
msbb_array19.covariates$Age[which(msbb_array19.covariates$Age=="89+")]=90 #For simplicity, all samples with AOD as 89+ were labeled 90
msbb_array19.covariates$Age=as.integer(msbb_array19.covariates$Age)
msbb_array19.covariates$pH=as.numeric(msbb_array19.covariates$pH)
msbb_array19.covariates$CDR=as.numeric(msbb_array19.covariates$CDR)
msbb_array19.covariates$PLQ_Mn=as.numeric(msbb_array19.covariates$PLQ_Mn)
msbb_array19.covariates$BrainBank=paste("X",msbb_array19.covariates$BrainBank,sep="")

msbb_array19.2=lapply(msbb_array19,function(x){rownames(x)<-x$ID;x})

msbb_array19.2_PLQ_Strat=msbb_array19.fingerprint_PLQ_Strat=msbb_array19_fingerprint.ZeroSumPathways_PLQ_Strat=msbb_array19_fingerprint.noZeromSumPathways_PLQ_Strat=msbb_array19_fingerprint.HighEntropy05_Indices_PLQ_Strat=msbb_array19_fingerprint.HighEntropy05_Fingerprint_PLQ_Strat=msbb_array19_fingerprint_HighEntropy05_Fingerprint.Core_PLQ_Strat=msbb_array19_fingerprint_HighEntropy05.CoreFingerprint_PLQ_Strat=vector(mode = "list",length = 3)
names(msbb_array19.2_PLQ_Strat)=names(msbb_array19.fingerprint_PLQ_Strat)=names(msbb_array19_fingerprint.ZeroSumPathways_PLQ_Strat)=names(msbb_array19_fingerprint.noZeromSumPathways_PLQ_Strat)=names(msbb_array19_fingerprint.HighEntropy05_Indices_PLQ_Strat)=names(msbb_array19_fingerprint.HighEntropy05_Fingerprint_PLQ_Strat)=names(msbb_array19_fingerprint_HighEntropy05_Fingerprint.Core_PLQ_Strat)=names(msbb_array19_fingerprint_HighEntropy05.CoreFingerprint_PLQ_Strat)=c(1:3)
msbb_array19.2_PLQ_Strat_SCE=msbb_array19.fingerprint_PLQ_Strat_SCE=msbb_array19_fingerprint.ZeroSumPathways_PLQ_Strat_SCE=msbb_array19_fingerprint.noZeromSumPathways_PLQ_Strat_SCE=msbb_array19_fingerprint.HighEntropy05_Indices_PLQ_Strat_SCE=msbb_array19_fingerprint.HighEntropy05_Fingerprint_PLQ_Strat_SCE=msbb_array19_fingerprint_HighEntropy05_Fingerprint.Core_PLQ_Strat_SCE=msbb_array19_fingerprint_HighEntropy05.CoreFingerprint_PLQ_Strat_SCE=vector(mode = "list",length = 3)
names(msbb_array19.2_PLQ_Strat_SCE)=names(msbb_array19.fingerprint_PLQ_Strat_SCE)=names(msbb_array19_fingerprint.ZeroSumPathways_PLQ_Strat_SCE)=names(msbb_array19_fingerprint.noZeromSumPathways_PLQ_Strat_SCE)=names(msbb_array19_fingerprint.HighEntropy05_Indices_PLQ_Strat_SCE)=names(msbb_array19_fingerprint.HighEntropy05_Fingerprint_PLQ_Strat_SCE)=names(msbb_array19_fingerprint_HighEntropy05_Fingerprint.Core_PLQ_Strat_SCE)=names(msbb_array19_fingerprint_HighEntropy05.CoreFingerprint_PLQ_Strat_SCE)=c(1:3)

msbb_array19.2_PLQ_Strat[[1]]=msbb_array19.covariates[which(msbb_array19.covariates$PLQ_Mn>=0&msbb_array19.covariates$PLQ_Mn<=10),1]
msbb_array19.2_PLQ_Strat[[2]]=msbb_array19.covariates[which(msbb_array19.covariates$PLQ_Mn>=11&msbb_array19.covariates$PLQ_Mn<=15),1]
msbb_array19.2_PLQ_Strat[[3]]=msbb_array19.covariates[which(msbb_array19.covariates$PLQ_Mn>=16),1]


multiID.u133a=msbb_array19.2[[1]]$ENTREZ_GENE_ID[grep(pattern="///",msbb_array19.2[[1]]$ENTREZ_GENE_ID)]
u133a_universe=c(msbb_array19.2[[1]]$ENTREZ_GENE_ID[-grep(pattern="///",msbb_array19.2[[1]]$ENTREZ_GENE_ID)],gsub(pattern=" ",replacement="",lapply(strsplit(x=multiID.u133a,split="///"),`[`,1)))
msbb_array19.2=lapply(msbb_array19.2,function(x){x$ENTREZ_GENE_ID<-u133a_universe;x})
msbb_array19.2.agg=lapply(msbb_array19.2,function(y)aggregate(x=y[,-c(1:4)],by=list(EntrezID=y$ENTREZ_GENE_ID),mean))
msbb_array19.2.agg2=lapply(msbb_array19.2.agg,function(x){rownames(x)<- x$EntrezID;x <- x[,-1];x})
for (p in 1:3){
  msbb_array19.fingerprint_PLQ_Strat[[p]]=lapply(lapply(msbb_array19.2,function(x)x[,which(colnames(x)%in%msbb_array19.2_PLQ_Strat[[p]])]),exprs2fingerprint,platform = "GPL96",species = "human",progressBar = T)
  #msbb_array19.fingerprint_PLQ_Strat_SCE[[p]]=lapply(lapply(msbb_array19.2.agg2,function(x)x[,which(colnames(x)%in%msbb_array19.2_PLQ_Strat[[p]])]),single.chip.enrichment,geneset = pathprint.Hs.gs,transformation="log.rank",statistic="mean",normalizedScore=F,progressBar=T)
  msbb_array19_fingerprint.ZeroSumPathways_PLQ_Strat[[p]]=lapply(msbb_array19.fingerprint_PLQ_Strat[[p]],function(x)which(rowSums(x)==0))
  #msbb_array19_fingerprint.ZeroSumPathways_PLQ_Strat_SCE[[p]]=lapply(msbb_array19.fingerprint_PLQ_Strat_SCE[[p]],function(x)which(rowSums(x)==0))
  msbb_array19_fingerprint.noZeromSumPathways_PLQ_Strat[[p]]=mapply(remove_ZeroSumPathways,msbb_array19.fingerprint_PLQ_Strat[[p]],msbb_array19_fingerprint.ZeroSumPathways_PLQ_Strat[[p]])
  msbb_array19_fingerprint.HighEntropy05_Indices_PLQ_Strat[[p]]=lapply(lapply(msbb_array19_fingerprint.noZeromSumPathways_PLQ_Strat[[p]],getEntropy,1),function(x)unname(which(x>0.5)))
  msbb_array19_fingerprint.HighEntropy05_Fingerprint_PLQ_Strat[[p]]=mapply(FUN = function(x,y)x[y,],msbb_array19_fingerprint.noZeromSumPathways_PLQ_Strat[[p]],msbb_array19_fingerprint.HighEntropy05_Indices_PLQ_Strat[[p]])
  msbb_array19_fingerprint_HighEntropy05_Fingerprint.Core_PLQ_Strat[[p]]=Reduce(intersect,lapply(msbb_array19_fingerprint.HighEntropy05_Fingerprint_PLQ_Strat[[p]],rownames))
  if(p==3){
    region=names(lapply(lapply(msbb_array19.2,function(x)x[,which(colnames(x)%in%msbb_array19.2_PLQ_Strat$`3`)]),function(x)dim(x)[2]))
    msbb_array19.fingerprint_PLQ_Strat[[p]]=lapply(lapply(msbb_array19.2[region],function(x)x[,which(colnames(x)%in%msbb_array19.2_PLQ_Strat[[p]])]),exprs2fingerprint,platform = "GPL96",species = "human",progressBar = T)
    msbb_array19.fingerprint_PLQ_Strat_SCE[[p]]=lapply(lapply(msbb_array19.2.agg2[region],function(x)x[,which(colnames(x)%in%msbb_array19.2_PLQ_Strat[[p]])]),single.chip.enrichment,geneset = pathprint.Hs.gs,transformation="log.rank",statistic="mean",normalizedScore=F,progressBar=T)
    msbb_array19_fingerprint.ZeroSumPathways_PLQ_Strat[[p]]=lapply(msbb_array19.fingerprint_PLQ_Strat[[p]],function(x)which(rowSums(x)==0))
    msbb_array19_fingerprint.noZeromSumPathways_PLQ_Strat[[p]]=mapply(remove_ZeroSumPathways,msbb_array19.fingerprint_PLQ_Strat[[p]],msbb_array19_fingerprint.ZeroSumPathways_PLQ_Strat[[p]])
    msbb_array19_fingerprint.HighEntropy05_Indices_PLQ_Strat[[p]]=lapply(lapply(msbb_array19_fingerprint.noZeromSumPathways_PLQ_Strat[[p]],getEntropy,1),function(x)unname(which(x>0.5)))
    msbb_array19_fingerprint.HighEntropy05_Fingerprint_PLQ_Strat[[p]]=mapply(FUN = function(x,y)x[y,],msbb_array19_fingerprint.noZeromSumPathways_PLQ_Strat[[p]],msbb_array19_fingerprint.HighEntropy05_Indices_PLQ_Strat[[p]])
    
  }
  
}
plq_strata_region_corefingerprint=vector(mode = "list",length = 17)
names(plq_strata_region_corefingerprint)=names(msbb_array19)
for(i in 1:17){
  plq_strata_region_corefingerprint[[i]]=Reduce(intersect,list(lapply(msbb_array19_fingerprint.HighEntropy05_Fingerprint_PLQ_Strat$`1`,rownames)[[i]],
                                                          lapply(msbb_array19_fingerprint.HighEntropy05_Fingerprint_PLQ_Strat$`2`,rownames)[[i]],
                                                          lapply(msbb_array19_fingerprint.HighEntropy05_Fingerprint_PLQ_Strat$`3`,rownames)[[i]]))
  
}

diffPathways_plq_strat=matrix(NA,nrow = 17,ncol = 9)
for(i in 1:17){
  plq_strata_region_corefingerprint=Reduce(intersect,list(lapply(msbb_array19_fingerprint.HighEntropy05_Fingerprint_PLQ_Strat$`1`,rownames)[[i]],
                                                          lapply(msbb_array19_fingerprint.HighEntropy05_Fingerprint_PLQ_Strat$`2`,rownames)[[i]],
                                                          lapply(msbb_array19_fingerprint.HighEntropy05_Fingerprint_PLQ_Strat$`3`,rownames)[[i]]))
  
  dfP_s12=diffPathways(fingerprints = cbind(msbb_array19.fingerprint_PLQ_Strat$`1`[[i]],msbb_array19.fingerprint_PLQ_Strat$`2`[[i]]),fac = c(rep("Strat1",length(colnames(msbb_array19.fingerprint_PLQ_Strat$`1`[[i]]))),rep("Strat2",length(colnames(msbb_array19.fingerprint_PLQ_Strat$`2`[[i]])))),threshold = 0.5)
  dfP_s23=diffPathways(fingerprints = cbind(msbb_array19.fingerprint_PLQ_Strat$`2`[[i]],msbb_array19.fingerprint_PLQ_Strat$`3`[[i]]),fac = c(rep("Strat2",length(colnames(msbb_array19.fingerprint_PLQ_Strat$`2`[[i]]))),rep("Strat3",length(colnames(msbb_array19.fingerprint_PLQ_Strat$`3`[[i]])))),threshold = 0.5)
  dfP_s13=diffPathways(fingerprints = cbind(msbb_array19.fingerprint_PLQ_Strat$`1`[[i]],msbb_array19.fingerprint_PLQ_Strat$`3`[[i]]),fac = c(rep("Strat1",length(colnames(msbb_array19.fingerprint_PLQ_Strat$`1`[[i]]))),rep("Strat3",length(colnames(msbb_array19.fingerprint_PLQ_Strat$`3`[[i]])))),threshold = 0.5)
  diffPathways_plq_strat[i,]=c(names(msbb_array19)[i],length(dfP_s12),paste(dfP_s12,collapse = ";"),length(dfP_s23),paste(dfP_s23,collapse = ";"),length(dfP_s13),paste(dfP_s13,collapse = ";"),length(plq_strata_region_corefingerprint),paste(plq_strata_region_corefingerprint,collapse = ";"))
  
}
colnames(diffPathways_plq_strat)=c("BrainRegion","#DiffPathways_Strata12","DiffPathways_Strata12","#DiffPathways_Strata23","DiffPathways_Strata23","#DiffPathways_Strata13","DiffPathways_Strata13","#CommonPathways_All_PLQ_Strata","CommonPathways_All_PLQ_Strata")
write.table(data.frame(diffPathways_plq_strat,stringsAsFactors = F),"MSBB_Diff_Common_PLQ_Strata_Pathprint_Summary.txt",sep = "\t",col.names = T,row.names = F,quote = F)

#Pan brain core pathways, by strata
CoreFingerprint.df=CoreFingerprint.colnames_length=CoreFingerprint.anno_df=vector(mode = "list",length = 3)
names(CoreFingerprint.df)=names(CoreFingerprint.colnames_length)=names(CoreFingerprint.anno_df)=c("Strata1","Strata2","Strata3")
for (p in c(1,3)){
  CoreFingerprint.df[[p]]=do.call("cbind",lapply(msbb_array19_fingerprint.HighEntropy05_Fingerprint_PLQ_Strat[[p]],function(x)x[which(rownames(x)%in%msbb_array19_fingerprint_HighEntropy05_Fingerprint.Core_PLQ_Strat[[p]]),]))  
  CoreFingerprint.colnames_length[[p]]=lapply(msbb_array19_fingerprint.HighEntropy05_Fingerprint_PLQ_Strat[[p]],function(x)length(colnames(x[which(rownames(x)%in%msbb_array19_fingerprint_HighEntropy05_Fingerprint.Core_PLQ_Strat[[p]]),])))
  CoreFingerprint.anno_df[[p]]=matrix(NA,nrow = length(colnames(CoreFingerprint.df[[p]])),ncol = 1)
  tmp1=c()
  for (i in 1:17){
    tmp=rep(names(CoreFingerprint.colnames_length[[p]])[i],length=CoreFingerprint.colnames_length[[p]][[i]])
    tmp1=append(tmp1,values = tmp)
    
  }
  CoreFingerprint.anno_df[[p]][,1]=tmp1
  rownames(CoreFingerprint.anno_df[[p]])=paste(colnames(CoreFingerprint.df[[p]]),tmp1,sep = ".")
  colnames(CoreFingerprint.anno_df[[p]])=c("Brain-region")
  colnames(CoreFingerprint.df[[p]])=paste(colnames(CoreFingerprint.df[[p]]),tmp1,sep = ".")
  CoreFingerprint.anno_df[[p]]=data.frame(CoreFingerprint.anno_df[[p]],stringsAsFactors = F)
  CoreFingerprint.df[[p]]=data.frame(CoreFingerprint.df[[p]],stringsAsFactors = F)
  
}

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
msbb_array19.NTr_PLQ_phenoData=msbb_array19.covariates[,c(1,6:7,10,12)]
msbb_array19.phenoDataIndices=lapply(lapply(msbb_array19.2,colnames),function(x)which(msbb_array19.NTr_PLQ_phenoData$BrainBank%in%x))
msbb_array19.NTr_PLQ_phenoData2=lapply(msbb_array19.phenoDataIndices,function(x)msbb_array19.NTr_PLQ_phenoData[x,])
msbb_array19.NTr_PLQ_phenoData2=lapply(msbb_array19.NTr_PLQ_phenoData2,function(x){rownames(x)<-x$BrainBank;x})
msbb_array19.NTr_PLQ_phenoData2=lapply(msbb_array19.NTr_PLQ_phenoData2,`[`,c(-1))

msbb_array19.eset=msbb_array19.DEG_PLQ=msbb_array19.DEG_NTR=list()
control.CDR=lapply(msbb_array19.NTr_PLQ_phenoData2,function(x)which(x$PLQ_Mn<=0))
disease.CDR=lapply(msbb_array19.NTr_PLQ_phenoData2,function(x)which(x$PLQ_Mn>=15))
for (i in 1:length(names(msbb_array19.NTr_PLQ_phenoData2))){
  design.df=matrix(NA,nrow=sum(length(control.CDR[[i]]),length(disease.CDR[[i]])),ncol=2)
  rownames(design.df)=c(rownames(msbb_array19.NTr_PLQ_phenoData2[[i]][control.CDR[[i]],]),rownames(msbb_array19.NTr_PLQ_phenoData2[[i]][disease.CDR[[i]],]))
  colnames(design.df)=c("control.CDR","disease.CDR")
  design.df[which(rownames(design.df)%in%rownames(msbb_array19.NTr_PLQ_phenoData2[[i]][control.CDR[[i]],])),1]=1
  design.df[-which(rownames(design.df)%in%rownames(msbb_array19.NTr_PLQ_phenoData2[[i]][control.CDR[[i]],])),1]=0
  design.df[which(rownames(design.df)%in%rownames(msbb_array19.NTr_PLQ_phenoData2[[i]][disease.CDR[[i]],])),2]=1
  design.df[-which(rownames(design.df)%in%rownames(msbb_array19.NTr_PLQ_phenoData2[[i]][disease.CDR[[i]],])),2]=0
  phenoData=AnnotatedDataFrame(data=msbb_array19.NTr_PLQ_phenoData2[[i]][which(rownames(msbb_array19.NTr_PLQ_phenoData2[[i]])%in%rownames(design.df)),])
  msbb_array19.eset[[i]]=ExpressionSet(assayData=as.matrix(msbb_array19.2[[i]][,which(colnames(msbb_array19.2[[i]])%in%rownames(design.df))]),phenoData=phenoData)
  fit=lmFit(msbb_array19.eset[[i]],design=design.df)
  fit=eBayes(fit)
  contMatrix=makeContrasts(CtrlvsDisease=disease.CDR-control.CDR,levels=design.df)
  fit2=contrasts.fit(fit,contMatrix)
  fit2=eBayes(fit2)
  msbb_array19.DEG_PLQ[[i]]=topTableF(fit2,adjust.method="BH",number=Inf)
}
names(msbb_array19.DEG_PLQ)=names(msbb_array19)

#Brain marker genes specific Pathprint
zhang_celltype_ADgenes=read.xls('../../../BrainExpression_Datasets/Zhang_19BrainRegions_Paper/Zhang_BrainCelltype_Markers.xlsx',skip=1,sheet=3,header=T,as.is=T)
zhang_celltype_ADgenes.list=zhang_celltype_PLQ_ADgenes.list=vector(mode = 'list',length = 5)
names(zhang_celltype_ADgenes.list)=names(zhang_celltype_PLQ_ADgenes.list)=sort(unique(zhang_celltype_ADgenes$Cell.type))
zhang_celltype_ADgenes.list$Astrocytes=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_ADgenes.list)[1],zhang_celltype_ADgenes$Cell.type)]
zhang_celltype_ADgenes.list$Endothelial=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_ADgenes.list)[2],zhang_celltype_ADgenes$Cell.type)]
zhang_celltype_ADgenes.list$Microglia=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_ADgenes.list)[3],zhang_celltype_ADgenes$Cell.type)]
zhang_celltype_ADgenes.list$Neurons=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_ADgenes.list)[4],zhang_celltype_ADgenes$Cell.type)]
zhang_celltype_ADgenes.list$Oligodendrocytes=zhang_celltype_ADgenes$Gene.symbol[grep(pattern = names(zhang_celltype_ADgenes.list)[5],zhang_celltype_ADgenes$Cell.type)]

hgu133a_probe_symbol.Map=select(x = hgu133a.db,keys = rownames(msbb_array19.2$AC),keytype = "PROBEID",columns = "SYMBOL")
zhang_celltype_ADgenes.probes=hgu133a_probe_symbol.Map$PROBEID[which(hgu133a_probe_symbol.Map$SYMBOL%in%zhang_celltype_ADgenes$Gene.symbol)]
msbb_array19.celltype=lapply(msbb_array19.2,function(x)x[which(rownames(x)%in%zhang_celltype_ADgenes.probes),])

msbb_array19_celltype.fingerprint_PLQ_Strat=msbb_array19_celltype_fingerprint.ZeroSumPathways_PLQ_Strat=msbb_array19_celltype_fingerprint.noZeromSumPathways_PLQ_Strat=msbb_array19_celltype_fingerprint.HighEntropy05_Indices_PLQ_Strat=msbb_array19_celltype_fingerprint.HighEntropy05_Fingerprint_PLQ_Strat=vector(mode = "list",length = 3)
names(msbb_array19_celltype.fingerprint_PLQ_Strat)=names(msbb_array19_celltype_fingerprint.ZeroSumPathways_PLQ_Strat)=names(msbb_array19_celltype_fingerprint.noZeromSumPathways_PLQ_Strat)=names(msbb_array19_celltype_fingerprint.HighEntropy05_Indices_PLQ_Strat)=names(msbb_array19_celltype_fingerprint.HighEntropy05_Fingerprint_PLQ_Strat)=c(1:3)
for(p in 1:3){
  msbb_array19_celltype.fingerprint_PLQ_Strat[[p]]=lapply(lapply(msbb_array19.celltype,function(x)x[,which(colnames(x)%in%msbb_array19.2_PLQ_Strat[[p]])]),exprs2fingerprint,platform = "GPL96",species = "human",progressBar = T)
  msbb_array19_celltype_fingerprint.ZeroSumPathways_PLQ_Strat[[p]]=lapply(msbb_array19_celltype.fingerprint_PLQ_Strat[[p]],function(x)which(rowSums(x)==0))
  msbb_array19_celltype_fingerprint.noZeromSumPathways_PLQ_Strat[[p]]=mapply(remove_ZeroSumPathways,msbb_array19_celltype.fingerprint_PLQ_Strat[[p]],msbb_array19_celltype_fingerprint.ZeroSumPathways_PLQ_Strat[[p]])
  msbb_array19_celltype_fingerprint.HighEntropy05_Indices_PLQ_Strat[[p]]=lapply(lapply(msbb_array19_celltype_fingerprint.noZeromSumPathways_PLQ_Strat[[p]],getEntropy,1),function(x)unname(which(x>0.5)))
  msbb_array19_celltype_fingerprint.HighEntropy05_Fingerprint_PLQ_Strat[[p]]=mapply(FUN = function(x,y)x[y,],msbb_array19_celltype_fingerprint.noZeromSumPathways_PLQ_Strat[[p]],msbb_array19_celltype_fingerprint.HighEntropy05_Indices_PLQ_Strat[[p]])
  
}
# 
# 
# names(msbb_array19.eset)=names(msbb_array19.DEG_PLQ)=names(msbb_array19)

#PathPrint on PLQ bins

