library(pathprint)
library(org.Hs.eg.db)
library(entropy)
library(limma)

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
msbb_array19.covariates2=msbb_array19.covariates[-grep(pattern = "X649|X740|X350|X708|X1426|X1009|X1117|X1287|X1203|X692",x = msbb_array19.covariates$BrainBank),]

msbb_array19.2=lapply(msbb_array19[-18],function(x){rownames(x)<-x$ID;x})
msbb_array19.2.agg=lapply(msbb_array19.2[-18],function(y)aggregate(x=y[,-c(1:4)],by=list(EntrezID=y$ENTREZ_GENE_ID),mean))
msbb_array19.2.agg2=lapply(msbb_array19.2.agg,function(x){rownames(x)<- x$EntrezID;x <- x[,-1];x})
msbb_array19_allSamples_Expression=do.call("cbind",msbb_array19.2.agg2)

msbb_array19_PLQ_Strat=msbb_array19_PLQ_Strat.SCE=vector(mode = "list",length = 5)
names(msbb_array19_PLQ_Strat)=names(msbb_array19_PLQ_Strat.SCE)=c("S1","S2","S3","S4","S5")
msbb_array19_PLQ_Strat$S1=msbb_array19.covariates2$BrainBank[msbb_array19.covariates2$PLQ_Mn>=1&msbb_array19.covariates2$PLQ_Mn<=10]
msbb_array19_PLQ_Strat$S2=msbb_array19.covariates2$BrainBank[msbb_array19.covariates2$PLQ_Mn>=11&msbb_array19.covariates2$PLQ_Mn<=15]
msbb_array19_PLQ_Strat$S3=msbb_array19.covariates2$BrainBank[msbb_array19.covariates2$PLQ_Mn>=16]
# msbb_array19_PLQ_Strat$S4=msbb_array19.covariates2$BrainBank[msbb_array19.covariates2$PLQ_Mn>=16&msbb_array19.covariates2$PLQ_Mn<=20]
# msbb_array19_PLQ_Strat$S5=msbb_array19.covariates2$BrainBank[msbb_array19.covariates2$PLQ_Mn>=21]

for(i in 3:1){
  msbb_array19_PLQ_Strat.SCE[[i]]=single.chip.enrichment2(exprs = msbb_array19_allSamples_Expression[,grep(pattern = paste(msbb_array19_PLQ_Strat[[i]],collapse = "|"),x = colnames(msbb_array19_allSamples_Expression))],geneset = pathprint.Hs.gs,statistic = "mean",progressBar = T)
}
#Precomputed on Sharc
msbb_array19_PLQ_Strat.SCE=readRDS("msbb_array19_PLQ_Strat_SCE.RDS")
msbb_array19_PLQ_Strat.diffPathways=matrix(NA,nrow=length(msbb_array19_PLQ_Strat),ncol=length(msbb_array19_PLQ_Strat))
dep_df1=vector(mode = "list",length = length(msbb_array19_PLQ_Strat))
for(t1 in 1:dim(msbb_array19_PLQ_Strat.diffPathways)[2]){
  for(t2 in 1:dim(msbb_array19_PLQ_Strat.diffPathways)[2]){
    group=factor(c(rep("S1",length=length(colnames(msbb_array19_PLQ_Strat.SCE[[t1]]))),(rep("S2",length=length(colnames(msbb_array19_PLQ_Strat.SCE[[t2]]))))))
    design_df=model.matrix(~0+group)
    colnames(design_df)=levels(group)
    contrasts_matrix=makeContrasts(S1-S2,levels = design_df)
    fit=lmFit(object = cbind(msbb_array19_PLQ_Strat.SCE[[t1]],msbb_array19_PLQ_Strat.SCE[[t2]]),design = design_df)
    fit2=contrasts.fit(fit,contrasts_matrix)
    fit2=eBayes(fit2)
    fit2=eBayes(fit2,trend = T)
    dep_df1[[t1]]=topTable(fit2,coef = 1,sort.by = "logFC",adjust.method = "BH",p.value = 0.05,number = 633)
    msbb_array19_PLQ_Strat.diffPathways[t1,t2]=length(rownames(dep_df1[[t1]][which(dep_df1[[t1]]$adj.P.Val<=0.05),]))
    
  }
  
}
msbb_array19_PLQ_Strat.diffPathways=round(msbb_array19_PLQ_Strat.diffPathways/633,digits = 4)
rownames(msbb_array19_PLQ_Strat.diffPathways)=colnames(msbb_array19_PLQ_Strat.diffPathways)=c("S1-PLQ<5","S2-PLQ>=6,<=10","S3-PLQ>=11,<=15","S4-PLQ>=16,<=20","S5-PLQ>=21")
msbb_array19_PLQ_Strat.anno_df=data.frame(Mean_Plaque=c(sum(msbb_array19.covariates2$PLQ_Mn[msbb_array19.covariates2$PLQ_Mn<=5]),
                                                        sum(msbb_array19.covariates2$PLQ_Mn[msbb_array19.covariates2$PLQ_Mn>=6&msbb_array19.covariates2$PLQ_Mn<=10]),
                                                        sum(msbb_array19.covariates2$PLQ_Mn[msbb_array19.covariates2$PLQ_Mn>=11&msbb_array19.covariates2$PLQ_Mn<=15]),
                                                        sum(msbb_array19.covariates2$PLQ_Mn[msbb_array19.covariates2$PLQ_Mn>=16&msbb_array19.covariates2$PLQ_Mn<=20]),
                                                        sum(msbb_array19.covariates2$PLQ_Mn[msbb_array19.covariates2$PLQ_Mn>=21])))

pheatmap(msbb_array19_PLQ_Strat.diffPathways,cluster_rows = T,cluster_cols = T,annotation = msbb_array19_PLQ_Strat.anno_df,cellwidth = 50,cellheight = 50,fontsize = 15)
