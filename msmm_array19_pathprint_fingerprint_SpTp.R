library(pathprint)
library(pheatmap)
library(org.Hs.eg.db)
library(entropy)

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

msbb_array19.2=lapply(msbb_array19[-18],function(x){rownames(x)<-x$ID;x})
msbb_array19.2.agg=lapply(msbb_array19.2[-18],function(y)aggregate(x=y[,-c(1:4)],by=list(EntrezID=y$ENTREZ_GENE_ID),mean))
msbb_array19.2.agg2=lapply(msbb_array19.2.agg,function(x){rownames(x)<- x$EntrezID;x <- x[,-1];x})
msbb_array19_allSamples_Expression=do.call("cbind",msbb_array19.2.agg2)
plq_mean_values=sort(unique(msbb_array19.covariates$PLQ_Mn))
names(plq_mean_values)=c(1:length(plq_mean_values))
msbb_array19_SpTp.SCE=vector(mode = "list",length = length(plq_mean_values))
msbb_array19_BrainRegion.SCE=vector(mode = "list",length = 17)
names(msbb_array19_BrainRegion.SCE)=names(msbb_array19)

msbb_array19_BrainRegion.SCE=lapply(msbb_array19.2.agg2,single.chip.enrichment2,geneset = pathprint.Hs.gs,statistic = "mean",normalizedScore = F)
msbb_array19_BrainRegion.diffPathways=vector(mode = "list",length = 633)
t_matrix=matrix(NA,nrow=633,ncol=17)
rownames(t_matrix)=rownames(msbb_array19_BrainRegion.SCE$AC)
colnames(t_matrix)=names(msbb_array19)
for(i in 1:633){
  for(j in 1:633){
    t_matrix[i,]=unlist(mapply(FUN = function(x,y)t.test(x,y,paired = F)$p.value,lapply(msbb_array19_BrainRegion.SCE,`[`,i,),lapply(msbb_array19_BrainRegion.SCE,`[`,j,)))
    msbb_array19_BrainRegion.diffPathways[[i]]=t_matrix
    
  }
}
#length(plq_mean_values)
ns=c(37,38,69,72,75,81,5,13,15,23,83,89,46,47,68,78,80,76)
msbb_array19_SpTp.SCE=vector(mode = "list",length = length(names(plq_mean_values[-ns])))
names(msbb_array19_SpTp.SCE)=names(plq_mean_values[-ns])
for(i in as.numeric(names(msbb_array19_SpTp.SCE))){
  search_samples=paste(msbb_array19.covariates$BrainBank[which(msbb_array19.covariates$PLQ_Mn%in%plq_mean_values[i])],collapse = "|")
  if(length(grep(pattern = search_samples,colnames(msbb_array19_allSamples_Expression)))==0){
    cat(paste(search_samples," not profiled and hence dropped!\n",sep = " "))
    msbb_array19_SpTp.SCE[[i]]=0
  }
  else if(length(grep(pattern = search_samples,colnames(msbb_array19_allSamples_Expression)))>0){
    cat(paste("Computing pathway fingerprint for ",length(grep(pattern = search_samples,colnames(msbb_array19_allSamples_Expression)))," samples & for mean plaque ",plq_mean_values[i],"\n",sep = ""))
    msbb_array19_SpTp.SCE[[i]]=single.chip.enrichment2(exprs = msbb_array19_allSamples_Expression[,grep(pattern = search_samples,colnames(msbb_array19_allSamples_Expression))],geneset = pathprint.Hs.gs,statistic = "mean",normalizedScore = F)  
  }
  
}






