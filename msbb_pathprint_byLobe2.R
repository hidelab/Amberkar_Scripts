library(pathprint)
library(org.Hs.eg.db)
library(dplyr)
library(magrittr)
library(data.table)
library(parallel)
library(magrittr)
library(enrichR)
library(GEOquery)

mapIds2<-function(IDs,IDFrom,IDTo){
  require(org.Hs.eg.db)
  idmap=mapIds(x = org.Hs.eg.db,keys = IDs,column = IDTo,keytype = IDFrom,multiVals = "first")
  na_vec=names(idmap[is.na(idmap)==T])
  idmap=idmap[is.na(idmap)==F]
  idmap_df=data.frame("From"=names(idmap),"To"=unlist(unname(idmap)),stringsAsFactors = F)
  return(list(map=idmap_df,noMap=na_vec))
}
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
dbs <- listEnrichrDbs()
kegg_dbs=grep(pattern = "KEGG",x = dbs$libraryName,value = T)[3]
ds = c("chipframe", "genesets","pathprint.Hs.gs" ,"platform.thresholds")
data(list = ds)

msbb_gse84422_series_matrix=getGEO(GEO = "GSE84422",GSEMatrix = T)
msbb_gse84422_series_matrix.GPL96=msbb_gse84422_series_matrix$`GSE84422-GPL96_series_matrix.txt.gz`
msbb_gse84422_series_matrix.GPL97=msbb_gse84422_series_matrix$`GSE84422-GPL97_series_matrix.txt.gz`
msbb_gse84422_series_matrix.GPL570=msbb_gse84422_series_matrix$`GSE84422-GPL570_series_matrix.txt.gz`

msbb_gse84422.fData=vector(mode = "list",length = 3)
names(msbb_gse84422.fData)=c("GPL96","GPL97","GPL570")
msbb_gse84422.fData$GPL96=fData(msbb_gse84422_series_matrix.GPL96)[,c(1:2,11:13)]
msbb_gse84422.fData$GPL97=fData(msbb_gse84422_series_matrix.GPL97)[,c(1:2,11:13)]
msbb_gse84422.fData$GPL570=fData(msbb_gse84422_series_matrix.GPL570)[,c(1:2,11:13)]


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



#Initiate list structures for storing results
msbb_gse84422.byLobe=msbb_gse84422_byLobe.SCE=vector(mode = "list",length = length(unique(lobe_bm_area.map$Lobe)))
names(msbb_gse84422.byLobe)=names(msbb_gse84422_byLobe.SCE)=unique(lobe_bm_area.map$Lobe)
msbb_gse84422_byLobe.SCE$Frontal_Lobe=msbb_gse84422_byLobe.SCE$Occipetal_Lobe=msbb_gse84422_byLobe.SCE$Temporal_Lobe=msbb_gse84422_byLobe.SCE$Parietal_Lobe=msbb_gse84422_byLobe.SCE$Dorsal_striatum=vector(mode = "list",length = 4)
names(msbb_gse84422_byLobe.SCE$Frontal_Lobe)=names(msbb_gse84422_byLobe.SCE$Occipetal_Lobe)=names(msbb_gse84422_byLobe.SCE$Temporal_Lobe)=names(msbb_gse84422_byLobe.SCE$Parietal_Lobe)=names(msbb_gse84422_byLobe.SCE$Dorsal_striatum)=c("S0","S1","S2","S3")



