library(pathprint)
library(org.Hs.eg.db)
library(dplyr)
library(magrittr)
library(data.table)
library(data.table)
library(parallel)
ncores=8

# Sys.setenv("plotly_username"="ssamberkar")
# Sys.setenv("plotly_api_key"="U5r0L12wRGNM98tLdTKz")
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
remove_ZeroSumPathways=function(x,y)list(x[-y,])

ds = c("chipframe", "genesets","pathprint.Hs.gs" ,"platform.thresholds")
data(list = ds)

msbb_array19.files=list.files(pattern = "*.tsv",full.names = T)
msbb_array19.files=msbb_array19.files[-c(18:19)]
msbb_array19=mapply(FUN = fread,msbb_array19.files,MoreArgs = list(header=T,sep="\t",data.table=F,showProgress=T))
names(msbb_array19)=c("AC","CN","DLPFC","FP","HP","IFG","ITG","MTG","OVC","PHG","PCC","PCG","PFC","PTMN","SPL","STG","TP")
msbb_array19.covariates=read.delim2("AMP-AD_MSBB_MSSM_covariates_mRNA_AffymetrixU133AB.txt",header = T,as.is = T)%>%
  mutate(Age=replace(Age,Age=="89+",90))%>%
  mutate(pH=as.numeric(pH))%>%
  mutate(CDR=as.numeric(CDR))%>%
  mutate(PLQ_Mn=as.numeric(PLQ_Mn))%>%
  mutate(BrainBank=paste("X",BrainBank,sep=""))
                        
msbb_array19.covariates2=msbb_array19.covariates[-grep(pattern = "X649|X740|X350|X708|X1426|X1009|X1117|X1287|X1203|X692",x = msbb_array19.covariates$BrainBank),]
ns=sort(c(37,38,69,72,75,81,13,83,89))
msbb_array19.covariates3=msbb_array19.covariates2[-which(msbb_array19.covariates2$BrainBank%in%msbb_array19.covariates2$BrainBank[ns]),]

#Collapsing Affy-IDs to Entrez-IDs

msbb_array19.2=lapply(msbb_array19[-18],function(x){rownames(x)<-x$ID;x})
multiID.u133a=msbb_array19.2[[1]]$ENTREZ_GENE_ID[grep(pattern="///",msbb_array19.2[[1]]$ENTREZ_GENE_ID)]
msbb_array19.2.agg=mclapply(msbb_array19.2[-18],function(y)aggregate(x=y[,-c(1:4)],by=list(EntrezID=y$ENTREZ_GENE_ID),mean),mc.cores = ncores)
msbb_array19.2.agg2=mclapply(msbb_array19.2.agg,function(x){rownames(x)<- x$EntrezID;x <- x[-grep(pattern = "///",x = rownames(x)),-1];x},mc.cores = ncores)

plq_mean_values=sort(unique(msbb_array19.covariates$PLQ_Mn))
names(plq_mean_values)=c(1:length(plq_mean_values))

#Regroup brain regions by lobes
lobe_bm_area.map=matrix(nrow=17,ncol=2)
lobe_bm_area.map[1,]=c("Frontal_Lobe","FP")
lobe_bm_area.map[2,]=c("Frontal_Lobe","AC")
lobe_bm_area.map[3,]=c("Frontal_Lobe","PFC")
lobe_bm_area.map[4,]=c("Occipetal_Lobe","OVC")
lobe_bm_area.map[5,]=c("Temporal_Lobe","ITG")
lobe_bm_area.map[6,]=c("Temporal_Lobe","MTG")
lobe_bm_area.map[7,]=c("Temporal_Lobe","STG")
lobe_bm_area.map[8,]=c("Parietal_Lobe","PCC")
lobe_bm_area.map[9,]=c("Temporal_Lobe","PHG")
lobe_bm_area.map[10,]=c("Temporal_Lobe","TP")
lobe_bm_area.map[11,]=c("Frontal_Lobe","PCG")
lobe_bm_area.map[12,]=c("Frontal_Lobe","IFG")
lobe_bm_area.map[13,]=c("Frontal_Lobe","DLPFC")
lobe_bm_area.map[14,]=c("Parietal_Lobe","SPL")
lobe_bm_area.map[16,]=c("Temporal_Lobe","HP")
lobe_bm_area.map[15,]=c("Dorsal_striatum","PTMN")
lobe_bm_area.map[17,]=c("Dorsal_striatum","CN")
lobe_bm_area.map=data.frame(Lobe=lobe_bm_area.map[,1],BrainRegion=lobe_bm_area.map[,2],stringsAsFactors = F)

msbb_array.byLobe=msbb_array_byLobe.SCE=vector(mode = "list",length = length(unique(lobe_bm_area.map$Lobe)))
names(msbb_array.byLobe)=names(msbb_array_byLobe.SCE)=unique(lobe_bm_area.map$Lobe)
msbb_array_byLobe.SCE$Frontal_Lobe=msbb_array_byLobe.SCE$Occipetal_Lobe=msbb_array_byLobe.SCE$Temporal_Lobe=msbb_array_byLobe.SCE$Parietal_Lobe=msbb_array_byLobe.SCE$Dorsal_striatum=vector(mode = "list",length = 4)
names(msbb_array_byLobe.SCE$Frontal_Lobe)=names(msbb_array_byLobe.SCE$Occipetal_Lobe)=names(msbb_array_byLobe.SCE$Temporal_Lobe)=names(msbb_array_byLobe.SCE$Parietal_Lobe)=names(msbb_array_byLobe.SCE$Dorsal_striatum)=c("S0","S1","S2","S3")

msbb_array19.PLQ_Strat=vector(mode = "list",length = 4)
names(msbb_array19.PLQ_Strat)=c("S0","S1","S2","S3")
msbb_array19.PLQ_Strat[[1]]=msbb_array19.covariates[which(msbb_array19.covariates$PLQ_Mn==0),1]
msbb_array19.PLQ_Strat[[2]]=msbb_array19.covariates[which(msbb_array19.covariates$PLQ_Mn>=1&msbb_array19.covariates$PLQ_Mn<=10),1]
msbb_array19.PLQ_Strat[[3]]=msbb_array19.covariates[which(msbb_array19.covariates$PLQ_Mn>=11&msbb_array19.covariates$PLQ_Mn<=20),1]
msbb_array19.PLQ_Strat[[4]]=msbb_array19.covariates[which(msbb_array19.covariates$PLQ_Mn>=21),1]


msbb_array.byLobe=foreach(i=1:length(names(msbb_array.byLobe)))%dopar%{
  array_byLobe=do.call("cbind",msbb_array19.2.agg2[which(names(msbb_array19.2.agg2)%in%lobe_bm_area.map$BrainRegion[lobe_bm_area.map$Lobe==names(msbb_array.byLobe)[i]])])  
} 
names(msbb_array.byLobe)=names(msbb_array_byLobe.SCE)=unique(lobe_bm_area.map$Lobe)
lobewise_samples=lapply(msbb_array.byLobe,colnames) 
plq_stratwise_samples=lapply(msbb_array19.PLQ_Strat,function(y)paste("\\b",gsub(pattern = "^X",replacement = "",x = y),"\\b",collapse = "|",sep=""))

for(n in 1:length(lobewise_samples)){
mapped_plq_samples=mclapply(plq_stratwise_samples,function(p)grep(pattern = p,x = lobewise_samples[[n]],value = T),mc.cores = ncores)
msbb_array_byLobe.SCE[[n]]=mclapply(mapped_plq_samples,function(e)single.chip.enrichment2(exprs = msbb_array.byLobe[[n]][,e],geneset = pathprint.Hs.gs,statistic = "mean",normalizedScore = F,progressBar = T),mc.cores = ncores)
}


msbb_array_byLobe.diffPathways=vector(mode = "list",length = length(msbb_array_byLobe.SCE))
names(msbb_array_byLobe.diffPathways)=names(msbb_array_byLobe.SCE)
msbb_array_byLobe.diffPathways$Frontal_Lobe=msbb_array_byLobe.diffPathways$Occipetal_Lobe=msbb_array_byLobe.diffPathways$Temporal_Lobe=msbb_array_byLobe.diffPathways$Parietal_Lobe=msbb_array_byLobe.diffPathways$Dorsal_striatum=vector(mode = "list",length = 3)
names(msbb_array_byLobe.diffPathways$Frontal_Lobe)=names(msbb_array_byLobe.diffPathways$Occipetal_Lobe)=names(msbb_array_byLobe.diffPathways$Temporal_Lobe)=names(msbb_array_byLobe.diffPathways$Parietal_Lobe)=names(msbb_array_byLobe.diffPathways$Dorsal_striatum)=names(msbb_array19.PLQ_Strat)[-1]
  
for(i in 1:length(msbb_array_byLobe.diffPathways)){
  
for(j in 2:4){
    if(is.null(dim(msbb_array_byLobe.SCE[[i]][[j]]))==F){
      if(j == 2){
      group=factor(c(rep(names(msbb_array_byLobe.SCE[[i]])[1],dim(msbb_array_byLobe.SCE[[i]][[1]])[2]),rep(names(msbb_array_byLobe.SCE[[i]])[j],dim(msbb_array_byLobe.SCE[[i]][[j]])[2])))
      design_df=model.matrix(~0+group)
      colnames(design_df)=levels(group)
      contrasts_matrix=makeContrasts(S0-S1,levels = design_df)
      fit=lmFit(object = cbind(data.frame(msbb_array_byLobe.SCE[[i]][[1]],stringsAsFactors = F),data.frame(msbb_array_byLobe.SCE[[i]][[j]],stringsAsFactors = F)),design = design_df)
      fit2=contrasts.fit(fit,contrasts_matrix)
      fit2=eBayes(fit2,trend = T)
      diff_pathways=rownames(topTable(fit2,coef = 1,sort.by = "logFC",adjust.method = "BH",p.value = 0.05,number = 633))
      if(length(diff_pathways)==0){
        msbb_array_byLobe.diffPathways[[i]][[j-1]]=0
      }
      msbb_array_byLobe.diffPathways[[i]][[j-1]]=sort(diff_pathways)
      }
      if(j == 3){
        group=factor(c(rep(names(msbb_array_byLobe.SCE[[i]])[1],dim(msbb_array_byLobe.SCE[[i]][[1]])[2]),rep(names(msbb_array_byLobe.SCE[[i]])[j],dim(msbb_array_byLobe.SCE[[i]][[j]])[2])))
        design_df=model.matrix(~0+group)
        colnames(design_df)=levels(group)
        contrasts_matrix=makeContrasts(S0-S2,levels = design_df)
        fit=lmFit(object = cbind(data.frame(msbb_array_byLobe.SCE[[i]][[1]],stringsAsFactors = F),data.frame(msbb_array_byLobe.SCE[[i]][[j]],stringsAsFactors = F)),design = design_df)
        fit2=contrasts.fit(fit,contrasts_matrix)
        fit2=eBayes(fit2,trend = T)
        diff_pathways=rownames(topTable(fit2,coef = 1,sort.by = "logFC",adjust.method = "BH",p.value = 0.05,number = 633))
        if(length(diff_pathways)==0){
          msbb_array_byLobe.diffPathways[[i]][[j-1]]=0
        }
        msbb_array_byLobe.diffPathways[[i]][[j-1]]=sort(diff_pathways)
      }
      if(j == 4){
        group=factor(c(rep(names(msbb_array_byLobe.SCE[[i]])[1],dim(msbb_array_byLobe.SCE[[i]][[1]])[2]),rep(names(msbb_array_byLobe.SCE[[i]])[j],dim(msbb_array_byLobe.SCE[[i]][[j]])[2])))
        design_df=model.matrix(~0+group)
        colnames(design_df)=levels(group)
        contrasts_matrix=makeContrasts(S0-S3,levels = design_df)
        fit=lmFit(object = cbind(data.frame(msbb_array_byLobe.SCE[[i]][[1]],stringsAsFactors = F),data.frame(msbb_array_byLobe.SCE[[i]][[j]],stringsAsFactors = F)),design = design_df)
        fit2=contrasts.fit(fit,contrasts_matrix)
        fit2=eBayes(fit2,trend = T)
        diff_pathways=rownames(topTable(fit2,coef = 1,sort.by = "logFC",adjust.method = "BH",p.value = 0.05,number = 633))
        if(length(diff_pathways)==0){
          msbb_array_byLobe.diffPathways[[i]][[j-1]]=0
        }
        msbb_array_byLobe.diffPathways[[i]][[j-1]]=sort(diff_pathways)
      }
    }
      
  }
}

##################################################################################################################
msbb_array_plq.byLobe=msbb_array_plq_corrGenes.byLobe=vector(mode = "list",length = 5)
names(msbb_array_plq.byLobe)=names(msbb_array_plq_corrGenes.byLobe)=names(msbb_array.byLobe)
for(i in 1:5){
  msbb_array_plq.byLobe[[i]]=unlist(lapply(unlist(lapply(strsplit(x = colnames(msbb_array.byLobe[[i]]),split = "\\."),`[[`,2)),function(x)msbb_array19.covariates2$PLQ_Mn[grep(pattern = paste("\\bX",x,"\\b",sep = ""),x=msbb_array19.covariates2$BrainBank)]))
  msbb_array_plq_corrGenes.byLobe[[i]]=foreach(g=1:dim(msbb_array.byLobe[[i]])[1]) %dopar% {data.frame(Gene=rownames(msbb_array.byLobe[[i]][g,]),
                                         Spearman.Rho=cor.test(x=unlist(unname(msbb_array.byLobe[[i]][g,])),y=msbb_array_plq.byLobe[[i]])[[4]],
                                         pval=cor.test(x=unlist(unname(msbb_array.byLobe[[i]][g,])),y=msbb_array_plq.byLobe[[i]])[[3]],
                                         stringsAsFactors=F)}
}


msbb_array_plq.byRegion=msbb_array_plq_corrGenes.byRegion=vector(mode = "list",length = 17)
names(msbb_array_plq.byRegion)=names(msbb_array_plq_corrGenes.byRegion)=names(msbb_array19)
for(i in 1:17){
  msbb_array_plq.byRegion[[i]]=msbb_array19.covariates2$PLQ_Mn[grep(pattern=paste('\\bX',colnames(msbb_array19.2.agg2[[i]]),'\\b',sep="",collapse='|'),x=msbb_array19.covariates2$BrainBank)]
  msbb_array_plq_corrGenes.byRegion[[i]]=foreach(g=1:dim(msbb_array19.2.agg2[[i]])[1]) %dopar% {data.frame(Gene=rownames(msbb_array19.2.agg2[[i]][g,]),
                                                                                                       Spearman.Rho=cor.test(x=unlist(unname(msbb_array19.2.agg2[[i]][g,])),y=msbb_array_plq.byRegion[[i]])[[4]],
                                                                                                       pval=cor.test(x=unlist(unname(msbb_array19.2.agg2[[i]][g,])),y=msbb_array_plq.byRegion[[i]])[[3]],
                                                                                                       stringsAsFactors=F)}
}
msbb_array_plq_corrGenes.byRegion2=lapply(lapply(msbb_array_plq_corrGenes.byRegion,function(a)data.frame(rbindlist(a),stringsAsFactors=F)),function(b){b<-b[-1,];b})
msbb_array_plq_corrGenes.byRegion2_adjp=lapply(msbb_array_plq_corrGenes.byRegion2,function(x){x$adj.p<-p.adjust(p=x$pval,method='fdr');x})

###################################################################################################
gse29378_series_matrix=fread("/Users/sandeepamberkar/Work/Data/AD_GSE29378/GSE29378_series_matrix.txt",sep = "\t",header = F,showProgress = T,data.table = F)
