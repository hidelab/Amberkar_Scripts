library(pathprint)
library(org.Hs.eg.db)
library(dplyr)
library(magrittr)
library(data.table)
library(parallel)
library(illuminaHumanv3.db)

library(clusterProfiler)
library(doParallel)

ncores=8

pathprint_membership.genes=lapply(pathprint.Hs.gs,function(x)unname(mapIds(x = org.Hs.eg.db,keys = as.character(x),keytype = "ENTREZID",column = "SYMBOL")))
names(pathprint_membership.genes)[grep(pattern = "/",names(pathprint_membership.genes))]=gsub(pattern = "/",replacement = " ",x = grep(pattern = "/",names(pathprint_membership.genes),value = T))
for(i in 1:633){
  write(pathprint_membership.genes[[i]],paste(names(pathprint_membership.genes)[i],"MemberGenes.txt",collapse = "_"),sep = "\n")
}
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

dbs <- listEnrichrDbs()
kegg_dbs=grep(pattern = "KEGG",x = dbs$libraryName,value = T)[3]
biocarta_dbs=grep(pattern = "BioCarta|Panther",x = dbs$libraryName,value = T)[4]
panther_dbs=grep(pattern = "BioCarta|Panther",x = dbs$libraryName,value = T)[5]


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
msbb_array19.2.agg=lapply(msbb_array19.2[-18],function(y)aggregate(x=y[,-c(1:4)],by=list(EntrezID=y$ENTREZ_GENE_ID),mean))
msbb_array19.2.agg2=lapply(msbb_array19.2.agg,function(x){rownames(x)<- x$EntrezID;x <- x[-grep(pattern = "///",x = rownames(x)),-1];x})

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

msbb_array19.CDR_Strat=vector(mode = "list",length = length(sort(unique(msbb_array19.covariates$CDR))))
names(msbb_array19.CDR_Strat)=sort(unique(msbb_array19.covariates$CDR))
for(i in 1:length(msbb_array19.CDR_Strat)){
  msbb_array19.CDR_Strat[[i]]=msbb_array19.covariates$BrainBank[msbb_array19.covariates$CDR==names(msbb_array19.CDR_Strat)[i]]  
}


msbb_array.byLobe=foreach(i=1:length(names(msbb_array.byLobe)))%dopar%{
  array_byLobe=do.call("cbind",msbb_array19.2.agg2[which(names(msbb_array19.2.agg2)%in%lobe_bm_area.map$BrainRegion[lobe_bm_area.map$Lobe==names(msbb_array.byLobe)[i]])])  
} 
msbb_array.byLobe2=foreach(i=1:length(names(msbb_array.byLobe2)))%dopar%{
  array_byLobe2=do.call("cbind",msbb_array19.3[which(names(msbb_array19.3)%in%lobe_bm_area.map$BrainRegion[lobe_bm_area.map$Lobe==names(msbb_array.byLobe2)[i]])])  
} 
names(msbb_array.byLobe)=names(msbb_array.byLobe2)=unique(lobe_bm_area.map$Lobe)


#msbb_array_byLobe.fingerprint=lapply(msbb_array.byLobe2,function(y)exprs2fingerprint(exprs = y,platform = "GPL570",species = "Homo sapiens",progressBar = T))
msbb_array_byLobe.fingerprint=readRDS("msbb_array_byLobe_fingerprint.RDS")
msbb_array_byLobe.fingerprint_Control=lapply(msbb_array_byLobe.fingerprint,function(x)x[,grep(pattern=Control_sample.vector,x=colnames(x))])
msbb_array_byLobe.fingerprint_AD=lapply(msbb_array_byLobe.fingerprint,function(x)x[,grep(pattern=AD_sample.vector,x=colnames(x))])
saveRDS(msbb_array_byLobe.fingerprint_Control,"msbb_array_byLobe_fingerprint_Control.RDS")
saveRDS(msbb_array_byLobe.fingerprint_AD,"msbb_array_byLobe_fingerprint_AD.RDS")
names(msbb_array.byLobe)=names(msbb_array_byLobe.SCE)=unique(lobe_bm_area.map$Lobe)
lobewise_samples=lapply(msbb_array.byLobe,colnames) 
plq_stratwise_samples=lapply(msbb_array19.PLQ_Strat,function(y)paste("\\b",gsub(pattern = "^X",replacement = "",x = y),"\\b",collapse = "|",sep=""))

for(n in 1:length(lobewise_samples)){
mapped_plq_samples=lapply(plq_stratwise_samples,function(p)grep(pattern = p,x = lobewise_samples[[n]],value = T))

msbb_array_byLobe.SCE[[n]]=lapply(mapped_plq_samples,function(e){if(length(e)>1){single.chip.enrichment2(exprs = msbb_array.byLobe[[n]][,e],geneset = pathprint.Hs.gs,statistic = "mean",normalizedScore = F,progressBar = T)}else cat(paste("\nNot enough samples!\n"))})
}

regionwise_samples=lapply(msbb_array19.2.agg2,colnames)
msbb_array_byRegion.SCE=vector(mode = "list",length = length(msbb_array19.2.agg2))
names(msbb_array_byRegion.SCE)=names(msbb_array19.2.agg2)
for(i in 1:length(msbb_array_byRegion.SCE)){
  msbb_array_byRegion.SCE[[i]]=vector(mode = "list",length = 4)
  names(msbb_array_byRegion.SCE)=names(msbb_array_byLobe.SCE$Frontal_Lobe)
}

#
for(n in 1:length(regionwise_samples)){
  mapped_plq_samples=lapply(plq_stratwise_samples,function(p)grep(pattern = p,x = regionwise_samples[[n]],value = T))
  msbb_array_byRegion.SCE[[n]]=mclapply(mapped_plq_samples,function(e)single.chip.enrichment2(exprs = msbb_array19.2.agg2[[n]][,e],geneset = pathprint.Hs.gs,statistic = "mean",normalizedScore = F,progressBar = T),mc.cores=2)
}
names(msbb_array_byRegion.SCE)=names(msbb_array19.2.agg2)

#Diff pathway analysis
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
      diff_pathways=topTable(fit2,coef = 1,sort.by = "logFC",adjust.method = "BH",p.value = 0.1,number = 633)%>%rownames_to_column("Pathway_Geneset")
      if(length(diff_pathways)==0){
        msbb_array_byLobe.diffPathways[[i]][[j-1]]=0
      }
      else
      msbb_array_byLobe.diffPathways[[i]][[j-1]]=diff_pathways
      }
      if(j == 3){
        group=factor(c(rep(names(msbb_array_byLobe.SCE[[i]])[1],dim(msbb_array_byLobe.SCE[[i]][[1]])[2]),rep(names(msbb_array_byLobe.SCE[[i]])[j],dim(msbb_array_byLobe.SCE[[i]][[j]])[2])))
        design_df=model.matrix(~0+group)
        colnames(design_df)=levels(group)
        contrasts_matrix=makeContrasts(S0-S2,levels = design_df)
        fit=lmFit(object = cbind(data.frame(msbb_array_byLobe.SCE[[i]][[1]],stringsAsFactors = F),data.frame(msbb_array_byLobe.SCE[[i]][[j]],stringsAsFactors = F)),design = design_df)
        fit2=contrasts.fit(fit,contrasts_matrix)
        fit2=eBayes(fit2,trend = T)
        diff_pathways=topTable(fit2,coef = 1,sort.by = "logFC",adjust.method = "BH",p.value = 0.1,number = 633)%>%rownames_to_column("Pathway_Geneset")
        if(length(diff_pathways)==0){
          msbb_array_byLobe.diffPathways[[i]][[j-1]]=0
        }
        else
        msbb_array_byLobe.diffPathways[[i]][[j-1]]=diff_pathways
      }
      if(j == 4){
        group=factor(c(rep(names(msbb_array_byLobe.SCE[[i]])[1],dim(msbb_array_byLobe.SCE[[i]][[1]])[2]),rep(names(msbb_array_byLobe.SCE[[i]])[j],dim(msbb_array_byLobe.SCE[[i]][[j]])[2])))
        design_df=model.matrix(~0+group)
        colnames(design_df)=levels(group)
        contrasts_matrix=makeContrasts(S0-S3,levels = design_df)
        fit=lmFit(object = cbind(data.frame(msbb_array_byLobe.SCE[[i]][[1]],stringsAsFactors = F),data.frame(msbb_array_byLobe.SCE[[i]][[j]],stringsAsFactors = F)),design = design_df)
        fit2=contrasts.fit(fit,contrasts_matrix)
        fit2=eBayes(fit2,trend = T)
        diff_pathways=topTable(fit2,coef = 1,sort.by = "logFC",adjust.method = "BH",p.value = 0.1,number = 633)%>%rownames_to_column("Pathway_Geneset")
        if(length(diff_pathways)==0){
          msbb_array_byLobe.diffPathways[[i]][[j-1]]=0
        }
        else
        msbb_array_byLobe.diffPathways[[i]][[j-1]]=diff_pathways
      }
    }
      
  }
}

msbb_array_byLobe.diffPathways$Frontal_Lobe=Filter(f = function(x)dim(x)[1]>0,x = msbb_array_byLobe.diffPathways$Frontal_Lobe)
msbb_array_byLobe.diffPathways$Occipetal_Lobe=Filter(f = function(x)dim(x)[1]>0,x = msbb_array_byLobe.diffPathways$Occipetal_Lobe)
msbb_array_byLobe.diffPathways$Temporal_Lobe=Filter(f = function(x)dim(x)[1]>0,x = msbb_array_byLobe.diffPathways$Temporal_Lobe)
msbb_array_byLobe.diffPathways$Parietal_Lobe=Filter(f = function(x)dim(x)[1]>0,x = msbb_array_byLobe.diffPathways$Parietal_Lobe)
msbb_array_byLobe.diffPathways$Dorsal_striatum=Filter(f = function(x)dim(x)[1]>0,x = msbb_array_byLobe.diffPathways$Dorsal_striatum)

rand_s2_DEPs=replicate(n = 10000,expr = Reduce(intersect,list(fl_s2_dep=sample(x = names(pathprint.Hs.gs),size = length(msbb_array_byLobe.diffPathways$Frontal_Lobe$S2$Pathway_Geneset),replace = F),
     tl_s2_dep=sample(x = names(pathprint.Hs.gs),size = length(msbb_array_byLobe.diffPathways$Temporal_Lobe$S2$Pathway_Geneset),replace = F),
     ol_s2_dep=sample(x = names(pathprint.Hs.gs),size = length(msbb_array_byLobe.diffPathways$Occipetal_Lobe$S2$Pathway_Geneset),replace = F))))


for(t in 1:length(msbb_array_byLobe.diffPathways)){
  for(n in 1:length(msbb_array_byLobe.diffPathways[[t]])){
    fwrite(msbb_array_byLobe.diffPathways[[t]][[n]],paste("MSBB_DEP_",names(msbb_array_byLobe.diffPathways[t]),"_",names(msbb_array_byLobe.diffPathways[[t]][n]),"_","table.txt",sep = ""),sep = "\t",col.names = T,row.names = F)
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
# Determine ABeta associated pathway signatures
msbb_array_plq_corrPathways.byLobe=vector(mode = "list",length = 5)
names(msbb_array_plq_corrPathways.byLobe)=names(msbb_array_byLobe.diffPathways)
for(i in 1:5){
  msbb_array_plq_corrPathways.byLobe[[i]]=vector(mode = "list",length = 3)
  names(msbb_array_plq_corrPathways.byLobe[[i]])=paste("S",1:3,sep="")
}
for(l in 1:5){
  cat(paste("Processing brain lobe",names(msbb_array_byLobe.SCE)[l],"\n",sep = " "))
  for(s in 1:3){
    cat(paste("for group",names(msbb_array_byLobe.SCE[[l]])[s+1],"\n",sep = " "))
    
      sample_reps=table(paste("X",unlist(lapply(strsplit(x = colnames(msbb_array_byLobe.SCE[[l]][[s+1]]),split = "\\."),`[[`,2)),sep = ""))
      plq_values.res=foreach(j=1:length(sample_reps))%do%{
        rep(msbb_array19.covariates2$PLQ_Mn[msbb_array19.covariates2$BrainBank%in%names(sample_reps[j])],sample_reps[j])
      }
      cat(paste("Computing pathway correlations\n",sep = " "))
      plq_corr_pathways=foreach(p=1:633)%do%{
        cor.test(x = msbb_array_byLobe.SCE[[l]][[s+1]][p,],y=unlist(plq_values.res),method = "spearman")[c(3:4)]
      }
      plq_corr_pathways.df=data.frame(rbindlist(plq_corr_pathways),stringsAsFactors = F)
      rownames(plq_corr_pathways.df)=names(pathprint.Hs.gs)
      plq_corr_pathways.df$adj.p=p.adjust(plq_corr_pathways.df$p.value,method = "BH")
      cat(paste("Filtering correlated pathways at adj.P <= 0.1\n",sep = " "))
      msbb_array_plq_corrPathways.byLobe[[l]][[s]]=rownames(plq_corr_pathways.df[plq_corr_pathways.df$adj.p<=0.1,])
    }
    
  
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

msbb_array_plq_corrGenes.byLobe2=lapply(lapply(msbb_array_plq_corrGenes.byLobe,function(x)data.frame(rbindlist(x),stringsAsFactors = F)%>%mutate(adjp=p.adjust(p=pval,method='fdr'))),function(y){y[-1,]})
msbb_array_plq_corrGenes.byLobe2=lapply(msbb_array_plq_corrGenes.byLobe2,function(x){x$Gene_HGNC=unname(mapIds(x = org.Hs.eg.db,keys = x$Gene,column = "SYMBOL",keytype = "ENTREZID"));x})
msbb_array_plq_corrGenes.byLobe2=Filter(f = function(y)dim(y)[1]>0,x = lapply(msbb_array_plq_corrGenes.byLobe2,function(x)x%>%filter(adjp<=0.1)%>%dplyr::select(c(Gene,Gene_HGNC,Spearman.Rho,Spearman.Rho,pval,adjp))))
msbb_array_plq_corrGenes.byLobe2=lapply(msbb_array_plq_corrGenes.byLobe2,function(x)x%>%filter(is.na(Gene_HGNC)==F))
#Only Frontal and Temporal Lobes have significant genes at adjp<0.05
msbb_array_plq_corrGenes_pFiltered_byLobe.Pathprint=Filter(f = function(a)length(a)>0,x = lapply(lapply(lapply(msbb_array_plq_corrGenes.byLobe2,function(x)x%>%filter(pval<=0.05)%>%dplyr::select(Gene)),`[[`,1),function(y)data.frame(enricher(gene = y,TERM2GENE = pathprint.gmt,pvalueCutoff = 0.05,pAdjustMethod = "BH"))[,1]))
msbb_array_plq_corrGenes_adjpFiltered_byLobe.Pathprint=Filter(f = function(a)length(a)>0,x=lapply(lapply(lapply(msbb_array_plq_corrGenes.byLobe2,function(x)x%>%filter(pval<=0.05)%>%dplyr::select(Gene))[c(1,3)],`[[`,1),function(y)data.frame(enricher(gene = y,TERM2GENE = pathprint.gmt,pvalueCutoff = 0.05,pAdjustMethod = "BH"))[,1]))
msbb_array_plq_corrGenes_adjpFiltered_byLobe.Pathprint_df=Filter(f = function(a)length(a)>0,x=lapply(lapply(lapply(msbb_array_plq_corrGenes.byLobe2,function(x)x%>%filter(pval<=0.05)%>%dplyr::select(Gene))[c(1,3)],`[[`,1),function(y)data.frame(enricher(gene = y,TERM2GENE = pathprint.gmt,pvalueCutoff = 0.05,pAdjustMethod = "BH"))))
msbb_array_plq_corrGenes_adjpFiltered_byLobe.KEGG=lapply(msbb_array_plq_corrGenes.byLobe2,function(x)enrichr(genes = x$Gene_HGNC,databases = kegg_dbs)[[1]]%>%filter(Adjusted.P.value<=0.1))

#Compare diffPatways across and within lobes
common_diffPathways_PLQstrata=matrix(NA,nrow=3,ncol=6)
lobe_diffPathways_jaccardMatrix=read.table("/Users/sandeepamberkar/Work/Data/MSMM-MSBB-HBTRC-Data/MSBB_Array19/Normalised_Data/Lobe_DiffPathways_JaccardMatrix.txt",sep="\t",header=T,row.names = 1,as.is = T)
colnames(lobe_diffPathways_jaccardMatrix)=c("Frontal.Lobe.S2","Occipital.Lobe.S2","Temporal.Lobe.S2","Frontal.Lobe.S3","Occipital.Lobe.S3","Temporal.Lobe.S3")
rownames(lobe_diffPathways_jaccardMatrix)=colnames(lobe_diffPathways_jaccardMatrix)
pheatmap(lobe_diffPathways_jaccardMatrix,cellwidth = 50,cellheight = 50,fontsize_row = 15,fontsize_col = 15,fontsize = 15,main = "DiffPathways comparison across plaque strata")


msbb_array_byLobe.DEG=msbb_array_byLobe_DEG.Pathprint=vector(mode = "list",length = 5)
names(msbb_array_byLobe.DEG)=names(msbb_array_byLobe_DEG.Pathprint)=names(msbb_array_byLobe.diffPathways)
msbb_array_byLobe.DEG$Frontal_Lobe=msbb_array_byLobe.DEG$Occipetal_Lobe=msbb_array_byLobe.DEG$Temporal_Lobe=msbb_array_byLobe.DEG$Parietal_Lobe=msbb_array_byLobe.DEG$Dorsal_striatum=vector(mode = "list",length = 3)
msbb_array_byLobe_DEG.Pathprint$Frontal_Lobe=msbb_array_byLobe_DEG.Pathprint$Occipetal_Lobe=msbb_array_byLobe_DEG.Pathprint$Temporal_Lobe=msbb_array_byLobe_DEG.Pathprint$Parietal_Lobe=msbb_array_byLobe_DEG.Pathprint$Dorsal_striatum=vector(mode = "list",length = 3)
names(msbb_array_byLobe.DEG$Frontal_Lobe)=names(msbb_array_byLobe.DEG$Occipetal_Lobe)=names(msbb_array_byLobe.DEG$Temporal_Lobe)=names(msbb_array_byLobe.DEG$Parietal_Lobe)=names(msbb_array_byLobe.DEG$Dorsal_striatum)=c("S1","S2","S3")
names(msbb_array_byLobe_DEG.Pathprint$Frontal_Lobe)=names(msbb_array_byLobe_DEG.Pathprint$Occipetal_Lobe)=names(msbb_array_byLobe_DEG.Pathprint$Temporal_Lobe)=names(msbb_array_byLobe_DEG.Pathprint$Parietal_Lobe)=names(msbb_array_byLobe_DEG.Pathprint$Dorsal_striatum)=c("S1","S2","S3")
#Lobewise, PLQ strat DEGs
rownames(msbb_array.byLobe$Frontal_Lobe)=rownames(msbb_array19.2.agg2$AC)
rownames(msbb_array.byLobe$Occipetal_Lobe)=rownames(msbb_array19.2.agg2$AC)
rownames(msbb_array.byLobe$Temporal_Lobe)=rownames(msbb_array19.2.agg2$AC)
rownames(msbb_array.byLobe$Parietal_Lobe)=rownames(msbb_array19.2.agg2$AC)
rownames(msbb_array.byLobe$Dorsal_striatum)=rownames(msbb_array19.2.agg2$AC)
for(t1 in 1:5){
  for(s in 2:length(names(msbb_array19.PLQ_Strat))){
    control_data=msbb_array.byLobe[[t1]][,grep(pattern = paste('\\b',gsub(pattern = "X",replacement = "",x = msbb_array19.PLQ_Strat$S0),'\\b',sep="",collapse = '|'),x = colnames(msbb_array.byLobe[[t1]]))]
    case_data=msbb_array.byLobe[[t1]][,grep(pattern = paste('\\b',gsub(pattern = "X",replacement = "",x = msbb_array19.PLQ_Strat[[s]]),'\\b',sep="",collapse = '|'),x = colnames(msbb_array.byLobe[[t1]]))]
    group=factor(c(rep(names(msbb_array19.PLQ_Strat)[1],length=length(colnames(control_data))),(rep(names(msbb_array19.PLQ_Strat)[s],length=length(colnames(case_data))))))
    design_df=model.matrix(~0+group)
    colnames(design_df)=levels(group)
    if(s==2){
      contrasts_matrix=limma::makeContrasts(S0-S1,levels = design_df)
      fit=limma::lmFit(object = cbind(control_data,case_data),design = design_df)
      fit2=limma::contrasts.fit(fit,contrasts_matrix)
      fit2=limma::eBayes(fit2)
      fit2=limma::eBayes(fit2,trend = T)
      msbb_array_byLobe.DEG[[t1]][[s-1]]=limma::topTable(fit2,coef = 1,sort.by = "logFC",adjust.method = "BH",p.value = 0.1,number = 18436)%>%rownames_to_column("Gene")
      if(dim(msbb_array_byLobe.DEG[[t1]][[s-1]])[1]==0){
        msbb_array_byLobe.DEG[[t1]][[s-1]]=0
      }
      else
        msbb_array_byLobe.DEG[[t1]][[s-1]]=limma::topTable(fit2,coef = 1,sort.by = "logFC",adjust.method = "BH",p.value = 0.1,number = 18436)%>%rownames_to_column("Gene")
        #msbb_array_byLobe.DEG[[t1]][[s-1]]=msbb_array_byLobe.DEG[[t1]][[s-1]]%>%mutate(Gene_HGNC=unname(mapIds(x = org.Hs.eg.db,keys = Gene,column = "SYMBOL",keytype = "ENTREZID")))%>%dplyr::select(c(Gene_HGNC,logFC,AveExpr,t,P.Value,adj.P.Val,B))
    }
    if(s==3){
      contrasts_matrix=limma::makeContrasts(S0-S2,levels = design_df)
      fit=limma::lmFit(object = cbind(control_data,case_data),design = design_df)
      fit2=limma::contrasts.fit(fit,contrasts_matrix)
      fit2=limma::eBayes(fit2)
      fit2=limma::eBayes(fit2,trend = T)
      msbb_array_byLobe.DEG[[t1]][[s-1]]=limma::topTable(fit2,coef = 1,sort.by = "logFC",adjust.method = "BH",p.value = 0.1,number = 18436)%>%rownames_to_column("Gene")
      if(dim(msbb_array_byLobe.DEG[[t1]][[s-1]])[1]==0){
        msbb_array_byLobe.DEG[[t1]][[s-1]]=0
      }
      else
        msbb_array_byLobe.DEG[[t1]][[s-1]]=limma::topTable(fit2,coef = 1,sort.by = "logFC",adjust.method = "BH",p.value = 0.1,number = 18436)%>%rownames_to_column("Gene")
        #msbb_array_byLobe.DEG[[t1]][[s-1]]=msbb_array_byLobe.DEG[[t1]][[s-1]]%>%mutate(Gene_HGNC=unname(mapIds(x = org.Hs.eg.db,keys = Gene,column = "SYMBOL",keytype = "ENTREZID")))%>%dplyr::select(c(Gene_HGNC,logFC,AveExpr,t,P.Value,adj.P.Val,B))
    }
    if(s==4){
      contrasts_matrix=limma::makeContrasts(S0-S3,levels = design_df)
      fit=limma::lmFit(object = cbind(control_data,case_data),design = design_df)
      fit2=limma::contrasts.fit(fit,contrasts_matrix)
      fit2=limma::eBayes(fit2)
      fit2=limma::eBayes(fit2,trend = T)
      msbb_array_byLobe.DEG[[t1]][[s-1]]=limma::topTable(fit2,coef = 1,sort.by = "logFC",adjust.method = "BH",p.value = 0.1,number = 18436)%>%rownames_to_column("Gene")
      if(dim(msbb_array_byLobe.DEG[[t1]][[s-1]])[1]==0){
        msbb_array_byLobe.DEG[[t1]][[s-1]]=0
      }
      else
        msbb_array_byLobe.DEG[[t1]][[s-1]]=limma::topTable(fit2,coef = 1,sort.by = "logFC",adjust.method = "BH",p.value = 0.1,number = 18436)%>%rownames_to_column("Gene")
        #msbb_array_byLobe.DEG[[t1]][[s-1]]=msbb_array_byLobe.DEG[[t1]][[s-1]]%>%mutate(Gene_HGNC=unname(mapIds(x = org.Hs.eg.db,keys = Gene,column = "SYMBOL",keytype = "ENTREZID")))%>%dplyr::select(c(Gene_HGNC,logFC,AveExpr,t,P.Value,adj.P.Val,B))
    }
    
    
  }
  
}
# msbb_array_byLobe.DEG=msbb_array_byLobe.DEG[-2]
# msbb_array_byLobe.DEG$Dorsal_striatum=msbb_array_byLobe.DEG$Dorsal_striatum[-1]
msbb_array_byLobe.DEG$Parietal_Lobe=msbb_array_byLobe.DEG$Parietal_Lobe[-c(1,3)]
msbb_array_byLobe.DEG$Occipetal_Lobe=msbb_array_byLobe.DEG$Occipetal_Lobe[-c(1,3)]
msbb_array_byLobe.DEG$Frontal_Lobe=lapply(msbb_array_byLobe.DEG$Frontal_Lobe,function(x)x%>%mutate(Gene_HGNC=paste(" ",unname(mapIds(x = org.Hs.eg.db,keys = Gene,column = "SYMBOL",keytype = "ENTREZID")),sep = ""))%>%dplyr::select(c(Gene_HGNC,Gene,logFC,AveExpr,t,P.Value,adj.P.Val,B)))
msbb_array_byLobe.DEG$Occipetal_Lobe=lapply(msbb_array_byLobe.DEG$Occipetal_Lobe,function(x)x%>%mutate(Gene_HGNC=paste(" ",unname(mapIds(x = org.Hs.eg.db,keys = Gene,column = "SYMBOL",keytype = "ENTREZID")),sep = ""))%>%dplyr::select(c(Gene_HGNC,Gene,logFC,AveExpr,t,P.Value,adj.P.Val,B)))
msbb_array_byLobe.DEG$Temporal_Lobe=lapply(msbb_array_byLobe.DEG$Temporal_Lobe,function(x)x%>%mutate(Gene_HGNC=paste(" ",unname(mapIds(x = org.Hs.eg.db,keys = Gene,column = "SYMBOL",keytype = "ENTREZID")),sep = ""))%>%dplyr::select(c(Gene_HGNC,Gene,logFC,AveExpr,t,P.Value,adj.P.Val,B)))
msbb_array_byLobe.DEG$Parietal_Lobe=lapply(msbb_array_byLobe.DEG$Parietal_Lobe,function(x)x%>%mutate(Gene_HGNC=paste(" ",unname(mapIds(x = org.Hs.eg.db,keys = Gene,column = "SYMBOL",keytype = "ENTREZID")),sep = ""))%>%dplyr::select(c(Gene_HGNC,Gene,logFC,AveExpr,t,P.Value,adj.P.Val,B)))
msbb_array_byLobe.DEG$Dorsal_striatum=lapply(msbb_array_byLobe.DEG$Dorsal_striatum,function(x)x%>%mutate(Gene_HGNC=paste(" ",unname(mapIds(x = org.Hs.eg.db,keys = Gene,column = "SYMBOL",keytype = "ENTREZID")),sep = ""))%>%dplyr::select(c(Gene_HGNC,Gene,logFC,AveExpr,t,P.Value,adj.P.Val,B)))




for(t in 1:length(msbb_array_byLobe.DEG)){
  for(n in 1:length(msbb_array_byLobe.DEG[[t]])){
    fwrite(msbb_array_byLobe.DEG[[t]][[n]],paste("MSBB_DEG_",names(msbb_array_byLobe.DEG[t]),"_",names(msbb_array_byLobe.DEG[[t]][n]),"_","table.txt",sep = ""),sep = "\t",col.names = T,row.names = F)
  }
}

frontal_lobe_common_DEGs.exprs=msbb_array.byLobe$Frontal_Lobe[which(rownames(msbb_array.byLobe$Frontal_Lobe)%in%Reduce(intersect,lapply(msbb_array_byLobe.DEG$Frontal_Lobe,function(x)x%>%pull(Gene)))),]

#Overlap comparison
u133_agg.Entrez=rownames(topTable(fit2,coef = 1,sort.by = "logFC",adjust.method = "BH",p.value = 1,number = 18436))
u133_agg.Symbol=unlist(unname(mapIds(x = org.Hs.eg.db,keytype = "ENTREZID",column = "SYMBOL",keys = u133_agg.Entrez)))

rand_DEG_overlaps.FrontalLobe=replicate(n = 1000,expr = intersect(intersect(sample(x = u133_agg.Symbol,size = length(msbb_array_byLobe.DEG$Frontal_Lobe$S1$Gene_HGNC),replace = F),
                                                                             sample(x = u133_agg.Symbol,size = length(msbb_array_byLobe.DEG$Frontal_Lobe$S2$Gene_HGNC),replace = F)),
                                                                   sample(x = u133_agg.Symbol,size = length(msbb_array_byLobe.DEG$Frontal_Lobe$S3$Gene_HGNC),replace = F)))

rand_DEG_overlaps.TemporalLobe=replicate(n = 1000,expr = intersect(intersect(sample(x = u133_agg.Symbol,size = length(msbb_array_byLobe.DEG$Temporal_Lobe$S1$Gene_HGNC),replace = F),
                                                            sample(x = u133_agg.Symbol,size = length(msbb_array_byLobe.DEG$Temporal_Lobe$S2$Gene_HGNC),replace = F)),
                                                  sample(x = u133_agg.Symbol,size = length(msbb_array_byLobe.DEG$Temporal_Lobe$S3$Gene_HGNC),replace = F)))




msbb_array_byLobe_DEG.Pathprint$Frontal_Lobe=lapply(msbb_array_byLobe.DEG$Frontal_Lobe,function(x)enricher(gene = rownames(x),TERM2GENE = pathprint.gmt)[,1])
msbb_array_byLobe_DEG.Pathprint$Occipetal_Lobe=lapply(msbb_array_byLobe.DEG$Occipetal_Lobe,function(x)enricher(gene = rownames(x),TERM2GENE = pathprint.gmt)[,1])
msbb_array_byLobe_DEG.Pathprint$Temporal_Lobe=lapply(msbb_array_byLobe.DEG$Temporal_Lobe,function(x)enricher(gene = rownames(x),TERM2GENE = pathprint.gmt)[,1])
msbb_array_byLobe_DEG.Pathprint$Parietal_Lobe=lapply(msbb_array_byLobe.DEG$Parietal_Lobe,function(x)enricher(gene = rownames(x),TERM2GENE = pathprint.gmt)[,1])
msbb_array_byLobe_DEG.Pathprint$Dorsal_striatum=lapply(msbb_array_byLobe.DEG$Dorsal_striatum,function(x)enricher(gene = rownames(x),TERM2GENE = pathprint.gmt)[,1])

msbb_array_byLobe_DEG.KEGG=list()
msbb_array_byLobe_DEG.KEGG$Frontal_Pole=lapply(msbb_array_byLobe.DEG$Frontal_Lobe,function(x)enrichr(genes = x$Gene_HGNC,databases = kegg_dbs)[[1]]%>%filter(Adjusted.P.value<=0.1))
msbb_array_byLobe_DEG.KEGG$Occipetal_Lobe=lapply(msbb_array_byLobe.DEG$Occipetal_Lobe,function(x)enrichr(genes = x$Gene_HGNC,databases = kegg_dbs)[[1]]%>%filter(Adjusted.P.value<=0.1))
msbb_array_byLobe_DEG.KEGG$Temporal_Lobe=lapply(msbb_array_byLobe.DEG$Temporal_Lobe,function(x)enrichr(genes = x$Gene_HGNC,databases = kegg_dbs)[[1]]%>%filter(Adjusted.P.value<=0.1))
msbb_array_byLobe_DEG.KEGG$Parietal_Lobe=lapply(msbb_array_byLobe.DEG$Parietal_Lobe,function(x)enrichr(genes = x$Gene_HGNC,databases = kegg_dbs)[[1]]%>%filter(Adjusted.P.value<=0.1))
msbb_array_byLobe_DEG.KEGG$Dorsal_striatum=lapply(msbb_array_byLobe.DEG$Dorsal_striatum,function(x)enrichr(genes = x$Gene_HGNC,databases = kegg_dbs)[[1]]%>%filter(Adjusted.P.value<=0.1))

#Identify commonalities between Frontal and Temporal Lobe
fl_tl_diffPathways_common_Gset=Reduce(intersect,lapply(msbb_array_byLobe.diffPathways[c(1,3)],`[[`,2))
fl_tl_corrPathways_common_Gset=Reduce(intersect,lapply(msbb_array_plq_corrPathways.byLobe[c(1,3)],`[[`,2))
fl_tl_DEG_common_Gset=Reduce(intersect,lapply(msbb_array_byLobe_DEG.Pathprint[c(1,3)],`[[`,2))

###################################################################################################
# Apply all the analyses on published dataset, PMID - 23705665
###################################################################################################
gse29378_series_matrix=getGEO(GEO = "GSE29378",GSEMatrix = T)
gse29378.exprs=data.frame(exprs(gse29378_series_matrix$GSE29378_series_matrix.txt.gz),stringsAsFactors = F)
gse29378.phenoData=pData(gse29378_series_matrix$GSE29378_series_matrix.txt.gz)
gse29378.featureData=fData(gse29378_series_matrix$GSE29378_series_matrix.txt.gz)

#gse29378_metadata_covariate.merged2=gse29378_metadata_covariate.merged[is.na(gse29378_metadata_covariate.merged$APOE.genotype)==F,][,-11]
#rownames(gse29378_metadata_covariate.merged2)=gse29378_metadata_covariate.merged2$GSM_ID

gse29378_PLQ_Strat=gse29378_PLQ_Strat.SCE=vector(mode = "list",length = 4)
names(gse29378_PLQ_Strat)=names(gse29378_PLQ_Strat.SCE)=c("S0","S1","S2","S3")
gse29378_PLQ_Strat$S0=gse29378.phenoData$geo_accession[gse29378.phenoData$`plaquescore(0-3):ch1`=="0"]
gse29378_PLQ_Strat$S1=gse29378.phenoData$geo_accession[gse29378.phenoData$`plaquescore(0-3):ch1`=="1"]
gse29378_PLQ_Strat$S2=gse29378.phenoData$geo_accession[gse29378.phenoData$`plaquescore(0-3):ch1`=="2"]
gse29378_PLQ_Strat$S3=gse29378.phenoData$geo_accession[gse29378.phenoData$`plaquescore(0-3):ch1`=="3"]


gse29378.exprs$EntrezID=as.character(gse29378.featureData$Entrez_Gene_ID)
gse29378.exprs$Symbol=gse29378.featureData$Symbol
gse29378_exprs.agg=aggregate(x = gse29378.exprs[,-which((colnames(gse29378.exprs)=="EntrezID")|(colnames(gse29378.exprs)=="Symbol"))],by=list(entrez=gse29378.exprs$EntrezID),mean)
gse29378_exprs.agg2=aggregate(x = gse29378.exprs[,-which((colnames(gse29378.exprs)=="EntrezID")|(colnames(gse29378.exprs)=="Symbol"))],by=list(symbol=gse29378.exprs$Symbol),mean)
rownames(gse29378_exprs.agg)=gse29378_exprs.agg$entrez
rownames(gse29378_exprs.agg2)=gse29378_exprs.agg2$symbol
gse29378_exprs.agg=gse29378_exprs.agg[,-1]
gse29378_exprs.agg2=gse29378_exprs.agg2[,-1]
gse29378_exprs.agg=log2(gse29378_exprs.agg)
gse29378_exprs.agg2=log2(gse29378_exprs.agg2)

gse29378_PLQ_Strat.SCE$S0=single.chip.enrichment2(exprs = gse29378_exprs.agg[,gse29378_PLQ_Strat$S0],geneset = pathprint.Hs.gs,statistic = "mean",normalizedScore = F,progressBar = T)
gse29378_PLQ_Strat.SCE$S1=single.chip.enrichment2(exprs = gse29378_exprs.agg[,gse29378_PLQ_Strat$S1],geneset = pathprint.Hs.gs,statistic = "mean",normalizedScore = F,progressBar = T)
gse29378_PLQ_Strat.SCE$S2=single.chip.enrichment2(exprs = gse29378_exprs.agg[,gse29378_PLQ_Strat$S2],geneset = pathprint.Hs.gs,statistic = "mean",normalizedScore = F,progressBar = T)
gse29378_PLQ_Strat.SCE$S3=single.chip.enrichment2(exprs = gse29378_exprs.agg[,gse29378_PLQ_Strat$S3],geneset = pathprint.Hs.gs,statistic = "mean",normalizedScore = F,progressBar = T)

gse29378_PLQ.diffPathways=vector(mode = "list",length = 4)
names(gse29378_PLQ.diffPathways)=names(gse29378_PLQ_Strat.SCE)
for(i in 2:4){
  if(i==2){
    group=factor(c(rep(names(gse29378_PLQ_Strat.SCE)[1],dim(gse29378_PLQ_Strat.SCE$S0)[2]),rep(names(gse29378_PLQ_Strat.SCE)[i],dim(gse29378_PLQ_Strat.SCE[[i]])[2])))
    design_df=model.matrix(~0+group)
    colnames(design_df)=levels(group)
    contrasts_matrix=makeContrasts(S0-S1,levels = design_df)
    fit=lmFit(object = cbind(data.frame(gse29378_PLQ_Strat.SCE$S0,stringsAsFactors = F),data.frame(gse29378_PLQ_Strat.SCE[[i]],stringsAsFactors = F)),design = design_df)
    fit2=contrasts.fit(fit,contrasts_matrix)
    fit2=eBayes(fit2,trend = T)
    diff_pathways=topTable(fit2,coef = 1,sort.by = "logFC",adjust.method = "BH",number =dim(gse29378_PLQ_Strat.SCE[[i]])[1])%>%rownames_to_column("Pathway_Geneset")
    if(length(diff_pathways)==0){
      gse29378_PLQ.diffPathways[[i]]=0
    }
    else gse29378_PLQ.diffPathways[[i]]=diff_pathways
  }
  if(i==3){
    group=factor(c(rep(names(gse29378_PLQ_Strat.SCE)[1],dim(gse29378_PLQ_Strat.SCE$S0)[2]),rep(names(gse29378_PLQ_Strat.SCE)[i],dim(gse29378_PLQ_Strat.SCE[[i]])[2])))
    design_df=model.matrix(~0+group)
    colnames(design_df)=levels(group)
    contrasts_matrix=makeContrasts(S0-S2,levels = design_df)
    fit=lmFit(object = cbind(data.frame(gse29378_PLQ_Strat.SCE$S0,stringsAsFactors = F),data.frame(gse29378_PLQ_Strat.SCE[[i]],stringsAsFactors = F)),design = design_df)
    fit2=contrasts.fit(fit,contrasts_matrix)
    fit2=eBayes(fit2,trend = T)
    diff_pathways=topTable(fit2,coef = 1,sort.by = "logFC",adjust.method = "BH",number =dim(gse29378_PLQ_Strat.SCE[[i]])[1])%>%rownames_to_column("Pathway_Geneset")
    if(length(diff_pathways)==0){
      gse29378_PLQ.diffPathways[[i]]=0
    }
    else gse29378_PLQ.diffPathways[[i]]=diff_pathways
  }
  if(i==4){
    group=factor(c(rep(names(gse29378_PLQ_Strat.SCE)[1],dim(gse29378_PLQ_Strat.SCE$S0)[2]),rep(names(gse29378_PLQ_Strat.SCE)[i],dim(gse29378_PLQ_Strat.SCE[[i]])[2])))
    design_df=model.matrix(~0+group)
    colnames(design_df)=levels(group)
    contrasts_matrix=makeContrasts(S0-S3,levels = design_df)
    fit=lmFit(object = cbind(data.frame(gse29378_PLQ_Strat.SCE$S0,stringsAsFactors = F),data.frame(gse29378_PLQ_Strat.SCE[[i]],stringsAsFactors = F)),design = design_df)
    fit2=contrasts.fit(fit,contrasts_matrix)
    fit2=eBayes(fit2,trend = T)
    diff_pathways=topTable(fit2,coef = 1,sort.by = "logFC",adjust.method = "BH",number =dim(gse29378_PLQ_Strat.SCE[[i]])[1])%>%rownames_to_column("Pathway_Geneset")
    if(length(diff_pathways)==0){
      gse29378_PLQ.diffPathways[[i]]=0
    }
    else gse29378_PLQ.diffPathways[[i]]=diff_pathways
  }
}
gse29378_PLQ.diffPathways=lapply(gse29378_PLQ.diffPathways[c(2:4)],function(x)x%>%filter(P.Value<=0.01))

gse29378_PLQ_Strat.exprs=gse29378_PLQ.DEGs=vector(mode = "list",length = 4)
names(gse29378_PLQ_Strat.exprs)=names(gse29378_PLQ.DEGs)=names(gse29378_PLQ_Strat.SCE)
gse29378_PLQ_Strat.exprs$S0=gse29378_exprs.agg2[,gse29378_PLQ_Strat$S0]
gse29378_PLQ_Strat.exprs$S1=gse29378_exprs.agg2[,gse29378_PLQ_Strat$S1]
gse29378_PLQ_Strat.exprs$S2=gse29378_exprs.agg2[,gse29378_PLQ_Strat$S2]
gse29378_PLQ_Strat.exprs$S3=gse29378_exprs.agg2[,gse29378_PLQ_Strat$S3]

for(i in 2:4){
  if(i==2){
    group=factor(c(rep(names(gse29378_PLQ_Strat.exprs)[1],dim(gse29378_PLQ_Strat.exprs$S0)[2]),rep(names(gse29378_PLQ_Strat.exprs)[i],dim(gse29378_PLQ_Strat.exprs[[i]])[2])))
    design_df=model.matrix(~0+group)
    colnames(design_df)=levels(group)
    contrasts_matrix=makeContrasts(S0-S1,levels = design_df)
    fit=lmFit(object = cbind(data.frame(gse29378_PLQ_Strat.exprs$S0,stringsAsFactors = F),data.frame(gse29378_PLQ_Strat.exprs[[i]],stringsAsFactors = F)),design = design_df)
    fit2=contrasts.fit(fit,contrasts_matrix)
    fit2=eBayes(fit2,trend = T)
    diff_genes=topTable(fit2,coef = 1,sort.by = "logFC",adjust.method = "BH",number = dim(gse29378_exprs.agg)[1])%>%rownames_to_column("Gene")
    #diff_genes$Gene_HGNC=mapIds(x = org.Hs.eg.db,keys = diff_genes$Gene,keytype = "ENTREZID",column = "SYMBOL")
    
    if(length(diff_genes)==0){
      gse29378_PLQ.DEGs[[i]]=0
    }
    else
      rownames(diff_genes)=rownames(gse29378_exprs.agg)
      gse29378_PLQ.DEGs[[i]]=diff_genes
  }
  if(i==3){
    group=factor(c(rep(names(gse29378_PLQ_Strat.exprs)[1],dim(gse29378_PLQ_Strat.exprs$S0)[2]),rep(names(gse29378_PLQ_Strat.exprs)[i],dim(gse29378_PLQ_Strat.exprs[[i]])[2])))
    design_df=model.matrix(~0+group)
    colnames(design_df)=levels(group)
    contrasts_matrix=makeContrasts(S0-S2,levels = design_df)
    fit=lmFit(object = cbind(data.frame(gse29378_PLQ_Strat.exprs$S0,stringsAsFactors = F),data.frame(gse29378_PLQ_Strat.exprs[[i]],stringsAsFactors = F)),design = design_df)
    fit2=contrasts.fit(fit,contrasts_matrix)
    fit2=eBayes(fit2,trend = T)
    diff_genes=topTable(fit2,coef = 1,sort.by = "logFC",adjust.method = "BH",number = dim(gse29378_exprs.agg)[1])%>%rownames_to_column("Gene")
    #diff_genes$Gene_HGNC=mapIds(x = org.Hs.eg.db,keys = diff_genes$Gene,keytype = "ENTREZID",column = "SYMBOL")
    if(length(diff_genes)==0){
      gse29378_PLQ.DEGs[[i]]=0
    }
    else
    rownames(diff_genes)=rownames(gse29378_exprs.agg)
    gse29378_PLQ.DEGs[[i]]=diff_genes
  }
  if(i==4){
    group=factor(c(rep(names(gse29378_PLQ_Strat.exprs)[1],dim(gse29378_PLQ_Strat.exprs$S0)[2]),rep(names(gse29378_PLQ_Strat.exprs)[i],dim(gse29378_PLQ_Strat.exprs[[i]])[2])))
    design_df=model.matrix(~0+group)
    colnames(design_df)=levels(group)
    contrasts_matrix=makeContrasts(S0-S3,levels = design_df)
    fit=lmFit(object = cbind(data.frame(gse29378_PLQ_Strat.exprs$S0,stringsAsFactors = F),data.frame(gse29378_PLQ_Strat.exprs[[i]],stringsAsFactors = F)),design = design_df)
    fit2=contrasts.fit(fit,contrasts_matrix)
    fit2=eBayes(fit2,trend = T)
    diff_genes=topTable(fit2,coef = 1,sort.by = "logFC",adjust.method = "BH",number = dim(gse29378_exprs.agg)[1])%>%rownames_to_column("Gene")
    #diff_genes$Gene_HGNC=mapIds(x = org.Hs.eg.db,keys = diff_genes$Gene,keytype = "ENTREZID",column = "SYMBOL")
    if(length(diff_genes)==0){
      gse29378_PLQ.DEGs[[i]]=0
    }
    else 
      rownames(diff_genes)=rownames(gse29378_exprs.agg)
    gse29378_PLQ.DEGs[[i]]=diff_genes
  }
}
