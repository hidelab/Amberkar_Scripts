library(org.Hs.eg.db)
library(dplyr)
library(magrittr)
library(data.table)
library(diffcoexp)
library(GEOquery)
jaccard=function(A,B){
  jc=set_cardinality(intersect(A,B))/set_cardinality(union(A,B))
  return(jc)
}
pathprint_membership.genes=lapply(pathprint.Hs.gs,function(x)unname(mapIds(x = org.Hs.eg.db,keys = as.character(x),keytype = "ENTREZID",column = "SYMBOL")))
names(pathprint_membership.genes)[grep(pattern = "/",names(pathprint_membership.genes))]=gsub(pattern = "/",replacement = " ",x = grep(pattern = "/",names(pathprint_membership.genes),value = T))

dbs <- listEnrichrDbs()
kegg_dbs=grep(pattern = "KEGG",x = dbs$libraryName,value = T)[3]
ds = c("chipframe", "genesets","pathprint.Hs.gs" ,"platform.thresholds")
data(list = ds)
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

msbb_gse84422_series_matrix=getGEO(GEO = "GSE84422",GSEMatrix = T)
msbb_gse84422_series_matrix.GPL96=msbb_gse84422_series_matrix$`GSE84422-GPL96_series_matrix.txt.gz`
msbb_gse84422_series_matrix.GPL97=msbb_gse84422_series_matrix$`GSE84422-GPL97_series_matrix.txt.gz`

msbb_gse84422.fData=vector(mode = "list",length = 3)
names(msbb_gse84422.fData)=c("GPL96","GPL97","GPL570")
msbb_gse84422.fData$GPL96=fData(msbb_gse84422_series_matrix.GPL96)[,c(1:2,11:13)]
msbb_gse84422.fData$GPL97=fData(msbb_gse84422_series_matrix.GPL97)[,c(1:2,11:13)]

msbb_gse84422.pData=vector(mode = "list",length = 2)
names(msbb_gse84422.pData)=c("GPL96","GPL97")
msbb_gse84422.pData$GPL96=pData(phenoData(msbb_gse84422_series_matrix.GPL96))
msbb_gse84422.pData$GPL97=pData(phenoData(msbb_gse84422_series_matrix.GPL97))
msbb_gse84422.pData$GPL96$pseudoSampleID=msbb_gse84422.pData$GPL97$pseudoSampleID=paste("pSample",1:dim(msbb_gse84422.pData$GPL96)[1],sep = "")
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$`brain region:ch1`=gsub(pattern = "^ ",replacement = "",x = x$`brain region:ch1`);x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$`clinical dementia rating:ch1`=as.numeric(gsub(pattern = "^ ",replacement = "",x = x$`clinical dementia rating:ch1`));x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$`braak neurofibrillary tangle score:ch1`=as.numeric(gsub(pattern = "^ ",replacement = "",x = x$`braak neurofibrillary tangle score:ch1`));x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$`average neuritic plaque density:ch1`=as.numeric(gsub(pattern = "^ ",replacement = "",x = x$`average neuritic plaque density:ch1`));x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$`neuropathological category:ch1`=gsub(pattern = "^ ",replacement = "",x = x$`neuropathological category:ch1`);x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$`SampleTypeCDR`="OTHER";x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$`SampleTypeBraak`="OTHER";x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$`SampleTypePLQ`="OTHER";x})

msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$SampleTypePLQ[which(x$`average neuritic plaque density:ch1`==0.00)]="PLQ_0";x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$SampleTypePLQ[which(x$`average neuritic plaque density:ch1`>=1.00&x$`average neuritic plaque density:ch1`<=10.50)]="PLQ_Low";x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$SampleTypePLQ[which(x$`average neuritic plaque density:ch1`>=10.60&x$`average neuritic plaque density:ch1`<=20.50)]="PLQ_Med";x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$SampleTypePLQ[which(x$`average neuritic plaque density:ch1`>=20.60)]="PLQ_High";x})

msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$SampleTypeBraak[which(x$`braak neurofibrillary tangle score:ch1`==0)]="Braak0";x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$SampleTypeBraak[which(x$`braak neurofibrillary tangle score:ch1`==1)]="Braak1";x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$SampleTypeBraak[which(x$`braak neurofibrillary tangle score:ch1`==2)]="Braak2";x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$SampleTypeBraak[which(x$`braak neurofibrillary tangle score:ch1`==3)]="Braak3";x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$SampleTypeBraak[which(x$`braak neurofibrillary tangle score:ch1`==4)]="Braak4";x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$SampleTypeBraak[which(x$`braak neurofibrillary tangle score:ch1`==5)]="Braak5";x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$SampleTypeBraak[which(x$`braak neurofibrillary tangle score:ch1`==6)]="Braak6";x})

msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$SampleTypeCDR[which(x$`clinical dementia rating:ch`==0)]="CDR0";x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$SampleTypeCDR[which(x$`clinical dementia rating:ch`==0.5)]="CDR05";x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$SampleTypeCDR[which(x$`clinical dementia rating:ch`==1)]="CDR1";x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$SampleTypeCDR[which(x$`clinical dementia rating:ch`==2)]="CDR2";x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$SampleTypeCDR[which(x$`clinical dementia rating:ch`==3)]="CDR3";x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$SampleTypeCDR[which(x$`clinical dementia rating:ch`==4)]="CDR4";x})
msbb_gse84422.pData=lapply(msbb_gse84422.pData,function(x){x$SampleTypeCDR[which(x$`clinical dementia rating:ch`==5)]="CDR5";x})



msbb_gse84422.exprs=vector(mode = "list",length = 2)
names(msbb_gse84422.exprs)=c("GPL96","GPL97")
msbb_gse84422.exprs$GPL96=exprs(msbb_gse84422_series_matrix.GPL96)
msbb_gse84422.exprs$GPL97=exprs(msbb_gse84422_series_matrix.GPL97)
colnames(msbb_gse84422.exprs$GPL96)=colnames(msbb_gse84422.exprs$GPL97)=msbb_gse84422.pData$GPL96$pseudoSampleID

msbb_gse84422_exprs.GPL96_97=rbind.data.frame(msbb_gse84422.exprs$GPL96,msbb_gse84422.exprs$GPL97)
msbb_gse84422_exprs.GPL96_97$GeneSymbol=c(msbb_gse84422.fData$GPL96$`Gene Symbol`,msbb_gse84422.fData$GPL97$`Gene Symbol`)
msbb_gse84422_exprs.GPL96_97$ENTREZ_GENE_ID=c(msbb_gse84422.fData$GPL96$ENTREZ_GENE_ID,msbb_gse84422.fData$GPL97$ENTREZ_GENE_ID)
msbb_gse84422_exprs_GPL96_97.agg=aggregate.data.frame(x=msbb_gse84422_exprs.GPL96_97[,-which(colnames(msbb_gse84422_exprs.GPL96_97)=="GeneSymbol"|colnames(msbb_gse84422_exprs.GPL96_97)=="ENTREZ_GENE_ID")],by=list(symbol=msbb_gse84422_exprs.GPL96_97$GeneSymbol),mean)
msbb_gse84422_exprs_GPL96_97.agg_Entrez=aggregate.data.frame(x=msbb_gse84422_exprs.GPL96_97[,-which(colnames(msbb_gse84422_exprs.GPL96_97)=="GeneSymbol"|colnames(msbb_gse84422_exprs.GPL96_97)=="ENTREZ_GENE_ID")],by=list(entrez=msbb_gse84422_exprs.GPL96_97$ENTREZ_GENE_ID),mean)
rownames(msbb_gse84422_exprs_GPL96_97.agg)=msbb_gse84422_exprs_GPL96_97.agg$symbol
rownames(msbb_gse84422_exprs_GPL96_97.agg_Entrez)=msbb_gse84422_exprs_GPL96_97.agg_Entrez$entrez
msbb_gse84422_exprs_GPL96_97.agg=msbb_gse84422_exprs_GPL96_97.agg[,-which(colnames(msbb_gse84422_exprs_GPL96_97.agg)=="symbol")]
msbb_gse84422_exprs_GPL96_97.agg_Entrez=msbb_gse84422_exprs_GPL96_97.agg_Entrez[,-which(colnames(msbb_gse84422_exprs_GPL96_97.agg_Entrez)=="entrez")]

msbb_gse84422_GPL96_97_byRegion.exprs=msbb_gse84422_GPL96_97_byRegion.exprs_Entrez=vector(mode = "list",length = 17)
names(msbb_gse84422_GPL96_97_byRegion.exprs)=names(msbb_gse84422_GPL96_97_byRegion.exprs_Entrez)=unique(msbb_gse84422.pData$GPL96$`brain region:ch1`)
msbb_gse84422_GPL96_97_samplesByRegion=lapply(unique(msbb_gse84422.pData$GPL96$`brain region:ch1`),function(y)msbb_gse84422.pData$GPL96%>%filter(`brain region:ch1`==y)%>%pull(pseudoSampleID))
msbb_gse84422_GPL96_97_byRegion.exprs=lapply(msbb_gse84422_GPL96_97_samplesByRegion, function(x)msbb_gse84422_exprs_GPL96_97.agg[,x])
msbb_gse84422_GPL96_97_byRegion.exprs_Entrez=lapply(msbb_gse84422_GPL96_97_samplesByRegion, function(x)msbb_gse84422_exprs_GPL96_97.agg_Entrez[,x])
names(msbb_gse84422_GPL96_97_byRegion.exprs)=names(msbb_gse84422_GPL96_97_byRegion.exprs_Entrez)=unique(msbb_gse84422.pData$GPL96$`brain region:ch1`)

#Regroup brain regions by lobes
lobe_bm_area.map=matrix(nrow=17,ncol=2)
lobe_bm_area.map[1,]=c("Frontal_Lobe","Frontal Pole")
lobe_bm_area.map[2,]=c("Frontal_Lobe","Anterior Cingulate")
lobe_bm_area.map[3,]=c("Frontal_Lobe","Prefrontal Cortex")
lobe_bm_area.map[4,]=c("Occipetal_Lobe","Occipital Visual Cortex")
lobe_bm_area.map[5,]=c("Temporal_Lobe","Inferior Temporal Gyrus")
lobe_bm_area.map[6,]=c("Temporal_Lobe","Middle Temporal Gyrus")
lobe_bm_area.map[7,]=c("Temporal_Lobe","Superior Temporal Gyrus")
lobe_bm_area.map[8,]=c("Parietal_Lobe","Posterior Cingulate Cortex")
lobe_bm_area.map[9,]=c("Temporal_Lobe","Parahippocampal Gyrus")
lobe_bm_area.map[10,]=c("Temporal_Lobe","Temporal Pole")
lobe_bm_area.map[11,]=c("Frontal_Lobe","Precentral Gyrus")
lobe_bm_area.map[12,]=c("Frontal_Lobe","Inferior Frontal Gyrus")
lobe_bm_area.map[13,]=c("Frontal_Lobe","Dorsolateral Prefrontal Cortex")
lobe_bm_area.map[14,]=c("Parietal_Lobe","Superior Parietal Lobule")
lobe_bm_area.map[16,]=c("Temporal_Lobe","Hippocampus")
lobe_bm_area.map[15,]=c("Dorsal_striatum","Putamen")
lobe_bm_area.map[17,]=c("Dorsal_striatum","Caudate Nucleus")
lobe_bm_area.map=data.frame(Lobe=lobe_bm_area.map[,1],BrainRegion=lobe_bm_area.map[,2],stringsAsFactors = F)

msbb_gse84422.byLobe=msbb_gse84422.byLobe_Entrez=vector(mode = "list",length = length(unique(lobe_bm_area.map$Lobe)))
names(msbb_gse84422.byLobe)=names(msbb_gse84422.byLobe_Entrez)=sort(unique(lobe_bm_area.map$Lobe))
msbb_gse84422.byLobe=foreach(i=names(msbb_gse84422.byLobe))%dopar%{
  array_byLobe=do.call("cbind",msbb_gse84422_GPL96_97_byRegion.exprs[lobe_bm_area.map$BrainRegion[lobe_bm_area.map$Lobe==i]])
}
msbb_gse84422.byLobe_Entrez=foreach(i=names(msbb_gse84422.byLobe_Entrez))%dopar%{
  array_byLobe=do.call("cbind",msbb_gse84422_GPL96_97_byRegion.exprs_Entrez[lobe_bm_area.map$BrainRegion[lobe_bm_area.map$Lobe==i]])
}
names(msbb_gse84422.byLobe)=names(msbb_gse84422.byLobe_Entrez)=sort(unique(lobe_bm_area.map$Lobe))

#Remove lobe names from sample name string
msbb_gse84422.byLobe=lapply(msbb_gse84422.byLobe,function(x){colnames(x)=unique(str_match(string=colnames(x),pattern = "pSample\\d+"));x})
msbb_gse84422.byLobe_Entrez=lapply(msbb_gse84422.byLobe_Entrez,function(x){colnames(x)=unique(str_match(string=colnames(x),pattern = "pSample\\d+"));x})

msbb_gse84422.CDR_Strat=vector(mode = "list",length = length(unique(lobe_bm_area.map$Lobe)))
names(msbb_gse84422.CDR_Strat)=sort(unique(lobe_bm_area.map$Lobe))
msbb_gse84422.CDR_Strat$Dorsal_striatum=msbb_gse84422.CDR_Strat$Frontal_Lobe=msbb_gse84422.CDR_Strat$Occipetal_Lobe=msbb_gse84422.CDR_Strat$Parietal_Lobe=msbb_gse84422.CDR_Strat$Temporal_Lobe=vector(mode = "list",length = length(unique(msbb_gse84422.pData$GPL96$SampleTypeCDR)))
names(msbb_gse84422.CDR_Strat$Dorsal_striatum)=names(msbb_gse84422.CDR_Strat$Frontal_Lobe)=names(msbb_gse84422.CDR_Strat$Occipetal_Lobe)=names(msbb_gse84422.CDR_Strat$Parietal_Lobe)=names(msbb_gse84422.CDR_Strat$Temporal_Lobe)=sort(unique(msbb_gse84422.pData$GPL96$SampleTypeCDR))
for(i in 1:length(msbb_gse84422.CDR_Strat)){
  for(j in 1:length(msbb_gse84422.CDR_Strat[[i]]))
  msbb_gse84422.CDR_Strat[[i]][[j]]=sort(msbb_gse84422.pData$GPL96%>%filter(SampleTypeCDR==names(msbb_gse84422.CDR_Strat[[i]])[[j]])%>%dplyr::select(c(`brain region:ch1`,pseudoSampleID))%>%dplyr::filter(`brain region:ch1`%in%lobe_bm_area.map$BrainRegion[lobe_bm_area.map$Lobe%in%names(msbb_gse84422.CDR_Strat)[i]])%>%pull(pseudoSampleID))
}
#CDR, grouped samples
msbb_gse84422.CDR_Samples=vector(mode = "list",length = length(sort(unique(msbb_gse84422.pData$GPL96$`clinical dementia rating:ch1`))))
names(msbb_gse84422.CDR_Samples)=sort(unique(msbb_gse84422.pData$GPL96$`clinical dementia rating:ch1`))
msbb_gse84422.CDR_Samples=lapply(names(msbb_gse84422.CDR_Samples),function(x)msbb_gse84422.pData$GPL96%>%filter(`clinical dementia rating:ch1`==x)%>%pull(pseudoSampleID))
names(msbb_gse84422.CDR_Samples)=sort(unique(msbb_gse84422.pData$GPL96$`clinical dementia rating:ch1`))

#CDR-wise SCE score
msbb_gse84422.CDR_SCE=vector(mode = "list",length = length(names(msbb_gse84422.byLobe_Entrez)))
names(msbb_gse84422.CDR_SCE)=names(msbb_gse84422.byLobe_Entrez)
msbb_gse84422.CDR_SCE$Dorsal_striatum=msbb_gse84422.CDR_SCE$Frontal_Lobe=msbb_gse84422.CDR_SCE$Occipetal_Lobe=msbb_gse84422.CDR_SCE$Parietal_Lobe=msbb_gse84422.CDR_SCE$Temporal_Lobe=vector(mode = "list",length = length(sort(unique(msbb_gse84422.pData$GPL96$SampleTypeCDR))))
names(msbb_gse84422.CDR_SCE$Dorsal_striatum)=names(msbb_gse84422.CDR_SCE$Frontal_Lobe)=names(msbb_gse84422.CDR_SCE$Occipetal_Lobe)=names(msbb_gse84422.CDR_SCE$Parietal_Lobe)=names(msbb_gse84422.CDR_SCE$Temporal_Lobe)=sort(unique(msbb_gse84422.pData$GPL96$`clinical dementia rating:ch1`))

for(l in 1:length(names(msbb_gse84422.CDR_SCE))){
  for(c in 1:length(msbb_gse84422.CDR_SCE[[l]])){
    msbb_gse84422.CDR_SCE[[l]][[c]]=single.chip.enrichment2(exprs = msbb_gse84422.byLobe_Entrez[[l]][,colnames(msbb_gse84422.byLobe_Entrez[[l]])%in%msbb_gse84422.CDR_Samples[[c]]],geneset = pathprint.Hs.gs,statistic = "mean",normalizedScore = F,progressBar = T)
  }
}

msbb_gse84422.CDR_SCE=lapply(msbb_gse84422.CDR_SCE,function(x){x=x[c('0','0.5','1','2','3','4','5')];x})


msbb_gse84422.CDR_DEPs=vector(mode = "list",length = length(names(msbb_gse84422.byLobe)))
names(msbb_gse84422.CDR_DEPs)=names(msbb_gse84422.byLobe)
msbb_gse84422.CDR_DEPs$Dorsal_striatum=msbb_gse84422.CDR_DEPs$Frontal_Lobe=msbb_gse84422.CDR_DEPs$Occipetal_Lobe=msbb_gse84422.CDR_DEPs$Parietal_Lobe=msbb_gse84422.CDR_DEPs$Temporal_Lobe=vector(mode = "list",length = length(sort(unique(msbb_gse84422.pData$GPL96$SampleTypeCDR))[-1]))
names(msbb_gse84422.CDR_DEPs$Dorsal_striatum)=names(msbb_gse84422.CDR_DEPs$Frontal_Lobe)=names(msbb_gse84422.CDR_DEPs$Occipetal_Lobe)=names(msbb_gse84422.CDR_DEPs$Parietal_Lobe)=names(msbb_gse84422.CDR_DEPs$Temporal_Lobe)=sort(unique(msbb_gse84422.pData$GPL96$`clinical dementia rating:ch1`))[-1]
for(l in 1:length(names(msbb_gse84422.CDR_SCE))){
  for(c in 2:length(names(msbb_gse84422.CDR_SCE[[l]]))){
    c_exprs=data.frame(msbb_gse84422.CDR_SCE[[l]]$`0`)
    d_exprs=data.frame(msbb_gse84422.CDR_SCE[[l]][[c]])
    group=factor(c(rep("CONTROL",dim(c_exprs)[2]),rep("AD",dim(d_exprs)[2])))
    design_df=model.matrix(~0+group)
    colnames(design_df)=levels(group)
    contrasts_matrix=makeContrasts(CONTROL-AD,levels = design_df)
    fit=lmFit(object = cbind.data.frame(c_exprs,d_exprs),design = design_df)
    fit2=contrasts.fit(fit,contrasts_matrix)
    fit2=eBayes(fit2,trend = T)
    msbb_gse84422.CDR_DEPs[[l]][[c-1]]=topTable(fit = fit2,coef = 1,number = dim(msbb_gse84422.CDR_SCE[[l]][[c]])[1],adjust.method = "BH",p.value = 0.2)%>%rownames_to_column("Geneset_Pathway")
    msbb_gse84422.CDR_DEPs[[l]][[c-1]]$Pathway_MemberGenes=unlist(lapply(msbb_gse84422.CDR_DEPs[[l]][[c-1]]$Geneset_Pathway,function(x)paste(unlist(unname(pathprint_membership.genes[names(pathprint_membership.genes)%in%x])),collapse = ",")))
    fwrite(msbb_gse84422.CDR_DEPs[[l]][[c-1]],paste(names(msbb_gse84422.CDR_DEPs)[[l]],names(msbb_gse84422.CDR_DEPs[[l]])[[c-1]],"CDR_DEPs.txt",sep = "_"),sep = "\t",col.names = T,row.names = F,quote = F)
    
  }
}

#NOSS DEPs
msbb_gse84422_CDR_DEP.Uniq=vector(mode = "list",length = length(names(msbb_gse84422.CDR_DEPs))-1)
names(msbb_gse84422_CDR_DEP.Uniq)=names(msbb_gse84422.CDR_DEPs)[-3]
msbb_gse84422_CDR_DEP.Uniq$Dorsal_striatum=vector(mode = "list",length = length(names(msbb_gse84422.CDR_DEPs$Dorsal_striatum)))
msbb_gse84422_CDR_DEP.Uniq$Frontal_Lobe=vector(mode = "list",length = length(names(msbb_gse84422.CDR_DEPs$Frontal_Lobe)))
msbb_gse84422_CDR_DEP.Uniq$Parietal_Lobe=vector(mode = "list",length = length(names(msbb_gse84422.CDR_DEPs$Parietal_Lobe)))
msbb_gse84422_CDR_DEP.Uniq$Temporal_Lobe=vector(mode = "list",length = length(names(msbb_gse84422.CDR_DEPs$Temporal_Lobe)))

names(msbb_gse84422_CDR_DEP.Uniq$Dorsal_striatum)=names(msbb_gse84422.CDR_DEPs$Dorsal_striatum)
names(msbb_gse84422_CDR_DEP.Uniq$Frontal_Lobe)=names(msbb_gse84422.CDR_DEPs$Frontal_Lobe)
names(msbb_gse84422_CDR_DEP.Uniq$Parietal_Lobe)=names(msbb_gse84422.CDR_DEPs$Parietal_Lobe)
names(msbb_gse84422_CDR_DEP.Uniq$Temporal_Lobe)=names(msbb_gse84422.CDR_DEPs$Temporal_Lobe)

for(b in 1:length(names(msbb_gse84422_CDR_DEP.Uniq$Dorsal_striatum))){
  msbb_gse84422_CDR_DEP.Uniq$Dorsal_striatum[[b]]=msbb_gse84422.CDR_DEPs$Dorsal_striatum[[b]][msbb_gse84422.CDR_DEPs$Dorsal_striatum[[b]]$Gene%in%setdiff(msbb_gse84422.CDR_DEPs$Dorsal_striatum[[b]]$Gene,c(Reduce(union,lapply(lapply(msbb_gse84422.CDR_DEPs$Dorsal_striatum,function(x)x%>%pull(Geneset_Pathway)),function(x)intersect(x,msbb_gse84422.CDR_DEPs$Dorsal_striatum[[b]]$Gene))[-b]),Reduce(intersect,lapply(msbb_gse84422.CDR_DEPs$Dorsal_striatum,function(x)x%>%pull(Geneset_Pathway))))),]
  fwrite(msbb_gse84422_CDR_DEP.Uniq$Dorsal_striatum[[b]],paste(names(msbb_gse84422.CDR_DEPs)[1],names(msbb_gse84422_CDR_DEP.Uniq$Dorsal_striatum)[b],"NOSS_CDR_DEPs.txt",sep = "_"),sep = "\t",col.names = T,row.names = F,quote = F)
}

for(b in 1:length(names(msbb_gse84422_CDR_DEP.Uniq$Frontal_Lobe))){
  msbb_gse84422_CDR_DEP.Uniq$Frontal_Lobe[[b]]=msbb_gse84422.CDR_DEPs$Frontal_Lobe[[b]][msbb_gse84422.CDR_DEPs$Frontal_Lobe[[b]]$Gene%in%setdiff(msbb_gse84422.CDR_DEPs$Frontal_Lobe[[b]]$Gene,c(Reduce(union,lapply(lapply(msbb_gse84422.CDR_DEPs$Frontal_Lobe,function(x)x%>%pull(Geneset_Pathway)),function(x)intersect(x,msbb_gse84422.CDR_DEPs$Frontal_Lobe[[b]]$Gene))[-b]),Reduce(intersect,lapply(msbb_gse84422.CDR_DEPs$Frontal_Lobe,function(x)x%>%pull(Geneset_Pathway))))),]
  fwrite(msbb_gse84422_CDR_DEP.Uniq$Frontal_Lobe[[b]],paste(names(msbb_gse84422.CDR_DEPs)[2],names(msbb_gse84422_CDR_DEP.Uniq$Frontal_Lobe)[b],"NOSS_CDR_DEPs.txt",sep = "_"),sep = "\t",col.names = T,row.names = F,quote = F)
}

for(b in 1:length(names(msbb_gse84422_CDR_DEP.Uniq$Parietal_Lobe))){
  msbb_gse84422_CDR_DEP.Uniq$Parietal_Lobe[[b]]=msbb_gse84422.CDR_DEPs$Parietal_Lobe[[b]][msbb_gse84422.CDR_DEPs$Parietal_Lobe[[b]]$Gene%in%setdiff(msbb_gse84422.CDR_DEPs$Parietal_Lobe[[b]]$Gene,c(Reduce(union,lapply(lapply(msbb_gse84422.CDR_DEPs$Parietal_Lobe,function(x)x%>%pull(Geneset_Pathway)),function(x)intersect(x,msbb_gse84422.CDR_DEPs$Parietal_Lobe[[b]]$Gene))[-b]),Reduce(intersect,lapply(msbb_gse84422.CDR_DEPs$Parietal_Lobe,function(x)x%>%pull(Geneset_Pathway))))),]
  fwrite(msbb_gse84422_CDR_DEP.Uniq$Parietal_Lobe[[b]],paste(names(msbb_gse84422.CDR_DEPs)[4],names(msbb_gse84422_CDR_DEP.Uniq$Parietal_Lobe)[b],"NOSS_CDR_DEPs.txt",sep = "_"),sep = "\t",col.names = T,row.names = F,quote = F)
}

for(b in 1:length(names(msbb_gse84422_CDR_DEP.Uniq$Temporal_Lobe))){
  msbb_gse84422_CDR_DEP.Uniq$Temporal_Lobe[[b]]=msbb_gse84422.CDR_DEPs$Temporal_Lobe[[b]][msbb_gse84422.CDR_DEPs$Temporal_Lobe[[b]]$Gene%in%setdiff(msbb_gse84422.CDR_DEPs$Temporal_Lobe[[b]]$Gene,c(Reduce(union,lapply(lapply(msbb_gse84422.CDR_DEPs$Temporal_Lobe,function(x)x%>%pull(Geneset_Pathway)),function(x)intersect(x,msbb_gse84422.CDR_DEPs$Temporal_Lobe[[b]]$Gene))[-b]),Reduce(intersect,lapply(msbb_gse84422.CDR_DEPs$Temporal_Lobe,function(x)x%>%pull(Geneset_Pathway))))),]    
  fwrite(msbb_gse84422_CDR_DEP.Uniq$Temporal_Lobe[[b]],paste(names(msbb_gse84422.CDR_DEPs)[5],names(msbb_gse84422_CDR_DEP.Uniq$Temporal_Lobe)[b],"NOSS_CDR_DEPs.txt",sep = "_"),sep = "\t",col.names = T,row.names = F,quote = F)  
}

msbb_gse84422_CDR_DEP.Uniq$Dorsal_striatum=Filter(f = function(x)dim(x)[1]>0,msbb_gse84422_CDR_DEP.Uniq$Dorsal_striatum)
msbb_gse84422_CDR_DEP.Uniq$Frontal_Lobe=Filter(f = function(x)dim(x)[1]>0,msbb_gse84422_CDR_DEP.Uniq$Frontal_Lobe)
msbb_gse84422_CDR_DEP.Uniq$Parietal_Lobe=Filter(f = function(x)dim(x)[1]>0,msbb_gse84422_CDR_DEP.Uniq$Parietal_Lobe)
msbb_gse84422_CDR_DEP.Uniq$Temporal_Lobe=Filter(f = function(x)dim(x)[1]>0,msbb_gse84422_CDR_DEP.Uniq$Temporal_Lobe)


#DEG analysis for CDR
msbb_gse84422.CDR_DEGs=vector(mode = "list",length = length(names(msbb_gse84422.byLobe)))
names(msbb_gse84422.CDR_DEGs)=names(msbb_gse84422.byLobe)
msbb_gse84422.CDR_DEGs$Dorsal_striatum=msbb_gse84422.CDR_DEGs$Frontal_Lobe=msbb_gse84422.CDR_DEGs$Occipetal_Lobe=msbb_gse84422.CDR_DEGs$Parietal_Lobe=msbb_gse84422.CDR_DEGs$Temporal_Lobe=vector(mode = "list",length = length(sort(unique(msbb_gse84422.pData$GPL96$SampleTypeCDR))[-1]))
names(msbb_gse84422.CDR_DEGs$Dorsal_striatum)=names(msbb_gse84422.CDR_DEGs$Frontal_Lobe)=names(msbb_gse84422.CDR_DEGs$Occipetal_Lobe)=names(msbb_gse84422.CDR_DEGs$Parietal_Lobe)=names(msbb_gse84422.CDR_DEGs$Temporal_Lobe)=sort(unique(msbb_gse84422.pData$GPL96$`clinical dementia rating:ch1`))[-1]

for(l in 1:length(names(msbb_gse84422.byLobe))){
  for(c in 2:length(msbb_gse84422.CDR_Strat[[l]])){
    c_exprs=msbb_gse84422.byLobe[[l]][,msbb_gse84422.CDR_Strat[[l]]$CDR0]
    d_exprs=msbb_gse84422.byLobe[[l]][,msbb_gse84422.CDR_Strat[[l]][[c]]]
    group=factor(c(rep("CONTROL",dim(c_exprs)[2]),rep("AD",dim(d_exprs)[2])))
    design_df=model.matrix(~0+group)
    colnames(design_df)=levels(group)
    contrasts_matrix=makeContrasts(CONTROL-AD,levels = design_df)
    fit=lmFit(object = cbind.data.frame(c_exprs,d_exprs),design = design_df)
    fit2=contrasts.fit(fit,contrasts_matrix)
    fit2=eBayes(fit2,trend = T)
    msbb_gse84422.CDR_DEGs[[l]][[c-1]]=topTable(fit = fit2,coef = 1,number = dim(msbb_gse84422.byLobe[[l]])[1],adjust.method = "BH",p.value = 0.1)%>%rownames_to_column("Gene")
    fwrite(msbb_gse84422.CDR_DEGs[[l]][[c-1]],paste(names(msbb_gse84422.CDR_DEGs)[[l]],names(msbb_gse84422.CDR_DEGs[[l]])[[c-1]],"CDR_DEGs.txt",sep = "_"),sep = "\t",col.names = T,row.names = F,quote = F)
  }
}

msbb_gse84422.CDR_DEGs$Dorsal_striatum=Filter(f = function(x)dim(x)[1]>0,x = msbb_gse84422.CDR_DEGs$Dorsal_striatum)
msbb_gse84422.CDR_DEGs$Frontal_Lobe=Filter(f = function(x)dim(x)[1]>0,x = msbb_gse84422.CDR_DEGs$Frontal_Lobe)
msbb_gse84422.CDR_DEGs$Occipetal_Lobe=Filter(f = function(x)dim(x)[1]>0,x = msbb_gse84422.CDR_DEGs$Occipetal_Lobe)
msbb_gse84422.CDR_DEGs$Parietal_Lobe=Filter(f = function(x)dim(x)[1]>0,x = msbb_gse84422.CDR_DEGs$Parietal_Lobe)
msbb_gse84422.CDR_DEGs$Temporal_Lobe=Filter(f = function(x)dim(x)[1]>0,x = msbb_gse84422.CDR_DEGs$Temporal_Lobe)

#Non-overlapping DEGs

msbb_gse84422_CDR_DEG.Uniq=vector(mode = "list",length = length(names(msbb_gse84422.CDR_DEGs))-1)
names(msbb_gse84422_CDR_DEG.Uniq)=names(msbb_gse84422.CDR_DEGs)[-3]
msbb_gse84422_CDR_DEG.Uniq$Dorsal_striatum=vector(mode = "list",length = length(names(msbb_gse84422.CDR_DEGs$Dorsal_striatum)))
msbb_gse84422_CDR_DEG.Uniq$Frontal_Lobe=vector(mode = "list",length = length(names(msbb_gse84422.CDR_DEGs$Frontal_Lobe)))
msbb_gse84422_CDR_DEG.Uniq$Parietal_Lobe=vector(mode = "list",length = length(names(msbb_gse84422.CDR_DEGs$Parietal_Lobe)))
msbb_gse84422_CDR_DEG.Uniq$Temporal_Lobe=vector(mode = "list",length = length(names(msbb_gse84422.CDR_DEGs$Temporal_Lobe)))
                                                
names(msbb_gse84422_CDR_DEG.Uniq$Dorsal_striatum)=names(msbb_gse84422.CDR_DEGs$Dorsal_striatum)
names(msbb_gse84422_CDR_DEG.Uniq$Frontal_Lobe)=names(msbb_gse84422.CDR_DEGs$Frontal_Lobe)
names(msbb_gse84422_CDR_DEG.Uniq$Parietal_Lobe)=names(msbb_gse84422.CDR_DEGs$Parietal_Lobe)
names(msbb_gse84422_CDR_DEG.Uniq$Temporal_Lobe)=names(msbb_gse84422.CDR_DEGs$Temporal_Lobe)

for(b in 1:length(names(msbb_gse84422_CDR_DEG.Uniq$Dorsal_striatum))){
  msbb_gse84422_CDR_DEG.Uniq$Dorsal_striatum[[b]]=msbb_gse84422.CDR_DEGs$Dorsal_striatum[[b]][msbb_gse84422.CDR_DEGs$Dorsal_striatum[[b]]$Gene%in%setdiff(msbb_gse84422.CDR_DEGs$Dorsal_striatum[[b]]$Gene,c(Reduce(union,lapply(lapply(msbb_gse84422.CDR_DEGs$Dorsal_striatum,function(x)x%>%pull(Gene)),function(x)intersect(x,msbb_gse84422.CDR_DEGs$Dorsal_striatum[[b]]$Gene))[-b]),Reduce(intersect,lapply(msbb_gse84422.CDR_DEGs$Dorsal_striatum,function(x)x%>%pull(Gene))))),]
  fwrite(msbb_gse84422_CDR_DEG.Uniq$Dorsal_striatum[[b]],paste(names(msbb_gse84422.CDR_DEGs)[1],names(msbb_gse84422_CDR_DEG.Uniq$Dorsal_striatum)[[b]],"NOSS_CDR_DEGs.txt",sep = "_"),sep = "\t",col.names = T,row.names = F,quote = F)
}

for(b in 1:length(names(msbb_gse84422_CDR_DEG.Uniq$Frontal_Lobe))){
  msbb_gse84422_CDR_DEG.Uniq$Frontal_Lobe[[b]]=msbb_gse84422.CDR_DEGs$Frontal_Lobe[[b]][msbb_gse84422.CDR_DEGs$Frontal_Lobe[[b]]$Gene%in%setdiff(msbb_gse84422.CDR_DEGs$Frontal_Lobe[[b]]$Gene,c(Reduce(union,lapply(lapply(msbb_gse84422.CDR_DEGs$Frontal_Lobe,function(x)x%>%pull(Gene)),function(x)intersect(x,msbb_gse84422.CDR_DEGs$Frontal_Lobe[[b]]$Gene))[-b]),Reduce(intersect,lapply(msbb_gse84422.CDR_DEGs$Frontal_Lobe,function(x)x%>%pull(Gene))))),]
  fwrite(msbb_gse84422_CDR_DEG.Uniq$Frontal_Lobe[[b]],paste(names(msbb_gse84422.CDR_DEGs)[2],names(msbb_gse84422_CDR_DEG.Uniq$Frontal_Lobe)[[b]],"NOSS_CDR_DEGs.txt",sep = "_"),sep = "\t",col.names = T,row.names = F,quote = F)
}

for(b in 1:length(names(msbb_gse84422_CDR_DEG.Uniq$Parietal_Lobe))){
  msbb_gse84422_CDR_DEG.Uniq$Parietal_Lobe[[b]]=msbb_gse84422.CDR_DEGs$Parietal_Lobe[[b]][msbb_gse84422.CDR_DEGs$Parietal_Lobe[[b]]$Gene%in%setdiff(msbb_gse84422.CDR_DEGs$Parietal_Lobe[[b]]$Gene,c(Reduce(union,lapply(lapply(msbb_gse84422.CDR_DEGs$Parietal_Lobe,function(x)x%>%pull(Gene)),function(x)intersect(x,msbb_gse84422.CDR_DEGs$Parietal_Lobe[[b]]$Gene))[-b]),Reduce(intersect,lapply(msbb_gse84422.CDR_DEGs$Parietal_Lobe,function(x)x%>%pull(Gene))))),]
    fwrite(msbb_gse84422_CDR_DEG.Uniq$Parietal_Lobe[[b]],paste(names(msbb_gse84422.CDR_DEGs)[4],names(msbb_gse84422_CDR_DEG.Uniq$Parietal_Lobe)[[b]],"NOSS_CDR_DEGs.txt",sep = "_"),sep = "\t",col.names = T,row.names = F,quote = F)
}

for(b in 1:length(names(msbb_gse84422_CDR_DEG.Uniq$Temporal_Lobe))){
  msbb_gse84422_CDR_DEG.Uniq$Temporal_Lobe[[b]]=msbb_gse84422.CDR_DEGs$Temporal_Lobe[[b]][msbb_gse84422.CDR_DEGs$Temporal_Lobe[[b]]$Gene%in%setdiff(msbb_gse84422.CDR_DEGs$Temporal_Lobe[[b]]$Gene,c(Reduce(union,lapply(lapply(msbb_gse84422.CDR_DEGs$Temporal_Lobe,function(x)x%>%pull(Gene)),function(x)intersect(x,msbb_gse84422.CDR_DEGs$Temporal_Lobe[[b]]$Gene))[-b]),Reduce(intersect,lapply(msbb_gse84422.CDR_DEGs$Temporal_Lobe,function(x)x%>%pull(Gene))))),]    
  fwrite(msbb_gse84422_CDR_DEG.Uniq$Temporal_Lobe[[b]],paste(names(msbb_gse84422.CDR_DEGs)[5],names(msbb_gse84422_CDR_DEG.Uniq$Temporal_Lobe)[[b]],"NOSS_CDR_DEGs.txt",sep = "_"),sep = "\t",col.names = T,row.names = F,quote = F)
}

msbb_gse84422_CDR_DEG.Uniq$Dorsal_striatum=Filter(f = function(x)dim(x)[1]>0,msbb_gse84422_CDR_DEG.Uniq$Dorsal_striatum)
msbb_gse84422_CDR_DEG.Uniq$Frontal_Lobe=Filter(f = function(x)dim(x)[1]>0,msbb_gse84422_CDR_DEG.Uniq$Frontal_Lobe)
msbb_gse84422_CDR_DEG.Uniq$Parietal_Lobe=Filter(f = function(x)dim(x)[1]>0,msbb_gse84422_CDR_DEG.Uniq$Parietal_Lobe)
msbb_gse84422_CDR_DEG.Uniq$Temporal_Lobe=Filter(f = function(x)dim(x)[1]>0,msbb_gse84422_CDR_DEG.Uniq$Temporal_Lobe)

# for(t in 1:length(msbb_gse84422.CDR_DEGs)){
#   for(b in 1:length(msbb_gse84422.CDR_DEGs[[t]])){
#     msbb_gse84422_CDR_DEG.Uniq[[t]][[b]]=msbb_gse84422.CDR_DEGs[[t]][[b]][msbb_gse84422.CDR_DEGs[[t]][[b]]$Gene%in%setdiff(msbb_gse84422.CDR_DEGs[[t]][[b]]$Gene,c(Reduce(union,lapply(lapply(msbb_gse84422.CDR_DEGs[[t]],function(x)x%>%pull(Gene)),function(x)intersect(x,msbb_gse84422.CDR_DEGs[[t]][[b]]$Gene))[-b]),Reduce(intersect,lapply(msbb_gse84422.CDR_DEGs[[t]],function(x)x%>%pull(Gene))))),]  
#     
#   }
#   
# }




msbb_gse84422_CDR_DEGs.KEGG=vector(mode = "list",length = 5)
names(msbb_gse84422_CDR_DEGs.KEGG)=names(msbb_gse84422.CDR_DEGs)
msbb_gse84422_CDR_DEGs.KEGG$Dorsal_striatum=msbb_gse84422_CDR_DEGs.KEGG$Frontal_Lobe=msbb_gse84422_CDR_DEGs.KEGG$Occipetal_Lobe=msbb_gse84422_CDR_DEGs.KEGG$Parietal_Lobe=msbb_gse84422_CDR_DEGs.KEGG$Temporal_Lobe=vector(mode = "list",length = length(msbb_gse84422.CDR_DEGs$Dorsal_striatum))
names(msbb_gse84422_CDR_DEGs.KEGG$Dorsal_striatum)=names(msbb_gse84422_CDR_DEGs.KEGG$Frontal_Lobe)=names(msbb_gse84422_CDR_DEGs.KEGG$Occipetal_Lobe)=names(msbb_gse84422_CDR_DEGs.KEGG$Parietal_Lobe)=names(msbb_gse84422_CDR_DEGs.KEGG$Temporal_Lobe)=names(msbb_gse84422.CDR_DEGs$Dorsal_striatum)

msbb_gse84422_CDR_DEGs.KEGG$Dorsal_striatum=lapply(msbb_gse84422.CDR_DEGs$Dorsal_striatum,function(x)enrichr(genes = x$Gene,databases = kegg_dbs)[[1]]%>%dplyr::filter(Adjusted.P.value<=0.2)%>%filter(grepl(pattern = ";",Genes)))
msbb_gse84422_CDR_DEGs.KEGG$Frontal_Lobe=lapply(msbb_gse84422.CDR_DEGs$Frontal_Lobe,function(x)enrichr(genes = x$Gene,databases = kegg_dbs)[[1]]%>%dplyr::filter(Adjusted.P.value<=0.2)%>%filter(grepl(pattern = ";",Genes)))
msbb_gse84422_CDR_DEGs.KEGG$Occipetal_Lobe=lapply(msbb_gse84422.CDR_DEGs$Occipetal_Lobe,function(x)enrichr(genes = x$Gene,databases = kegg_dbs)[[1]]%>%dplyr::filter(Adjusted.P.value<=0.2)%>%filter(grepl(pattern = ";",Genes)))
msbb_gse84422_CDR_DEGs.KEGG$Parietal_Lobe=lapply(msbb_gse84422.CDR_DEGs$Parietal_Lobe,function(x)enrichr(genes = x$Gene,databases = kegg_dbs)[[1]]%>%dplyr::filter(Adjusted.P.value<=0.2)%>%filter(grepl(pattern = ";",Genes)))
msbb_gse84422_CDR_DEGs.KEGG$Temporal_Lobe=lapply(msbb_gse84422.CDR_DEGs$Temporal_Lobe,function(x)enrichr(genes = x$Gene,databases = kegg_dbs)[[1]]%>%dplyr::filter(Adjusted.P.value<=0.2)%>%filter(grepl(pattern = ";",Genes)))

msbb_gse84422_CDR_DEG_Uniq.KEGG=vector(mode = "list",length = 4)
names(msbb_gse84422_CDR_DEG_Uniq.KEGG)=names(msbb_gse84422_CDR_DEG.Uniq)
msbb_gse84422_CDR_DEG_Uniq.KEGG$Dorsal_striatum=msbb_gse84422_CDR_DEG_Uniq.KEGG$Frontal_Lobe=msbb_gse84422_CDR_DEG_Uniq.KEGG$Parietal_Lobe=msbb_gse84422_CDR_DEG_Uniq.KEGG$Temporal_Lobe=vector(mode = "list",length = length(msbb_gse84422_CDR_DEG.Uniq$Dorsal_striatum))
names(msbb_gse84422_CDR_DEG_Uniq.KEGG$Dorsal_striatum)=names(msbb_gse84422_CDR_DEG_Uniq.KEGG$Frontal_Lobe)=names(msbb_gse84422_CDR_DEG_Uniq.KEGG$Parietal_Lobe)=names(msbb_gse84422_CDR_DEG_Uniq.KEGG$Temporal_Lobe)=names(msbb_gse84422_CDR_DEG.Uniq$Dorsal_striatum)

msbb_gse84422_CDR_DEG_Uniq.KEGG$Dorsal_striatum=lapply(msbb_gse84422_CDR_DEG.Uniq$Dorsal_striatum,function(x)enrichr(genes = x$Gene,databases = kegg_dbs)[[1]]%>%dplyr::filter(Adjusted.P.value<=0.2)%>%filter(grepl(pattern = ";",Genes)))
msbb_gse84422_CDR_DEG_Uniq.KEGG$Frontal_Lobe=lapply(msbb_gse84422_CDR_DEG.Uniq$Frontal_Lobe,function(x)enrichr(genes = x$Gene,databases = kegg_dbs)[[1]]%>%dplyr::filter(Adjusted.P.value<=0.2)%>%filter(grepl(pattern = ";",Genes)))
msbb_gse84422_CDR_DEG_Uniq.KEGG$Parietal_Lobe=lapply(msbb_gse84422_CDR_DEG.Uniq$Parietal_Lobe,function(x)enrichr(genes = x$Gene,databases = kegg_dbs)[[1]]%>%dplyr::filter(Adjusted.P.value<=0.2)%>%filter(grepl(pattern = ";",Genes)))
msbb_gse84422_CDR_DEG_Uniq.KEGG$Temporal_Lobe=lapply(msbb_gse84422_CDR_DEG.Uniq$Temporal_Lobe,function(x)enrichr(genes = x$Gene,databases = kegg_dbs)[[1]]%>%dplyr::filter(Adjusted.P.value<=0.2)%>%filter(grepl(pattern = ";",Genes)))

msbb_gse84422_CDR_DEG_Uniq.KEGG$Dorsal_striatum=Filter(f = function(x)dim(x)[1]>0,msbb_gse84422_CDR_DEG_Uniq.KEGG$Dorsal_striatum)
msbb_gse84422_CDR_DEG_Uniq.KEGG$Frontal_Lobe=Filter(f = function(x)dim(x)[1]>0,msbb_gse84422_CDR_DEG_Uniq.KEGG$Frontal_Lobe)
msbb_gse84422_CDR_DEG_Uniq.KEGG$Parietal_Lobe=Filter(f = function(x)dim(x)[1]>0,msbb_gse84422_CDR_DEG_Uniq.KEGG$Parietal_Lobe)
msbb_gse84422_CDR_DEG_Uniq.KEGG$Temporal_Lobe=Filter(f = function(x)dim(x)[1]>0,msbb_gse84422_CDR_DEG_Uniq.KEGG$Temporal_Lobe)


#Braak strat
msbb_gse84422.Braak_Strat=vector(mode = "list",length = length(unique(lobe_bm_area.map$Lobe)))
names(msbb_gse84422.Braak_Strat)=sort(unique(lobe_bm_area.map$Lobe))
msbb_gse84422.Braak_Strat$Dorsal_striatum=msbb_gse84422.Braak_Strat$Frontal_Lobe=msbb_gse84422.Braak_Strat$Occipetal_Lobe=msbb_gse84422.Braak_Strat$Parietal_Lobe=msbb_gse84422.Braak_Strat$Temporal_Lobe=vector(mode = "list",length = length(unique(msbb_gse84422.pData$GPL96$SampleTypeBraak)))
names(msbb_gse84422.Braak_Strat$Dorsal_striatum)=names(msbb_gse84422.Braak_Strat$Frontal_Lobe)=names(msbb_gse84422.Braak_Strat$Occipetal_Lobe)=names(msbb_gse84422.Braak_Strat$Parietal_Lobe)=names(msbb_gse84422.Braak_Strat$Temporal_Lobe)=sort(unique(msbb_gse84422.pData$GPL96$SampleTypeBraak))
for(i in 1:length(msbb_gse84422.Braak_Strat)){
  for(j in 1:length(msbb_gse84422.Braak_Strat[[i]]))
    msbb_gse84422.Braak_Strat[[i]][[j]]=sort(msbb_gse84422.pData$GPL96%>%filter(SampleTypeBraak==names(msbb_gse84422.Braak_Strat[[i]])[[j]])%>%dplyr::select(c(`brain region:ch1`,pseudoSampleID))%>%dplyr::filter(`brain region:ch1`%in%lobe_bm_area.map$BrainRegion[lobe_bm_area.map$Lobe%in%names(msbb_gse84422.Braak_Strat)[i]])%>%pull(pseudoSampleID))
}
msbb_gse84422.Braak_Strat=Filter(f = function(x)length(x$Braak0)>1,x = msbb_gse84422.Braak_Strat)
msbb_gse84422.Braak_DEGs=vector(mode = "list",length = length(names(msbb_gse84422.Braak_Strat)))
names(msbb_gse84422.Braak_DEGs)=names(msbb_gse84422.Braak_Strat)
msbb_gse84422.Braak_DEGs$Frontal_Lobe=msbb_gse84422.Braak_DEGs$Temporal_Lobe=vector(mode = "list",length = length(sort(unique(msbb_gse84422.pData$GPL96$SampleTypeBraak))[-1]))
names(msbb_gse84422.Braak_DEGs$Frontal_Lobe)=names(msbb_gse84422.Braak_DEGs$Temporal_Lobe)=sort(unique(msbb_gse84422.pData$GPL96$`braak neurofibrillary tangle score:ch1`))[-1]

for(l in 1:length(names(msbb_gse84422.Braak_Strat))){
  for(c in 2:length(msbb_gse84422.Braak_Strat[[l]])){
    c_exprs=msbb_gse84422.byLobe[[names(msbb_gse84422.Braak_Strat)[l]]][,msbb_gse84422.Braak_Strat[[l]]$Braak0]
    d_exprs=msbb_gse84422.byLobe[[names(msbb_gse84422.Braak_Strat)[l]]][,msbb_gse84422.Braak_Strat[[l]][[c]]]
    group=factor(c(rep("CONTROL",dim(c_exprs)[2]),rep("AD",dim(d_exprs)[2])))
    design_df=model.matrix(~0+group)
    colnames(design_df)=levels(group)
    contrasts_matrix=makeContrasts(CONTROL-AD,levels = design_df)
    fit=lmFit(object = cbind.data.frame(c_exprs,d_exprs),design = design_df)
    fit2=contrasts.fit(fit,contrasts_matrix)
    fit2=eBayes(fit2,trend = T)
    msbb_gse84422.Braak_DEGs[[l]][[c-1]]=topTable(fit = fit2,coef = 1,number = dim(msbb_gse84422.byLobe[[l]])[1],adjust.method = "BH",p.value = 0.1)%>%rownames_to_column("Gene")
    fwrite(msbb_gse84422.Braak_DEGs[[l]][[c-1]],paste(names(msbb_gse84422.Braak_DEGs)[[l]],names(msbb_gse84422.Braak_DEGs[[l]])[[c-1]],"Braak_DEGs.txt",sep = "_"),sep = "\t",col.names = T,row.names = F,quote = F)
  }
}

msbb_gse84422_Braak_DEG.Uniq=vector(mode = "list",length = length(names(msbb_gse84422.Braak_DEGs)))
names(msbb_gse84422_Braak_DEG.Uniq)=names(msbb_gse84422.Braak_DEGs)
msbb_gse84422_Braak_DEG.Uniq$Frontal_Lobe=msbb_gse84422_Braak_DEG.Uniq$Temporal_Lobe=vector(mode = "list",length = 6)
names(msbb_gse84422_Braak_DEG.Uniq$Frontal_Lobe)=names(msbb_gse84422_Braak_DEG.Uniq$Temporal_Lobe)=names(msbb_gse84422.Braak_DEGs$Frontal_Lobe)

for(t in 1:length(msbb_gse84422.Braak_DEGs)){
  for(b in 1:length(msbb_gse84422.Braak_DEGs[[t]])){
    msbb_gse84422_Braak_DEG.Uniq[[t]][[b]]=setdiff(msbb_gse84422.Braak_DEGs[[t]][[b]]$Gene,c(Reduce(union,lapply(lapply(msbb_gse84422.Braak_DEGs[[t]],function(x)x%>%pull(Gene)),function(x)intersect(x,msbb_gse84422.Braak_DEGs[[t]][[b]]$Gene))[-b]),Reduce(intersect,lapply(msbb_gse84422.Braak_DEGs[[t]],function(x)x%>%pull(Gene)))))  
  }
  
}
  

msbb_gse84422_Braak_DEG.KEGG=vector(mode = "list",length = 2)
names(msbb_gse84422_Braak_DEG.KEGG)=names(msbb_gse84422.Braak_DEGs)
msbb_gse84422_Braak_DEG.KEGG$Frontal_Lobe=lapply(msbb_gse84422.Braak_DEGs$Frontal_Lobe,function(x)enrichr(genes = x$Gene,databases = kegg_dbs)[[1]]%>%dplyr::filter(Adjusted.P.value<=0.1))
msbb_gse84422_Braak_DEG.KEGG$Temporal_Lobe=lapply(msbb_gse84422.Braak_DEGs$Temporal_Lobe,function(x)enrichr(genes = x$Gene,databases = kegg_dbs)[[1]]%>%dplyr::filter(Adjusted.P.value<=0.1))

msbb_gse84422.Braak_SCE=vector(mode = "list",length = length(msbb_gse84422.Braak_Strat))
names(msbb_gse84422.Braak_SCE)=names(msbb_gse84422.Braak_Strat)
msbb_gse84422.Braak_SCE$Frontal_Lobe=msbb_gse84422.Braak_SCE$Temporal_Lobe=vector(mode = "list",length = 7)
names(msbb_gse84422.Braak_SCE$Frontal_Lobe)=names(msbb_gse84422.Braak_SCE$Temporal_Lobe)=paste("Braak",0:6,sep = "")
for(c in 1:length(msbb_gse84422.Braak_Strat$Frontal_Lobe)){
  msbb_gse84422.Braak_SCE$Frontal_Lobe[[c]]=single.chip.enrichment2(exprs = msbb_gse84422.byLobe_Entrez$Frontal_Lobe[,colnames(msbb_gse84422.byLobe_Entrez$Frontal_Lobe)%in%msbb_gse84422.Braak_Strat$Frontal_Lobe[[c]]],geneset = pathprint.Hs.gs,statistic = "mean",normalizedScore = F,progressBar = T)
  msbb_gse84422.Braak_SCE$Temporal_Lobe[[c]]=single.chip.enrichment2(exprs = msbb_gse84422.byLobe_Entrez$Temporal_Lobe[,colnames(msbb_gse84422.byLobe_Entrez$Temporal_Lobe)%in%msbb_gse84422.Braak_Strat$Temporal_Lobe[[c]]],geneset = pathprint.Hs.gs,statistic = "mean",normalizedScore = F,progressBar = T)
}

#Braak based DEPs
msbb_gse84422.Braak_DEPs=vector(mode = "list",length = length(names(msbb_gse84422.Braak_SCE)))
names(msbb_gse84422.Braak_DEPs)=names(msbb_gse84422.Braak_SCE)
msbb_gse84422.Braak_DEPs$Frontal_Lobe=msbb_gse84422.Braak_DEPs$Temporal_Lobe=vector(mode = "list",length = length(sort(unique(msbb_gse84422.pData$GPL96$SampleTypeBraak))[-1]))
names(msbb_gse84422.Braak_DEPs$Frontal_Lobe)=names(msbb_gse84422.Braak_DEPs$Temporal_Lobe)=names(msbb_gse84422.Braak_SCE$Frontal_Lobe)[-1]

for(l in 1:length(names(msbb_gse84422.Braak_SCE))){
  for(c in 2:length(names(msbb_gse84422.Braak_SCE[[l]]))){
    c_exprs=data.frame(msbb_gse84422.Braak_SCE[[l]]$Braak0)
    d_exprs=data.frame(msbb_gse84422.Braak_SCE[[l]][[c]])
    group=factor(c(rep("CONTROL",dim(c_exprs)[2]),rep("AD",dim(d_exprs)[2])))
    design_df=model.matrix(~0+group)
    colnames(design_df)=levels(group)
    contrasts_matrix=makeContrasts(CONTROL-AD,levels = design_df)
    fit=lmFit(object = cbind.data.frame(c_exprs,d_exprs),design = design_df)
    fit2=contrasts.fit(fit,contrasts_matrix)
    fit2=eBayes(fit2,trend = T)
    msbb_gse84422.Braak_DEPs[[l]][[c-1]]=topTable(fit = fit2,coef = 1,number = dim(msbb_gse84422.Braak_SCE[[l]][[c]])[1],adjust.method = "BH",p.value = 0.2)%>%rownames_to_column("Geneset_Pathway")
    msbb_gse84422.Braak_DEPs[[l]][[c-1]]$Pathway_MemberGenes=unlist(lapply(msbb_gse84422.Braak_DEPs[[l]][[c-1]]$Geneset_Pathway,function(x)paste(unlist(unname(pathprint_membership.genes[names(pathprint_membership.genes)%in%x])),collapse = ",")))
    fwrite(msbb_gse84422.Braak_DEPs[[l]][[c-1]],paste(names(msbb_gse84422.Braak_DEPs)[[l]],names(msbb_gse84422.Braak_DEPs[[l]])[[c-1]],"Braak_DEPs.txt",sep = "_"),sep = "\t",col.names = T,row.names = F,quote = F)
  }
}


#PLQ based DEGs
msbb_gse84422.PLQ_Strat=vector(mode = "list",length = length(unique(lobe_bm_area.map$Lobe)))
names(msbb_gse84422.PLQ_Strat)=sort(unique(lobe_bm_area.map$Lobe))
msbb_gse84422.PLQ_Strat$Dorsal_striatum=msbb_gse84422.PLQ_Strat$Frontal_Lobe=msbb_gse84422.PLQ_Strat$Occipetal_Lobe=msbb_gse84422.PLQ_Strat$Parietal_Lobe=msbb_gse84422.PLQ_Strat$Temporal_Lobe=vector(mode = "list",length = length(unique(msbb_gse84422.pData$GPL96$SampleTypePLQ)))
names(msbb_gse84422.PLQ_Strat$Dorsal_striatum)=names(msbb_gse84422.PLQ_Strat$Frontal_Lobe)=names(msbb_gse84422.PLQ_Strat$Occipetal_Lobe)=names(msbb_gse84422.PLQ_Strat$Parietal_Lobe)=names(msbb_gse84422.PLQ_Strat$Temporal_Lobe)=sort(unique(msbb_gse84422.pData$GPL96$SampleTypePLQ))
for(i in 1:length(msbb_gse84422.PLQ_Strat)){
  for(j in 1:length(msbb_gse84422.PLQ_Strat[[i]]))
    msbb_gse84422.PLQ_Strat[[i]][[j]]=(msbb_gse84422.pData$GPL96%>%filter(SampleTypePLQ==names(msbb_gse84422.PLQ_Strat[[i]])[[j]])%>%dplyr::select(c(`brain region:ch1`,pseudoSampleID))%>%dplyr::filter(`brain region:ch1`%in%lobe_bm_area.map$BrainRegion[lobe_bm_area.map$Lobe%in%names(msbb_gse84422.PLQ_Strat)[i]])%>%pull(pseudoSampleID))
}
msbb_gse84422.PLQ_Strat=lapply(msbb_gse84422.PLQ_Strat,function(x){x=x[c('PLQ_0','PLQ_Low','PLQ_Med','PLQ_High')];x})

msbb_gse84422.PLQ_DEGs=vector(mode = "list",length = length(names(msbb_gse84422.byLobe)))
names(msbb_gse84422.PLQ_DEGs)=names(msbb_gse84422.byLobe)
msbb_gse84422.PLQ_DEGs$Dorsal_striatum=msbb_gse84422.PLQ_DEGs$Frontal_Lobe=msbb_gse84422.PLQ_DEGs$Occipetal_Lobe=msbb_gse84422.PLQ_DEGs$Parietal_Lobe=msbb_gse84422.PLQ_DEGs$Temporal_Lobe=vector(mode = "list",length = length(sort(unique(msbb_gse84422.pData$GPL96$SampleTypePLQ))[-1]))
names(msbb_gse84422.PLQ_DEGs$Dorsal_striatum)=names(msbb_gse84422.PLQ_DEGs$Frontal_Lobe)=names(msbb_gse84422.PLQ_DEGs$Occipetal_Lobe)=names(msbb_gse84422.PLQ_DEGs$Parietal_Lobe)=names(msbb_gse84422.PLQ_DEGs$Temporal_Lobe)=sort(unique(msbb_gse84422.pData$GPL96$SampleTypePLQ))[-1]

msbb_gse84422.PLQ_DEGs=lapply(msbb_gse84422.PLQ_DEGs,function(x){x=x[c('PLQ_Low','PLQ_Med','PLQ_High')];x})
for(l in 1:length(names(msbb_gse84422.byLobe))){
  for(c in 2:length(msbb_gse84422.PLQ_Strat[[l]])){
    c_exprs=msbb_gse84422.byLobe[[names(msbb_gse84422.PLQ_Strat)[l]]][,msbb_gse84422.PLQ_Strat[[l]]$PLQ_0]
    d_exprs=msbb_gse84422.byLobe[[names(msbb_gse84422.PLQ_Strat)[l]]][,msbb_gse84422.PLQ_Strat[[l]][[c]]]
    group=factor(c(rep("CONTROL",dim(c_exprs)[2]),rep("AD",dim(d_exprs)[2])))
    design_df=model.matrix(~0+group)
    colnames(design_df)=levels(group)
    contrasts_matrix=makeContrasts(CONTROL-AD,levels = design_df)
    fit=lmFit(object = cbind.data.frame(c_exprs,d_exprs),design = design_df)
    fit2=contrasts.fit(fit,contrasts_matrix)
    fit2=eBayes(fit2,trend = T)
    msbb_gse84422.PLQ_DEGs[[l]][[c-1]]=topTable(fit = fit2,coef = 1,number = dim(msbb_gse84422.byLobe[[l]])[1],adjust.method = "BH",p.value = 0.1)%>%rownames_to_column("Gene")
    fwrite(msbb_gse84422.PLQ_DEGs[[l]][[c-1]],paste(names(msbb_gse84422.PLQ_DEGs)[[l]],names(msbb_gse84422.PLQ_DEGs[[l]])[[c-1]],"PLQ_DEGs.txt",sep = "_"),sep = "\t",col.names = T,row.names = F,quote = F)
  }
}
msbb_gse84422.PLQ_DEGs$Dorsal_striatum=Filter(f = function(x)dim(x)[1]>0,x = msbb_gse84422.PLQ_DEGs$Dorsal_striatum)
msbb_gse84422.PLQ_DEGs$Frontal_Lobe=Filter(f = function(x)dim(x)[1]>0,x = msbb_gse84422.PLQ_DEGs$Frontal_Lobe)
msbb_gse84422.PLQ_DEGs$Occipetal_Lobe=Filter(f = function(x)dim(x)[1]>0,x = msbb_gse84422.PLQ_DEGs$Occipetal_Lobe)
msbb_gse84422.PLQ_DEGs$Parietal_Lobe=Filter(f = function(x)dim(x)[1]>0,x = msbb_gse84422.PLQ_DEGs$Parietal_Lobe)
msbb_gse84422.PLQ_DEGs$Temporal_Lobe=Filter(f = function(x)dim(x)[1]>0,x = msbb_gse84422.PLQ_DEGs$Temporal_Lobe)

#Non overlapping PLQ DEGs

msbb_gse84422_PLQ_DEG.Uniq=vector(mode = "list",length = length(names(msbb_gse84422.PLQ_DEGs)))
names(msbb_gse84422_PLQ_DEG.Uniq)=names(msbb_gse84422.PLQ_DEGs)
msbb_gse84422_PLQ_DEG.Uniq$Dorsal_striatum=msbb_gse84422_PLQ_DEG.Uniq$Frontal_Lobe=msbb_gse84422_PLQ_DEG.Uniq$Occipetal_Lobe=msbb_gse84422_PLQ_DEG.Uniq$Parietal_Lobe=msbb_gse84422_PLQ_DEG.Uniq$Temporal_Lobe=vector(mode = "list",length = length(msbb_gse84422.PLQ_DEGs$Dorsal_striatum))
names(msbb_gse84422_PLQ_DEG.Uniq$Dorsal_striatum)=names(msbb_gse84422_PLQ_DEG.Uniq$Frontal_Lobe)=names(msbb_gse84422_PLQ_DEG.Uniq$Occipetal_Lobe)=names(msbb_gse84422_PLQ_DEG.Uniq$Parietal_Lobe)=names(msbb_gse84422_PLQ_DEG.Uniq$Temporal_Lobe)=names(msbb_gse84422.PLQ_DEGs$Dorsal_striatum)

for(t in 1:length(msbb_gse84422.PLQ_DEGs)){
  for(b in 1:length(msbb_gse84422.PLQ_DEGs[[t]])){
    msbb_gse84422_PLQ_DEG.Uniq[[t]][[b]]=msbb_gse84422.PLQ_DEGs[[t]][[b]][msbb_gse84422.PLQ_DEGs[[t]][[b]]$Gene%in%setdiff(msbb_gse84422.PLQ_DEGs[[t]][[b]]$Gene,c(Reduce(union,lapply(lapply(msbb_gse84422.PLQ_DEGs[[t]],function(x)x%>%pull(Gene)),function(x)intersect(x,msbb_gse84422.PLQ_DEGs[[t]][[b]]$Gene))[-b]),Reduce(intersect,lapply(msbb_gse84422.PLQ_DEGs[[t]],function(x)x%>%pull(Gene))))),]  
    fwrite(msbb_gse84422_PLQ_DEG.Uniq[[t]][[b]],paste(names(msbb_gse84422.PLQ_DEGs)[[t]],names(msbb_gse84422.PLQ_DEGs[[t]])[[b]],"NOSS_PLQ_DEGs.txt",sep = "_"),sep = "\t",col.names = T,row.names = F,quote = F)
  }
  
}


msbb_gse84422_PLQ_DEGs.KEGG=vector(mode = "list",length = 5)
names(msbb_gse84422_PLQ_DEGs.KEGG)=names(msbb_gse84422.PLQ_DEGs)
msbb_gse84422_PLQ_DEGs.KEGG$Dorsal_striatum=vector(mode = "list",length = length(msbb_gse84422.PLQ_DEGs$Dorsal_striatum))
msbb_gse84422_PLQ_DEGs.KEGG$Frontal_Lobe=vector(mode = "list",length = length(msbb_gse84422.PLQ_DEGs$Frontal_Lobe))
msbb_gse84422_PLQ_DEGs.KEGG$Occipetal_Lobe=vector(mode = "list",length = length(msbb_gse84422.PLQ_DEGs$Occipetal_Lobe))
msbb_gse84422_PLQ_DEGs.KEGG$Parietal_Lobe=vector(mode = "list",length = length(msbb_gse84422.PLQ_DEGs$Parietal_Lobe))
msbb_gse84422_PLQ_DEGs.KEGG$Temporal_Lobe=vector(mode = "list",length = length(msbb_gse84422.PLQ_DEGs$Temporal_Lobe))

names(msbb_gse84422_PLQ_DEGs.KEGG$Dorsal_striatum)=names(msbb_gse84422.PLQ_DEGs$Dorsal_striatum)
names(msbb_gse84422_PLQ_DEGs.KEGG$Frontal_Lobe)=names(msbb_gse84422.PLQ_DEGs$Frontal_Lobe)
names(msbb_gse84422_PLQ_DEGs.KEGG$Occipetal_Lobe)=names(msbb_gse84422.PLQ_DEGs$Occipetal_Lobe)
names(msbb_gse84422_PLQ_DEGs.KEGG$Parietal_Lobe)=names(msbb_gse84422.PLQ_DEGs$Parietal_Lobe)
names(msbb_gse84422_PLQ_DEGs.KEGG$Temporal_Lobe)=names(msbb_gse84422.PLQ_DEGs$Temporal_Lobe)

msbb_gse84422_PLQ_DEGs.KEGG$Dorsal_striatum=lapply(msbb_gse84422.PLQ_DEGs$Dorsal_striatum,function(x)enrichr(genes = x$Gene,databases = kegg_dbs)[[1]]%>%dplyr::filter(Adjusted.P.value<0.1))
msbb_gse84422_PLQ_DEGs.KEGG$Frontal_Lobe=lapply(msbb_gse84422.PLQ_DEGs$Frontal_Lobe,function(x)enrichr(genes = x$Gene,databases = kegg_dbs)[[1]]%>%dplyr::filter(Adjusted.P.value<0.1))
msbb_gse84422_PLQ_DEGs.KEGG$Occipetal_Lobe=lapply(msbb_gse84422.PLQ_DEGs$Occipetal_Lobe,function(x)enrichr(genes = x$Gene,databases = kegg_dbs)[[1]]%>%dplyr::filter(Adjusted.P.value<0.1))
msbb_gse84422_PLQ_DEGs.KEGG$Parietal_Lobe=lapply(msbb_gse84422.PLQ_DEGs$Parietal_Lobe,function(x)enrichr(genes = x$Gene,databases = kegg_dbs)[[1]]%>%dplyr::filter(Adjusted.P.value<0.1))
msbb_gse84422_PLQ_DEGs.KEGG$Temporal_Lobe=lapply(msbb_gse84422.PLQ_DEGs$Temporal_Lobe,function(x)enrichr(genes = x$Gene,databases = kegg_dbs)[[1]]%>%dplyr::filter(Adjusted.P.value<0.1))

msbb_gse84422_PLQ_DEG_Uniq.KEGG=vector(mode = "list",length = 5)
names(msbb_gse84422_PLQ_DEG_Uniq.KEGG)=names(msbb_gse84422_PLQ_DEG.Uniq)
msbb_gse84422_PLQ_DEG_Uniq.KEGG$Dorsal_striatum=msbb_gse84422_PLQ_DEG_Uniq.KEGG$Frontal_Lobe=msbb_gse84422_PLQ_DEG_Uniq.KEGG$Occipetal_Lobe=msbb_gse84422_PLQ_DEG_Uniq.KEGG$Parietal_Lobe=msbb_gse84422_PLQ_DEG_Uniq.KEGG$Temporal_Lobe=vector(mode = "list",length = 3)
names(msbb_gse84422_PLQ_DEG_Uniq.KEGG$Dorsal_striatum)=names(msbb_gse84422_PLQ_DEG.Uniq$Dorsal_striatum)
names(msbb_gse84422_PLQ_DEG_Uniq.KEGG$Frontal_Lobe)=names(msbb_gse84422_PLQ_DEG.Uniq$Frontal_Lobe)
names(msbb_gse84422_PLQ_DEG_Uniq.KEGG$Occipetal_Lobe)=names(msbb_gse84422_PLQ_DEG.Uniq$Occipetal_Lobe)
names(msbb_gse84422_PLQ_DEG_Uniq.KEGG$Parietal_Lobe)=names(msbb_gse84422_PLQ_DEG.Uniq$Parietal_Lobe)
names(msbb_gse84422_PLQ_DEG_Uniq.KEGG$Temporal_Lobe)=names(msbb_gse84422_PLQ_DEG.Uniq$Temporal_Lobe)

msbb_gse84422_PLQ_DEG_Uniq.KEGG$Dorsal_striatum=lapply(msbb_gse84422_PLQ_DEG.Uniq$Dorsal_striatum,function(x)enrichr(genes = x$Gene,databases = kegg_dbs)[[1]]%>%dplyr::filter(Adjusted.P.value<0.1))
msbb_gse84422_PLQ_DEG_Uniq.KEGG$Frontal_Lobe=lapply(msbb_gse84422_PLQ_DEG.Uniq$Frontal_Lobe,function(x)enrichr(genes = x$Gene,databases = kegg_dbs)[[1]]%>%dplyr::filter(Adjusted.P.value<0.1))
msbb_gse84422_PLQ_DEG_Uniq.KEGG$Occipetal_Lobe=lapply(msbb_gse84422_PLQ_DEG.Uniq$Occipetal_Lobe,function(x)enrichr(genes = x$Gene,databases = kegg_dbs)[[1]]%>%dplyr::filter(Adjusted.P.value<0.1))
msbb_gse84422_PLQ_DEG_Uniq.KEGG$Parietal_Lobe=lapply(msbb_gse84422_PLQ_DEG.Uniq$Parietal_Lobe,function(x)enrichr(genes = x$Gene,databases = kegg_dbs)[[1]]%>%dplyr::filter(Adjusted.P.value<0.1))
msbb_gse84422_PLQ_DEG_Uniq.KEGG$Temporal_Lobe=lapply(msbb_gse84422_PLQ_DEG.Uniq$Temporal_Lobe,function(x)enrichr(genes = x$Gene,databases = kegg_dbs)[[1]]%>%dplyr::filter(Adjusted.P.value<0.1))


#DEPs for PLQ strata
msbb_gse84422.PLQ_SCE=vector(mode = "list",length = length(names(msbb_gse84422.byLobe_Entrez)))
names(msbb_gse84422.PLQ_SCE)=names(msbb_gse84422.byLobe_Entrez)
msbb_gse84422.PLQ_SCE$Dorsal_striatum=msbb_gse84422.PLQ_SCE$Frontal_Lobe=msbb_gse84422.PLQ_SCE$Occipetal_Lobe=msbb_gse84422.PLQ_SCE$Parietal_Lobe=msbb_gse84422.PLQ_SCE$Temporal_Lobe=vector(mode = "list",length = length(sort(unique(msbb_gse84422.pData$GPL96$SampleTypePLQ))))
names(msbb_gse84422.PLQ_SCE$Dorsal_striatum)=names(msbb_gse84422.PLQ_SCE$Frontal_Lobe)=names(msbb_gse84422.PLQ_SCE$Occipetal_Lobe)=names(msbb_gse84422.PLQ_SCE$Parietal_Lobe)=names(msbb_gse84422.PLQ_SCE$Temporal_Lobe)=c('PLQ_0','PLQ_Low','PLQ_Med','PLQ_High')

for(l in 1:length(names(msbb_gse84422.PLQ_SCE))){
  for(c in 1:length(msbb_gse84422.PLQ_SCE[[l]])){
    msbb_gse84422.PLQ_SCE[[l]][[c]]=single.chip.enrichment2(exprs = msbb_gse84422.byLobe_Entrez[[l]][,colnames(msbb_gse84422.byLobe_Entrez[[l]])%in%msbb_gse84422.PLQ_Strat[[l]][[c]]],geneset = pathprint.Hs.gs,statistic = "mean",normalizedScore = F,progressBar = T)
  }
}

msbb_gse84422.PLQ_DEPs=vector(mode = "list",length = length(names(msbb_gse84422.byLobe)))
names(msbb_gse84422.PLQ_DEPs)=names(msbb_gse84422.byLobe)
msbb_gse84422.PLQ_DEPs$Dorsal_striatum=msbb_gse84422.PLQ_DEPs$Frontal_Lobe=msbb_gse84422.PLQ_DEPs$Occipetal_Lobe=msbb_gse84422.PLQ_DEPs$Parietal_Lobe=msbb_gse84422.PLQ_DEPs$Temporal_Lobe=vector(mode = "list",length = length(sort(unique(msbb_gse84422.pData$GPL96$SampleTypePLQ))[-1]))
names(msbb_gse84422.PLQ_DEPs$Dorsal_striatum)=names(msbb_gse84422.PLQ_DEPs$Frontal_Lobe)=names(msbb_gse84422.PLQ_DEPs$Occipetal_Lobe)=names(msbb_gse84422.PLQ_DEPs$Parietal_Lobe)=names(msbb_gse84422.PLQ_DEPs$Temporal_Lobe)=names(msbb_gse84422.PLQ_SCE$Dorsal_striatum)[-1]
for(l in 1:length(names(msbb_gse84422.PLQ_SCE))){
  for(c in 2:length(names(msbb_gse84422.PLQ_SCE[[l]]))){
    c_exprs=data.frame(msbb_gse84422.PLQ_SCE[[l]]$PLQ_0)
    d_exprs=data.frame(msbb_gse84422.PLQ_SCE[[l]][[c]])
    group=factor(c(rep("CONTROL",dim(c_exprs)[2]),rep("AD",dim(d_exprs)[2])))
    design_df=model.matrix(~0+group)
    colnames(design_df)=levels(group)
    contrasts_matrix=makeContrasts(CONTROL-AD,levels = design_df)
    fit=lmFit(object = cbind.data.frame(c_exprs,d_exprs),design = design_df)
    fit2=contrasts.fit(fit,contrasts_matrix)
    fit2=eBayes(fit2,trend = T)
    msbb_gse84422.PLQ_DEPs[[l]][[c-1]]=topTable(fit = fit2,coef = 1,number = dim(msbb_gse84422.PLQ_SCE[[l]][[c]])[1],adjust.method = "BH",p.value = 0.2)%>%rownames_to_column("Geneset_Pathway")
    msbb_gse84422.PLQ_DEPs[[l]][[c-1]]$Pathway_MemberGenes=unlist(lapply(msbb_gse84422.PLQ_DEPs[[l]][[c-1]]$Geneset_Pathway,function(x)paste(unlist(unname(pathprint_membership.genes[names(pathprint_membership.genes)%in%x])),collapse = ",")))
    fwrite(msbb_gse84422.PLQ_DEPs[[l]][[c-1]],paste(names(msbb_gse84422.PLQ_DEPs)[[l]],names(msbb_gse84422.PLQ_DEPs[[l]])[[c-1]],"PLQ_DEPs.txt",sep = "_"),sep = "\t",col.names = T,row.names = F,quote = F)
  }
}

#Include pathway membership genes, map to gene symbols and remove empty elements in list
msbb_gse84422.PLQ_DEPs$Dorsal_striatum=Filter(f = function(a)dim(a)>1|is.null(a)!=T,x = msbb_gse84422.PLQ_DEPs$Dorsal_striatum)
msbb_gse84422.PLQ_DEPs$Frontal_Lobe=Filter(f = function(a)dim(a)>1|is.null(a)!=T,x = msbb_gse84422.PLQ_DEPs$Frontal_Lobe)
msbb_gse84422.PLQ_DEPs$Occipetal_Lobe=Filter(f = function(a)dim(a)>1|is.null(a)!=T,x = msbb_gse84422.PLQ_DEPs$Occipetal_Lobe)
msbb_gse84422.PLQ_DEPs$Parietal_Lobe=Filter(f = function(a)dim(a)>1|is.null(a)!=T,x = msbb_gse84422.PLQ_DEPs$Parietal_Lobe)
msbb_gse84422.PLQ_DEPs$Temporal_Lobe=Filter(f = function(a)dim(a)>1|is.null(a)!=T,x = msbb_gse84422.PLQ_DEPs$Temporal_Lobe)



msbb_gse84422_plq_corrGenes.byLobe=vector(mode = "list",length = 5)
names(msbb_gse84422_plq_corrGenes.byLobe)=names(msbb_gse84422.byLobe)
msbb_gse84422_plq_values=msbb_gse84422.pData$GPL96$`average neuritic plaque density:ch1`
names(msbb_gse84422_plq_values)=msbb_gse84422.pData$GPL96$pseudoSampleID
#Computation should be done on HPC
for(i in 1:length(msbb_gse84422.byLobe)){
  
  msbb_gse84422_plq_corrGenes.byLobe[[i]]=foreach(g=1:dim(msbb_gse84422.byLobe[[i]])[1]) %dopar% {data.frame(Gene=rownames(msbb_gse84422.byLobe[[i]][g,]),
                                                                                                       Spearman.Rho=cor.test(x=unlist(unname(msbb_gse84422.byLobe[[i]][g,])),y=unname(msbb_gse84422_plq_values[colnames(msbb_gse84422.byLobe[[i]])]),method = "spearman")[['estimate']],
                                                                                                       pval=cor.test(x=unlist(unname(msbb_gse84422.byLobe[[i]][g,])),y=unname(msbb_gse84422_plq_values[colnames(msbb_gse84422.byLobe[[i]])]),method = "spearman")[['p.value']],
                                                                                                       stringsAsFactors=F)}
}

#Read results from HPC computations
msbb_gse84422_plq_corrGenes.byLobe=readRDS("msbb_gse84422_PLQ_CorrGenes_byLobe.RDS")
msbb_gse84422_plq_corrGenes_byLobe.df=lapply(msbb_gse84422_plq_corrGenes.byLobe, function(x)data.frame(rbindlist(x),stringsAsFactors = F))
msbb_gse84422_plq_corrGenes_byLobe.df=lapply(msbb_gse84422_plq_corrGenes_byLobe.df,function(x)x%>%mutate(adj.p=p.adjust(p = pval,method = "fdr")))
msbb_gse84422_plq_corrGenes_byLobe.df=lapply(msbb_gse84422_plq_corrGenes_byLobe.df,function(x)x%>%dplyr::filter(adj.p<=0.2))
msbb_gse84422_plq_corrGenes_byLobe.df=Filter(f = function(x)dim(x)[1]>0,x = msbb_gse84422_plq_corrGenes_byLobe.df)
for(t in 1:length(msbb_gse84422_plq_corrGenes_byLobe.df)){
  fwrite(msbb_gse84422_plq_corrGenes_byLobe.df[[t]],paste(names(msbb_gse84422_plq_corrGenes_byLobe.df)[t],"PAGs.txt",sep = "_"),sep = "\t",col.names = T)
}
msbb_gse84422_plq_corrGenes_byLobe.df
msbb_gse84422_plq_corrGenes_byLobe.KEGG=vector(mode = "list",length = 4)
names(msbb_gse84422_plq_corrGenes_byLobe.KEGG)=names(msbb_gse84422_plq_corrGenes_byLobe.df)
msbb_gse84422_plq_corrGenes_byLobe.KEGG=lapply(msbb_gse84422_plq_corrGenes_byLobe.df,function(x)enrichr(genes = x$Gene,databases = kegg_dbs)[[1]]%>%dplyr::filter(Adjusted.P.value<0.1))



gene_expr=msbb_gse84422.byLobe$Dorsal_striatum[1,1:10]
plq_values=msbb_gse84422.pData$GPL96$`average neuritic plaque density:ch1`[1:10]

#Comparisons between DEGs for each covariate

msbb_gse84422_CDR_DEG_compMatrix.list=vector(mode = "list",length = 5)
msbb_gse84422_PLQ_DEG_compMatrix.list=vector(mode = "list",length = 5)
msbb_gse84422_CDR_DEG_compMatrix_InterLobular.list=vector(mode = "list",length = 5)
names(msbb_gse84422_CDR_DEG_compMatrix.list)=names(msbb_gse84422.CDR_DEGs)
names(msbb_gse84422_PLQ_DEG_compMatrix.list)=names(msbb_gse84422.PLQ_DEGs)

for(i in 1:length(msbb_gse84422.CDR_DEGs)){
  msbb_gse84422_CDR_DEG_compMatrix.list[[i]]=matrix(NA,nrow=length(msbb_gse84422.CDR_DEGs$Dorsal_striatum),ncol = length(msbb_gse84422.CDR_DEGs$Dorsal_striatum))
  rownames(msbb_gse84422_CDR_DEG_compMatrix.list[[i]])=colnames(msbb_gse84422_CDR_DEG_compMatrix.list[[i]])=names(msbb_gse84422.CDR_DEGs[[i]])
}
for(i in 1:length(msbb_gse84422.PLQ_DEGs)){
  msbb_gse84422_PLQ_DEG_compMatrix.list[[i]]=matrix(NA,nrow=length(msbb_gse84422.PLQ_DEGs$Dorsal_striatum),ncol = length(msbb_gse84422.PLQ_DEGs$Dorsal_striatum))
  rownames(msbb_gse84422_PLQ_DEG_compMatrix.list[[i]])=colnames(msbb_gse84422_PLQ_DEG_compMatrix.list[[i]])=names(msbb_gse84422.PLQ_DEGs[[i]])
}

for(l in 1:5){
  for(i in 1:6){
    msbb_gse84422_CDR_DEG_compMatrix.list[[l]][i,]=unlist(lapply(lapply(msbb_gse84422.CDR_DEGs[[l]],function(x)x%>%pull(Gene)),function(y)length(intersect(y,msbb_gse84422.CDR_DEGs[[l]][[i]]$Gene))/length(union(y,msbb_gse84422.CDR_DEGs[[l]][[i]]$Gene))))
  }
}

for(l in 1:5){
  for(i in 1:3){
    msbb_gse84422_PLQ_DEG_compMatrix.list[[l]][i,]=unlist(lapply(lapply(msbb_gse84422.PLQ_DEGs[[l]],function(x)x%>%pull(Gene)),function(y)length(intersect(y,msbb_gse84422.PLQ_DEGs[[l]][[i]]$Gene))/length(union(y,msbb_gse84422.PLQ_DEGs[[l]][[i]]$Gene))))
  }
}  

Reduce(intersect,list(B1=Reduce(intersect,lapply(msbb_gse84422.Braak_DEGs,function(x)x$`1`%>%pull(Gene))),
     B2=Reduce(intersect,lapply(msbb_gse84422.Braak_DEGs,function(x)x$`2`%>%pull(Gene))),
     B3=Reduce(intersect,lapply(msbb_gse84422.Braak_DEGs,function(x)x$`3`%>%pull(Gene))),
     B4=Reduce(intersect,lapply(msbb_gse84422.Braak_DEGs,function(x)x$`4`%>%pull(Gene))),
     B5=Reduce(intersect,lapply(msbb_gse84422.Braak_DEGs,function(x)x$`5`%>%pull(Gene))),
     B6=Reduce(intersect,lapply(msbb_gse84422.Braak_DEGs,function(x)x$`6`%>%pull(Gene)))))

Reduce(intersect,list(FL=Reduce(intersect,lapply(msbb_gse84422.Braak_DEGs$Frontal_Lobe,function(x)x%>%pull(Gene))),
TL=Reduce(intersect,lapply(msbb_gse84422.Braak_DEGs$Temporal_Lobe,function(x)x%>%pull(Gene)))))

msbb_gse84422_Interlobular_SuperDEGs.PLQ=Reduce(intersect,list(Low=Reduce(intersect,lapply(msbb_gse84422.PLQ_DEGs,function(x)x$`PLQ_Low`%>%pull(Gene))),
                      Med=Reduce(intersect,lapply(msbb_gse84422.PLQ_DEGs,function(x)x$`PLQ_Med`%>%pull(Gene))),
                      High=Reduce(intersect,lapply(msbb_gse84422.PLQ_DEGs,function(x)x$`PLQ_High`%>%pull(Gene)))))

msbb_gse84422_Intralobular_SuperDEGs.PLQ=Reduce(intersect,list(DorStr=Reduce(intersect,(lapply(msbb_gse84422.PLQ_DEGs$Dorsal_striatum,function(x)x%>%pull(Gene)))),
     FrtLobe=Reduce(intersect,(lapply(msbb_gse84422.PLQ_DEGs$Dorsal_striatum,function(x)x%>%pull(Gene)))),
     OccLobe=Reduce(intersect,(lapply(msbb_gse84422.PLQ_DEGs$Dorsal_striatum,function(x)x%>%pull(Gene)))),
     ParLobe=Reduce(intersect,(lapply(msbb_gse84422.PLQ_DEGs$Dorsal_striatum,function(x)x%>%pull(Gene)))),
     TempLobe=Reduce(intersect,(lapply(msbb_gse84422.PLQ_DEGs$Dorsal_striatum,function(x)x%>%pull(Gene))))))


# Plot DEGs
msbb_gse84422_PLQ_DEG_Summary=read.table("MSBB_GSE84422_PLQ_DEG_Summary_for_ggplot2.txt",sep = "\t",header = T,as.is = T)
msbb_gse84422_PLQ_DEG_Summary$PLQ_Load[msbb_gse84422_PLQ_DEG_Summary$PLQ_Load=="Low"]="1Low"
msbb_gse84422_PLQ_DEG_Summary$PLQ_Load[msbb_gse84422_PLQ_DEG_Summary$PLQ_Load=="Med"]="2Med"
msbb_gse84422_PLQ_DEG_Summary$PLQ_Load[msbb_gse84422_PLQ_DEG_Summary$PLQ_Load=="High"]="3High"
plot_plq_DEGs=ggplot(data = msbb_gse84422_PLQ_DEG_Summary,aes(x=PLQ_Load,y=DEGs,group=Brain_Lobe))+geom_line(aes(color=Brain_Lobe),size=1)+geom_point(aes(color=Brain_Lobe))+theme(axis.title.x = element_text(face = "bold",size = 12,colour = "black"),axis.text.x = element_text(size = 14,colour = "black"),axis.title.y = element_text(face = "bold",size = 12,colour = "black"),legend.text = element_text(face = "bold",size = 12),axis.text.y = element_text(size = 14,colour = "black"))+theme_light(base_size = 14)+ggtitle("Brain lobe-wise DEGs - PLQ",subtitle = "Generic DEGs")+theme(plot.title = element_text(size = 16))
plot_plq_NOSS_DEGs=ggplot(data = msbb_gse84422_PLQ_DEG_Summary,aes(x=PLQ_Load,y=NOSS_DEGs,group=Brain_Lobe))+geom_line(aes(color=Brain_Lobe),size=1)+geom_point(aes(color=Brain_Lobe))+theme(axis.title.x = element_text(face = "bold",size = 12,colour = "black"),axis.text.x = element_text(size = 14,colour = "black"),axis.title.y = element_text(face = "bold",size = 12,colour = "black"),legend.text = element_text(face = "bold",size = 12),axis.text.y = element_text(size = 14,colour = "black"))+theme_light(base_size = 14)+ggtitle("Brain lobe-wise DEGs - PLQ",subtitle = "Non-overlapping State specific DEGs")+theme(plot.title = element_text(size = 16))

msbb_gse84422_PLQ_DEP_Summary=read.table("../GSE84422/Lobe_Analysis/MSBB_GSE84422_PLQ_DEP_Summary_for_ggplot2.txt",sep = "\t",header = T,as.is = T)
msbb_gse84422_PLQ_DEP_Summary$PLQ_Load[msbb_gse84422_PLQ_DEP_Summary$PLQ_Load=="Low"]="1Low"
msbb_gse84422_PLQ_DEP_Summary$PLQ_Load[msbb_gse84422_PLQ_DEP_Summary$PLQ_Load=="Medium"]="2Med"
msbb_gse84422_PLQ_DEP_Summary$PLQ_Load[msbb_gse84422_PLQ_DEP_Summary$PLQ_Load=="High"]="3High"
plot_plq_DEPs=ggplot(data = msbb_gse84422_PLQ_DEP_Summary,aes(x=PLQ_Load,y=DEPs,group=BrainLobe))+geom_line(aes(color=BrainLobe),size=1)+geom_point(aes(color=BrainLobe))+theme(axis.title.x = element_text(face = "bold",size = 12,colour = "black"),axis.text.x = element_text(size = 14,colour = "black"),axis.title.y = element_text(face = "bold",size = 12,colour = "black"),legend.text = element_text(face = "bold",size = 16),axis.text.y = element_text(size = 16,colour = "black"))+theme_light(base_size = 14)+ggtitle("Brain lobe-wise DEPs - PLQ")+theme(plot.title = element_text(size = 16))



msbb_gse84422_CDR_DEG_Summary=read.table("MSBB_GSE84422_CDR_DEG_Summary_for_ggplot2.txt",sep = "\t",header = T,as.is = T)

msbb_gse84422_CDR_DEG_Summary$CDR=as.character(msbb_gse84422_CDR_DEG_Summary$CDR)
msbb_gse84422_CDR_DEG_Summary$DEGs=log10(msbb_gse84422_CDR_DEG_Summary$DEGs)
msbb_gse84422_CDR_DEG_Summary$NOSS_DEGs=log10(msbb_gse84422_CDR_DEG_Summary$NOSS_DEGs)
plot_cdr_DEGs=ggplot(data = msbb_gse84422_CDR_DEG_Summary,aes(x=CDR,y=DEGs,group=`Brain.Lobe`))+geom_line(aes(color=`Brain.Lobe`),size=1)+geom_point(aes(color=`Brain.Lobe`))+theme(axis.title.x = element_text(face = "bold",size = 12,colour = "black"),axis.text.x = element_text(size = 12,colour = "black"),axis.title.y = element_text(face = "bold",size = 12,colour = "black"),legend.text = element_text(face = "bold",size = 12),axis.text.y = element_text(size = 12,colour = "black"))+theme_light(base_size = 14)+ggtitle("Brain lobe-wise DEGs - CDR",subtitle = "Generic DEGs")+theme(plot.title = element_text(size = 16))
plot_cdr_NOSS_DEGs=ggplot(data = msbb_gse84422_CDR_DEG_Summary,aes(x=CDR,y=NOSS_DEGs,group=`Brain.Lobe`))+geom_line(aes(color=`Brain.Lobe`),size=1)+geom_point(aes(color=`Brain.Lobe`))+theme(axis.title.x = element_text(face = "bold",size = 12,colour = "black"),axis.text.x = element_text(size = 12,colour = "black"),axis.title.y = element_text(face = "bold",size = 12,colour = "black"),legend.text = element_text(face = "bold",size = 12),axis.text.y = element_text(size = 12,colour = "black"))+theme_light(base_size = 14)+ggtitle("Brain lobe-wise DEGs - CDR",subtitle = "Non-overlapping State specific DEGs")+theme(plot.title = element_text(size = 16))

msbb_gse84422_CDR_DEP_Summary=read.table("../GSE84422/Lobe_Analysis/MSBB_GSE84422_CDR_DEP_Summary_for_ggplot2.txt",sep = "\t",header = T,as.is = T)
msbb_gse84422_CDR_DEP_Summary=msbb_gse84422_CDR_DEP_Summary[,c(1:3)]
msbb_gse84422_CDR_DEP_Summary$CDR=as.character(msbb_gse84422_CDR_DEP_Summary$CDR)
plot_cdr_DEPs=ggplot(data = msbb_gse84422_CDR_DEP_Summary,aes(x=CDR,y=DEPs,group=BrainLobe))+geom_line(aes(color=BrainLobe),size=1)+geom_point(aes(color=BrainLobe))+theme(axis.title.x = element_text(face = "bold",size = 12,colour = "black"),axis.text.x = element_text(size = 12,colour = "black"),axis.title.y = element_text(face = "bold",size = 12,colour = "black"),legend.text = element_text(face = "bold",size = 12),axis.text.y = element_text(size = 12,colour = "black"))+theme_light(base_size = 14)+ggtitle("Brain lobe-wise DEPs - CDR")+theme(plot.title = element_text(size = 16))


msbb_gse84422_Braak_DEG_Summary=read.table("MSBB_GSE84422_Braak_DEG_Summary_for_ggplot2.txt",sep = "\t",header = T,as.is = T)


plot_braak_DEGs=ggplot(data = msbb_gse84422_Braak_DEG_Summary,aes(x=Braak,y=DEGs,group=`Brain_Lobe`))+geom_line(aes(color=`Brain_Lobe`),size=1)+geom_point(aes(color=`Brain_Lobe`))+theme(axis.title.x = element_text(face = "bold",size = 18,colour = "black"),axis.text.x = element_text(size = 18,colour = "black"),axis.title.y = element_text(face = "bold",size = 18,colour = "black"),legend.text = element_text(face = "bold",size = 18),axis.text.y = element_text(size = 18,colour = "black"))+theme_light(base_size = 14)+ggtitle("Brain lobe-wise DEGs - Braak",subtitle = "Generic DEGs")+theme(plot.title = element_text(size = 16))
plot_braak_NOSS_DEGs=ggplot(data = msbb_gse84422_Braak_DEG_Summary,aes(x=Braak,y=NOSS_DEGs,group=`Brain_Lobe`))+geom_line(aes(color=`Brain_Lobe`),size=1)+geom_point(aes(color=`Brain_Lobe`))+theme(axis.title.x = element_text(face = "bold",size = 18,colour = "black"),axis.text.x = element_text(size = 18,colour = "black"),axis.title.y = element_text(face = "bold",size = 18,colour = "black"),legend.text = element_text(face = "bold",size = 18),axis.text.y = element_text(size = 18,colour = "black"))+theme_light(base_size = 14)+ggtitle("Brain lobe-wise DEGs - Braak",subtitle = "Non-overlapping State specific DEGs")+theme(plot.title = element_text(size = 16))

# Validation in Miller et al.

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
gse29378_exprs.agg=aggregate(x = gse29378.exprs[,-which(colnames(gse29378.exprs)=="EntrezID")],by=list(entrez=gse29378.exprs$EntrezID),mean)
rownames(gse29378_exprs.agg)=gse29378_exprs.agg$entrez
gse29378_exprs.agg=gse29378_exprs.agg[,-1]
gse29378_exprs.agg=log2(gse29378_exprs.agg)


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
gse29378_PLQ.diffPathways=lapply(gse29378_PLQ.diffPathways[-1],function(x)x%>%filter(adj.P.Val<=0.2))

gse29378_PLQ_Strat.exprs=gse29378_PLQ.DEGs=vector(mode = "list",length = 4)
names(gse29378_PLQ_Strat.exprs)=names(gse29378_PLQ.DEGs)=names(gse29378_PLQ_Strat.SCE)
gse29378_PLQ_Strat.exprs$S0=gse29378_exprs.agg[,gse29378_PLQ_Strat$S0]
gse29378_PLQ_Strat.exprs$S1=gse29378_exprs.agg[,gse29378_PLQ_Strat$S1]
gse29378_PLQ_Strat.exprs$S2=gse29378_exprs.agg[,gse29378_PLQ_Strat$S2]
gse29378_PLQ_Strat.exprs$S3=gse29378_exprs.agg[,gse29378_PLQ_Strat$S3]

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
    diff_genes$Gene_HGNC=mapIds(x = org.Hs.eg.db,keys = diff_genes$Gene,keytype = "ENTREZID",column = "SYMBOL")
    
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
    diff_genes$Gene_HGNC=mapIds(x = org.Hs.eg.db,keys = diff_genes$Gene,keytype = "ENTREZID",column = "SYMBOL")
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
    diff_genes$Gene_HGNC=mapIds(x = org.Hs.eg.db,keys = diff_genes$Gene,keytype = "ENTREZID",column = "SYMBOL")
    if(length(diff_genes)==0){
      gse29378_PLQ.DEGs[[i]]=0
    }
    else 
      rownames(diff_genes)=rownames(gse29378_exprs.agg)
    gse29378_PLQ.DEGs[[i]]=diff_genes
  }
}


