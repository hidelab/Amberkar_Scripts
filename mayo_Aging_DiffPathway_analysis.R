library(synapseClient)
library(dplyr)
library(magrittr)
library(data.table)
library(org.Hs.eg.db)
library(pathprint)
library(limma)
synapseLogin()
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


ds = c("chipframe", "genesets","pathprint.Hs.gs" ,"platform.thresholds")
data(list = ds)
mayo_PA_norm_counts.synapse=synGet(id='syn7440549')
mayo_PA_covariates.synapse=synGet(id='syn7437038')
mayo_AD_norm_counts.synapse=synGet(id='syn4650265')
mayo_AD_covariates.synapse=synGet(id='syn8466814')
mayo_AD_norm_counts=fread(mayo_AD_norm_counts.synapse@filePath,sep = "\t",header = T,stringsAsFactors = F,data.table = F)
mayo_PA_norm_counts=fread(mayo_PA_norm_counts.synapse@filePath,sep = "\t",header = T,stringsAsFactors = F,data.table = F)
mayo_PA_covariates=fread(mayo_PA_covariates.synapse@filePath,sep = "\t",header = T,stringsAsFactors = F,data.table = F)
mayo_AD_covariates=fread(mayo_AD_covariates.synapse@filePath,sep = "\t",header = T,stringsAsFactors = F,data.table = F)
mayo_PA_samples=mayo_PA_covariates%>%filter(AgeAtDeath>=90)%>%pull(SubjectID)
mayo_allPA_samples=mayo_PA_covariates$SubjectID
mayo_Control_samples=mayo_AD_covariates%>%filter(Tissue.SourceDiagnosis=="TCX.CONTROL"&AgeAtDeath>=90)%>%pull(SampleID)
mayo_yControl_samples=mayo_AD_covariates%>%filter(Tissue.SourceDiagnosis=="TCX.CONTROL"&AgeAtDeath<=75)%>%pull(SampleID)
mayo_AD_samples=mayo_AD_covariates%>%filter(Tissue.SourceDiagnosis=="TCX.AD"&AgeAtDeath>=90)%>%pull(SampleID)
mayo_AD_samples2=mayo_AD_covariates%>%filter(Tissue.SourceDiagnosis=="TCX.AD")%>%pull(SampleID)
mayo_PA_AD_samples_To_analyse=c(mayo_PA_samples,mayo_AD_samples)
mayo_Control_AD_samples_To_analyse=c(mayo_Control_samples,mayo_AD_samples2)
mayo_yControl_PA_samples_To_analyse=c(mayo_yControl_samples,mayo_PA_samples)

mayo_AD_PA_norm_counts=cbind.data.frame(mayo_PA_norm_counts,mayo_AD_norm_counts)
mayo_yControl_PA_norm_counts=cbind.data.frame(mayo_AD_norm_counts[,colnames(mayo_AD_norm_counts)%in%mayo_yControl_samples],mayo_PA_norm_counts)
mayo_AD_PA_norm_counts=mayo_AD_PA_norm_counts[which((rowSums(mayo_AD_PA_norm_counts>0)>=ncol(mayo_AD_PA_norm_counts)/3)=="TRUE"),]
mayo_AD_norm_counts=mayo_AD_norm_counts[which((rowSums(mayo_AD_norm_counts>0)>=ncol(mayo_AD_norm_counts)/3)=="TRUE"),]
mayo_yControl_PA_norm_counts=mayo_yControl_PA_norm_counts[which((rowSums(mayo_yControl_PA_norm_counts>0)>=ncol(mayo_yControl_PA_norm_counts)/3)=="TRUE"),]

mayo_AD_PA_norm_counts$EntrezID=unname(mapIds(x = org.Hs.eg.db,keys = mayo_AD_PA_norm_counts$ensembl_id,column = "ENTREZID",keytype = "ENSEMBL",multiVals = "first"))
mayo_AD_norm_counts$EntrezID=unname(mapIds(x = org.Hs.eg.db,keys = mayo_AD_norm_counts$ensembl_id,column = "ENTREZID",keytype = "ENSEMBL",multiVals = "first"))
mayo_yControl_PA_norm_counts$EntrezID=unname(mapIds(x = org.Hs.eg.db,keys = mayo_yControl_PA_norm_counts$ensembl_id,column = "ENTREZID",keytype = "ENSEMBL",multiVals = "first"))

mayo_AD_PA_norm_counts.agg=aggregate.data.frame(x = mayo_AD_PA_norm_counts[,-c(1,43,320)],by = list(entrez=mayo_AD_PA_norm_counts$EntrezID),mean)
mayo_AD_norm_counts.agg=aggregate.data.frame(x = mayo_AD_norm_counts[,-c(1,278)],by = list(entrez=mayo_AD_norm_counts$EntrezID),mean)
mayo_yControl_PA_norm_counts.agg=aggregate.data.frame(x = mayo_yControl_PA_norm_counts[,-c(15,57)],by = list(entrez=mayo_yControl_PA_norm_counts$EntrezID),mean)
rownames(mayo_AD_PA_norm_counts.agg)=mayo_AD_PA_norm_counts.agg$entrez
rownames(mayo_AD_norm_counts.agg)=mayo_AD_norm_counts.agg$entrez
rownames(mayo_yControl_PA_norm_counts.agg)=mayo_yControl_PA_norm_counts.agg$entrez
mayo_AD_PA_norm_counts.agg=mayo_AD_PA_norm_counts.agg[,-which(colnames(mayo_AD_PA_norm_counts.agg)=="entrez")]
mayo_AD_norm_counts.agg=mayo_AD_norm_counts.agg[,-which(colnames(mayo_AD_norm_counts.agg)=="entrez")]
mayo_yControl_PA_norm_counts.agg=mayo_yControl_PA_norm_counts.agg[,-which(colnames(mayo_yControl_PA_norm_counts.agg)=="entrez")]

mayo_AD_PA_norm_counts.SCE=single.chip.enrichment2(exprs = mayo_AD_PA_norm_counts.agg[,colnames(mayo_AD_PA_norm_counts.agg)%in%mayo_PA_AD_samples_To_analyse],geneset = pathprint.Hs.gs,statistic = "mean",progressBar = T)
mayo_AD_norm_counts.SCE=single.chip.enrichment2(exprs = mayo_AD_norm_counts.agg[,colnames(mayo_AD_norm_counts.agg)%in%mayo_Control_AD_samples_To_analyse],geneset = pathprint.Hs.gs,statistic = "mean",progressBar = T)
mayo_yControl_PA_norm_counts.SCE=single.chip.enrichment2(exprs = mayo_yControl_PA_norm_counts.agg,geneset = pathprint.Hs.gs,statistic = "mean",progressBar = T)

mayo_cohort.Diffpathways=vector(mode = "list",length = 3)
names(mayo_cohort.Diffpathways)=c("PA-AD","Control-AD","yControl-PA")
for(i in 1:3){
  if(i==1){
    group=factor(c(rep("PA",length(mayo_PA_samples)),rep("AD",length(mayo_AD_samples))))
    design_df=model.matrix(~0+group)
    colnames(design_df)=levels(group)
    contrasts_matrix=makeContrasts(PA-AD,levels = design_df)
    fit=lmFit(object = mayo_AD_PA_norm_counts.SCE,design = design_df)
    fit2=contrasts.fit(fit,contrasts_matrix)
    fit2=eBayes(fit2,trend = T)
    mayo_cohort.Diffpathways[[i]]=topTable(fit2,coef = 1,number = 633,p.value = 1,adjust.method = "BH")%>%rownames_to_column("Pathway/GS")
  }
  if(i==2){
    group=factor(c(rep("Control",length(mayo_Control_samples)),rep("AD",length(mayo_AD_samples2))))
    design_df=model.matrix(~0+group)
    colnames(design_df)=levels(group)
    contrasts_matrix=makeContrasts(Control-AD,levels = design_df)
    fit=lmFit(object = mayo_AD_norm_counts.SCE,design = design_df)
    fit2=contrasts.fit(fit,contrasts_matrix)
    fit2=eBayes(fit2,trend = T)
    mayo_cohort.Diffpathways[[i]]=topTable(fit2,coef = 1,number = 633,p.value = 1,adjust.method = "BH")%>%rownames_to_column("Pathway/GS")
  }
  if(i==3){
    group=factor(c(rep("yControl",length(mayo_yControl_samples)),rep("PA",length(mayo_allPA_samples))))
    design_df=model.matrix(~0+group)
    colnames(design_df)=levels(group)
    contrasts_matrix=makeContrasts(yControl-PA,levels = design_df)
    fit=lmFit(object = mayo_yControl_PA_norm_counts.SCE,design = design_df)
    fit2=contrasts.fit(fit,contrasts_matrix)
    fit2=eBayes(fit2,trend = T)
    mayo_cohort.Diffpathways[[i]]=topTable(fit2,coef = 1,number = 633,p.value = 1,adjust.method = "BH")%>%rownames_to_column("Pathway/GS")
  }
}

#For boxplots
mayo_AD_PA_diffPathways.df=data.frame(mayo_AD_PA_norm_counts.SCE[which(rownames(mayo_AD_PA_norm_counts.SCE)%in%mayo_cohort.Diffpathways$`PA-AD`$`Pathway/GS`),])
mayo_AD_PA_diffPathways.df$SampleType="SAMPLE"
mayo_AD_PA_diffPathways.df$SampleType[1:21]="Path.Aging"
mayo_AD_PA_diffPathways.df$SampleType[22:40]="AD"
mayo_AD_PA_diffPathways.df=mayo_AD_PA_diffPathways.df%>%rownames_to_column("Pathway_Geneset")


#Read doo data
setwd("/Users/sandeepamberkar/Dropbox/MGH 3D AD RNAseq analysis JLa5/Complete datasets/")
doo_AD_3D_data.counts=fread("JLa5_complete_CPM_table.txt",sep = "\t",header = T,data.table = F)
doo_AD_3D_data.counts=doo_AD_3D_data.counts[which((rowSums(doo_AD_3D_data.counts>1)>=ncol(doo_AD_3D_data.counts)/3)=="TRUE"),]
doo_AD_3D_data.counts$EntrezID=unname(mapIds(x = EnsDb.Hsapiens.v79,keys = doo_AD_3D_data.counts$ID,column = "ENTREZID",keytype = "GENEID"))
doo_AD_3D_data_counts.agg=aggregate.data.frame(x = doo_AD_3D_data.counts[,-c(1:3,28)],by=list(entrez=doo_AD_3D_data.counts$EntrezID),mean)
rownames(doo_AD_3D_data_counts.agg)=doo_AD_3D_data_counts.agg$entrez
doo_AD_3D_data_counts.agg=doo_AD_3D_data_counts.agg[,-1]
#Columns 18:24 have the BACEi data
doo_BACEi_LY28BSI.SCE=single.chip.enrichment2(exprs = doo_AD_3D_data_counts.agg[,c(18:24)],geneset = pathprint.Hs.gs,statistic = "mean",progressBar = T)
#Columns 10:17 have DMSO data
doo_6wk_DMSO.SCE=single.chip.enrichment2(exprs = doo_AD_3D_data_counts.agg[,c(17:10)],geneset = pathprint.Hs.gs,statistic = "mean",progressBar = T)

doo_A5_ADdrug_DMSO.SCE=single.chip.enrichment2(exprs = doo_AD_3D_data_counts.agg[,c(4:6,10:13)],geneset = pathprint.Hs.gs,statistic = "mean",progressBar = T)

doo_AD_3D.diffPathways=list()

group=factor(c(rep("CONTROL",3),rep("AD",4)))
design_df=model.matrix(~0+group)
colnames(design_df)=levels(group)
contrasts_matrix=makeContrasts(CONTROL-AD,levels = design_df)
fit=lmFit(object = doo_BACEi_LY28BSI.SCE,design = design_df)
fit2=contrasts.fit(fit,contrasts_matrix)
fit2=eBayes(fit2,trend = T)
doo_AD_3D.diffPathways[[1]]=topTable(fit2,coef = 1,number = 633,p.value = 0.2,adjust.method = "BH")%>%rownames_to_column("Pathway/Geneset")

group=factor(c(rep("CONTROL_DMSO",4),rep("AD_DMSO",4)))
design_df=model.matrix(~0+group)
colnames(design_df)=levels(group)
contrasts_matrix=makeContrasts(CONTROL_DMSO-AD_DMSO,levels = design_df)
fit=lmFit(object = doo_6wk_DMSO.SCE,design = design_df)
fit2=contrasts.fit(fit,contrasts_matrix)
fit2=eBayes(fit2,trend = T)
doo_AD_3D.diffPathways[[2]]=topTable(fit2,coef = 1,number = 633,p.value = 0.2,adjust.method = "BH")%>%rownames_to_column("Pathway/Geneset")

group=factor(c(rep("CONTROL",3),rep("AD",4)))
design_df=model.matrix(~0+group)
colnames(design_df)=levels(group)
contrasts_matrix=makeContrasts(CONTROL-AD,levels = design_df)
fit=lmFit(object = doo_A5_ADdrug_DMSO.SCE,design = design_df)
fit2=contrasts.fit(fit,contrasts_matrix)
fit2=eBayes(fit2,trend = T)
doo_AD_3D.diffPathways[[3]]=topTable(fit2,coef = 1,number = 633,p.value = 0.2,adjust.method = "BH")%>%rownames_to_column("Pathway/Geneset")
names(doo_AD_3D.diffPathways)=c("BACEi_DMSO","6wk_DMSO","A5_ADdrug_DMSO")
write.table(doo_AD_3D.diffPathways$BACEi_DMSO,"/Users/sandeepamberkar/Work/Collaborations/Resiliome/Doo_BACEi_DMSO_diffPathways.txt",sep = "\t",col.names = T,row.names = F)
write.table(doo_AD_3D.diffPathways$`6wk_DMSO`,"/Users/sandeepamberkar/Work/Collaborations/Resiliome/Doo_6wk_AD_DMSO_diffPathways.txt",sep = "\t",col.names = T,row.names = F)
write.table(doo_AD_3D.diffPathways$A5_ADdrug_DMSO,"/Users/sandeepamberkar/Work/Collaborations/Resiliome/Doo_A5_ADdrug_DMSO_DiffPathways.txt",sep = "\t",col.names = T,row.names = F)

#DiffPathways Berchtold dataset
berchtold_series_matrix=getGEO(GEO = "GSE48350",GSEMatrix = T)
berchtold.exprs=data.frame(exprs(berchtold_series_matrix$GSE48350_series_matrix.txt.gz),stringsAsFactors = F)
berchtold.fData=fData(berchtold_series_matrix$GSE48350_series_matrix.txt.gz)
berchtold.pData=pData(phenoData(berchtold_series_matrix$GSE48350_series_matrix.txt.gz))
berchtold.exprs$EntrezID=berchtold.fData$ENTREZ_GENE_ID
berchtold_exprs.agg=aggregate.data.frame(x = berchtold.exprs[,-c(254)],by=list(entrez=berchtold.exprs$EntrezID),mean)
rownames(berchtold_exprs.agg)=berchtold_exprs.agg$entrez
berchtold_exprs.agg=berchtold_exprs.agg[,-which(colnames(berchtold_exprs.agg)=="entrez")]

berchtold_CONTROL_samples=berchtold.pData$geo_accession[grep(pattern = "AD",x = berchtold.pData$title[which(as.integer(berchtold.pData$`age (yrs):ch1`)>90)],invert = T)]
berchtold_AD_samples=berchtold.pData$geo_accession[grep(pattern = "AD",x = berchtold.pData$title[which(as.integer(berchtold.pData$`age (yrs):ch1`)>90)])]
berchtold_brain_regions=unique(berchtold.pData$`brain region:ch1`)[-5]
berchtold_brain_region.exprs=vector(mode = "list",length = 4)
names(berchtold_brain_region.exprs)=berchtold_brain_regions
for(i in 1:4){
  berchtold_brain_region.exprs[[i]]=berchtold_exprs.agg[,which(colnames(berchtold_exprs.agg)%in%c(intersect(berchtold_CONTROL_samples,berchtold.pData$geo_accession[which(berchtold.pData$`brain region:ch1`==berchtold_brain_regions[i])]),
                                                                                                  intersect(berchtold_AD_samples,berchtold.pData$geo_accession[which(berchtold.pData$`brain region:ch1`==berchtold_brain_regions[i])])))]
}
#Compute pathway expression by brain region
berchtold_brain_region.SCE=lapply(berchtold_brain_region.exprs,function(x)single.chip.enrichment2(exprs =x,geneset = pathprint.Hs.gs,statistic = "mean",progressBar = T))

#Compute differential pathways
berchtold_brain_region.diffPathways=vector(mode = "list",length = 4)
names(berchtold_brain_region.diffPathways)=names(berchtold_brain_region.SCE)
for(d in 1:4){
  group=factor(c(rep("CONTROL",length(intersect(berchtold_CONTROL_samples,colnames(berchtold_brain_region.SCE[[d]])))),rep("AD",length(intersect(berchtold_AD_samples,colnames(berchtold_brain_region.SCE[[d]]))))))
  design_df=model.matrix(~0+group)
  colnames(design_df)=levels(group)
  contrasts_matrix=makeContrasts(CONTROL-AD,levels = design_df)
  fit=lmFit(object = berchtold_brain_region.SCE[[d]],design = design_df)
  fit2=contrasts.fit(fit,contrasts_matrix)
  fit2=eBayes(fit2,trend = T)
  berchtold_brain_region.diffPathways[[d]]=topTable(fit2,coef = 1,number = 633,p.value = 1,adjust.method = "BH")%>%rownames_to_column("Pathway/Geneset")%>%dplyr::filter(P.Value<=0.05)
  write.table(berchtold_brain_region.diffPathways[[d]],paste("Berchtold",names(berchtold_brain_region.diffPathways)[d],"diffPathways.txt",sep = "_"),sep = "\t",col.names = T,row.names = F,quote=F)
}



mayo_pathAD_CRP_pathways=data.frame(mayo_AD_PA_norm_counts.SCE[which(rownames(mayo_AD_PA_norm_counts.SCE)=="Cytoplasmic Ribosomal Proteins (Wikipathways)"),])
colnames(mayo_pathAD_CRP_pathways)=c("PathExpr")
mayo_pathAD_CRP_pathways$Sampletype=c(rep("PathAging",21),rep("AD",19))
mayo_pathAD_CRP_pathways$PathExpr=mayo_pathAD_CRP_pathways$PathExpr
ggplot(data = mayo_pathAD_CRP_pathways,aes(x=Sampletype,y = PathExpr,fill=Sampletype))+geom_boxplot(outlier.colour = "red")+theme_light(base_size = 12)+ggtitle(label = "Cytoplasmic Ribosomal Proteins (Wikipathways) -- Pathway expression",subtitle = "Mayo Pathological aged vs. AD")

doo_6wk_DMSO_CRP_pathways=data.frame(doo_6wk_DMSO.SCE[which(rownames(doo_6wk_DMSO.SCE)=="Cytoplasmic Ribosomal Proteins (Wikipathways)"),])
doo_6wk_DMSO_CRP_pathways$Sampletype=c(rep("CONTROL",4),rep("AD",4))
colnames(doo_6wk_DMSO_CRP_pathways)=c("PathExpr","Sampletype")
ggplot(data = doo_6wk_DMSO_CRP_pathways,aes(x=Sampletype,y = PathExpr,fill=Sampletype))+geom_boxplot(outlier.color = "red")+theme_light(base_size = 12)+ggtitle(label = "Cytoplasmic Ribosomal Proteins (Wikipathways) -- Pathway expression",subtitle = "Doo 6 week DMSO, Control vs. AD")




