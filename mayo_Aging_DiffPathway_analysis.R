library(synapser)
library(dplyr)
library(magrittr)
library(data.table)
library(org.Hs.eg.db)
library(pathprint)
library(limma)
synLogin()
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
pathprint_membership.genes=lapply(pathprint.Hs.gs,function(x)unname(mapIds(x = org.Hs.eg.db,keys = as.character(x),keytype = "ENTREZID",column = "SYMBOL")))
names(pathprint_membership.genes)[grep(pattern = "/",names(pathprint_membership.genes))]=gsub(pattern = "/",replacement = " ",x = grep(pattern = "/",names(pathprint_membership.genes),value = T))


ds = c("chipframe", "genesets","pathprint.Hs.gs" ,"platform.thresholds")
data(list = ds)
mayo_PA_norm_counts.synapse=synGet('syn7440549')
mayo_PA_covariates.synapse=synGet('syn7437038')
mayo_AD_norm_counts.synapse=synGet('syn4650265')
mayo_AD_covariates.synapse=synGet('syn8466814')
mayo_AD_norm_counts=fread(mayo_AD_norm_counts.synapse$path,sep = "\t",header = T,stringsAsFactors = F,data.table = F)
mayo_PA_norm_counts=fread(mayo_PA_norm_counts.synapse$path,sep = "\t",header = T,stringsAsFactors = F,data.table = F)
mayo_PA_covariates=fread(mayo_PA_covariates.synapse$path,sep = "\t",header = T,stringsAsFactors = F,data.table = F)
mayo_AD_covariates=fread(mayo_AD_covariates.synapse$path,sep = "\t",header = T,stringsAsFactors = F,data.table = F)
mayo_PA_samples=mayo_PA_covariates%>%filter(AgeAtDeath>=90)%>%pull(SubjectID)
mayo_allPA_samples=mayo_PA_covariates$SubjectID
mayo_Control_samples=mayo_AD_covariates%>%filter(Tissue.SourceDiagnosis=="TCX.CONTROL"&AgeAtDeath>=90)%>%pull(SampleID)
mayo_yControl_samples=mayo_AD_covariates%>%filter(Tissue.SourceDiagnosis=="TCX.CONTROL"&AgeAtDeath<=75)%>%pull(SampleID)
mayo_oControl_samples=mayo_AD_covariates%>%dplyr::filter(Tissue.SourceDiagnosis=="TCX.CONTROL"&AgeAtDeath>=90)%>%pull(SampleID)
mayo_AD_samples=mayo_AD_covariates%>%filter(Tissue.SourceDiagnosis=="TCX.AD"&AgeAtDeath>=90)%>%pull(SampleID)
mayo_AD_samples2=mayo_AD_covariates%>%filter(Tissue.SourceDiagnosis=="TCX.AD")%>%pull(SampleID)
mayo_PA_AD_samples_To_analyse=c(mayo_PA_samples,mayo_AD_samples)
mayo_Control_AD_samples_To_analyse=c(mayo_Control_samples,mayo_AD_samples2)
mayo_yControl_PA_samples_To_analyse=c(mayo_yControl_samples,mayo_PA_samples)

mayo_AD_PA_norm_counts=cbind.data.frame(mayo_PA_norm_counts,mayo_AD_norm_counts)
mayo_yControl_Control_norm_counts=cbind.data.frame(mayo_AD_norm_counts[,colnames(mayo_AD_norm_counts)%in%mayo_yControl_samples],mayo_AD_norm_counts[,colnames(mayo_AD_norm_counts)%in%mayo_oControl_samples])
mayo_yControl_Control_norm_counts$ensembl_id=mayo_AD_norm_counts$ensembl_id
mayo_AD_PA_norm_counts=mayo_AD_PA_norm_counts[which((rowSums(mayo_AD_PA_norm_counts>0)>=ncol(mayo_AD_PA_norm_counts)/3)=="TRUE"),]
mayo_AD_norm_counts=mayo_AD_norm_counts[which((rowSums(mayo_AD_norm_counts>0)>=ncol(mayo_AD_norm_counts)/3)=="TRUE"),]
mayo_yControl_PA_norm_counts=mayo_yControl_PA_norm_counts[which((rowSums(mayo_yControl_PA_norm_counts>0)>=ncol(mayo_yControl_PA_norm_counts)/3)=="TRUE"),]
mayo_yControl_Control_norm_counts=mayo_yControl_oControl_norm_counts[which((rowSums(mayo_yControl_oControl_norm_counts>0)>=ncol(mayo_yControl_oControl_norm_counts)/3)=="TRUE"),]

mayo_AD_PA_norm_counts$EntrezID=unname(mapIds(x = org.Hs.eg.db,keys = mayo_AD_PA_norm_counts$ensembl_id,column = "ENTREZID",keytype = "ENSEMBL",multiVals = "first"))
mayo_AD_norm_counts$EntrezID=unname(mapIds(x = org.Hs.eg.db,keys = mayo_AD_norm_counts$ensembl_id,column = "ENTREZID",keytype = "ENSEMBL",multiVals = "first"))
mayo_yControl_PA_norm_counts$EntrezID=unname(mapIds(x = org.Hs.eg.db,keys = mayo_yControl_PA_norm_counts$ensembl_id,column = "ENTREZID",keytype = "ENSEMBL",multiVals = "first"))
mayo_yControl_oControl_norm_counts$EntrezID=unname(mapIds(x = EnsDb.Hsapiens.v79,keys = mayo_yControl_oControl_norm_counts$ensembl_id,column = "ENTREZID",keytype = "GENEID",multiVals = "first"))


mayo_AD_PA_norm_counts.agg=aggregate.data.frame(x = mayo_AD_PA_norm_counts[,-c(1,43,320)],by = list(entrez=mayo_AD_PA_norm_counts$EntrezID),mean)
mayo_AD_norm_counts.agg=aggregate.data.frame(x = mayo_AD_norm_counts[,-c(1,278)],by = list(entrez=mayo_AD_norm_counts$EntrezID),mean)
mayo_yControl_PA_norm_counts.agg=aggregate.data.frame(x = mayo_yControl_PA_norm_counts[,-c(15,57)],by = list(entrez=mayo_yControl_PA_norm_counts$EntrezID),mean)
mayo_yControl_oControl_norm_counts.agg=aggregate.data.frame(x = mayo_yControl_PA_norm_counts[,-c(15,57)],by = list(entrez=mayo_yControl_PA_norm_counts$EntrezID),mean)
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
    mayo_cohort.Diffpathways[[i]]=topTable(fit2,coef = 1,number = 633,p.value = 1,adjust.method = "BH")%>%rownames_to_column("Geneset_Pathway")
    mayo_cohort.Diffpathways[[i]]$Pathway_MemberGenes=unlist(lapply(mayo_cohort.Diffpathways[[i]]$Geneset_Pathway,function(x)paste(unlist(unname(pathprint_membership.genes[names(pathprint_membership.genes)%in%x])),collapse = ",")))
  }
  if(i==2){
    group=factor(c(rep("Control",length(mayo_Control_samples)),rep("AD",length(mayo_AD_samples2))))
    design_df=model.matrix(~0+group)
    colnames(design_df)=levels(group)
    contrasts_matrix=makeContrasts(Control-AD,levels = design_df)
    fit=lmFit(object = mayo_AD_norm_counts.SCE,design = design_df)
    fit2=contrasts.fit(fit,contrasts_matrix)
    fit2=eBayes(fit2,trend = T)
    mayo_cohort.Diffpathways[[i]]=topTable(fit2,coef = 1,number = 633,p.value = 1,adjust.method = "BH")%>%rownames_to_column("Geneset_Pathway")
    mayo_cohort.Diffpathways[[i]]$Pathway_MemberGenes=unlist(lapply(mayo_cohort.Diffpathways[[i]]$Geneset_Pathway,function(x)paste(unlist(unname(pathprint_membership.genes[names(pathprint_membership.genes)%in%x])),collapse = ",")))
  }
  if(i==3){
    group=factor(c(rep("yControl",length(mayo_yControl_samples)),rep("PA",length(mayo_allPA_samples))))
    design_df=model.matrix(~0+group)
    colnames(design_df)=levels(group)
    contrasts_matrix=makeContrasts(yControl-PA,levels = design_df)
    fit=lmFit(object = mayo_yControl_PA_norm_counts.SCE,design = design_df)
    fit2=contrasts.fit(fit,contrasts_matrix)
    fit2=eBayes(fit2,trend = T)
    mayo_cohort.Diffpathways[[i]]=topTable(fit2,coef = 1,number = 633,p.value = 1,adjust.method = "BH")%>%rownames_to_column("Geneset_Pathway")
    mayo_cohort.Diffpathways[[i]]$Pathway_MemberGenes=unlist(lapply(mayo_cohort.Diffpathways[[i]]$Geneset_Pathway,function(x)paste(unlist(unname(pathprint_membership.genes[names(pathprint_membership.genes)%in%x])),collapse = ",")))
  }
}
mayo_cohort.Diffpathways=lapply(mayo_cohort.Diffpathways,function(x)x%>%mutate(log_FoldChange=log2(abs(logFC)))%>%dplyr::filter(`adj.P.Val`<=0.1)%>%dplyr::select(c(Geneset_Pathway,log_FoldChange,AveExpr,t,`P.Value`,`adj.P.Val`,B,Pathway_MemberGenes)))
mayo_cohort.Diffpathways=lapply(mayo_cohort.Diffpathways,function(x){x$t=round(x$t,digits = 3);x$`adj.P.Val`=round(x$`adj.P.Val`,digits = 3);x$B=round(x$B,digits = 3);x$log_FoldChange=round(x$log_FoldChange,digits = 3);x$AveExpr=round(x$AveExpr,digits = 3);x})
fwrite(mayo_cohort.Diffpathways$`PA-AD`,"/Users/sandeepamberkar/Work/Collaborations/Resiliome/Mayo/Mayo_PathAging_AD_DiffExprPathways.txt",sep = "\t",col.names = T,quote = F)
fwrite(mayo_cohort.Diffpathways$`Control-AD`,"/Users/sandeepamberkar/Work/Collaborations/Resiliome/Mayo/Mayo_Control_AD_DiffExprPathways.txt",sep = "\t",col.names = T,quote = F)
fwrite(mayo_cohort.Diffpathways$`yControl-PA`,"/Users/sandeepamberkar/Work/Collaborations/Resiliome/Mayo/Mayo_youngControl_PA_DiffExprPathways.txt",sep = "\t",col.names = T,quote = F)


mayo_PA_AD_yControl_PA.CommonDEPs=intersect(mayo_cohort.Diffpathways$`PA-AD`%>%filter(adj.P.Val<=0.1)%>%pull('Geneset_Pathway'),mayo_cohort.Diffpathways$`yControl-PA`%>%filter(adj.P.Val<=0.1)%>%pull('Geneset_Pathway'))
mayo_Ctrl_AD_yControl_PA.CommonDEPs=intersect(mayo_cohort.Diffpathways$`Control-AD`%>%filter(adj.P.Val<=0.1)%>%pull('Geneset_Pathway'),mayo_cohort.Diffpathways$`yControl-PA`%>%filter(adj.P.Val<=0.1)%>%pull('Geneset_Pathway'))



#For boxplots
mayo_AD_PA_diffPathways.df=data.frame(mayo_AD_PA_norm_counts.SCE[which(rownames(mayo_AD_PA_norm_counts.SCE)%in%mayo_cohort.Diffpathways$`PA-AD`$`Pathway/GS`),])
mayo_AD_PA_diffPathways.df$SampleType="SAMPLE"
mayo_AD_PA_diffPathways.df$SampleType[1:21]="Path.Aging"
mayo_AD_PA_diffPathways.df$SampleType[22:40]="AD"
mayo_AD_PA_diffPathways.df=mayo_AD_PA_diffPathways.df%>%rownames_to_column("Pathway_Geneset")


#Read doo data - Drug treated and mod pathology
setwd("/Users/sandeepamberkar/Dropbox/MGH 3D AD RNAseq analysis JLa5/Complete datasets/")
doo_AD_3D_data.counts=fread("JLa5_complete_CPM_table.txt",sep = "\t",header = T,data.table = F)
doo_AD_3D_data.counts=doo_AD_3D_data.counts[which((rowSums(doo_AD_3D_data.counts>1)>=ncol(doo_AD_3D_data.counts)/3)=="TRUE"),]
doo_AD_3D_data.counts$EntrezID=unname(mapIds(x = EnsDb.Hsapiens.v79,keys = doo_AD_3D_data.counts$ID,column = "ENTREZID",keytype = "GENEID"))
doo_AD_3D_data_counts.agg=aggregate.data.frame(x = doo_AD_3D_data.counts[,-c(1:3,28)],by=list(entrez=doo_AD_3D_data.counts$EntrezID),mean)
rownames(doo_AD_3D_data_counts.agg)=doo_AD_3D_data_counts.agg$entrez
doo_AD_3D_data_counts.agg=doo_AD_3D_data_counts.agg[,-1]
#Read doo data - Doo high pathology
setwd("/Users/sandeepamberkar/Work/Collaborations/Resiliome/Doo_3D_AD/")
doo_AD_3D_high_path_data.counts=fread("040418_part1_complete_RPKM_table.txt",sep = "\t",header = T,data.table = F)
doo_AD_3D_Low_High_Path.counts=doo_AD_3D_high_path_data.counts[,c(1:2,3:5,)]
doo_AD_3D_high_path_data.counts=doo_AD_3D_high_path_data.counts[which((rowSums(doo_AD_3D_high_path_data.counts>1)>=ncol(doo_AD_3D_high_path_data.counts)/3)=="TRUE"),]
doo_AD_3D_high_path_data.counts$EntrezID=unname(mapIds(x = EnsDb.Hsapiens.v79,keys = doo_AD_3D_high_path_data.counts$ID,column = "ENTREZID",keytype = "GENEID"))
doo_AD_3D_high_path_data_counts.agg=aggregate.data.frame(x = doo_AD_3D_high_path_data.counts[-c(1:3,19)],by=list(entrez=doo_AD_3D_high_path_data.counts$EntrezID),mean)
rownames(doo_AD_3D_high_path_data_counts.agg)=doo_AD_3D_high_path_data_counts.agg$entrez
doo_AD_3D_high_path_data_counts.agg=doo_AD_3D_high_path_data_counts.agg[,-1]

#Columns 7:9 are D4, 10:12 H10
doo_AD_3D_high_path_D4.SCE=single.chip.enrichment2(exprs = doo_AD_3D_high_path_data_counts.agg[,c(1:3,7:9)],geneset = pathprint.Hs.gs,statistic = "mean",progressBar = T)

doo_AD_3D_high_path_H10.SCE=single.chip.enrichment2(exprs = doo_AD_3D_high_path_data_counts.agg[,c(1:3,10:12)],geneset = pathprint.Hs.gs,statistic = "mean",progressBar = T)

#Columns 18:24 have the BACEi data
doo_BACEi_LY28BSI.SCE=single.chip.enrichment2(exprs = doo_AD_3D_data_counts.agg[,c(18:24)],geneset = pathprint.Hs.gs,statistic = "mean",progressBar = T)
#Columns 10:17 have DMSO data
doo_6wk_DMSO.SCE=single.chip.enrichment2(exprs = doo_AD_3D_data_counts.agg[,c(17:10)],geneset = pathprint.Hs.gs,statistic = "mean",progressBar = T)

doo_A5_ADdrug_DMSO.SCE=single.chip.enrichment2(exprs = doo_AD_3D_data_counts.agg[,c(4:6,10:13)],geneset = pathprint.Hs.gs,statistic = "mean",progressBar = T)

doo_Lo_High_Pathology.SCE=single.chip.enrichment2(exprs = doo_AD_3D_data_counts.agg[,c(4:6,10:13)],geneset = pathprint.Hs.gs,statistic = "mean",progressBar = T)
doo_AD_3D.diffPathways=list()

group=factor(c(rep("CONTROL",3),rep("AD",4)))
design_df=model.matrix(~0+group)
colnames(design_df)=levels(group)
contrasts_matrix=makeContrasts(CONTROL-AD,levels = design_df)
fit=lmFit(object = doo_BACEi_LY28BSI.SCE,design = design_df)
fit2=contrasts.fit(fit,contrasts_matrix)
fit2=eBayes(fit2,trend = T)
doo_AD_3D.diffPathways[[1]]=topTable(fit2,coef = 1,number = 633,p.value = 1,adjust.method = "BH")%>%rownames_to_column("Geneset_Pathway")

group=factor(c(rep("CONTROL_DMSO",4),rep("AD_DMSO",4)))
design_df=model.matrix(~0+group)
colnames(design_df)=levels(group)
contrasts_matrix=makeContrasts(CONTROL_DMSO-AD_DMSO,levels = design_df)
fit=lmFit(object = doo_6wk_DMSO.SCE,design = design_df)
fit2=contrasts.fit(fit,contrasts_matrix)
fit2=eBayes(fit2,trend = T)
doo_AD_3D.diffPathways[[2]]=topTable(fit2,coef = 1,number = 633,p.value = 1,adjust.method = "BH")%>%rownames_to_column("Geneset_Pathway")

group=factor(c(rep("CONTROL",3),rep("AD",4)))
design_df=model.matrix(~0+group)
colnames(design_df)=levels(group)
contrasts_matrix=makeContrasts(CONTROL-AD,levels = design_df)
fit=lmFit(object = doo_A5_ADdrug_DMSO.SCE,design = design_df)
fit2=contrasts.fit(fit,contrasts_matrix)
fit2=eBayes(fit2,trend = T)
doo_AD_3D.diffPathways[[3]]=topTable(fit2,coef = 1,number = 633,p.value = 1,adjust.method = "BH")%>%rownames_to_column("Geneset_Pathway")

group=factor(c(rep("CONTROL",3),rep("AD",3)))
design_df=model.matrix(~0+group)
colnames(design_df)=levels(group)
contrasts_matrix=makeContrasts(CONTROL-AD,levels = design_df)
fit=lmFit(object = doo_AD_3D_high_path_D4.SCE,design = design_df)
fit2=contrasts.fit(fit,contrasts_matrix)
fit2=eBayes(fit2,trend = T)
doo_AD_3D.diffPathways[[4]]=topTable(fit2,coef = 1,number = 633,p.value = 1,adjust.method = "BH")%>%rownames_to_column("Geneset_Pathway")

group=factor(c(rep("CONTROL",3),rep("AD",3)))
design_df=model.matrix(~0+group)
colnames(design_df)=levels(group)
contrasts_matrix=makeContrasts(CONTROL-AD,levels = design_df)
fit=lmFit(object = doo_AD_3D_high_path_H10.SCE,design = design_df)
fit2=contrasts.fit(fit,contrasts_matrix)
fit2=eBayes(fit2,trend = T)
doo_AD_3D.diffPathways[[5]]=topTable(fit2,coef = 1,number = 633,p.value = 1,adjust.method = "BH")%>%rownames_to_column("Geneset_Pathway")

doo_AD_3D.diffPathways[[1]]$Pathway_MemberGenes=unlist(lapply(doo_AD_3D.diffPathways[[1]]$Geneset_Pathway,function(x)paste(unlist(unname(pathprint_membership.genes[names(pathprint_membership.genes)%in%x])),collapse = ",")))
doo_AD_3D.diffPathways[[2]]$Pathway_MemberGenes=unlist(lapply(doo_AD_3D.diffPathways[[2]]$Geneset_Pathway,function(x)paste(unlist(unname(pathprint_membership.genes[names(pathprint_membership.genes)%in%x])),collapse = ",")))
doo_AD_3D.diffPathways[[3]]$Pathway_MemberGenes=unlist(lapply(doo_AD_3D.diffPathways[[3]]$Geneset_Pathway,function(x)paste(unlist(unname(pathprint_membership.genes[names(pathprint_membership.genes)%in%x])),collapse = ",")))
doo_AD_3D.diffPathways[[4]]$Pathway_MemberGenes=unlist(lapply(doo_AD_3D.diffPathways[[4]]$Geneset_Pathway,function(x)paste(unlist(unname(pathprint_membership.genes[names(pathprint_membership.genes)%in%x])),collapse = ",")))
doo_AD_3D.diffPathways[[5]]$Pathway_MemberGenes=unlist(lapply(doo_AD_3D.diffPathways[[5]]$Geneset_Pathway,function(x)paste(unlist(unname(pathprint_membership.genes[names(pathprint_membership.genes)%in%x])),collapse = ",")))
doo_AD_3D.diffPathways=lapply(doo_AD_3D.diffPathways,function(x)x%>%mutate(log2_FC=log2(abs(logFC)))%>%dplyr::select(c(Geneset_Pathway,log2_FC,AveExpr,t,`P.Value`,`adj.P.Val`,B,Pathway_MemberGenes)))
names(doo_AD_3D.diffPathways)=c("BACEi_Control","6wk_DMSO","A5_ADdrug_DMSO","D4_Control","H10_Control")

write.table(doo_AD_3D.diffPathways$BACEi_Control,"/Users/sandeepamberkar/Work/Collaborations/Resiliome/Doo_BACEi_Control_diffPathways.txt",sep = "\t",col.names = T,row.names = F)
write.table(doo_AD_3D.diffPathways$`6wk_DMSO`,"/Users/sandeepamberkar/Work/Collaborations/Resiliome/Doo_6wk_AD_DMSO_diffPathways.txt",sep = "\t",col.names = T,row.names = F)
write.table(doo_AD_3D.diffPathways$A5_ADdrug_DMSO,"/Users/sandeepamberkar/Work/Collaborations/Resiliome/Doo_A5_ADdrug_DMSO_DiffPathways.txt",sep = "\t",col.names = T,row.names = F)
write.table(doo_AD_3D.diffPathways$D4_Control,"/Users/sandeepamberkar/Work/Collaborations/Resiliome/Doo_D4_Control_DiffPathways.txt",sep = "\t",col.names = T,row.names = F)
write.table(doo_AD_3D.diffPathways$H10_Control,"/Users/sandeepamberkar/Work/Collaborations/Resiliome/Doo_H10_Control_DiffPathways.txt",sep = "\t",col.names = T,row.names = F)

#Intersect Mayo and Doo High pathology DEPs

Cognition_HiPathology_DEPs.list=lapply(mayo_cohort.Diffpathways,function(x)intersect(x%>%dplyr::filter(adj.P.Val<=0.1)%>%pull(Geneset_Pathway),union(doo_AD_3D.diffPathways$H10_Control%>%dplyr::filter(adj.P.Val<=0.1)%>%pull(Geneset_Pathway),doo_AD_3D.diffPathways$D4_Control%>%dplyr::filter(adj.P.Val<=0.1)%>%pull(Geneset_Pathway))))
Cognition_LoPathology_DEPs.list=lapply(mayo_cohort.Diffpathways,function(x)intersect(x%>%dplyr::filter(adj.P.Val<=0.1)%>%pull(Geneset_Pathway),union(doo_AD_3D.diffPathways$`6wk_DMSO`%>%dplyr::filter(adj.P.Val<=0.2)%>%pull(Geneset_Pathway),doo_AD_3D.diffPathways$A5_ADdrug_DMSO%>%dplyr::filter(adj.P.Val<=0.2)%>%pull(Geneset_Pathway))))
write(Cognition_HiPathology_DEPs.list$`PA-AD`,'../Hi-CDP-DEPs.txt',sep = "\n")
write(Cognition_HiPathology_DEPs.list$`yControl-PA`,'../Hi-AP-DEPs.txt',sep = "\n")
write(Cognition_LoPathology_DEPs.list$`PA-AD`,'../Lo-CDP-DEPs.txt',sep = "\n")
write(Cognition_LoPathology_DEPs.list$`yControl-PA`,'../Lo-AP-DEPs.txt',sep = "\n")





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
  berchtold_brain_region.diffPathways[[d]]=topTable(fit2,coef = 1,number = 633,p.value = 1,adjust.method = "BH")%>%rownames_to_column("Geneset_Pathway")%>%dplyr::filter(P.Value<=0.05)
  berchtold_brain_region.diffPathways[[d]]$Pathway_MemberGenes=unlist(lapply(berchtold_brain_region.diffPathways[[d]]$Geneset_Pathway,function(x)paste(unlist(unname(pathprint_membership.genes[names(pathprint_membership.genes)%in%x])),collapse = ",")))
  write.table(berchtold_brain_region.diffPathways[[d]],paste("Berchtold",names(berchtold_brain_region.diffPathways)[d],"diffPathways.txt",sep = "_"),sep = "\t",col.names = T,row.names = F,quote=F)
}



mayo_pathAD_CRP_pathways=data.frame(mayo_AD_PA_norm_counts.SCE[which(rownames(mayo_AD_PA_norm_counts.SCE)=="Cytoplasmic Ribosomal Proteins (Wikipathways)"),])
colnames(mayo_pathAD_CRP_pathways)=c("PathExpr")
mayo_pathAD_CRP_pathways$Sampletype=c(rep("Path_Aging",21),rep("Mayo-AD",19))
mayo_pathAD_CRP_pathways$PathExpr=mayo_pathAD_CRP_pathways$PathExpr
ggplot(data = mayo_pathAD_CRP_pathways,aes(x=Sampletype,y = PathExpr,fill=Sampletype))+geom_boxplot(outlier.colour = "red")+theme_light(base_size = 12)+ggtitle(label = "Cytoplasmic Ribosomal Proteins (Wikipathways) -- Pathway expression",subtitle = "Mayo Pathological aged vs. AD")

mayo_yControl_PA_CRP_pathways.df=data.frame(PathExpr=mayo_yControl_PA_norm_counts.SCE[which(rownames(mayo_yControl_PA_norm_counts.SCE)=='Translation Factors (Wikipathways)'),])
mayo_yControl_PA_CRP_pathways.df$Sampletype=c(rep('young_Controls',14),rep('Path_Aging',41))
colnames(mayo_yControl_PA_CRP_pathways.df)=c("PathExpr","Sampletype")

doo_6wk_DMSO_CRP_pathways=data.frame(doo_6wk_DMSO.SCE[which(rownames(doo_6wk_DMSO.SCE)=="Cytoplasmic Ribosomal Proteins (Wikipathways)"),])
doo_6wk_DMSO_CRP_pathways$Sampletype=c(rep("Doo-Control-DMSO",4),rep("A5-DMSO",4))
colnames(doo_6wk_DMSO_CRP_pathways)=c("PathExpr","Sampletype")
ggplot(data = doo_6wk_DMSO_CRP_pathways,aes(x=Sampletype,y = PathExpr,fill=Sampletype))+geom_boxplot(outlier.color = "red")+theme_light(base_size = 12)+ggtitle(label = "Cytoplasmic Ribosomal Proteins (Wikipathways) -- Pathway expression",subtitle = "Doo 6 week DMSO, Control vs. AD")


# doo_6wk_DMSO_Ribosome_pathway=data.frame(doo_6wk_DMSO.SCE[which(rownames(doo_6wk_DMSO.SCE)=="Ribosome (KEGG)"),])
# doo_6wk_DMSO_Ribosome_pathway$SampleType=c(rep("Control-DMSO",4),rep("A5-DMSO",4))
# colnames(doo_6wk_DMSO_Ribosome_pathway)=c("PathExpr","Sampletype")
# ggplot(data = doo_6wk_DMSO_Ribosome_pathway,aes(x=Sampletype,y = PathExpr,fill=Sampletype))+geom_boxplot(outlier.color = "red")+theme_light(base_size = 12)+theme(axis.text.x = element_text(size = 12))+ggtitle(label = "Ribosome (KEGG) -- Pathway expression",subtitle = "Doo 6 week DMSO, Control vs. AD")

doo_D4_CRP_pathways.df=data.frame(doo_AD_3D_high_path_D4.SCE[which(rownames(doo_AD_3D_high_path_D4.SCE)=='Cytoplasmic Ribosomal Proteins (Wikipathways)'),])
doo_D4_CRP_pathways.df$SampleType=c(rep("Doo-Control2",3),rep("D4",3))
colnames(doo_D4_CRP_pathways.df)=c("PathExpr","Sampletype")

doo_H10_CRP_pathways.df=data.frame(doo_AD_3D_high_path_H10.SCE[which(rownames(doo_AD_3D_high_path_H10.SCE)=='Cytoplasmic Ribosomal Proteins (Wikipathways)'),])
doo_H10_CRP_pathways.df$SampleType=c(rep("Doo-Control2",3),rep("H10",3))
colnames(doo_H10_CRP_pathways.df)=c("PathExpr","Sampletype")

ggplot(data = rbind.data.frame(mayo_pathAD_CRP_pathways,mayo_yControl_PA_CRP_pathways.df)%>%mutate(PathExpr=log10(PathExpr)),aes(x=Sampletype,y = PathExpr,fill=Sampletype))+geom_boxplot(outlier.color = "red")+theme_light(base_size = 14)+theme(axis.text.x = element_text(size = 14))+ggtitle(label = "Cytoplasmic Ribosomal Proteins (Wikipathways)",subtitle = "All datasets")+theme(legend.text = element_text(size = 11))+ylab("log(Pathway.Expression)")

#Plot Translation Factors (Wikipathways), common among yControls, PA and Doo 6wk AD
mayo_yControl_PA_TranslFactors_WikiPathways.df=data.frame(PathExpr=mayo_yControl_PA_norm_counts.SCE[which(rownames(mayo_yControl_PA_norm_counts.SCE)=='Translation Factors (Wikipathways)'),])
mayo_yControl_PA_TranslFactors_WikiPathways.df$Sampletype=c(rep('young_Controls',14),rep('Path_Aging',41))
ggplot(data = mayo_yControl_PA_TranslFactors_WikiPathways.df,aes(x=Sampletype,y = PathExpr,fill=Sampletype))+geom_boxplot(outlier.color = "red")+theme_light(base_size = 12)+theme(axis.text.x = element_text(size = 12))+ggtitle(label = "Translation Factors (Wikipathways)",subtitle = "Mayo young Controls (<75) ~ Path.Aged (>=90), no AD")

# mayo_PA_AD_TranslFactors_WikiPathways.df=data.frame(PathExpr=mayo_AD_PA_norm_counts.SCE[which(rownames(mayo_yControl_PA_norm_counts.SCE)=='Translation Factors (Wikipathways)'),])
# mayo_PA_AD_TranslFactors_WikiPathways.df$Sampletype=c(rep('young_Controls',14),rep('Path_Aging',41))
# ggplot(data = mayo_yControl_PA_TranslFactors_WikiPathways.df,aes(x=Sampletype,y = PathExpr,fill=Sampletype))+geom_boxplot(outlier.color = "red")+theme_light(base_size = 12)+theme(axis.text.x = element_text(size = 12))+ggtitle(label = "Translation Factors (Wikipathways)",subtitle = "Mayo young Controls (<75) ~ Path.Aged (>=90), no AD")

mayo_PA_AD_TranslFactors_WikiPathways.df=data.frame(mayo_AD_PA_norm_counts.SCE[which(rownames(mayo_AD_PA_norm_counts.SCE)=='Translation Factors (Wikipathways)'),])
colnames(mayo_PA_AD_TranslFactors_WikiPathways.df)=c("PathExpr")
mayo_PA_AD_TranslFactors_WikiPathways.df$Sampletype=c(rep("Path_Aging",21),rep("Mayo_AD",19))
mayo_PA_AD_TranslFactors_WikiPathways.df$PathExpr=mayo_PA_AD_TranslFactors_WikiPathways.df$PathExpr
ggplot(data = mayo_PA_AD_TranslFactors_WikiPathways.df,aes(x=Sampletype,y = PathExpr,fill=Sampletype))+geom_boxplot(outlier.colour = "red")+theme_light(base_size = 12)+ggtitle(label = "Translation Factors (Wikipathways) -- Pathway expression",subtitle = "Mayo Path.Aged (>=90) vs. AD (>=90)")

doo_6wk_DMSO_TranslFactors_WikiPathways.df=data.frame(doo_6wk_DMSO.SCE[which(rownames(doo_6wk_DMSO.SCE)=='Translation Factors (Wikipathways)'),])
doo_6wk_DMSO_TranslFactors_WikiPathways.df$SampleType=c(rep("Doo-Control-DMSO",4),rep("A5-DMSO",4))
colnames(doo_6wk_DMSO_TranslFactors_WikiPathways.df)=c("PathExpr","Sampletype")
ggplot(data = doo_6wk_DMSO_TranslFactors_WikiPathways.df,aes(x=Sampletype,y = PathExpr,fill=Sampletype))+geom_boxplot(outlier.color = "red")+theme_light(base_size = 12)+theme(axis.text.x = element_text(size = 12))+ggtitle(label = "Translation Factors (Wikipathways)",subtitle = "Doo 6 week DMSO, Control vs. AD")
ggplot(data = rbind.data.frame(mayo_yControl_PA_TranslFactors_WikiPathways.df,mayo_PA_AD_TranslFactors_WikiPathways.df,doo_6wk_DMSO_TranslFactors_WikiPathways.df),aes(x=Sampletype,y = PathExpr,fill=Sampletype))+geom_boxplot(outlier.color = "red")+theme_light(base_size = 12)+theme(axis.text.x = element_text(size = 12))+ggtitle(label = "Translation Factors (Wikipathways)",subtitle = "All datasets")

doo_D4_TranslFactors_WikiPathways.df=data.frame(doo_AD_3D_high_path_D4.SCE[which(rownames(doo_AD_3D_high_path_D4.SCE)=='Translation Factors (Wikipathways)'),])
doo_D4_TranslFactors_WikiPathways.df$SampleType=c(rep("Doo-Control2",3),rep("D4",3))
colnames(doo_D4_TranslFactors_WikiPathways.df)=c("PathExpr","Sampletype")

doo_H10_TranslFactors_WikiPathways.df=data.frame(doo_AD_3D_high_path_H10.SCE[which(rownames(doo_AD_3D_high_path_H10.SCE)=='Translation Factors (Wikipathways)'),])
doo_H10_TranslFactors_WikiPathways.df$SampleType=c(rep("Doo-Control2",3),rep("H10",3))
colnames(doo_H10_TranslFactors_WikiPathways.df)=c("PathExpr","Sampletype")

ggplot(data = rbind.data.frame(mayo_yControl_PA_TranslFactors_WikiPathways.df,mayo_PA_AD_TranslFactors_WikiPathways.df,doo_6wk_DMSO_TranslFactors_WikiPathways.df,doo_H10_TranslFactors_WikiPathways.df,doo_D4_TranslFactors_WikiPathways.df)%>%mutate(PathExpr=log10(PathExpr)),aes(x=Sampletype,y = PathExpr,fill=Sampletype))+geom_boxplot(outlier.color = "red")+theme_light(base_size = 12)+theme(axis.text.x = element_text(size = 12))+ggtitle(label = "Translation Factors (Wikipathways)",subtitle = "All datasets")+theme(legend.text = element_text(size = 11))+ylab("log(Pathway.Expression)")

#Compare Mayo Control-AD, Doo A5-Drug-AD DEPs using permutation test
random_mayo_DEPs.PA_AD=replicate(n = 10000,expr = sample(x = names(pathprint.Hs.gs),size = length(mayo_cohort.Diffpathways$`PA-AD`%>%dplyr::filter(adj.P.Val<=0.1)%>%pull(Geneset_Pathway)),replace = F))
random_mayo_DEPs.yCtrl_PA=replicate(n = 10000,expr = sample(x = names(pathprint.Hs.gs),size = length(mayo_cohort.Diffpathways$`yControl-PA`%>%dplyr::filter(adj.P.Val<=0.1)%>%pull(Geneset_Pathway)),replace = F))
random_mayo_DEPs.Ctrl_AD=replicate(n = 10000,expr = sample(x = names(pathprint.Hs.gs),size = length(mayo_cohort.Diffpathways$`Control-AD`%>%dplyr::filter(adj.P.Val<=0.1)%>%pull(Geneset_Pathway)),replace = F))
random_doo_A5_ADdrug_DMSO_DEPs=replicate(n = 1000,expr = sample(x = names(pathprint.Hs.gs),size = length(doo_AD_3D.diffPathways$A5_ADdrug_DMSO$Geneset_Pathway),replace = F))
random_DEP_overlaps.list=list()
for(l in 1:10000){
  random_DEP_overlaps.list[[l]]=length(intersect(random_mayo_DEPs.PA_AD[,l],random_mayo_DEPs.Ctrl_AD[,l]))
}

#Compute genetic enrichments in DEPs

tanzi_ranked_genes=read.table("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Gene_Lists/All_genes_ranked.txt",header = T,sep = "\t",stringsAsFactors = F)
tanzi_ranked_genes.great=read.table("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Gene_Lists/Ranked_genes_corrected_allele_top_4000_SNPs_genes_from_great.txt",header = T,sep = "\t",stringsAsFactors = F)
tanzi_ranked_genes.great$SNP_Rank=as.numeric(gsub(pattern = "SNP_rank_",replacement = "",x = tanzi_ranked_genes.great$SNP_Rank))
tanzi_ranked_genes.AD=unique(tanzi_ranked_genes.great%>%mutate(SNP_Rank=as.integer(gsub(pattern = "SNP_rank_",replacement = "",x = SNP_Rank)))%>%dplyr::filter(Corrected_association=="Case"&SNP_Rank<=1000)%>%pull(Gene))
tanzi_ranked_genes.Control=unique(tanzi_ranked_genes.great%>%mutate(SNP_Rank=as.integer(gsub(pattern = "SNP_rank_",replacement = "",x = SNP_Rank)))%>%dplyr::filter(Corrected_association=="Control"&SNP_Rank<=1000)%>%pull(Gene))
tanzi_damaging_snps=scan("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Gene_Lists/Damaging_exonic_GeneSymbol.txt",sep = "\n",what = "char")
tanzi_protective_snps=scan("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Gene_Lists/Protective_exonic_GeneSymbol.txt",sep = "\n",what = "char")
Control_SNP_genes=scan("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Gene_Lists/Control_genes_unfiltered.txt",sep = "\n",what = "char",skip = 1)
pathprint=read.gmt(gmtfile = "/Users/sandeepamberkar/Work/Software/pathprint.gmt")

tanzi_damaging.pathprint=enricher(gene = unlist(unname(mapIds(x = org.Hs.eg.db,keys = tanzi_damaging_snps,column = "ENTREZID",keytype = "SYMBOL"))),TERM2GENE = pathprint,pvalueCutoff = 0.1,pAdjustMethod = "BH")
tanzi_protective.pathprint=enricher(gene = unlist(unname(mapIds(x = org.Hs.eg.db,keys = tanzi_protective_snps,column = "ENTREZID",keytype = "SYMBOL"))),TERM2GENE = pathprint,pvalueCutoff = 0.2,pAdjustMethod = "BH")
tanzi_SAGs.pathprint=enricher(gene = unlist(unname(mapIds(x = org.Hs.eg.db,keys = union(tanzi_ranked_genes.Control,tanzi_ranked_genes.AD),column = "ENTREZID",keytype = "SYMBOL"))),TERM2GENE = pathprint,pvalueCutoff = 0.1,pAdjustMethod = "BH")
Control_SNP_genes.pathprint=enricher(gene = unlist(unname(mapIds(x = org.Hs.eg.db,keys = Control_SNP_genes,column = "ENTREZID",keytype = "SYMBOL"))),TERM2GENE = pathprint,pvalueCutoff = 0.2,pAdjustMethod = "BH")

#PCA on Control and PA samples
mayo_AD_PA_counts=cbind.data.frame(mayo_AD_norm_counts[,c('ensembl_id',mayo_Control_samples)],mayo_PA_norm_counts[,mayo_PA_samples])
mayo_AD_PA_counts=mayo_AD_PA_counts[which((rowSums(mayo_AD_PA_counts>0)>=ncol(mayo_AD_PA_counts)/3)=="TRUE"),]
mayo_AD_PA_counts$GeneSymbol=unname(mapIds(x = EnsDb.Hsapiens.v79,keys = mayo_AD_PA_counts$ensembl_id,keytype = "GENEID",column = "SYMBOL"))
mayo_AD_PA_counts.agg=aggregate.data.frame(x = mayo_AD_PA_counts[,-c(1,42)],by = list(symbol=mayo_AD_PA_counts$GeneSymbol),mean)
rownames(mayo_AD_PA_counts.agg)=mayo_AD_PA_counts.agg$symbol
mayo_AD_PA_counts.agg=mayo_AD_PA_counts.agg[,-1]
colnames(mayo_AD_PA_counts.agg)=colnames(mayo_AD_PA_counts)
mayo_AD_PA_counts.agg2=data.frame(t(mayo_AD_PA_counts.agg))
mayo_AD_PA_counts.agg2=cbind.data.frame(SampleType=c(rep("Mayo_Control",19),rep("Mayo_PA",21)),mayo_AD_PA_counts.agg2)
mayo_AD_PA_counts.pca=prcomp(mayo_AD_PA_counts.agg2[,-1])
scores=as.data.frame(mayo_AD_PA_counts.pca$x)
