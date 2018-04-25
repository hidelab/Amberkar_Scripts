library(pathprint)
library(org.Hs.eg.db)
library(entropy)
Sys.setenv("plotly_username"="ssamberkar")
Sys.setenv("plotly_api_key"="U5r0L12wRGNM98tLdTKz")

setwd("/shared/hidelab2/user/md4zsa/Work/Data/MSMM_Array/MSBB_Array19/Normalised_Data")
mapIds2<-function(IDs,IDFrom,IDTo){
  require(org.Hs.eg.db)
  require(illuminaHumanv3)
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
plq_mean_values=sort(unique(msbb_array19.covariates$PLQ_Mn))
names(plq_mean_values)=c(1:length(plq_mean_values))
msbb_array19_SpTp.SCE=vector(mode = "list",length = length(plq_mean_values))
msbb_array19_BrainRegion.SCE=vector(mode = "list",length = 17)
names(msbb_array19_BrainRegion.SCE)=names(msbb_array19)

msbb_array19_BrainRegion.SCE=lapply(msbb_array19.2.agg2,single.chip.enrichment2,geneset = pathprint.Hs.gs,statistic = "mean",normalizedScore = F)
msbb_array19_BrainRegion.diffPathways=vector(mode = "list",length = 633)

saveRDS(msbb_array19_BrainRegion.SCE,"msbb_array19_BrainRegion_SCE.RDS")

saveRDS(msbb_array19_BrainRegion.diffPathways,"msbb_array19_BrainRegion_diffPathways.RDS")
#length(plq_mean_values)
ns=sort(c(37,38,69,72,75,81,13,83,89))
msbb_array19.covariates3=msbb_array19.covariates2[-which(msbb_array19.covariates2$BrainBank%in%msbb_array19.covariates2$BrainBank[ns]),]
plq_mean_values2=plq_mean_values[-ns]
msbb_array19_SpTp.SCE=vector(mode = "list",length = length(names(plq_mean_values2)))
names(msbb_array19_SpTp.SCE)=names(plq_mean_values2)
for(i in as.numeric(names(msbb_array19_SpTp.SCE_Random))){
  search_samples=paste(msbb_array19.covariates3$BrainBank[which(msbb_array19.covariates3$PLQ_Mn%in%plq_mean_values2[i])],collapse = "|")
  # if(match(i,ns)==T){
  #   cat(paste("Dropping single sample ",search_samples,"...\n",sep = ""))
  #   next
  # }
  if(length(grep(pattern = "X",search_samples))>=1){
    cat(paste("Computing pathway fingerprint for ",length(grep(pattern = search_samples,colnames(msbb_array19_allSamples_Expression)))," samples & for mean plaque ",plq_mean_values[i],"\n",sep = ""))
    msbb_array19_SpTp.SCE_Random[[i]]=single.chip.enrichment2(exprs = msbb_array19_allSamples_Expression[,grep(pattern = search_samples,colnames(msbb_array19_allSamples_Expression))],geneset = pathprint.Hs.gs.random,statistic = "mean",normalizedScore = F)
  }
  
}
# for(i in as.numeric(names(msbb_array19_SpTp.SCE))){
#   search_samples=paste(msbb_array19.covariates3$BrainBank[which(msbb_array19.covariates3$PLQ_Mn%in%plq_mean_values2[i])],collapse = "|")
#   if(length(grep(pattern = "X",search_samples))==0){
#     
#     cat(paste("Dropped sample ",search_samples," at index ",i,"...\n",sep=""))
#     }
# }
saveRDS(msbb_array19_SpTp.SCE,"msbb_array19_SpTp_SCE.RDS")
msbb_array19_SpTp.SCE=readRDS("msbb_array19_SpTp_SCE.RDS")
msbb_array19_SpTp.SCE2=msbb_array19_SpTp.SCE[lapply(msbb_array19_SpTp.SCE,length)>0]
msbb_array19_SpTp.SCE2_Random=msbb_array19_SpTp_SCE.Random[lapply(msbb_array19_SpTp_SCE.Random,length)>0]
#Diffpathways by limma
dep_df1=vector(mode = "list",length = length(msbb_array19_SpTp.SCE2))
dep_df1.random=vector(mode = "list",length = length(msbb_array19_SpTp.SCE2_Random))
dep_df2=vector(mode = "list",length = length(msbb_array19))
names(dep_df2)=names(msbb_array19)
for(t1 in 1:dim(msbb_array19_BrainRegion.diffPathways)[2]){
  for(t2 in 1:dim(msbb_array19_BrainRegion.diffPathways)[2]){
    group=factor(c(rep("R1",length=length(colnames(msbb_array19_BrainRegion.SCE[[t1]]))),(rep("R2",length=length(colnames(msbb_array19_BrainRegion.SCE[[t2]]))))))
    design_df=model.matrix(~0+group)
    colnames(design_df)=levels(group)
    contrasts_matrix=makeContrasts(R1-R2,levels = design_df)
    fit=lmFit(object = cbind(msbb_array19_BrainRegion.SCE[[t1]],msbb_array19_BrainRegion.SCE[[t2]]),design = design_df)
    fit2=contrasts.fit(fit,contrasts_matrix)
    fit2=eBayes(fit2)
    fit2=eBayes(fit2,trend = T)
    dep_df2[[t1]]=topTable(fit2,coef = 1,sort.by = "logFC",adjust.method = "BH",p.value = 0.01,number = 633)
    msbb_array19_BrainRegion.diffPathways[t1,t2]=length(rownames(dep_df[which(dep_df$adj.P.Val<=0.05),]))
    
  }
  
}

msbb_array19_SpTp.diffPathways=matrix(NA,nrow=length(msbb_array19_SpTp.SCE2),ncol=3)
for(t1 in 1:dim(msbb_array19_SpTp.diffPathways)[1]){
    
         group=factor(c(rep("lowPLQ",length=length(colnames(msbb_array19_SpTp.SCE2[[1]]))),(rep("highPLQ",length=length(colnames(msbb_array19_SpTp.SCE2[[t1]]))))))
         design_df=model.matrix(~0+group)
         colnames(design_df)=levels(group)
         contrasts_matrix=makeContrasts(lowPLQ-highPLQ,levels = design_df)
         fit=lmFit(object = cbind(msbb_array19_SpTp.SCE2[[1]],msbb_array19_SpTp.SCE2[[t1]]),design = design_df)
         fit2=contrasts.fit(fit,contrasts_matrix)
         fit2=eBayes(fit2,trend = T)
         dep_df1=topTable(fit2,coef = 1,sort.by = "logFC",adjust.method = "BH",p.value = 0.05,number = 633)
         if(dim(dep_df1)[1]==0){
           msbb_array19_SpTp.diffPathways[t1,]=c(unname(plq_mean_values2[names(msbb_array19_SpTp.SCE2)[t1]]),0,0)
         }
         else
         msbb_array19_SpTp.diffPathways[t1,]=rbind(Mean_PLQ_Value=unname(plq_mean_values2[names(msbb_array19_SpTp.SCE2)[t1]]),
                                                   Nr.DiffPathways=length(rownames(dep_df1[which(dep_df1$adj.P.Val<=0.05),])),
                                                   DiffPathways=paste(rownames(dep_df1[which(dep_df1$adj.P.Val<=0.05),]),collapse = ";"))
           
}
msbb_array19_SpTp.diffPathways=data.frame(msbb_array19_SpTp.diffPathways,stringsAsFactors = F)
colnames(msbb_array19_SpTp.diffPathways)=c("Patient_MeanPLQ","diffPathways","Pathways")
msbb_array19_SpTp.diffPathways$Patient_MeanPLQ=as.numeric(msbb_array19_SpTp.diffPathways$Patient_MeanPLQ)
msbb_array19_SpTp.diffPathways$diffPathways=as.integer(msbb_array19_SpTp.diffPathways$diffPathways)

filt1=grep(pattern = "\\b0\\b",x = msbb_array19_SpTp.diffPathways[,2])
msbb_array19_SpTp.diffPathways2=msbb_array19_SpTp.diffPathways[-filt1,]
msbb_array19_SpTp.diffPathways2_ScatterPlot.loess=loessFit(y = log(msbb_array19_SpTp.diffPathways2$diffPathways),x = msbb_array19_SpTp.diffPathways2$Patient_MeanPLQ,span = 0.8)
ft <- list(family = "sans serif",size = 12,color = 'black')
msbb_array19_SpTp.diffPathways2_ScatterPlot=msbb_array19_SpTp.diffPathways2%>%plot_ly(x=~Patient_MeanPLQ)%>%add_markers(y=~log(diffPathways))%>%add_lines(x=~Patient_MeanPLQ,y = fitted(msbb_array19_SpTp.diffPathways2_ScatterPlot.loess))%>%layout(font=ft,title="#Differential pathways across all patients' mean PLQ",showlegend=F)
plotly_IMAGE(msbb_array19_SpTp.diffPathways2_ScatterPlot,width = 1100,height = 1100,format = "png",out_file = "msbb_array19_SpTp.diffPathways2_ScatterPlot.png")


msbb_array19_SpTp.top_diffPathways_list=vector(mode = "list",length = 5)

for (i in 1:5){
  msbb_array19_SpTp.top_diffPathways_list[[i]]=sort(unlist(strsplit(msbb_array19_SpTp.diffPathways[order(msbb_array19_SpTp.diffPathways$`#diffPathways`,decreasing = T)[1:5][i],3],split = ";")))
}
names(msbb_array19_SpTp.top_diffPathways_list)=msbb_array19_SpTp.diffPathways[order(msbb_array19_SpTp.diffPathways$`#diffPathways`,decreasing = T)[1:5],1]
msbb_array19_SpTp.top_diffPathways_Conserved=Reduce(intersect,msbb_array19_SpTp.top_diffPathways_list)
findSamples=msbb_array19.covariates3$BrainBank[which(msbb_array19.covariates3$PLQ_Mn%in%as.numeric(names(msbb_array19_SpTp.top_diffPathways_list)))]
findSamples=findSamples[-grep(pattern = "956",x =  colnames(SpTp_3D_mat))]
findSamples.list=lapply(msbb_array19_SpTp.SCE2,function(x)grep(pattern = paste("\\b",findSamples,"\\b",collapse = "|",sep = ""),x = colnames(x),value = T))[lapply(lapply(msbb_array19_SpTp.SCE2,function(x)grep(pattern = paste("\\b",findSamples,"\\b",collapse = "|",sep = ""),x = colnames(x),value = T)),length)>0]
SpTp_3D_mat=cbind(msbb_array19_SpTp.SCE2$`12`[which(rownames(msbb_array19_SpTp.SCE2$`12`)%in%msbb_array19_SpTp.top_diffPathways_Conserved),findSamples.list$`12`],
msbb_array19_SpTp.SCE2$`18`[which(rownames(msbb_array19_SpTp.SCE2$`18`)%in%msbb_array19_SpTp.top_diffPathways_Conserved),findSamples.list$`18`],
msbb_array19_SpTp.SCE2$`39`[which(rownames(msbb_array19_SpTp.SCE2$`39`)%in%msbb_array19_SpTp.top_diffPathways_Conserved),findSamples.list$`39`],
msbb_array19_SpTp.SCE2$`60`[which(rownames(msbb_array19_SpTp.SCE2$`60`)%in%msbb_array19_SpTp.top_diffPathways_Conserved),findSamples.list$`60`],
msbb_array19_SpTp.SCE2$`88`[which(rownames(msbb_array19_SpTp.SCE2$`88`)%in%msbb_array19_SpTp.top_diffPathways_Conserved),findSamples.list$`88`])
f1=list(faimly="Arial",size=12,color="black")
ax=list(title="diffPathways",titlefont=f1,showtickLabels=T,ticks="outside",dtick=0.25,ticklen=5,tickwidth=2)
ay=list(title="Brain region",titlefont=f1,showtickLabels=T,ticks="outside",dtick=0.25,ticklen=5,tickwidth=2)
brain_region_vector=unlist(lapply(strsplit(x = colnames(SpTp_3D_mat[,1:61]),split = "\\."),`[[`,1))
axis <- list(
  titlefont = f1,
  tickfont = f1,
  showgrid = T
)
scene = list(
  xaxis = c(axis,title="Conserved diffPathways"),
  yaxis = c(axis,title="Brain regions"),
  zaxis = c(axis,title="Pathway activity"),
  camera = list(eye = list(x = 2, y = 2, z = 2)))
v=plot_ly(x=~rownames(SpTp_3D_mat),y=~brain_region_vector,z=~SpTp_3D_mat,colors = "YlOrRd")%>%add_surface()%>%layout(xaxis=ax,yaxis=ay,scene=scene)
names(colnames(SpTp_3D_mat))=NA
names(colnames(SpTp_3D_mat))[19:49]="12.16"
names(colnames(SpTp_3D_mat))[15:16]="5.230769"
names(colnames(SpTp_3D_mat))[17:18]="8.16"
names(colnames(SpTp_3D_mat))[50:61]="35.28"
names(colnames(SpTp_3D_mat))[1:14]="3.92"
SpTp_3D_mat.anno_df=data.frame(Mean_PLQ=factor(names(colnames(SpTp_3D_mat))))
rownames(SpTp_3D_mat.anno_df)=colnames(SpTp_3D_mat)
pheatmap(SpTp_3D_mat,cluster_cols = T,annotation = SpTp_3D_mat.anno_df)

msbb_array19_SpTp.diffPathways_random=matrix(NA,nrow=60,ncol=60)
for(t1 in 1:dim(msbb_array19_SpTp.diffPathways_random)[2]){
  for(t2 in 1:dim(msbb_array19_SpTp.diffPathways_random)[2]){
    group=factor(c(rep("lowPLQ",length=length(colnames(msbb_array19_SpTp.SCE2_Random[[t1]]))),(rep("highPLQ",length=length(colnames(msbb_array19_SpTp.SCE2_Random[[t2]]))))))
    design_df=model.matrix(~0+group)
    colnames(design_df)=levels(group)
    contrasts_matrix=makeContrasts(lowPLQ-highPLQ,levels = design_df)
    fit=lmFit(object = cbind(msbb_array19_SpTp.SCE2_Random[[t1]],msbb_array19_SpTp.SCE2_Random[[t2]]),design = design_df)
    fit2=contrasts.fit(fit,contrasts_matrix)
    fit2=eBayes(fit2)
    fit2=eBayes(fit2,trend = T)
    dep_df1.random[[t1]]=topTable(fit2,coef = 1,sort.by = "logFC",adjust.method = "BH",p.value = 0.01,number = 633)
    msbb_array19_SpTp.diffPathways_random[t1,t2]=length(rownames(dep_df1.random[[t1]][which(dep_df1.random[[t1]]$adj.P.Val<=0.01),]))
    
  }
  
}

dep_df1=readRDS("msbb_array19_SpTp_dep_df1.RDS")
#Comparison of diff. pathways ~ mean plaque values at various top thresholds
dep_df1_common_pathways_range=seq(15,150,by = 5)
dep_df1_common_pathways_rangeMatrix=matrix(NA,nrow=length(dep_df1_common_pathways_range),ncol = 3)
for(i in 1:length(dep_df1_common_pathways_range)){
  dep_df1_common_pathways=Reduce(intersect,lapply(dep_df1[lapply(dep_df1,function(x)dim(x)[1])>dep_df1_common_pathways_range[i]],rownames))
  dep_df1_common_pathways_indices=names(lapply(dep_df1,function(x)which(rownames(x)%in%dep_df1_common_pathways))[lapply(lapply(dep_df1,function(x)which(rownames(x)%in%dep_df1_common_pathways)),length)==length(dep_df1_common_pathways)])
  dep_df1_common_pathways_rangeMatrix[i,]=c(dep_df1_common_pathways_range[i],length(Reduce(intersect,lapply(dep_df1[lapply(dep_df1,function(x)dim(x)[1])>dep_df1_common_pathways_range[i]],rownames))),length(dep_df1_common_pathways_indices))
}
rownames(dep_df1_common_pathways_rangeMatrix)=dep_df1_common_pathways_rangeMatrix[,1]
dep_df1_common_pathways_rangeMatrix=data.frame(dep_df1_common_pathways_rangeMatrix[,-1],stringsAsFactors = F)
#colnames(dep_df1_common_pathways_rangeMatrix)=c("#Common differentially expressed pathways","#Mean PLQ values")
write.table(dep_df1_common_pathways_rangeMatrix,col.names = T,"msbb_array19_SpTp_Common_diffPathways_RangeMatrix.txt",sep = "\t",row.names = T)

dep_df1_common_pathways_expr=data.frame(t(do.call("cbind",lapply(dep_df1[dep_df1_common_pathways_indices],function(x)x[which(rownames(x)%in%dep_df1_common_pathways),2]))))
dep_df1_common_pathways_expr$CommonPathways=rep(dep_df1_common_pathways,length=38)
f <- list(family = "Courier New, monospace",size = 18,color = "#7f7f7f")
x <- list(title = "#Top differentially active pathways",titlefont = f)
p=plot_ly(data = dep_df1_common_pathways_rangeMatrix,x=~as.numeric(rownames(dep_df1_common_pathways_rangeMatrix)),y=~X1,name="Nr.of.Common.diffPathways",mode = "lines",type="scatter")%>%add_trace(y=~X2,name='Mean-PLQ values(Patients)')%>%layout(xaxis=x)
export(p,file = "Common_diffPathways_Range.png")


msbb_array19_SpTp.diffPathways=read.table("msbb_array19_SpTp_diffPathwaysMatrix.txt",sep = "\t",header = T,as.is = T)
msbb_array19_SpTp.diffPathways_p001=read.table("msbb_array19_SpTp_diffPathwaysMatrix_p001.txt",sep = "\t",header = T,as.is = T)
msbb_array19_SpTp.diffPathways=msbb_array19_SpTp.diffPathways/633
msbb_array19_SpTp.diffPathways_p001=msbb_array19_SpTp.diffPathways_p001/633
colnames(msbb_array19_SpTp.diffPathways)=plq_mean_values[as.numeric(names(msbb_array19_SpTp.SCE2))]
colnames(msbb_array19_SpTp.diffPathways_p001)=plq_mean_values[as.numeric(names(msbb_array19_SpTp.SCE2))]
rownames(msbb_array19_SpTp.diffPathways)=plq_mean_values[as.numeric(names(msbb_array19_SpTp.SCE2))]
rownames(msbb_array19_SpTp.diffPathways_p001)=plq_mean_values[as.numeric(names(msbb_array19_SpTp.SCE2))]
msbb_array19_SpTp_diffPathways.anno_df=data.frame(PLQ_Mean=plq_mean_values[as.numeric(names(msbb_array19_SpTp.SCE2))],stringsAsFactors = F)
rownames(msbb_array19_SpTp_diffPathways.anno_df)=colnames(msbb_array19_SpTp.diffPathways)

pheatmap(msbb_array19_SpTp.diffPathways,annotation_col = msbb_array19_SpTp_diffPathways.anno_df,annotation_row = msbb_array19_SpTp_diffPathways.anno_df,show_rownames = F,show_colnames = F,cellwidth = 10,cellheight = 10,main="Temporal (mean plaque) pathway landscape - pairwise differential pathway activation")

msbb_array19_BrainRegion.diffPathways=matrix(NA,nrow=17,ncol=17)
for(t1 in 1:dim(msbb_array19_BrainRegion.diffPathways)[2]){
  for(t2 in 1:dim(msbb_array19_BrainRegion.diffPathways)[2]){
    group=factor(c(rep("R1",length=length(colnames(msbb_array19_BrainRegion.SCE[[t1]]))),(rep("R2",length=length(colnames(msbb_array19_BrainRegion.SCE[[t2]]))))))
    design_df=model.matrix(~0+group)
    colnames(design_df)=levels(group)
    contrasts_matrix=makeContrasts(R1-R2,levels = design_df)
    fit=lmFit(object = cbind(msbb_array19_BrainRegion.SCE[[t1]],msbb_array19_BrainRegion.SCE[[t2]]),design = design_df)
    fit2=contrasts.fit(fit,contrasts_matrix)
    fit2=eBayes(fit2)
    fit2=eBayes(fit2,trend = T)
    dep_df=topTable(fit2,coef = 1,sort.by = "logFC",adjust.method = "BH",p.value = 0.01,number = 633)
    msbb_array19_BrainRegion.diffPathways[t1,t2]=length(rownames(dep_df[which(dep_df$adj.P.Val<=0.01),]))
    
  }
  
}
msbb_array19_BrainRegion.diffPathways=msbb_array19_BrainRegion.diffPathways/633
colnames(msbb_array19_BrainRegion.diffPathways)=rownames(msbb_array19_BrainRegion.diffPathways)=names(msbb_array19)
pheatmap(msbb_array19_BrainRegion.diffPathways,main="Spatial (brain region) pathway landscape - pairwise differential pathway activation")


#


#DEG genes, based on mean plq
msbb_array19.eset=msbb_array19.DEG_PLQ=vector(mode = "list",length = 4)
names(msbb_array19.DEG_PLQ)=c("S0","S1","S2","S3")
msbb_array19.DEG_PLQ$S0=msbb_array19.DEG_PLQ$S1=msbb_array19.DEG_PLQ$S2=msbb_array19.DEG_PLQ$S3=vector(mode = "list",length = 17)
msbb_array19.eset$S0=msbb_array19.eset$S1=msbb_array19.eset$S2=msbb_array19.eset$S3=vector(mode = "list",length = 17)
control.plq=lapply(msbb_array19.NTr_PLQ_phenoData2,function(x)which(x$PLQ_Mn==0))
disease.plq=lapply(msbb_array19.NTr_PLQ_phenoData2,function(x)which(x$PLQ_Mn>=21))
# control.ntr=lapply(msbb_array19.NTr_PLQ_phenoData2,function(x)which(x$NTrSum<5))
# disease.ntr=lapply(msbb_array19.NTr_PLQ_phenoData2,function(x)which(x$NTrSum>10))
for (i in 1:length(names(msbb_array19.NTr_PLQ_phenoData2))){
  design.df=matrix(NA,nrow=sum(length(control.plq[[i]]),length(disease.plq[[i]])),ncol=2)
  rownames(design.df)=c(rownames(msbb_array19.NTr_PLQ_phenoData2[[i]][control.plq[[i]],]),rownames(msbb_array19.NTr_PLQ_phenoData2[[i]][disease.plq[[i]],]))
  colnames(design.df)=c("control.plq","disease.plq")
  design.df[which(rownames(design.df)%in%rownames(msbb_array19.NTr_PLQ_phenoData2[[i]][control.plq[[i]],])),1]=1
  design.df[-which(rownames(design.df)%in%rownames(msbb_array19.NTr_PLQ_phenoData2[[i]][control.plq[[i]],])),1]=0
  design.df[which(rownames(design.df)%in%rownames(msbb_array19.NTr_PLQ_phenoData2[[i]][disease.plq[[i]],])),2]=1
  design.df[-which(rownames(design.df)%in%rownames(msbb_array19.NTr_PLQ_phenoData2[[i]][disease.plq[[i]],])),2]=0
  phenoData=AnnotatedDataFrame(data=msbb_array19.NTr_PLQ_phenoData2[[i]][which(rownames(msbb_array19.NTr_PLQ_phenoData2[[i]])%in%rownames(design.df)),])
  msbb_array19.eset$S3[[i]]=ExpressionSet(assayData=as.matrix(msbb_array19.agg2[[i]][,which(colnames(msbb_array19.agg2[[i]])%in%rownames(design.df))]),phenoData=phenoData)  
  fit=lmFit(msbb_array19.eset$S3[[i]],design=design.df)
  fit=eBayes(fit)
  contMatrix=makeContrasts(CtrlvsDisease=disease.plq-control.plq,levels=design.df)
  fit2=contrasts.fit(fit,contMatrix)
  fit2=eBayes(fit2)
  msbb_array19.DEG_PLQ$S3[[i]]=topTableF(fit2,adjust.method="BH",number=Inf)
}


names(msbb_array19.eset)=names(msbb_array19.DEG_NTR)=names(msbb_array19.DEG_PLQ)=names(msbb_array19)

msbb_array19.SCE2=mapply(single.chip.enrichment2,msbb_array19.agg,MoreArgs = list(geneset = pathprint.Hs.gs,"mean",T))

msbb_array19.2_PLQ_Strat=msbb_array19.fingerprint_PLQ_Strat=msbb_array19_fingerprint.ZeroSumPathways_PLQ_Strat=msbb_array19_fingerprint.noZeromSumPathways_PLQ_Strat=msbb_array19_fingerprint.HighEntropy05_Indices_PLQ_Strat=msbb_array19_fingerprint.HighEntropy05_Fingerprint_PLQ_Strat=msbb_array19_fingerprint_HighEntropy05_Fingerprint.Core_PLQ_Strat=msbb_array19_fingerprint_HighEntropy05.CoreFingerprint_PLQ_Strat=vector(mode = "list",length = 3)
names(msbb_array19.2_PLQ_Strat)=names(msbb_array19.fingerprint_PLQ_Strat)=names(msbb_array19_fingerprint.ZeroSumPathways_PLQ_Strat)=names(msbb_array19_fingerprint.noZeromSumPathways_PLQ_Strat)=names(msbb_array19_fingerprint.HighEntropy05_Indices_PLQ_Strat)=names(msbb_array19_fingerprint.HighEntropy05_Fingerprint_PLQ_Strat)=names(msbb_array19_fingerprint_HighEntropy05_Fingerprint.Core_PLQ_Strat)=names(msbb_array19_fingerprint_HighEntropy05.CoreFingerprint_PLQ_Strat)=c(1:3)
msbb_array19.2_PLQ_Strat_SCE=msbb_array19.fingerprint_PLQ_Strat_SCE=msbb_array19_fingerprint.ZeroSumPathways_PLQ_Strat_SCE=msbb_array19_fingerprint.noZeromSumPathways_PLQ_Strat_SCE=msbb_array19_fingerprint.HighEntropy05_Indices_PLQ_Strat_SCE=msbb_array19_fingerprint.HighEntropy05_Fingerprint_PLQ_Strat_SCE=msbb_array19_fingerprint_HighEntropy05_Fingerprint.Core_PLQ_Strat_SCE=msbb_array19_fingerprint_HighEntropy05.CoreFingerprint_PLQ_Strat_SCE=vector(mode = "list",length = 3)
names(msbb_array19.2_PLQ_Strat_SCE)=names(msbb_array19.fingerprint_PLQ_Strat_SCE)=names(msbb_array19_fingerprint.ZeroSumPathways_PLQ_Strat_SCE)=names(msbb_array19_fingerprint.noZeromSumPathways_PLQ_Strat_SCE)=names(msbb_array19_fingerprint.HighEntropy05_Indices_PLQ_Strat_SCE)=names(msbb_array19_fingerprint.HighEntropy05_Fingerprint_PLQ_Strat_SCE)=names(msbb_array19_fingerprint_HighEntropy05_Fingerprint.Core_PLQ_Strat_SCE)=names(msbb_array19_fingerprint_HighEntropy05.CoreFingerprint_PLQ_Strat_SCE)=c(1:3)

msbb_array19.2_PLQ_Strat=vector(mode = "list",length = 4)
names(msbb_array19.2_PLQ_Strat)=c("S0","S1","S2","S3")
msbb_array19.2_PLQ_Strat[[1]]=msbb_array19.covariates[which(msbb_array19.covariates$PLQ_Mn==0),1]
msbb_array19.2_PLQ_Strat[[2]]=msbb_array19.covariates[which(msbb_array19.covariates$PLQ_Mn>=1&msbb_array19.covariates$PLQ_Mn<=10),1]
msbb_array19.2_PLQ_Strat[[3]]=msbb_array19.covariates[which(msbb_array19.covariates$PLQ_Mn>=11&msbb_array19.covariates$PLQ_Mn<=20),1]
msbb_array19.2_PLQ_Strat[[4]]=msbb_array19.covariates[which(msbb_array19.covariates$PLQ_Mn>=21),1]

names(msbb_array19.fingerprint_PLQ_Strat)=names(msbb_array19.2_PLQ_Strat)
msbb_array19.fingerprint_PLQ_Strat$S0=single.chip.enrichment2(exprs = msbb_array19_allSamples_Expression[,grep(pattern = paste("\\b",msbb_array19.2_PLQ_Strat$S0,"\\b",sep = "",collapse = "|"),x = colnames(msbb_array19_allSamples_Expression))],geneset = pathprint.Hs.gs,statistic = "mean",progressBar = T)
msbb_array19.fingerprint_PLQ_Strat$S1=single.chip.enrichment2(exprs = msbb_array19_allSamples_Expression[,grep(pattern = paste("\\b",msbb_array19.2_PLQ_Strat$S1,"\\b",sep = "",collapse = "|"),x = colnames(msbb_array19_allSamples_Expression))],geneset = pathprint.Hs.gs,statistic = "mean",progressBar = T)
msbb_array19.fingerprint_PLQ_Strat$S2=single.chip.enrichment2(exprs = msbb_array19_allSamples_Expression[,grep(pattern = paste("\\b",msbb_array19.2_PLQ_Strat$S2,"\\b",sep = "",collapse = "|"),x = colnames(msbb_array19_allSamples_Expression))],geneset = pathprint.Hs.gs,statistic = "mean",progressBar = T)
msbb_array19.fingerprint_PLQ_Strat$S3=single.chip.enrichment2(exprs = msbb_array19_allSamples_Expression[,grep(pattern = paste("\\b",msbb_array19.2_PLQ_Strat$S3,"\\b",sep = "",collapse = "|"),x = colnames(msbb_array19_allSamples_Expression))],geneset = pathprint.Hs.gs,statistic = "mean",progressBar = T)

msbb_array19_PLQ_Strat.diffPathways=matrix(NA,ncol=3,nrow=3)
for(i in 2:4){
  group=factor(c(rep("S0",length=length(colnames(msbb_array19.fingerprint_PLQ_Strat$S0))),(rep("S3",length=length(colnames(msbb_array19.fingerprint_PLQ_Strat$S3))))))
  design_df=model.matrix(~0+group)
  colnames(design_df)=levels(group)
  contrasts_matrix=makeContrasts(S0-S3,levels = design_df)
  fit=lmFit(object = cbind(msbb_array19.fingerprint_PLQ_Strat$S0,msbb_array19.fingerprint_PLQ_Strat$S3),design = design_df)
  fit2=contrasts.fit(fit,contrasts_matrix)
  fit2=eBayes(fit2,trend = T)
  dep_df1=topTable(fit2,coef = 1,sort.by = "logFC",adjust.method = "BH",p.value = 0.05,number = 633)
  msbb_array19_PLQ_Strat.diffPathways[3,]=rbind(PLQ_Strata="S3",
                                              Nr.DiffPathways=length(rownames(dep_df1[which(dep_df1$adj.P.Val<=0.05),])),
                                              DiffPathways=paste(rownames(dep_df1[which(dep_df1$adj.P.Val<=0.05),]),collapse = ";"))
  
}


diffPathways_S0_Core=msbb_array19.fingerprint_PLQ_Strat$S0[which(rownames(msbb_array19.fingerprint_PLQ_Strat$S0)%in%msbb_array19_PLQ_Strat.diffPathways_Core),]
diffPathways_S1_Core=msbb_array19.fingerprint_PLQ_Strat$S1[which(rownames(msbb_array19.fingerprint_PLQ_Strat$S1)%in%msbb_array19_PLQ_Strat.diffPathways_Core),]
diffPathways_S2_Core=msbb_array19.fingerprint_PLQ_Strat$S2[which(rownames(msbb_array19.fingerprint_PLQ_Strat$S2)%in%msbb_array19_PLQ_Strat.diffPathways_Core),]
diffPathways_S3_Core=msbb_array19.fingerprint_PLQ_Strat$S3[which(rownames(msbb_array19.fingerprint_PLQ_Strat$S3)%in%msbb_array19_PLQ_Strat.diffPathways_Core),]
diffPathways_StratCore.df=cbind(diffPathways_S0_Core,diffPathways_S1_Core,diffPathways_S2_Core,diffPathways_S3_Core)
names(colnames(diffPathways_StratCore.df))=c(rep("S0",164),rep("S1",503),rep("S2",168),rep("S3",61))
diffPathways_StratCore.anno_df=data.frame(Strata=factor(names(colnames(diffPathways_StratCore.df))),Brain_region=factor(colnames(diffPathways_StratCore.df)))
rownames(diffPathways_StratCore.anno_df)=colnames(diffPathways_StratCore.df)
#colnames(diffPathways_StratCore.df)=unlist(lapply(strsplit(colnames(diffPathways_StratCore.df),split = "\\."),`[[`,1))
diffPathways_StratCore.anno_df$Brain_region=unlist(lapply(strsplit(colnames(diffPathways_StratCore.df),split = "\\."),`[[`,1))
pheatmap(diffPathways_StratCore.df,annotation = diffPathways_StratCore.anno_df,cluster_cols = F,labels_col = "",main = "Conserved diffPathways across Strata",fontsize = 10,breaks = c(3,3.5,4),color = c("blue","white","red"))

fp_plq_deg_apoptosis_expr.df=rbind(as.data.frame(t(msbb_array19.agg2$FP[which(msbb_array19.agg2$FP$Gene.Symbol%in%fp_plq_deg_apoptosis_enriched_genes),lapply(msbb_array19.NTr_PLQ_phenoData2,function(x)which(x$PLQ_Mn>=1&x$PLQ_Mn<=10))$`FP`][,-1])),as.data.frame(t(x = msbb_array19.agg2$FP[which(msbb_array19.agg2$FP$Gene.Symbol%in%fp_plq_deg_apoptosis_enriched_genes),control.plq$FP])))
colnames(fp_plq_deg_apoptosis_expr.df)=msbb_array19.agg2$FP[which(msbb_array19.agg2$FP$Gene.Symbol%in%fp_plq_deg_apoptosis_enriched_genes),lapply(msbb_array19.NTr_PLQ_phenoData2,function(x)which(x$PLQ_Mn>=1&x$PLQ_Mn<=10))$`FP`][,1]
boxplot(fp_plq_deg_apoptosis_expr.df[1:12,-8][,1],fp_plq_deg_apoptosis_expr.df[13:43,-8][,1],
        fp_plq_deg_apoptosis_expr.df[1:12,-8][,2],fp_plq_deg_apoptosis_expr.df[13:43,-8][,2],
        fp_plq_deg_apoptosis_expr.df[1:12,-8][,3],fp_plq_deg_apoptosis_expr.df[13:43,-8][,3],
        fp_plq_deg_apoptosis_expr.df[1:12,-8][,4],fp_plq_deg_apoptosis_expr.df[13:43,-8][,4],
        fp_plq_deg_apoptosis_expr.df[1:12,-8][,5],fp_plq_deg_apoptosis_expr.df[13:43,-8][,5],
        fp_plq_deg_apoptosis_expr.df[1:12,-8][,6],fp_plq_deg_apoptosis_expr.df[13:43,-8][,6],
        fp_plq_deg_apoptosis_expr.df[1:12,-8][,7],fp_plq_deg_apoptosis_expr.df[13:43,-8][,7],
        col=c("green1","red1"),
        names=c(paste(colnames(fp_plq_deg_apoptosis_expr.df[1:12,-8])[1]),paste(colnames(fp_plq_deg_apoptosis_expr.df[1:12,-8])[1]),
                paste(colnames(fp_plq_deg_apoptosis_expr.df[1:12,-8])[2]),paste(colnames(fp_plq_deg_apoptosis_expr.df[1:12,-8])[2]),
                paste(colnames(fp_plq_deg_apoptosis_expr.df[1:12,-8])[3]),paste(colnames(fp_plq_deg_apoptosis_expr.df[1:12,-8])[3]),
                paste(colnames(fp_plq_deg_apoptosis_expr.df[1:12,-8])[4]),paste(colnames(fp_plq_deg_apoptosis_expr.df[1:12,-8])[4]),
                paste(colnames(fp_plq_deg_apoptosis_expr.df[1:12,-8])[5]),paste(colnames(fp_plq_deg_apoptosis_expr.df[1:12,-8])[5]),
                paste(colnames(fp_plq_deg_apoptosis_expr.df[1:12,-8])[6]),paste(colnames(fp_plq_deg_apoptosis_expr.df[1:12,-8])[6]),
                paste(colnames(fp_plq_deg_apoptosis_expr.df[1:12,-8])[7]),paste(colnames(fp_plq_deg_apoptosis_expr.df[1:12,-8])[7])),
        main="Frontal Pole, S0~S1 DEGs enriched in KEGG - Apoptosis")
legend(legend = c("0 PLQ","S1 PLQ"),fill = c("green1","red1"),x = 12.5,y = 6.5)

msbb_array19_BrainRegion.diffPathways=vector(mode = "list",length = 3)
names(msbb_array19_BrainRegion.diffPathways)=c("S1","S2","S3")
msbb_array19_BrainRegion.diffPathways$S1=msbb_array19_BrainRegion.diffPathways$S2=msbb_array19_BrainRegion.diffPathways$S3=vector(mode = "list",length = 17)
names(msbb_array19_BrainRegion.diffPathways$S1)=names(msbb_array19_BrainRegion.diffPathways$S2)=names(msbb_array19_BrainRegion.diffPathways$S3)=names(msbb_array19)
plq_strat_values=list(1:10,11:20,21:max(msbb_array19.covariates3$PLQ_Mn))
for (s in 1:3){
  for(i in c(1:8,10:17)){
    brain_region_SCE=data.frame(msbb_array19.SCE2[[i]],stringsAsFactors = F)
    brain_region_SCE[634,]=msbb_array19.covariates$PLQ_Mn[which(msbb_array19.covariates$BrainBank%in%colnames(brain_region_SCE))]
    brain_region_SCE.Control=brain_region_SCE[-634,which(brain_region_SCE[634,]==0)]
    brain_region_SCE.AD_PLQ=brain_region_SCE[-634,which(brain_region_SCE[634,]<=max(plq_strat_values[[s]])&brain_region_SCE[634,]>=min(plq_strat_values[[s]]))]
    group=factor(c(rep("lowPLQ",length(brain_region_SCE.Control)),rep("highPLQ",length(brain_region_SCE.AD_PLQ))))
    design_df=model.matrix(~0+group)
    colnames(design_df)=levels(group)
    contrasts_matrix=makeContrasts(lowPLQ-highPLQ,levels = design_df)
    fit=lmFit(object = cbind(brain_region_SCE.Control,brain_region_SCE.AD_PLQ),design = design_df)
    fit2=contrasts.fit(fit,contrasts_matrix)
    fit2=eBayes(fit2,trend = T)
    msbb_array19_BrainRegion.diffPathways[[s]][[i]]=sort(rownames(topTable(fit2,coef = 1,sort.by = "logFC",adjust.method = "BH",p.value = 0.20,number = 633)))
    
  }
  msbb_array19_BrainRegion.diffPathways[[s]]=msbb_array19_BrainRegion.diffPathways[[s]][lapply(msbb_array19_BrainRegion.diffPathways[[s]],length)>0]
}

gse29378.SCE=single.chip.enrichment2(exprs = log(gse29378_matrix.agg),geneset = pathprint.Hs.gs,statistic = "mean",progressBar = T)
gse29378_plq=sort(unique(gse29378_metadata_covariate.merged$Plaques))
gse29378.diffPathways=vector(mode = "list",length = 3)
names(gse29378.diffPathways)=c("S1","S2","S3")
for(i in 2:4){
  gse29378.Control=gse29378_metadata_covariate.merged$GSM_ID[gse29378_metadata_covariate.merged$Plaques==0][is.na(gse29378_metadata_covariate.merged$GSM_ID[gse29378_metadata_covariate.merged$Plaques==0])!=T]
  gse29378.AD=gse29378_metadata_covariate.merged$GSM_ID[gse29378_metadata_covariate.merged$Plaques==gse29378_plq[i]][is.na(gse29378_metadata_covariate.merged$GSM_ID[gse29378_metadata_covariate.merged$Plaques==gse29378_plq[i]])!=T]
  group=factor(c(rep("lowPLQ",length(gse29378.Control)),rep("highPLQ",length(gse29378.AD))))
  design_df=model.matrix(~0+group)
  colnames(design_df)=levels(group)
  contrasts_matrix=makeContrasts(lowPLQ-highPLQ,levels = design_df)
  fit=lmFit(object = cbind(gse29378.SCE[,gse29378.Control],gse29378.SCE[,gse29378.AD]),design = design_df)
  fit2=contrasts.fit(fit,contrasts_matrix)
  fit2=eBayes(fit2,trend = T)
  gse29378.diffPathways[[i-1]]=rownames(topTable(fit2,coef = 1,sort.by = "logFC",adjust.method = "BH",p.value = 0.05,number = 633))
}
gse29378.diffPathways$S1=gse29378.diffPathways$S2=gse29378.diffPathways$S3=vector(mode = "list",length = 2)
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
lobe_bm_area.map[15,]=c("Frontal_Lobe","PFC")
lobe_bm_area.map[16,]=c("Temporal_Lobe","HP")






