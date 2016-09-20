library(igraph)
library(rsgcc)
library(org.Hs.eg.db)
library(WGCNA)
library(DESeq2)
library(cocor)
setwd("/Users/sandeepamberkar/Work/Data/MSMM-MSBB-HBTRC-Data/MSMM_RNAseq/MSMM_RNAseq_FinalRelease/")
cat(paste("Reading data ...\n"))

msmm_rnaseq_data=read.table("AMP-AD_MSBB_MSSM_IlluminaHiSeq2500_final.TMM_normalized.race_age_RIN_PMI_batch_corrected.tsv",sep = "\t",header = T,as.is = T)
msmm_rnaseq_raw_data=read.table("AMP-AD_MSBB_MSSM_IlluminaHiSeq2500_raw_counts_final.tsv",sep = "\t",header = T,as.is = T)
msmm_rnaseq_clinical=read.table("MSBB_clinical.csv",header = T,sep = ",",as.is = T)
msmm_rnaseq_metaSample=read.table("AMP-AD_MSBB_MSSM_IlluminaHiSeq2500.sample.meta.tsv",header = T,sep = "\t",as.is = T)
msmm_rnaseq_covariates=read.table("AMP-AD_MSBB_MSSM_meta.traits.tsv",header = T,sep = "\t",as.is = T)
msmm_rnaseq_covariates2=msmm_rnaseq_covariates[which(msmm_rnaseq_covariates$Sample.ID%in%msmm_rnaseq_metaSample$LibID[grep(pattern = "^BM_",msmm_rnaseq_metaSample$LibID)]),]
msmm_rnaseq_covariates2$AOD[which(msmm_rnaseq_covariates2$AOD=="90+")]=90
msmm_rnaseq_covariates2=msmm_rnaseq_covariates2[-which(is.na(msmm_rnaseq_covariates2$bbscore)==T),]

cat(paste("Averaging mean expression for multiple Ensembl probeIDs to Gene Symbols ...\n"))
msmm_rnaseq_data2=msmm_rnaseq_data[,c(1,7,which(colnames(msmm_rnaseq_data)%in%msmm_rnaseq_metaSample$NewBarcode[grep(pattern = "^BM_",msmm_rnaseq_metaSample$LibID)]))]
msmm_rnaseq_raw_data2=msmm_rnaseq_raw_data[,c(1,7,which(colnames(msmm_rnaseq_data)%in%msmm_rnaseq_metaSample$NewBarcode[grep(pattern = "^BM_",msmm_rnaseq_metaSample$LibID)]))]
colnames(msmm_rnaseq_data2)[-c(1:2)]=msmm_rnaseq_metaSample$LibID[1:469]
colnames(msmm_rnaseq_raw_data2)[-c(1:2)]=msmm_rnaseq_metaSample$LibID[1:469]
msmm_rnaseq_data2.agg=aggregate(x = msmm_rnaseq_data2[,-c(1:2)],by=list(geneSymbol=msmm_rnaseq_data2$geneSymbol),mean)
msmm_rnaseq_raw_data2.agg=aggregate(x = msmm_rnaseq_raw_data2[,-c(1:2)],by=list(geneSymbol=msmm_rnaseq_raw_data2$geneSymbol),mean)

msmm_rnaseq=msmm_rnaseq_alog=msmm_rnaseq_raw=vector(mode = "list",length = 3)
names(msmm_rnaseq)=names(msmm_rnaseq_alog)=names(msmm_rnaseq_raw)=c("BM_10","BM_22","BM_36")
for (l in 1:length(names(msmm_rnaseq))){
  msmm_rnaseq[[l]]=msmm_rnaseq_data2.agg[,c(1,grep(pattern = names(msmm_rnaseq)[l],colnames(msmm_rnaseq_data2.agg)))]  
  rownames(msmm_rnaseq[[l]])=msmm_rnaseq[[l]]$geneSymbol
  msmm_rnaseq[[l]]=msmm_rnaseq[[l]][,-1]
}
for (l in 1:length(names(msmm_rnaseq_alog))){
  msmm_rnaseq_alog[[l]]=data.frame(lapply(msmm_rnaseq[[l]],function(x)2^x))
  rownames(msmm_rnaseq_alog[[l]])=rownames(msmm_rnaseq[[l]])
}
for (l in 1:length(names(msmm_rnaseq))){
  msmm_rnaseq_raw[[l]]=msmm_rnaseq_raw_data2.agg[,c(1,grep(pattern = names(msmm_rnaseq)[l],colnames(msmm_rnaseq_raw_data2.agg)))]  
  rownames(msmm_rnaseq_raw[[l]])=msmm_rnaseq_raw[[l]]$geneSymbol
  msmm_rnaseq_raw[[l]]=msmm_rnaseq_raw[[l]][,-1]
  msmm_rnaseq_raw[[l]]=apply(msmm_rnaseq_raw[[l]],2,round,digits=0)
}

msmm_rnaseq_Plaque.Corr=vector(mode = "list",length = 3)
names(msmm_rnaseq_Plaque.Corr)=names(msmm_rnaseq)
# 
cat(paste("Grouping samples by brain region ...\n"))
low_Plaque_Samples=high_Plaque_Samples=vector(mode="list",length = 3)
names(low_Plaque_Samples)=names(high_Plaque_Samples)=names(msmm_rnaseq)

low_Plaque_Samples[[1]]=grep(pattern = "^BM_10",msmm_rnaseq_covariates2$Sample.ID[which(msmm_rnaseq_covariates2$PlaqueMean<=1)],value = T)
low_Plaque_Samples[[2]]=grep(pattern = "^BM_22",msmm_rnaseq_covariates2$Sample.ID[which(msmm_rnaseq_covariates2$PlaqueMean<=1)],value = T)
low_Plaque_Samples[[3]]=grep(pattern = "^BM_36",msmm_rnaseq_covariates2$Sample.ID[which(msmm_rnaseq_covariates2$PlaqueMean<=1)],value = T)

high_Plaque_Samples[[1]]=grep(pattern = "^BM_10",msmm_rnaseq_covariates2$Sample.ID[which(msmm_rnaseq_covariates2$PlaqueMean>=15)],value = T)
high_Plaque_Samples[[2]]=grep(pattern = "^BM_22",msmm_rnaseq_covariates2$Sample.ID[which(msmm_rnaseq_covariates2$PlaqueMean>=15)],value = T)
high_Plaque_Samples[[3]]=grep(pattern = "^BM_36",msmm_rnaseq_covariates2$Sample.ID[which(msmm_rnaseq_covariates2$PlaqueMean>=15)],value = T)

cat(paste(rep("#",100),collapse = ""),"\n")
cat(paste("Beginning correlation computations for Plaque pathology ...\n"))
cat(paste(rep("#",100),collapse = ""),"\n")
tmp1=matrix(NA,nrow=dim(msmm_rnaseq[[1]])[1],ncol=3)
#tmp4=matrix(NA,nrow=dim(msmm_rnaseq[[1]])[1],ncol=3)
for (i in 1:length(names(msmm_rnaseq))){
  
  for (j in 1:dim(msmm_rnaseq[[i]])[1]){
    cat(paste("Processing gene - ",rownames(msmm_rnaseq[[i]])[j],"for brain region",names(msmm_rnaseq)[i],"& PLQ \n",sep=" "))
    x1=unname(unlist(msmm_rnaseq[[i]][j,which(colnames(msmm_rnaseq[[i]])%in%msmm_rnaseq_covariates2$Sample.ID)]))
    y1=msmm_rnaseq_covariates2$PlaqueMean[which(msmm_rnaseq_covariates2$Sample.ID%in%colnames(msmm_rnaseq[[i]]))]
    rho=cor.test(x1,y1,method = "spearman",alternative = "two.sided",exact = T)
    rho.sp=unname(rho$estimate)
    rho.p=rho$p.value
    genes=rownames(msmm_rnaseq[[i]][,which(colnames(msmm_rnaseq[[i]])%in%msmm_rnaseq_covariates2$Sample.ID)])[j]
    tmp1[j,1:3]=cbind(genes,rho.sp,rho.p)
    msmm_rnaseq_Plaque.Corr[[i]]=data.frame(Genes=tmp1[,1],Rho.Spearman=as.numeric(tmp1[,2]),Pval=as.numeric(tmp1[,3]),stringsAsFactors = F)
    
    # cat(paste("Processing gene - ",rownames(msmm_rnaseq[[i]][,high_Plaque_Samples[[i]]])[j],"for brain region",names(msmm_rnaseq)[i],"& PLQ \n",sep=" "))
    # x2=unname(unlist(msmm_rnaseq[[i]][j,which(colnames(msmm_rnaseq[[i]][,high_Plaque_Samples[[i]]])%in%msmm_rnaseq_covariates2$Sample.ID)]))
    # y2=msmm_rnaseq_covariates2$PlaqueMean[which(msmm_rnaseq_covariates2$Sample.ID%in%colnames(msmm_rnaseq[[i]][,high_Plaque_Samples[[i]]]))]
    # rho.sp.high=unname(cor.test(x2,y2,method = "spearman")$estimate)
    # rho.p.high=cor.test(x2,y2)$p.value
    # genes.high=rownames(msmm_rnaseq[[i]][,which(colnames(msmm_rnaseq[[i]][,high_Plaque_Samples[[i]]])%in%msmm_rnaseq_covariates2$Sample.ID)])[j]
    # tmp4[j,1:3]=cbind(genes.high,rho.sp.high,rho.p.high)
    # msmm_rnaseq_Plaque.Corr[[i]][[2]]=data.frame(Genes=tmp4[,1],Rho.Spearman=as.numeric(tmp4[,2]),Pval=as.numeric(tmp4[,3]),stringsAsFactors = F)
  }
}

cat(paste("Filtering genes at pval < 5% ...\n"))
names(msmm_rnaseq_Plaque.Corr)=names(msmm_rnaseq)   
msmm_rnaseq_Plaque.Corr005.Genes=lapply(msmm_rnaseq_Plaque.Corr,function(x)x[which(x$Pval<=0.01),1])
# msmm_rnaseq_Plaque.Corr005.LowGenes=lapply(lapply(msmm_rnaseq_Plaque.Corr,`[[`,1),function(x)x[which(x$Pval<=0.05),1])
# msmm_rnaseq_Plaque.Corr005.HighGenes=lapply(lapply(msmm_rnaseq_Plaque.Corr,`[[`,2),function(x)x[which(x$Pval<=0.05),1])

msmm_rnaseq_Plaque.Corr_Genes005.ExprData=vector(mode = "list",length = 3)
names(msmm_rnaseq_Plaque.Corr_Genes005.ExprData)=names(msmm_rnaseq)
# msmm_rnaseq_Plaque.Corr_Genes005.ExprData[[1]]=msmm_rnaseq_Plaque.Corr_Genes005.ExprData[[2]]=msmm_rnaseq_Plaque.Corr_Genes005.ExprData[[3]]=vector(mode = "list",length = 2)
# names(msmm_rnaseq_Plaque.Corr_Genes005.ExprData[[1]])=names(msmm_rnaseq_Plaque.Corr_Genes005.ExprData[[2]])=names(msmm_rnaseq_Plaque.Corr_Genes005.ExprData[[3]])=c("low","high")

cat(paste("Filtering whole dataset for significantly correlated plaque-pathology genes ...\n"))
for (i in 1:length(names(msmm_rnaseq))){
  msmm_rnaseq_Plaque.Corr_Genes005.ExprData[[i]]=lapply(msmm_rnaseq,function(x)x[which(rownames(x)%in%msmm_rnaseq_Plaque.Corr005.Genes[[i]]),])[[i]]
  #msmm_rnaseq_Plaque.Corr_Genes005.ExprData[[i]][[2]]=lapply(msmm_rnaseq,function(x)x[which(rownames(x)%in%msmm_rnaseq_Plaque.Corr005.HighGenes[[i]]),])[[i]]
}

cat(paste(rep("#",100),"\n",sep = ""))
cat(paste("Beginning DEG analysis ...\n"))
cat(paste(rep("#",100),"\n",sep = ""))
#DEG Analysis
msmm_rnaseq.colData=msmm_rnaseq.Plaque_DEG=vector(mode = "list",length=3)
names(msmm_rnaseq.colData)=names(msmm_rnaseq.Plaque_DEG)=c("BM_10","BM_22","BM_36")
for (i in 1:length(names(msmm_rnaseq.colData))){
  x=matrix(NA,nrow=sum(length(low_Plaque_Samples[[i]]),length(high_Plaque_Samples[[i]])),ncol=2)
  rownames(x)=c(low_Plaque_Samples[[i]],high_Plaque_Samples[[i]])
  x[,1]=c(rep("low",length(low_Plaque_Samples[[i]])),rep("high",length(high_Plaque_Samples[[i]])))
  x[,2]=c(rep("paired-end",sum(length(low_Plaque_Samples[[i]]),length(high_Plaque_Samples[[i]]))))
  msmm_rnaseq.colData[[i]]=data.frame(x,stringsAsFactors = F)
  #rownames(msmm_rnaseq.colData[[i]])=sort(c(low_Tangle_Samples[grep(pattern = names(msmm_rnaseq.colData)[i],low_Tangle_Samples)],high_Tangle_Samples[grep(pattern = names(msmm_rnaseq.colData)[i],high_Tangle_Samples)]))
  colnames(msmm_rnaseq.colData[[i]])=c("condition","type")
  rm(x)
}
# 
for(j in 1:length(names(msmm_rnaseq))){
   dds=DESeqDataSetFromMatrix(countData = msmm_rnaseq_raw[[j]][,rownames(msmm_rnaseq.colData[[j]])],colData = msmm_rnaseq.colData[[j]],design=~condition)
   msmm_rnaseq.Plaque_DEG[[j]]=dds
}

msmm_rnaseq.Plaque_DEG=lapply(msmm_rnaseq.Plaque_DEG,DESeq)
msmm_rnaseq.Plaque_DEG=lapply(msmm_rnaseq.Plaque_DEG,varianceStabilizingTransformation)
lowCounts=lapply(lapply(msmm_rnaseq.Plaque_DEG,counts),function(x)(rowSums(x)>1))
lowCounts_genes=lapply(lowCounts,function(x)which(x==F))
msmm_rnaseq.Plaque_DEG.Genes=lapply(lapply(lapply(msmm_rnaseq.Plaque_DEG,results),as.data.frame),function(x)rownames(x[which((x$padj<=0.05)&(abs(x$log2FoldChange)>=log2(2))),]))
########################################################################################################################

msmm_rnaseq_Plaque.Corr_Genes_PLQ.corrMatrix=vector(mode = "list",length = 3)
names(msmm_rnaseq_Plaque.Corr_Genes_PLQ.corrMatrix)=names(msmm_rnaseq)
msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.corrMatrix=msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.corrMatrix=vector(mode = "list",length = 3)
names(msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.corrMatrix)=names(msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.corrMatrix)=names(msmm_rnaseq)
for (i in 1:length(names(msmm_rnaseq))){
  #msmm_rnaseq_Plaque.Corr_Genes_PLQ.corrMatrix[[i]]=apply(cor(x=t(msmm_rnaseq_Plaque.Corr_Genes005.ExprData[[i]]),y=t(msmm_rnaseq_Plaque.Corr_Genes005.ExprData[[i]]),method = "spearman"),1,round,digits=3)
  msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.corrMatrix[[i]]=apply(cor(x=t(msmm_rnaseq_data2.agg[,low_Plaque_Samples[[i]]]),y=t(msmm_rnaseq_data2.agg[,low_Plaque_Samples[[i]]]),method = "spearman"),1,round,digits=3)
  msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.corrMatrix[[i]]=apply(cor(x=t(msmm_rnaseq_data2.agg[,high_Plaque_Samples[[i]]]),y=t(msmm_rnaseq_data2.agg[,high_Plaque_Samples[[i]]]),method = "spearman"),1,round,digits=3)
}
msmm_rnaseq_Plaque.Corr_Genes_PLQ.adjMatrix=vector(mode = "list",length = 3)
names(msmm_rnaseq_Plaque.Corr_Genes_PLQ.adjMatrix)=names(msmm_rnaseq)
msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.adjMatrix=msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.adjMatrix=vector(mode = "list",length = 3)
names(msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.adjMatrix)=names(msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.adjMatrix)=names(msmm_rnaseq)

for (t in 1:length(names(msmm_rnaseq))){
  tmp_matrix=matrix(NA,nrow=dim(msmm_rnaseq_Plaque.Corr_Genes_PLQ.corrMatrix[[t]])[1],ncol = dim(msmm_rnaseq_Plaque.Corr_Genes_PLQ.corrMatrix[[t]])[2])
  colnames(tmp_matrix)=rownames(tmp_matrix)=rownames(msmm_rnaseq_Plaque.Corr_Genes_PLQ.corrMatrix[[t]])
  for (g1 in 1:dim(msmm_rnaseq_Plaque.Corr_Genes_PLQ.corrMatrix[[t]])[1]){
    tmp1=msmm_rnaseq_Plaque.Corr_Genes_PLQ.corrMatrix[[t]][g1,]
    tmp1[unname(which(abs(tmp1)>=0.5))]=1
    tmp1[-unname(which(abs(tmp1)>=0.5))]=0
    tmp_matrix[g1,]=tmp1
    
  }
  msmm_rnaseq_Plaque.Corr_Genes_PLQ.adjMatrix[[t]]=tmp_matrix
}

for (t in 1:length(names(msmm_rnaseq))){
  tmp_matrix=matrix(NA,nrow=dim(msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.corrMatrix[[t]])[1],ncol = dim(msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.corrMatrix[[t]])[2])
  colnames(tmp_matrix)=rownames(tmp_matrix)=rownames(msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.corrMatrix[[t]])
  for (g1 in 1:dim(msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.corrMatrix[[t]])[1]){
    tmp1=msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.corrMatrix[[t]][g1,]
    tmp1[unname(which(abs(tmp1)>=0.25))]=1
    tmp1[-unname(which(abs(tmp1)>=0.25))]=0
    tmp_matrix[g1,]=tmp1
    
  }
  msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.adjMatrix[[t]]=tmp_matrix
}

for (t in 1:length(names(msmm_rnaseq))){
  tmp_matrix=matrix(NA,nrow=dim(msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.corrMatrix[[t]])[1],ncol = dim(msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.corrMatrix[[t]])[2])
  colnames(tmp_matrix)=rownames(tmp_matrix)=rownames(msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.corrMatrix[[t]])
  for (g1 in 1:dim(msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.corrMatrix[[t]])[1]){
    tmp1=msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.corrMatrix[[t]][g1,]
    tmp1[unname(which(abs(tmp1)>=0.25))]=1
    tmp1[-unname(which(abs(tmp1)>=0.25))]=0
    tmp_matrix[g1,]=tmp1
    
  }
  msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.adjMatrix[[t]]=tmp_matrix
}
####################################################################################################
msmm_rnaseq_Plaque.Corr_Genes_PLQ.Coexpp=vector(mode = "list",length = 3)
names(msmm_rnaseq_Plaque.Corr_Genes_PLQ.Coexpp)=names(msmm_rnaseq)
msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.Coexpp=msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.Coexpp=vector(mode = "list",length = 3)
names(msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.Coexpp)=names(msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.Coexpp)=names(msmm_rnaseq)
common_genes.coexpp=common_genes.coexppNet=vector(mode = "list",length = 3)
names(common_genes.coexpp)=names(common_genes.coexppNet)=names(msmm_rnaseq)
common_genes.coexppNet[[1]]=common_genes.coexppNet[[2]]=common_genes.coexppNet[[3]]=vector(mode = "list",length = 2)
names(common_genes.coexppNet[[1]])=names(common_genes.coexppNet[[3]])=names(common_genes.coexppNet[[3]])=c("low","high")

# for (i in 1:length(names(msmm_rnaseq))){
#   
#   msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.Coexpp[[i]]=graph.adjacency(adjmatrix = msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.adjMatrix[[i]],mode = "undirected",weighted = T,diag = F)
#   msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.Coexpp[[i]]=graph.adjacency(adjmatrix = msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.adjMatrix[[i]],mode = "undirected",weighted = T,diag = F)
#   V(msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.Coexpp[[i]])$deg=degree(graph = msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.Coexpp[[i]])
#   cat(paste("Computing Betweenness of",names(msmm_rnaseq)[i],"Low PLQ coexpression network for",length(V(msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.Coexpp[[i]])),"genes ...\n",sep = " "))
#   V(msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.Coexpp[[i]])$bwn=betweenness(graph = msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.Coexpp[[i]],directed = F)
#   cat(paste("Computing edge Betweenness of",names(msmm_rnaseq)[i],"Low PLQ coexpression network for",length(E(msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.Coexpp[[i]])),"edges ...\n",sep = " "))
#   E(msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.Coexpp[[i]])$bwn=edge_betweenness(graph = msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.Coexpp[[i]],directed = F)
#   cat(paste("Computing Closenness of",names(msmm_rnaseq)[i],"Low PLQ coexpression network for",length(V(msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.Coexpp[[i]])),"genes ...\n",sep = " "))
#   V(msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.Coexpp[[i]])$cls=closeness(graph = msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.Coexpp[[i]])
#   cat(paste("Computing Page Rank of",names(msmm_rnaseq)[i],"Low PLQ coexpression network for",length(V(msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.Coexpp[[i]])),"genes ...\n",sep = " "))
#   V(msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.Coexpp[[i]])$pgr=page_rank(graph = msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.Coexpp[[i]],directed = F)$vector
#   cat(paste("Computing clustering coefficient of",names(msmm_rnaseq)[i],"Low PLQ coexpression network for",length(V(msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.Coexpp[[i]])),"genes ...\n",sep = " "))
#   V(msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.Coexpp[[i]])$cc=transitivity(graph = msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.Coexpp[[i]],type="local")
#   cat(paste("Colouring nodes of",names(msmm_rnaseq)[i],"Low PLQ coexpression network for",length(V(msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.Coexpp[[i]])),"genes ...\n",sep = " "))
#   V(msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.Coexpp[[i]])$colour=c("250,250,250")
#   V(msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.Coexpp[[i]])$deg=degree(graph = msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.Coexpp[[i]])
#   cat(paste("Computing Betweenness of",names(msmm_rnaseq)[i],"High PLQ coexpression network for",length(V(msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.Coexpp[[i]])),"genes ...\n",sep = " "))
#   V(msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.Coexpp[[i]])$bwn=betweenness(graph = msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.Coexpp[[i]],directed = F)
#   cat(paste("Computing edge Betweenness of",names(msmm_rnaseq)[i],"High PLQ coexpression network for",length(E(msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.Coexpp[[i]])),"edges ...\n",sep = " "))
#   E(msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.Coexpp[[i]])$bwn=edge_betweenness(graph = msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.Coexpp[[i]],directed = F)
#   cat(paste("Computing Closenness of",names(msmm_rnaseq)[i],"High PLQ coexpression network for",length(V(msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.Coexpp[[i]])),"genes ...\n",sep = " "))
#   V(msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.Coexpp[[i]])$cls=closeness(graph = msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.Coexpp[[i]])
#   cat(paste("Computing Page Rank of",names(msmm_rnaseq)[i],"High PLQ coexpression network for",length(V(msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.Coexpp[[i]])),"genes ...\n",sep = " "))
#   V(msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.Coexpp[[i]])$pgr=page_rank(graph = msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.Coexpp[[i]],directed = F)$vector
#   cat(paste("Computing clustering coefficient of",names(msmm_rnaseq)[i],"High PLQ coexpression network for",length(V(msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.Coexpp[[i]])),"genes ...\n",sep = " "))
#   V(msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.Coexpp[[i]])$cc=transitivity(graph = msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.Coexpp[[i]],type="local")
#   cat(paste("Colouring nodes of",names(msmm_rnaseq)[i],"High PLQ coexpression network for",length(V(msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.Coexpp[[i]])),"genes ...\n",sep = " "))
#   V(msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.Coexpp[[i]])$colour=c("250,250,250")
#   cat(paste(rep("#",100),collapse = ""))
#   cat(paste("\n"))
#   cat(paste("Finding common genes between low and high ...\n"))
#   common_genes.coexpp[[i]]=intersect(V(msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.Coexpp[[i]])$name,V(msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.Coexpp[[i]])$name)
#   common_genes.coexppNet[[i]][[1]]=induced_subgraph(graph = msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.Coexpp[[i]],vids = common_genes.coexpp[[i]])
#   common_genes.coexppNet[[i]][[2]]=induced_subgraph(graph = msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.Coexpp[[i]],vids = common_genes.coexpp[[i]])
#   cat(paste("Colouring common genes with a different colour ...\n"))
#   V(msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.Coexpp[[i]])$colour[which(V(msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.Coexpp[[i]])$name%in%common_genes.coexpp[[i]])]=c("238,106,80")
#   V(msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.Coexpp[[i]])$colour[which(V(msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.Coexpp[[i]])$name%in%common_genes.coexpp[[i]])]=c("238,106,80")
#   write.graph(graph = msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.Coexpp[[i]],file = paste(names(msmm_rnaseq)[i],"_LowCoexppGenes.gml",sep=""),format = "gml")
#   write.graph(graph = msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.Coexpp[[i]],file = paste(names(msmm_rnaseq)[i],"_HighCoexppGenes.gml",sep=""),format = "gml")
#   cat(paste("Writing graph files ...\n",sep = " "))
#   write.graph(common_genes.coexppNet[[i]][[1]],file = paste(names(msmm_rnaseq)[i],"_Common_LowCoexppGenes.gml",sep = ""),format = "gml")
#   write.graph(common_genes.coexppNet[[i]][[2]],file = paste(names(msmm_rnaseq)[i],"_Common_HighCoexppGenes.gml",sep = ""),format = "gml")
#   cat(paste("Done!\n"))
# }


# topology_CM=vector(mode = "list",length = 3)
# names(topology_CM)=names(msmm_rnaseq)
# 
# topology=c("Degree(Connectivity)","Betweenness","Closeness","PageRank")
# # for (m in 1:dim(topology_compare_matrix)[1]){
# for (t in 1:length(names(msmm_rnaseq))){
#     topology_compare_matrix=matrix(NA,nrow = 5,ncol = 4)
#     rownames(topology_compare_matrix)=c("Degree","Betweenness","Closeness","PGR","Clustering coefficient")
#     colnames(topology_compare_matrix)=c("#Low_Coexpp","#High_Coexpp","P.value","CommonGenes")
#     
#     common_genes=intersect(V(msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.Coexpp[[t]])$name,V(msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.Coexpp[[t]])$name)
#     topology_compare_matrix[1,]=c(length(V(msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.Coexpp[[t]])$name),
#                                   length(V(msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.Coexpp[[t]])$name),
#                                   wilcox.test(V(msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.Coexpp[[t]])$deg,V(msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.Coexpp[[t]])$deg,alternative = "two.sided")$p.value,
#                                   paste(common_genes,collapse = ","))
#     topology_compare_matrix[2,]=c(length(V(msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.Coexpp[[t]])$name),
#                                   length(V(msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.Coexpp[[t]])$name),
#                                   wilcox.test(V(msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.Coexpp[[t]])$bwn,V(msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.Coexpp[[t]])$bwn,alternative = "two.sided")$p.value,
#                                   paste(common_genes,collapse = ","))
#     topology_compare_matrix[3,]=c(length(V(msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.Coexpp[[t]])$name),
#                                   length(V(msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.Coexpp[[t]])$name),
#                                   wilcox.test(V(msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.Coexpp[[t]])$cls,V(msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.Coexpp[[t]])$cls,alternative = "two.sided")$p.value,
#                                   paste(common_genes,collapse = ","))
#     topology_compare_matrix[4,]=c(length(V(msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.Coexpp[[t]])$name),
#                                   length(V(msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.Coexpp[[t]])$name),
#                                   wilcox.test(V(msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.Coexpp[[t]])$pgr,V(msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.Coexpp[[t]])$pgr,alternative = "two.sided")$p.value,
#                                   paste(common_genes,collapse = ","))
#     topology_compare_matrix[5,]=c(length(V(msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.Coexpp[[t]])$name),
#                                   length(V(msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.Coexpp[[t]])$name),
#                                   wilcox.test(V(msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.Coexpp[[t]])$cc,V(msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.Coexpp[[t]])$cc,alternative = "two.sided")$p.value,
#                                   paste(common_genes,collapse = ","))
#     topology_CM[[t]]=data.frame(topology_compare_matrix,stringsAsFactors = F)
# }


# for (t in 1:length(msmm_rnaseq)){
#   V(msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.Coexpp[[t]])$colour[which(V(msmm_rnaseq_Plaque.Corr_Genes_LowPLQ.Coexpp[[t]])$name%in%common_genes.coexpp[[t]])]=c("238,106,80")
#   V(msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.Coexpp[[t]])$colour[which(V(msmm_rnaseq_Plaque.Corr_Genes_HighPLQ.Coexpp[[t]])$name%in%common_genes.coexpp[[t]])]=c("238,106,80")
# }


#save.image(file = "msmm_rnaseq_FinalDataset_contScore_Coexpp.RData")
#Comparing correlations using the Cocor package
msmm_rnaseq.cocor=vector(mode = "list",length = 3)
msmm_rnaseq.cocor[[1]]=msmm_rnaseq.cocor[[2]]=msmm_rnaseq.cocor[[3]]=vector(mode = "list",length = 2)
names(msmm_rnaseq.cocor)=names(msmm_rnaseq)
names(msmm_rnaseq.cocor[[1]])=names(msmm_rnaseq.cocor[[2]])=names(msmm_rnaseq.cocor[[3]])=c("pval","zstat")
for (t in 1:length(names(msmm_rnaseq))){
  tmp5=matrix(data = NA,nrow = length(msmm_rnaseq_Plaque.Corr005.Genes[[t]][1:20]),ncol = length(msmm_rnaseq_Plaque.Corr005.Genes[[t]][1:20]))
  rownames(tmp5)=colnames(tmp5)=msmm_rnaseq_Plaque.Corr005.Genes[[t]][1:20]
  tmp6=matrix(data = NA,nrow = length(msmm_rnaseq_Plaque.Corr005.Genes[[t]][1:20]),ncol = length(msmm_rnaseq_Plaque.Corr005.Genes[[t]][1:20]))
  rownames(tmp6)=colnames(tmp6)=msmm_rnaseq_Plaque.Corr005.Genes[[t]][1:20]
  for (a in 1:dim(tmp5)[1]){
    cat(paste("Filtering gene ",rownames(msmm_rnaseq[[t]])[a],sep = ""))
    for (b in 1:dim(tmp5)[2]){
      cat(paste("And ",rownames(msmm_rnaseq[[t]])[b],"for brain region",names(msmm_rnaseq)[t],"\n",sep = ""))
      x1=msmm_rnaseq[[t]][which(rownames(msmm_rnaseq[[t]])%in%msmm_rnaseq_Plaque.Corr005.Genes[[t]][a]),low_Plaque_Samples[[t]]]
      y1=msmm_rnaseq[[t]][which(rownames(msmm_rnaseq[[t]])%in%msmm_rnaseq_Plaque.Corr005.Genes[[t]][b]),low_Plaque_Samples[[t]]]
      x2=msmm_rnaseq[[t]][which(rownames(msmm_rnaseq[[t]])%in%msmm_rnaseq_Plaque.Corr005.Genes[[t]][a]),high_Plaque_Samples[[t]]]
      y2=msmm_rnaseq[[t]][which(rownames(msmm_rnaseq[[t]])%in%msmm_rnaseq_Plaque.Corr005.Genes[[t]][b]),high_Plaque_Samples[[t]]]
      cat(paste("Computing cocor for ",rownames(msmm_rnaseq[[t]])[a],"and",rownames(msmm_rnaseq[[t]])[b],"\n",sep = ""))
      cocorCalc=cocor.indep.groups(r1.jk = as.numeric(cor(x = t(x1),y = t(y1),method = "spearman")),r2.hm = as.numeric(cor(x = t(x2),y = t(y2),method = "spearman")),n1 = length(low_Plaque_Samples[[t]]),n2 = length(high_Plaque_Samples[[t]]),test = "fisher1925")
      tmp5[a,b]=cocorCalc@fisher1925$p.value
      msmm_rnaseq.cocor[[t]][[1]]=tmp5
      tmp6[a,b]=cocorCalc@fisher1925$statistic  
      msmm_rnaseq.cocor[[t]][[2]]=tmp6
    }
  }
  cat(paste("Done!\n"))
}

sig_edges=apply(msmm_rnaseq.cocor[[1]][[1]],1,function(x)which(x<=0.05))
sig_edges_indices=lapply(sig_edges,length)

########################################################################
test1=msmm_rnaseq.cocor_filtered$p001$BM_36[order(msmm_rnaseq.cocor_filtered$p001$BM_36$abs.corr.change,decreasing = T),][1:5,]
gene1=msmm_rnaseq$BM_36[which(rownames(msmm_rnaseq$BM_36)%in%union(test1$Gene.A,test1$Gene.B)),which(colnames(msmm_rnaseq$BM_36)%in%high_Plaque_Samples$BM_36)][3,]
gene2=msmm_rnaseq$BM_36[which(rownames(msmm_rnaseq$BM_36)%in%union(test1$Gene.A,test1$Gene.B)),which(colnames(msmm_rnaseq$BM_36)%in%high_Plaque_Samples$BM_36)][4,]
plq=msmm_rnaseq_covariates$PlaqueMean[which(msmm_rnaseq_covariates$Sample.ID%in%colnames(msmm_rnaseq$BM_36))]
df_anno=data.frame(PMI=msmm_rnaseq_covariates$PMI[which(msmm_rnaseq_covariates$Sample.ID%in%high_Plaque_Samples$BM_36)],
                   CERAD=msmm_rnaseq_covariates$CERAD[which(msmm_rnaseq_covariates$Sample.ID%in%high_Plaque_Samples$BM_36)],
                   Braak=msmm_rnaseq_covariates$bbscore[which(msmm_rnaseq_covariates$Sample.ID%in%high_Plaque_Samples$BM_36)],
                   AOD=msmm_rnaseq_covariates$AOD[which(msmm_rnaseq_covariates$Sample.ID%in%high_Plaque_Samples$BM_36)])
rownames(df_anno)=colnames(gene1)
pheatmap(as.matrix(msmm_rnaseq$BM_36[which(rownames(msmm_rnaseq$BM_36)%in%union(test1$Gene.A,test1$Gene.B)),which(colnames(msmm_rnaseq$BM_36)%in%high_Plaque_Samples$BM_36)]),annotation = df_anno)
