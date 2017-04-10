library(DESeq2)
library(org.Hs.eg.db)
library(clusterProfiler)

msbb_rnaseq2016_rawData=read.table("AMP-AD_MSBB_MSSM_IlluminaHiSeq2500_raw_counts_September_2016.txt",sep = "\t",header = T,as.is = T)
msbb_rnaseq_covariates=read.csv("MSBB_RNAseq_covariates.csv",header = T,as.is = T)
msbb_rnaseq_clinical_covariates=read.csv("MSBB_clinical.csv",header = T,as.is = T)
colnames(msbb_rnaseq2016_rawData)=unlist(lapply(strsplit(x = colnames(msbb_rnaseq2016_rawData),split = "X"),`[[`,2))
msbb_ensembl_symbol2=data.frame(Ensembl=names(mapIds(x = org.Hs.eg.db,keys=rownames(msbb_rnaseq2016_rawData),column = "SYMBOL",keytype = "ENSEMBL",multiVals = "first")),Symbol=mapIds(x = org.Hs.eg.db,keys=rownames(msbb_rnaseq2016_rawData),column = "SYMBOL",keytype = "ENSEMBL",multiVals = "first"),stringsAsFactors = F)
msbb_rnaseq_covariates.merged=merge(x = msbb_rnaseq_clinical_covariates,y=msbb_rnaseq_covariates,by=c("individualIdentifier","individualIdentifier"))
msbb_rnaseq_covariates.merged2=msbb_rnaseq_covariates.merged[grep(pattern = "unmapped",x = msbb_rnaseq_covariates.merged$fileName,invert = T),]

msbb_rnaseq2016_rawData2=msbb_rnaseq2016_rawData[-which(rownames(msbb_rnaseq2016_rawData)%in%msbb_ensembl_symbol2$Ensembl[which(is.na(msbb_ensembl_symbol2$Symbol)==T)]),]
msbb_rnaseq2016_rawData2=data.frame(cbind(Symbol=msbb_ensembl_symbol2$Symbol[-which(is.na(msbb_ensembl_symbol2$Symbol)==T)],msbb_rnaseq2016_rawData2),stringsAsFactors = F)
msbb_rnaseq2016_rawData2.agg=aggregate(x=msbb_rnaseq2016_rawData2[,-1],by=list(geneSymbol=msbb_rnaseq2016_rawData2$Symbol),mean)
colnames(msbb_rnaseq2016_rawData2.agg)=gsub(pattern = "X",replacement = "",x = colnames(msbb_rnaseq2016_rawData2.agg))
rownames(msbb_rnaseq2016_rawData2.agg)=msbb_rnaseq2016_rawData2.agg$geneSymbol
lowPlaque_samples=highPlaque_samples=msbb_rnaseq.colData=msbb_rnaseq_covariates.merged_final=msbb_rnaseq2016_byRegion_rawData=vector(mode = "list",length = 4)
names(lowPlaque_samples)=names(highPlaque_samples)=names(msbb_rnaseq.colData)=names(msbb_rnaseq2016_byRegion_rawData)=names(msbb_rnaseq_covariates.merged_final)=c("FP","IFG","PHG","STG")

# msbb_rnaseq_covariates.merged=merge(x = msbb_rnaseq_clinical_covariates,y=msbb_rnaseq_covariates,by=c("individualIdentifier","individualIdentifier"))
#msbb_rnaseq_covariates.merged2=msbb_rnaseq_covariates.merged[grep(pattern = "unmapped|resequenced",x = msbb_rnaseq_covariates.merged$fileName,invert = T),]
#Read region-wise covariates independently from text files to avoid whitespace issues
msbb_rnaseq_covariates.merged_final$FP=read.table("MSBB_RNAseq2016_FP_covariates.txt",sep = "\t",header = T,as.is = T)
msbb_rnaseq_covariates.merged_final$IFG=read.table("MSBB_RNAseq2016_IFG_covariates.txt",sep = "\t",header = T,as.is = T)
msbb_rnaseq_covariates.merged_final$PHG=read.table("MSBB_RNAseq2016_PHG_covariates.txt",sep = "\t",header = T,as.is = T)
msbb_rnaseq_covariates.merged_final$STG=read.table("MSBB_RNAseq2016_STG_covariates.txt",sep = "\t",header = T,as.is = T)

lowPlaque_samples=lapply(msbb_rnaseq_covariates.merged_final,function(x)x$sampleIdentifier[which(x$PlaqueMean<=1)])
highPlaque_samples=lapply(msbb_rnaseq_covariates.merged_final,function(x)x$sampleIdentifier[which(x$PlaqueMean>=15)])
#DEG  analysis by DESeq2
msbb_rnaseq2016_byRegion_rawData$FP=msbb_rnaseq2016_rawData2[,which(colnames(msbb_rnaseq2016_rawData2)%in%unlist(lapply(strsplit(x = msbb_rnaseq_covariates.merged2$sampleIdentifier[which(msbb_rnaseq_covariates.merged2$BrodmannArea=="BM10")],split = "_"),`[[`,3)))]
msbb_rnaseq2016_byRegion_rawData$IFG=msbb_rnaseq2016_rawData2[,which(colnames(msbb_rnaseq2016_rawData2)%in%unlist(lapply(strsplit(x = msbb_rnaseq_covariates.merged2$sampleIdentifier[which(msbb_rnaseq_covariates.merged2$BrodmannArea=="BM44")],split = "_"),`[[`,3)))]
msbb_rnaseq2016_byRegion_rawData$PHG=msbb_rnaseq2016_rawData2[,which(colnames(msbb_rnaseq2016_rawData2)%in%unlist(lapply(strsplit(x = msbb_rnaseq_covariates.merged2$sampleIdentifier[which(msbb_rnaseq_covariates.merged2$BrodmannArea=="BM36")],split = "_"),`[[`,3)))]
msbb_rnaseq2016_byRegion_rawData$STG=msbb_rnaseq2016_rawData2[,which(colnames(msbb_rnaseq2016_rawData2)%in%unlist(lapply(strsplit(x = msbb_rnaseq_covariates.merged2$sampleIdentifier[which(msbb_rnaseq_covariates.merged2$BrodmannArea=="BM22")],split = "_"),`[[`,3)))]

msbb_rnaseq.colData$FP=data.frame(condition=c(rep("low",length(unlist(lapply(strsplit(lowPlaque_samples$FP,split = "_"),`[[`,3)))),rep("high",length(unlist(lapply(strsplit(highPlaque_samples$FP,split = "_"),`[[`,3))))),type="paired-end")
msbb_rnaseq.colData$IFG=data.frame(condition=c(rep("low",length(unlist(lapply(strsplit(lowPlaque_samples$IFG,split = "_"),`[[`,3)))),rep("high",length(unlist(lapply(strsplit(highPlaque_samples$IFG,split = "_"),`[[`,3))))),type="paired-end")
msbb_rnaseq.colData$PHG=data.frame(condition=c(rep("low",length(unlist(lapply(strsplit(lowPlaque_samples$PHG,split = "_"),`[[`,3)))),rep("high",length(unlist(lapply(strsplit(highPlaque_samples$PHG,split = "_"),`[[`,3))))),type="paired-end")
msbb_rnaseq.colData$STG=data.frame(condition=c(rep("low",length(unlist(lapply(strsplit(lowPlaque_samples$STG,split = "_"),`[[`,3)))),rep("high",length(unlist(lapply(strsplit(highPlaque_samples$STG,split = "_"),`[[`,3))))),type="paired-end")
rownames(msbb_rnaseq.colData$FP)=c(unlist(lapply(strsplit(lowPlaque_samples$FP,split = "_"),`[[`,3)),unlist(lapply(strsplit(highPlaque_samples$FP,split = "_"),`[[`,3)))
rownames(msbb_rnaseq.colData$IFG)=c(unlist(lapply(strsplit(lowPlaque_samples$IFG,split = "_"),`[[`,3)),unlist(lapply(strsplit(highPlaque_samples$IFG,split = "_"),`[[`,3)))
rownames(msbb_rnaseq.colData$PHG)=c(unlist(lapply(strsplit(lowPlaque_samples$PHG,split = "_"),`[[`,3)),unlist(lapply(strsplit(highPlaque_samples$PHG,split = "_"),`[[`,3)))
rownames(msbb_rnaseq.colData$STG)=c(unlist(lapply(strsplit(lowPlaque_samples$STG,split = "_"),`[[`,3)),unlist(lapply(strsplit(highPlaque_samples$STG,split = "_"),`[[`,3)))
msbb_rnaseq_LowHigh_samples=msbb_rnaseq_LowHigh.dds=vector(mode = "list",length = 4)
names(msbb_rnaseq_LowHigh_samples)=names(msbb_rnaseq_LowHigh.dds)=names(msbb_rnaseq.colData)
msbb_rnaseq_LowHigh_samples$FP=msbb_rnaseq2016_rawData2.agg[,rownames(msbb_rnaseq.colData$FP)]
msbb_rnaseq_LowHigh_samples$IFG=msbb_rnaseq2016_rawData2.agg[,rownames(msbb_rnaseq.colData$IFG)]
msbb_rnaseq_LowHigh_samples$PHG=msbb_rnaseq2016_rawData2.agg[,rownames(msbb_rnaseq.colData$PHG)]
msbb_rnaseq_LowHigh_samples$STG=msbb_rnaseq2016_rawData2.agg[,rownames(msbb_rnaseq.colData$STG)]

msbb_rnaseq_LowHigh.dds$FP=DESeqDataSetFromMatrix(countData = t(data.frame(apply(msbb_rnaseq_LowHigh_samples$FP,1,as.integer))),colData = msbb_rnaseq.colData$FP,design = ~condition)
msbb_rnaseq_LowHigh.dds$IFG=DESeqDataSetFromMatrix(countData = t(data.frame(apply(msbb_rnaseq_LowHigh_samples$IFG,1,as.integer))),colData = msbb_rnaseq.colData$IFG,design = ~condition)
msbb_rnaseq_LowHigh.dds$PHG=DESeqDataSetFromMatrix(countData = t(data.frame(apply(msbb_rnaseq_LowHigh_samples$PHG,1,as.integer))),colData = msbb_rnaseq.colData$PHG,design = ~condition)
msbb_rnaseq_LowHigh.dds$STG=DESeqDataSetFromMatrix(countData = t(data.frame(apply(msbb_rnaseq_LowHigh_samples$STG,1,as.integer))),colData = msbb_rnaseq.colData$STG,design = ~condition)

msbb_rnaseq_LowHigh.dds=lapply(msbb_rnaseq_LowHigh.dds,DESeq)
msbb_rnaseq_LowHigh.results=lapply(lapply(msbb_rnaseq_LowHigh.dds,results,tidy=T),data.frame,stringsAsFactors = F)
msbb_rnaseq_LowHigh.DEG_Genes=lapply(msbb_rnaseq_LowHigh.results,function(x)rownames(x[which(abs(x$log2FoldChange)>1&x$padj<0.05),]))

write.table(msbb_rnaseq_LowHigh.results$FP,"MSBB_RNAseq2016_FP_DEG.txt",sep = "\t",col.names = T,row.names = F,quote = F)
write.table(msbb_rnaseq_LowHigh.results$IFG,"MSBB_RNAseq2016_IFG_DEG.txt",sep = "\t",col.names = T,row.names = F,quote = F)
write.table(msbb_rnaseq_LowHigh.results$PHG,"MSBB_RNAseq2016_PHG_DEG.txt",sep = "\t",col.names = T,row.names = F,quote = F)
write.table(msbb_rnaseq_LowHigh.results$STG,"MSBB_RNAseq2016_STG_DEG.txt",sep = "\t",col.names = T,row.names = F,quote = F)













