library(foreach)
library(doParallel)

msbb_rnaseq2016_data=read.table("AMP-AD_MSBB_MSSM_IlluminaHiSeq2500_normalized_counts_September_2016.txt",sep="\t",header = T,as.is = T)
msbb_rnaseq_covariates=read.csv("MSBB_RNAseq_covariates.csv",header = T,as.is = T)
msbb_rnaseq_clinical_covariates=read.csv("MSBB_clinical.csv",header = T,as.is = T)
colnames(msbb_rnaseq2016_data)=unlist(lapply(strsplit(x = colnames(msbb_rnaseq2016_data),split = "X"),`[[`,2))
msbb_ensembl_symbol=data.frame(Ensembl=names(mapIds(x = org.Hs.eg.db,keys=rownames(msbb_rnaseq2016_data),column = "SYMBOL",keytype = "ENSEMBL",multiVals = "first")),Symbol=mapIds(x = org.Hs.eg.db,keys=rownames(msbb_rnaseq2016_data),column = "SYMBOL",keytype = "ENSEMBL",multiVals = "first"),stringsAsFactors = F)
msbb_rnaseq_covariates.merged=merge(x = msbb_rnaseq_clinical_covariates,y=msbb_rnaseq_covariates,by=c("individualIdentifier","individualIdentifier"))
msbb_rnaseq_covariates.merged2=msbb_rnaseq_covariates.merged[grep(pattern = "unmapped|resequenced",x = msbb_rnaseq_covariates.merged$fileName,invert = T),]


msbb_rnaseq2016_byRegion=msbb_rnaseq_covariates.merged_final=msbb_rnaseq2016_PLQGenes=vector(mode = "list",length = 4)
names(msbb_rnaseq2016_byRegion)=names(msbb_rnaseq_covariates.merged_final)=names(msbb_rnaseq2016_PLQGenes)=c("FP","IFG","PHG","STG") 
msbb_rnaseq2016_data2=msbb_rnaseq2016_data[-which(rownames(msbb_rnaseq2016_data)%in%msbb_ensembl_symbol$Ensembl[which(is.na(msbb_ensembl_symbol$Symbol)==T)]),]
msbb_rnaseq2016_data2=data.frame(cbind(Symbol=msbb_ensembl_symbol$Symbol[-which(is.na(msbb_ensembl_symbol$Symbol)==T)],msbb_rnaseq2016_data2),stringsAsFactors = F)
msbb_rnaseq2016_data2.agg=aggregate(x=msbb_rnaseq2016_data2[,-1],by=list(geneSymbol=msbb_rnaseq2016_data2$Symbol),mean)
#rm_genes=grep(pattern = "_|\\.|-",msbb_rnaseq2016_data2.agg$geneSymbol,invert = F)
rownames(msbb_rnaseq2016_data2.agg)=msbb_rnaseq2016_data2.agg$geneSymbol
msbb_rnaseq2016_data2.agg=msbb_rnaseq2016_data2.agg[,-1]
colnames(msbb_rnaseq2016_data2.agg)=gsub(pattern = "X",replacement = "",x = colnames(msbb_rnaseq2016_data2.agg))

msbb_rnaseq_covariates.merged_final$FP=msbb_rnaseq_covariates.merged2[which(unlist(lapply(strsplit(x=msbb_rnaseq_covariates.merged2$sampleIdentifier,split = "_"),`[[`,3))%in%colnames(msbb_rnaseq2016_byRegion$FP)),c(1:11,13,15:16)]
msbb_rnaseq_covariates.merged_final$IFG=msbb_rnaseq_covariates.merged2[which(unlist(lapply(strsplit(x=msbb_rnaseq_covariates.merged2$sampleIdentifier,split = "_"),`[[`,3))%in%colnames(msbb_rnaseq2016_byRegion$IFG)),c(1:11,13,15:16)]
msbb_rnaseq_covariates.merged_final$PHG=msbb_rnaseq_covariates.merged2[which(unlist(lapply(strsplit(x=msbb_rnaseq_covariates.merged2$sampleIdentifier,split = "_"),`[[`,3))%in%colnames(msbb_rnaseq2016_byRegion$PHG)),c(1:11,13,15:16)]
msbb_rnaseq_covariates.merged_final$STG=msbb_rnaseq_covariates.merged2[which(unlist(lapply(strsplit(x=msbb_rnaseq_covariates.merged2$sampleIdentifier,split = "_"),`[[`,3))%in%colnames(msbb_rnaseq2016_byRegion$STG)),c(1:11,13,15:16)]

msbb_rnaseq2016_byRegion$FP=msbb_rnaseq2016_data2.agg[,which(colnames(msbb_rnaseq2016_data2.agg)%in%unlist(lapply(strsplit(x = msbb_rnaseq_covariates.merged2$sampleIdentifier[which(msbb_rnaseq_covariates.merged2$BrodmannArea=="BM10")],split = "_"),`[[`,3)))]
msbb_rnaseq2016_byRegion$IFG=msbb_rnaseq2016_data2.agg[,which(colnames(msbb_rnaseq2016_data2.agg)%in%unlist(lapply(strsplit(x = msbb_rnaseq_covariates.merged2$sampleIdentifier[which(msbb_rnaseq_covariates.merged2$BrodmannArea=="BM44")],split = "_"),`[[`,3)))]
msbb_rnaseq2016_byRegion$PHG=msbb_rnaseq2016_data2.agg[,which(colnames(msbb_rnaseq2016_data2.agg)%in%unlist(lapply(strsplit(x = msbb_rnaseq_covariates.merged2$sampleIdentifier[which(msbb_rnaseq_covariates.merged2$BrodmannArea=="BM36")],split = "_"),`[[`,3)))]
msbb_rnaseq2016_byRegion$STG=msbb_rnaseq2016_data2.agg[,which(colnames(msbb_rnaseq2016_data2.agg)%in%unlist(lapply(strsplit(x = msbb_rnaseq_covariates.merged2$sampleIdentifier[which(msbb_rnaseq_covariates.merged2$BrodmannArea=="BM22")],split = "_"),`[[`,3)))]

lowPlaque_samples=highPlaque_samples=vector(mode = "list",length = 4)
names(lowPlaque_samples)=names(highPlaque_samples)=names(msbb_rnaseq2016_byRegion)
lowPlaque_samples=lapply(msbb_rnaseq_covariates.merged_final,function(x)x$sampleIdentifier[which(x$PlaqueMean<=1)])
highPlaque_samples=lapply(msbb_rnaseq_covariates.merged_final,function(x)x$sampleIdentifier[which(x$PlaqueMean>=15)])

i=0
blocksize=300
blocks=round(length(rownames(msbb_rnaseq2016_byRegion$FP))/blocksize,digits = 0)#length(rownames(msbb_rnaseq2016_byRegion[[1]]))/5
start<-i*blocksize+1
end<-min((i+1)*blocksize, length(rownames(msbb_rnaseq2016_byRegion$FP)))
broadman_area=c("BM10","BM44","BM36","BM22")
for (r in 1:length(msbb_rnaseq2016_byRegion)){
  while(start < end){
  
    input<-start:end
    pb = txtProgressBar(min=0,max=length(input),style=3,initial=0)
    cat("\n")
    plaque=msbb_rnaseq_covariates.merged_final[[r]]$PlaqueMean[which(msbb_rnaseq_covariates.merged_final[[r]]$BrodmannArea==broadman_area[[r]])]
    geneCounts=list()
    for (l in 1:length(input)){
      geneCounts[[l]]=unlist(unname(msbb_rnaseq2016_byRegion[[r]][l,which(colnames(msbb_rnaseq2016_byRegion[[r]])%in%unlist(lapply(strsplit(x = msbb_rnaseq_covariates.merged_final[[r]]$sampleIdentifier,split = "_"),`[[`,3)))]))  
    }
    names(geneCounts)=rownames(msbb_rnaseq2016_byRegion[[r]][input,which(colnames(msbb_rnaseq2016_byRegion[[r]])%in%unlist(lapply(strsplit(x = msbb_rnaseq_covariates.merged_final[[r]]$sampleIdentifier,split = "_"),`[[`,3)))])
    res = mclapply(geneCounts,cor.test,plaque,method = "spearman",mc.cores = nc)
    close(pb)
    rho.p=mclapply(res,function(x)x$p.value,mc.cores = nc)
    rho=mclapply(res,function(x)unname(x$estimate),mc.cores = nc)
    result=data.frame(Genes=names(geneCounts),Rho=unlist(rho),Rho.p=unlist(rho.p),stringsAsFactors = F)
    write.table(result, file=paste(names(msbb_rnaseq2016_byRegion)[r],"_",i, ".txt",sep = ""), sep="\t",col.names = T, row.names=FALSE, quote = FALSE)
    #write.table(result, file=paste0("/shared/hidelab2/user/md4zsa/Work/Data/MSMM_RNAseq/MSMM_RNAseq_FinalRelease2/",names(msbb_rnaseq2016_byRegion)[r], i, ".txt"), sep="\t",col.names = T, row.names=FALSE, quote = FALSE)
    i<-i+1
    start<-i*blocksize+1
    end<-min((i+1)*blocksize, length(rownames(msbb_rnaseq2016_byRegion$FP)))  
  }
  
}
cat(paste("Done!\n"))

