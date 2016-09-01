library(pathprint)
library(pheatmap)
library(foreach)
library(entropy)
data("pathprint.Hs.gs")
library(metaArray)
library(doMC)

#Read in RNAseq normalised count data, as obtained directly from Synapse
msbb_data=read.table("AMP-AD_MSBB_MSSM_IlluminaHiSeq2500_mRNA_normalized-sex-race-age-RIN-PMI-batch-site.corrected.csv",sep=",",header=T,as.is=T,row.names=1)
#Read in Ensembl-Entrez mapping, as generated from Biomart
msbb_ensg_entrez_mapped=read.table("msbb_ensg_mapped2.txt",sep="\t",header=T,as.is=T)
#Read in metadata that includes controls and other phenotype infp
msbb_metadata=read.table("AMP-AD_MSBB_MSSM_covariates_mRNA_IlluminaHiSeq2500.tsv",sep="\t",header=T,as.is=T)
#Make controls from Covariates for MSMM based on CERJ scores, all 1's are classified as normal brain;CDR=0 as No AD
msbb_data_meta_controls=msbb_metadata[which((msbb_metadata$CERJ=="1A"|msbb_metadata$CERJ=="1B"|msbb_metadata$CERJ=="1C"|msbb_metadata$CERJ=="1") & msbb_metadata$CDR==0),]
#Select those rows in the data, that have corresponding Entrez identifiers. Since
#several Ensembl IDs map to a single Entrez, only those with a 1:1 mapping are preserved.
#Also not all Ensembl IDs map to an Entrez ID. Hence, the 1:1 mapping.
#This results in 18031 rows in the data from 24413 in the original dataset.
index=which(unique(rownames(msbb_data))%in%msbb_ensg_entrez_mapped$Ensembl.Gene.ID)
msbb_data_entrez=msbb_ensg_entrez_mapped$EntrezGene.ID[which(rownames(msbb_data[index,])%in%msbb_ensg_entrez_mapped$Ensembl.Gene.ID)]
#Aggregating all Entrez and their corresponsing data points
msbb.entrez.agg=aggregate(msbb_data[index,],by=list(EntrezID=msbb_data_entrez),FUN=max)
rownames(msbb.entrez.agg)=msbb.entrez.agg[,1]
msbb.entrez.agg=msbb.entrez.agg[,-1]
#Compute SCE for the data
msbb.SCE=single.chip.enrichment(exprs=msbb.entrez.agg,geneset=pathprint.Hs.gs,transformation="squared.rank",statistic="mean",normalizedScore=F,progressBar=T)
save(msbb.SCE,file="MSBB.SCE.Rdata")

#Build the POE
pb <- txtProgressBar(min = 0,max = length(pathprint.Hs.gs),style = 3)
POE.msbb=foreach(i=names(pathprint.Hs.gs)) %do%{
  sample=msbb.SCE[i,]
  sample.POE=sample
  sample.POE[]=NA
  sample.valid=sample[!(is.na(sample))]
  fit=fit.em(sample.valid,cl=rep(0,length(sample.valid)))
  setTxtProgressBar(pb,match(i,names(pathprint.Hs.gs)))
  sample.POE[!(is.na(sample))]=fit$expr
}
POE.msbb.matrix=as.matrix(t(as.data.frame(POE.msbb)))
rownames(POE.msbb.matrix)=names(pathprint.Hs.gs)
colnames(POE.msbb.matrix)=colnames(msbb_data)
save(POE.msbb.matrix,file="MSBB.POE.matrix.RData")
pheatmap(POE.msbb.matrix[150:250,150:250],color=c("blue","white","red"),
         cluster_rows=T,cluster_cols=T,
         show_colnames=F,clustering_distance_rows="manhattan",
         clustering_distance_cols="manhattan",
         fontsize=9,fontsize_row=9,cell_width=8,cell_height=8,
         width=100,legend_breaks=c(-1,0,1),clustering_method="average",
         legend=T,annotation=msbb_annotation,
         main="PathPrint Mt.Sinai 448 samples")
# 
# pdf(file="POE.Mayo.Pathprint.HMap.pdf",pointsize=9)
# #HEatmap for Mayo dataset
# 
#          #cellwidth=5,cellheight=5
#          #)
# #HEatmap for Mt Sinai dataset
# pheatmap(POE.mayo_msbb[1:100,1:100],color=c("blue","white","red"),
#         cluster_rows=T,cluster_cols=T,
#         show_colnames=F,clustering_distance_rows="manhattan",
#         clustering_distance_cols="manhattan",
#         fontsize=6,fontsize_row=6,
#         height=100,
#         legend_breaks=c(-1,0,1),clustering_method="average",legend=T,
#         main="PathPrint Mayo/Mt.Sinai 726 samples"
#         )
# dev.off()
# msbb_annotation=data.frame(Var1=factor(1:dim(POE.matrix)[2],labels="Test"),Var2=factor(1:dim(POE.matrix)[2],labels="Tissue"))
# rownames(msbb_annotation)=colnames(POE.matrix)
# levels(msbb_annotation$Var1)[-msbb_controls_annotationIndex]="Test"
# levels(msbb_annotation$Var1)[msbb_controls_annotationIndex]="Control"
# levels(msbb_annotation$Var2)[grep("BM_36",rownames(msbb_annotation$Var2))]="BM_36"
# levels(msbb_annotation$Var2)[grep("BM_22",rownames(msbb_annotation$Var2))]="BM_22"
# levels(msbb_annotation$Var2)[grep("BM_10",rownames(msbb_annotation$Var2))]="BM_10"
# 
# POE.mayo_msbb=foreach (i = names(pathprint.Hs.gs)) %do%{
# sample = SCE_Mayo_MSBB[i,]
# sample.POE.mayo_msbb=sample
# sample.POE.mayo_msbb[]=NA
# sample.valid=sample[!(is.na(sample))]
# fit=fit.em(sample.valid,cl=rep(0,length(sample.valid)))
# setTxtProgressBar(pb, match(i, names(pathprint.Hs.gs)))
# sample.POE.mayo_msbb[!(is.na(sample))]=fit$expr
# }
# pheatmap(POE.msbb.matrix[500:600,350:440],color=c("blue","white","red"),
#          cluster_rows=T,cluster_cols=T,
#          show_colnames=F,clustering_distance_rows="manhattan",
#          clustering_distance_cols="manhattan",
#          fontsize=9,fontsize_row=9,cell_width=8,cell_height=8,
#          legend_breaks=c(-1,0,1),clustering_method="average",
#          legend=T,
#          main="PathPrint Mt.Sinai 448 samples")
# ###########################################################
# # Annotate for Combined dataset
# levels(comm_mayo_msbb_annotation$Dataset)[which(rownames(comm_mayo_msbb_annotation)%in%rownames(msbb_data_controls))]="MSBB_Controls"
# levels(comm_mayo_msbb_annotation$Dataset)[279:726][-which(rownames(comm_mayo_msbb_annotation)[279:726]%in%rownames(msbb_data_controls))]="MSBB_Test"
# comm_mayo_msbb_annotation$Dataset[which(rownames(comm_mayo_msbb_annotation)%in%rownames(msbb_data_controls))]="MSBB_Controls"
# comm_mayo_msbb_annotation$Dataset[279:726][-which(rownames(comm_mayo_msbb_annotation)[279:726]%in%rownames(msbb_data_controls))]="MSBB_Test"
# comm_mayo_msbb_annotation$Dataset[279:726][-which(rownames(comm_mayo_msbb_annotation)[279:726]%in%rownames(msbb_data_controls))]="MSBB_Test"
# mayo_data_controls2=paste("X",mayo_data_controls2,sep="")
# levels(comm_mayo_msbb_annotation$Dataset)[which(rownames(comm_mayo_msbb_annotation)[1:278]%in%mayo_data_controls2)]="Mayo_Controls"
# 
# 
# 
# msbb_coexpp_enrichedReactomePathways=enrichMap(gene=c(msbb_coexpp[[1]]$ENTREZ_GENE_ID,msbb_coexpp[[2]]$ENTREZ_GENE_ID,msbb_coexpp[[3]]$ENTREZ_GENE_ID),organism="human",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.2,universe=synsy, minGSSize = 5,readable = FALSE)