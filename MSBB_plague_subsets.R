setwd("/home/md1wwxx/sharedmd1wwxx/MSBBCorrChange/data")
counts<-readRDS(file="syn5898491_counts.RDS")
covariates<-readRDS(file="syn6100548_covariates.RDS")
MSBB_clinical<-readRDS(file="syn6101474_MSBB_clinical.RDS")
setdiff(colnames(counts), covariates$sampleIdentifier)
setdiff(covariates$sampleIdentifier, colnames(counts))
covariates$sampleIdentifier<-sub("BM_22", "hB_RNA", covariates$sampleIdentifier)
covariates$sampleIdentifier<-sub("BM_10", "hB_RNA", covariates$sampleIdentifier)
covariates$sampleIdentifier<-sub("BM_36", "hB_RNA", covariates$sampleIdentifier)
covariates$sampleIdentifier<-sub("_resequenced", "", covariates$sampleIdentifier)
setdiff(colnames(counts), covariates$sampleIdentifier)
setdiff(covariates$sampleIdentifier, colnames(counts))
sampleInfo<-merge(covariates, MSBB_clinical, by.x="individualIdentifier", by.y="individualIdentifier")
unique(sampleInfo$BrodmannArea)
sampleInfo<-subset(sampleInfo, subset=(sampleIdentifier %in% colnames(counts)))
sampleInfo<-sampleInfo[!duplicated(sampleInfo$sampleIdentifier),]
nrow(sampleInfo)
table(sampleInfo$BrodmannArea)
BM10lowPlaque<-subset(sampleInfo, subset=(BrodmannArea=="BM10" & PlaqueMean<=1))
nrow(BM10lowPlaque)
length(intersect(BM10lowPlaque$sampleIdentifier, colnames(counts)))
BM10highPlaque<-subset(sampleInfo, subset=(BrodmannArea=="BM10" & PlaqueMean>=15))
nrow(BM10highPlaque)
length(intersect(BM10highPlaque$sampleIdentifier, colnames(counts)))

BM22lowPlaque<-subset(sampleInfo, subset=(BrodmannArea=="BM22" & PlaqueMean<=1))
nrow(BM22lowPlaque)
length(intersect(BM22lowPlaque$sampleIdentifier, colnames(counts)))
BM22highPlaque<-subset(sampleInfo, subset=(BrodmannArea=="BM22" & PlaqueMean>=15))
nrow(BM22highPlaque)
length(intersect(BM22highPlaque$sampleIdentifier, colnames(counts)))

BM36lowPlaque<-subset(sampleInfo, subset=(BrodmannArea=="BM36" & PlaqueMean<=1))
nrow(BM36lowPlaque)
length(intersect(BM36lowPlaque$sampleIdentifier, colnames(counts)))
BM36highPlaque<-subset(sampleInfo, subset=(BrodmannArea=="BM36" & PlaqueMean>=15))
nrow(BM36highPlaque)
length(intersect(BM36highPlaque$sampleIdentifier, colnames(counts)))

BM44lowPlaque<-subset(sampleInfo, subset=(BrodmannArea=="BM44" & PlaqueMean<=1))
nrow(BM44lowPlaque)
length(intersect(BM44lowPlaque$sampleIdentifier, colnames(counts)))
BM44highPlaque<-subset(sampleInfo, subset=(BrodmannArea=="BM44" & PlaqueMean>=15))
nrow(BM44highPlaque)
length(intersect(BM44highPlaque$sampleIdentifier, colnames(counts)))










