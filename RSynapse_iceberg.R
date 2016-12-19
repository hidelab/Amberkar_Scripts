#module load apps/R
#source('http://depot.sagebase.org/CRAN.R')
#pkgInstall("synapseClient")
#pkgInstall("Rsftp")
#pkgInstall("RCurl")
#pkgInstall("rjson")
#pkgInstall("digest")
#pkgInstall("RUnit")

require(synapseClient)
synapseLogin('W.Wei@Sheffield.ac.uk', 'Synapse2491753')
#syn5898491
target<-"syn5898491"
synGet(id= target, version=NULL, downloadFile=T, downloadLocation='/home/md1wwxx/sharedmd1wwxx/MSBBCorrChange/data', ifcollision="keep.both", load=F)
setwd("/home/md1wwxx/sharedmd1wwxx/MSBBCorrChange/data")
dir()
counts<-read.table("AMP-AD_MSBB_MSSM_IlluminaHiSeq2500_final.TMM_normalized.race_age_RIN_PMI_batch_corrected.tsv", header = TRUE,skip=0, row.names=NULL, sep="\t", stringsAsFactors=FALSE,na.strings=c("", " "), quote="",  comment.char ="")
saveRDS(counts, file="syn5898491_counts.RDS")

#syn6100548
target<-"syn6100548"
synGet(id= target, version=NULL, downloadFile=T, downloadLocation='/home/md1wwxx/sharedmd1wwxx/MSBBCorrChange/data', ifcollision="keep.both", load=F)
dir()
covariates<-read.table("MSBB_RNAseq_covariates.csv", header = TRUE,skip=0, row.names=NULL, sep=",", stringsAsFactors=FALSE,na.strings=c("", " "), quote="",  comment.char ="")
saveRDS(covariates, file="syn6100548_covariates.RDS")

target<-"syn6101474"
synGet(id= target, version=NULL, downloadFile=T, downloadLocation='/home/md1wwxx/sharedmd1wwxx/MSBBCorrChange/data', ifcollision="keep.both", load=F)
dir()
MSBB_clinical<-read.table("MSBB_clinical.csv", header = TRUE,skip=0, row.names=NULL, sep=",", stringsAsFactors=FALSE,na.strings=c("", " "), quote="",  comment.char ="")
saveRDS(MSBB_clinical, file="syn6101474_MSBB_clinical.RDS")







