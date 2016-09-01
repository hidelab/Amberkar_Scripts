msmm_rnaseq_b1_Braak=msmm_rnaseq_b1_AOD=msmm_rnaseq_b1_Braak.raw=msmm_rnaseq_b1_AOD.raw=list()
msmm_rnaseq_b1_raw=read.table("msmm_rnaseq_B1_raw_counts.txt",header = T,sep = "\t",as.is = T)
rownames(msmm_rnaseq_b1_raw)=msmm_rnaseq_b1_raw$GeneId
msmm_rnaseq_b1_raw=msmm_rnaseq_b1_raw[,-1]
msmm_rnaseq_b1=read.csv("msmm_covar_corrected_B1.csv",header = T,row.names = 1,as.is = T)
msmm_rnaseq_covariates=read.delim2("msmm_rnaseq_B1_covariates.tsv",header = T,sep = "\t",row.names = 1,as.is = T)
msmm_rnaseq_covariates$AOD[which(msmm_rnaseq_covariates$AOD=="89+")]=90
msmm_rnaseq_covariates$AOD=as.integer(msmm_rnaseq_covariates$AOD)

msmm_rnaseq_b1_AOD[[1]]=msmm_rnaseq_b1[,which(colnames(msmm_rnaseq_b1)%in%rownames(msmm_rnaseq_covariates)[which(msmm_rnaseq_covariates$AOD>=61 & msmm_rnaseq_covariates$AOD<=70)])]
msmm_rnaseq_b1_AOD[[2]]=msmm_rnaseq_b1[,which(colnames(msmm_rnaseq_b1)%in%rownames(msmm_rnaseq_covariates)[which(msmm_rnaseq_covariates$AOD>=71 & msmm_rnaseq_covariates$AOD<=80)])]
msmm_rnaseq_b1_AOD[[3]]=msmm_rnaseq_b1[,which(colnames(msmm_rnaseq_b1)%in%rownames(msmm_rnaseq_covariates)[which(msmm_rnaseq_covariates$AOD>=81 & msmm_rnaseq_covariates$AOD<=90)])]

msmm_rnaseq_b1_AOD.raw[[1]]=msmm_rnaseq_b1_raw[,which(colnames(msmm_rnaseq_b1_raw)%in%rownames(msmm_rnaseq_covariates)[which(msmm_rnaseq_covariates$AOD>=61 & msmm_rnaseq_covariates$AOD<=70)])]
msmm_rnaseq_b1_AOD.raw[[2]]=msmm_rnaseq_b1_raw[,which(colnames(msmm_rnaseq_b1_raw)%in%rownames(msmm_rnaseq_covariates)[which(msmm_rnaseq_covariates$AOD>=71 & msmm_rnaseq_covariates$AOD<=80)])]
msmm_rnaseq_b1_AOD.raw[[3]]=msmm_rnaseq_b1_raw[,which(colnames(msmm_rnaseq_b1_raw)%in%rownames(msmm_rnaseq_covariates)[which(msmm_rnaseq_covariates$AOD>=81 & msmm_rnaseq_covariates$AOD<=90)])]

msmm_rnaseq_b1_Braak.raw[[1]]=msmm_rnaseq_b1_raw[,which(colnames(msmm_rnaseq_b1_raw)%in%rownames(msmm_rnaseq_covariates)[which(msmm_rnaseq_covariates$bbscore==0)])]
msmm_rnaseq_b1_Braak.raw[[3]]=msmm_rnaseq_b1_raw[,which(colnames(msmm_rnaseq_b1_raw)%in%rownames(msmm_rnaseq_covariates)[which(msmm_rnaseq_covariates$bbscore>=5)])]
msmm_rnaseq_b1_Braak.raw[[2]]=msmm_rnaseq_b1_raw[,which(colnames(msmm_rnaseq_b1_raw)%in%rownames(msmm_rnaseq_covariates)[which(msmm_rnaseq_covariates$bbscore>=1 & msmm_rnaseq_covariates$bbscore<=4)])]

msmm_rnaseq_b1_Braak[[1]]=msmm_rnaseq_b1[,which(colnames(msmm_rnaseq_b1)%in%rownames(msmm_rnaseq_covariates)[which(msmm_rnaseq_covariates$bbscore==0)])]
msmm_rnaseq_b1_Braak[[3]]=msmm_rnaseq_b1[,which(colnames(msmm_rnaseq_b1)%in%rownames(msmm_rnaseq_covariates)[which(msmm_rnaseq_covariates$bbscore>=5)])]
msmm_rnaseq_b1_Braak[[2]]=msmm_rnaseq_b1[,which(colnames(msmm_rnaseq_b1)%in%rownames(msmm_rnaseq_covariates)[which(msmm_rnaseq_covariates$bbscore>=1 & msmm_rnaseq_covariates$bbscore<=4)])]

names(msmm_rnaseq_b1_AOD)=names(msmm_rnaseq_b1_AOD.raw)=c("msmm_rnaseq_b1_AOD1","msmm_rnaseq_b1_AOD2","msmm_rnaseq_b1_AOD3")
names(msmm_rnaseq_b1_Braak)=names(msmm_rnaseq_b1_Braak.raw)=c("msmm_rnaseq_b1_noBraak","msmm_rnaseq_b1_modBraak","msmm_rnaseq_b1_hiBraak")
for (i in 1:3)
  write.table(msmm_rnaseq_b1_AOD[[i]],paste("MSMM_RNAseq_AOD/",names(msmm_rnaseq_b1_AOD)[i],".txt",sep=""),sep="\t",col.names = F,row.names = F,quote = F)
  write.table(msmm_rnaseq_b1_AOD.raw[[i]],paste("MSMM_RNAseq_AOD/",names(msmm_rnaseq_b1_AOD.raw)[i],"_raw.txt",sep=""),sep="\t",col.names = F,row.names = F,quote = F)
  write.table(msmm_rnaseq_b1_Braak[[i]],paste("MSMM_RNAseq_Braak/",names(msmm_rnaseq_b1_Braak)[i],".txt",sep=""),sep="\t",col.names = F,row.names = F,quote = F)
  write.table(msmm_rnaseq_b1_Braak.raw[[i]],paste("MSMM_RNAseq_Braak/",names(msmm_rnaseq_b1_Braak.raw)[i],"_raw.txt",sep=""),sep="\t",col.names = F,row.names = F,quote = F)

for (j in 1:3)
  write.table(msmm_rnaseq_b1_AOD[[j]][1:5000,],paste("MSMM_RNAseq_AOD/aod_data5k/",names(msmm_rnaseq_b1_AOD)[j],".txt",sep=""),sep="\t",col.names = F,row.names = F,quote = F)
  write.table(msmm_rnaseq_b1_AOD[[j]][1:5000,],paste("MSMM_RNAseq_AOD/aod_data5k/",names(msmm_rnaseq_b1_AOD)[j],".txt",sep=""),sep="\t",col.names = F,row.names = F,quote = F)
  write.table(msmm_rnaseq_b1_Braak[[j]][1:5000,],paste("MSMM_RNAseq_Braak/braak_data5k/",names(msmm_rnaseq_b1_AOD)[j],".txt",sep=""),sep="\t",col.names = F,row.names = F,quote = F)

write(rownames(msmm_rnaseq_b1_AOD[[1]]),"MSMM_RNAseq_AOD/aod_data5k/geneid.txt",sep="\n")
write(rownames(msmm_rnaseq_b1_Braak[[1]]),"MSMM_RNAseq_Braak/braak_data5k/geneid.txt",sep="\n")
write(rownames(msmm_rnaseq_b1_AOD.raw[[1]]),"MSMM_RNAseq_AOD/aod_data5k_raw/geneid.txt",sep="\n")
write(rownames(msmm_rnaseq_b1_Braak.raw[[1]]),"MSMM_RNAseq_Braak/braak_data5k_raw/geneid.txt",sep="\n")

system("sed -n '1,5001p' MSMM_RNAseq_AOD/msmm_rnaseq_b1_AOD1.txt > MSMM_RNAseq_AOD/aod_data5k/msmm_rnaseq_b1_AOD1_5kgenes.txt")
system("sed -n '1,5001p' MSMM_RNAseq_AOD/msmm_rnaseq_b1_AOD2.txt > MSMM_RNAseq_AOD/aod_data5k/msmm_rnaseq_b1_AOD2_5kgenes.txt")
system("sed -n '1,5001p' MSMM_RNAseq_AOD/msmm_rnaseq_b1_AOD3.txt > MSMM_RNAseq_AOD/aod_data5k/msmm_rnaseq_b1_AOD3_5kgenes.txt")

write(system("basename $(find MSMM_RNAseq_AOD/aod_data5k -type f -name '*genes.txt'|grep '5k')",intern = T),"MSMM_RNAseq_AOD/aod_data5k/datasets_lists.txt",sep="\n")
system("more MSMM_RNAseq_AOD/aod_data5k/datasets_lists.txt")





