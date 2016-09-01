#!/bin/bash
#$ -q long   # <- the name of the Q you want to submit to
#$ -rmem 32G # <- Real memory request
#$ -mem 64G # <- Virtual memory request
#$ -binding linear:16 #  <- seek 16 cores
#$ -S /bin/bash   # <- run the job under bash
#$ -N TDP43_RNAseq # <- name of the job in the qstat output
#$ -o TDP43_RNAseq.out # <- name of the output file.
#$ -e TDP43_RNAseq.stderr # <- name of the stderr file.
library("doParallel")
library("Rsubread")
library("foreach")


cl=makeCluster(16)
registerDoParallel(cl)
registerDoParallel(cores=16)
data_dir="/fastdata/md4zsa/TDP_Omics_Study/Data/"
out_dir="/fastdata/md4zsa/TDP_Omics_Study/Results"
index_path="/data/md4zsa/TDP_Omics_Study/Results/Rsubread_indices/hs_hg38"
chrs=c(seq(1,22,by=1),"MT","X","Y")
ens_genes=read.table(paste(data_dir,"Ensembl82_GeneAnnotation.txt",sep="/"),sep="\t",header = T,as.is=T)
ens_transcript=read.table(paste(data_dir,"Ensembl82_TranscriptAnnotation.txt",sep="/"),sep="\t",header = T,as.is=T)
ens_genes=ens_genes[which(ens_genes$Chr%in%chrs),]
ens_transcript=ens_transcript[which(ens_transcript$Chr%in%chrs),]
ens_genes=ens_genes[,-6]
ens_transcript=ens_transcript[,-6]
colnames(ens_genes)=c("GeneID","Chr","Start","End","Strand")
colnames(ens_transcript)=c("GeneID","Chr","Start","End","Strand")

#Check if reading from readfiles of the same sample
fname=list()
fc_TDP=list()
R1_files=list.files(path=paste(data_dir,"R1",sep="/"),full.names=T)
R2_files=list.files(path=paste(data_dir,"R2",sep="/"),full.names=T)
R1R2_files=cbind(R1_files,R2_files)
parts=list(1:10,11:20,21:30,31:40)
#Check if reading from readfiles of the same sample
for (r in 1:dim(R1R2_files)[1]){
  cat(paste("Checking file integrity for ",R1R2_files[r],"...\n",sep=""))
  if(strsplit(basename(R1R2_files[r,]),split="_")[[1]][1]==strsplit(basename(R1R2_files[r,]),split="_")[[2]][1]){
    fname[[r]]=strsplit(basename(R1R2_files[r,]),split="_")[[1]][1]
  }
}

foreach(p=parts) %dopar% {
  align(index=index_path,
        readfile1=R1_files[p],
        readfile2=R2_files[p],
        nthreads=32,
        output_file=paste(out_dir,"/",fname[p],".","bam",sep=""),
        phredOffset=64)
}

#Move .indel files generated from Rsubread to separate directory

results_bam=system(paste("find ",out_dir,"/*.bam",sep=""),intern = T)
for (i in 1:length(results_bam)) {
     paste("Counting for sample ",results_bam[i],"...",sep="")
     fc_TDP_genes[[i]]=featureCounts(results_bam[i],annot.ext=hg38_cds,nthreads=32,isPairedEnd=T,countMultiMappingReads = T,requireBothEndsMapped = T)
     write.table(fc_TDP[[i]]$counts,paste(out_dir,"/GeneCounts/",fname[i],".counts",sep=""),sep="\t",col.names=T,row.names=T,quote=F)
     write.table(fc_TDP[[i]]$annotation,paste(out_dir,"/Annotation/",fname[i],".annotation",sep=""),sep="\t",col.names=T,row.names=T,quote=F)
   }
paste("Alignment and featureCounts calculation completed in proc.time()[3]/60 mins.")

