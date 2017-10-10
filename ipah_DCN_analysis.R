library(data.table)
library(cocor)
library(parallel)
library(org.Hs.eg.db)
ncore = detectCores()
mapIds2<-function(IDs,IDFrom,IDTo){
  idmap=mapIds(x = org.Hs.eg.db,keys = IDs,column = IDTo,keytype = IDFrom,multiVals = "first")
  na_vec=names(idmap[is.na(idmap)==T])
  idmap=idmap[is.na(idmap)==F]
  idmap_df=data.frame("From"=names(idmap),"To"=unlist(unname(idmap)),stringsAsFactors = F)
  return(list(map=idmap_df,noMap=na_vec))
}
ProcessElement <- function(ic){
  A = ceiling((sqrt(8*(ic+1)-7)+1)/2)
  B = ic-choose(floor(1/2+sqrt(2*ic)),2)
  
  c_A = as.numeric(c_counts[A,])
  c_B = as.numeric(c_counts[B,])
  
  t_A = as.numeric(t_counts[A,])
  t_B = as.numeric(t_counts[B,])
  
  # get correlation between the summaries for the unique genes
  tmp = data.frame(IC= ic, Gene.A=gene.names[A], Gene.B=gene.names[B], r.c=NA, p.c=NA, n.c=NA, r.t=NA, p.t=NA, n.t=NA, p.cocor=NA)
  #if( (var(c_A) * var(c_B))!=0)
  tmp$n.c<-sum(!is.na(c_A + c_B))
  if (tmp$n.c >=10)
  {
    c_cortest<-cor.test(c_A, c_B, method="spearman")
    tmp$r.c<-c_cortest$estimate
    tmp$p.c<-c_cortest$p.value
  }
  
  #if( (var(t_A) * var(t_B))!=0)
  tmp$n.t<-sum(!is.na(t_A + t_B))
  if(tmp$n.t >=10)
  {
    t_cortest<-cor.test(t_A, t_B, method="spearman")
    tmp$r.t<-t_cortest$estimate
    tmp$p.t<-t_cortest$p.value
  }
  
  if ( (!is.na(tmp$r.c)) && (!is.na(tmp$r.t)) )
  {
    diffcor<-cocor.indep.groups(tmp$r.c, tmp$r.t, tmp$n.c, tmp$n.t)
    tmp$p.cocor<-diffcor@fisher1925$p.value
  }
  
  #setTxtProgressBar(pb,ic %% n_part)
  setTxtProgressBar(pb,ic %% blocksize)
  return(tmp)
}

setwd("/Users/sandeepamberkar/Work/Data/IPAH/")
ipah_metadata=readRDS("/Users/sandeepamberkar/Work/Data/IPAH/lawrie_sample_group.RDS")
count_files=system("find /Users/sandeepamberkar/Work/Data/IPAH/result -name 'quant.sf'",intern=T)
ipah_curated=read.csv("PAH Targets Bertero 2014 and AL.csv",header = T,as.is = T)
ipah_metadata$External.ID[which(ipah_metadata$group=="HV")]=gsub(pattern="_v1",replacement="",x=ipah_metadata$External.ID[which(ipah_metadata$group=="HV")])
healthy_samples=grep(paste(ipah_metadata$External.ID[which(ipah_metadata$group=="HV")],collapse="|"),count_files,value=T)
ipah_samples=grep(paste(ipah_metadata$External.ID[which(ipah_metadata$group=="IPAH")],collapse="|"),count_files,value=T)
healthy_counts.df=data.frame(lapply(lapply(healthy_samples,function(x)fread(input=x,sep="\t",header=T,stringsAsFactors=F,showProgress=T,data.table=F)),`[[`,4),stringsAsFactors=F)
ipah_counts.df=data.frame(lapply(lapply(ipah_samples,function(x)fread(input=x,sep="\t",header=T,stringsAsFactors=F,showProgress=T,data.table=F)),`[[`,4),stringsAsFactors=F)
colnames(healthy_counts.df)=unlist(lapply(strsplit(x=healthy_samples,split="/"),`[[`,8))
colnames(ipah_counts.df)=unlist(lapply(strsplit(x=ipah_samples,split="/"),`[[`,8))
f1=fread(count_files[1],sep="\t",header=T,stringsAsFactors=F,data.table=F)
rownames(ipah_counts.df)=f1$Name
rownames(healthy_counts.df)=f1$Name

ipah_colData=matrix(NA,sum(length(healthy_samples),length(ipah_samples)),ncol=2)
ipah_colData[,1]=c(rep("HV",length(healthy_samples)),rep("IPAH",length(ipah_samples)))
ipah_colData[,2]=rep("paired-end",sum(length(healthy_samples),length(ipah_samples)))
rownames(ipah_colData)=c(colnames(healthy_counts.df),colnames(ipah_counts.df))
ipah_colData=data.frame(condition=ipah_colData[,1],type=ipah_colData[,2],stringsAsFactors = F)
ipah_all_samples.df=data.frame(healthy_counts.df,ipah_counts.df,stringsAsFactors = F)
ipah_all_samples.df=data.frame(apply(ipah_all_samples.df,2,round,digits=0))
colnames(ipah_all_samples.df)=gsub(pattern = "^X",replacement = "",x = colnames(ipah_all_samples.df))
rownames(ipah_all_samples.df)=f1$Name

ipah_all_samples.dds=DESeqDataSetFromMatrix(countData = ipah_all_samples.df,colData = ipah_colData,design = ~condition)
ipah_all_samples.dds=ipah_all_samples.dds[rowSums(counts(ipah_all_samples.dds))>1,]
ipah_all_samples.dds=DESeq(ipah_all_samples.dds)
ipah_all_samples.results=data.frame(results(ipah_all_samples.dds)[which((results(ipah_all_samples.dds)[,6]<=0.05)&(abs(results(ipah_all_samples.dds)[,2]>log2(2)))),],stringsAsFactors = F)
ipah_all_samples.results=ipah_all_samples.results[order(ipah_all_samples.results$log2FoldChange,decreasing = T),]
ipah_DETs=unlist(lapply(strsplit(rownames(ipah_all_samples.results[which(((abs(ipah_all_samples.results$log2FoldChange)>log(2))==T)&(ipah_all_samples.results$padj<=0.05)),]),split = "\\."),`[[`,1))
ipah_all_samples.dds=DESeqDataSetFromMatrix(countData = ipah_all_samples.df,colData = ipah_colData,design = ~condition)
ipah_all_samples.dds=ipah_all_samples.dds[rowSums(counts(ipah_all_samples.dds))>1,]
ipah_all_samples.dds=DESeq(ipah_all_samples.dds)
ipah_all_samples.results=data.frame(results(ipah_all_samples.dds)[which((results(ipah_all_samples.dds)[,6]<=0.05)&(abs(results(ipah_all_samples.dds)[,2]>log2(1.5)))),],stringsAsFactors = F)
ipah_all_samples.results[,c(1:6)]=round(ipah_all_samples.results[,c(1:6)],digits = 3)
ipah_all_samples.results=ipah_all_samples.results[order(ipah_all_samples.results$log2FoldChange,decreasing = T),]
ipah_DETs=unlist(lapply(strsplit(rownames(ipah_all_samples.results[which(((abs(ipah_all_samples.results$log2FoldChange)>log(1.5))==T)&(ipah_all_samples.results$padj<=0.05)),]),split = "\\."),`[[`,1))
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ipah_all_samples.results$external_transcript_name=getBM(attributes=c("external_transcript_name","ensembl_transcript_id","ensembl_gene_id"),filters = c("ensembl_transcript_id"),values = unlist(lapply(strsplit(rownames(ipah_all_samples.results),split = "\\."),`[[`,1)),mart = ensembl)[,1]
ipah_all_samples.results$Symbol=getBM(attributes=c("external_transcript_name","ensembl_transcript_id","ensembl_gene_id","hgnc_symbol","description"),filters = c("ensembl_transcript_id"),values = unlist(lapply(strsplit(rownames(ipah_all_samples.results),split = "\\."),`[[`,1)),mart = ensembl)[,4]
ipah_all_samples.results$Desc=getBM(attributes=c("external_transcript_name","ensembl_transcript_id","ensembl_gene_id","hgnc_symbol","description"),filters = c("ensembl_transcript_id"),values = unlist(lapply(strsplit(rownames(ipah_all_samples.results),split = "\\."),`[[`,1)),mart = ensembl)[,5]


