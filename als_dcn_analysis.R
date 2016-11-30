library(DESeq2)
library(org.Hs.eg.db)
library(igraph)
library(cocor)
library(clusterProfiler)
library(data.table)
library(parallel)

setwd("/Users/sandeepamberkar/Work/Data/ALS_Claire_Datasets/2016-03-07_tdp43_project/")
fcx_als_gene_counts_rawData=read.table("annotated_combined.counts",header = T,sep = "\t",row.names = 1,as.is = F)
fcx_als_gene_tpm_rawData=read.table("combined.gene.sf.tpm",header = T,sep = "\t",row.names = 1,as.is = F)

fcx_als_gene_counts=data.frame(fcx_als_gene_counts_rawData[,c(28,grep(pattern = "C9",x = colnames(fcx_als_gene_counts_rawData)),grep(pattern = "control",x = colnames(fcx_als_gene_counts_rawData)))])
fcx_als_gene_counts.Control=grep(pattern = "control",x = colnames(fcx_als_gene_counts_rawData),value = T)
fcx_als_gene_counts.C9=grep(pattern = "C9",x = colnames(fcx_als_gene_counts_rawData),value = T)
fcx_als_gene_counts.colData=data.frame(condition=c(rep("C9",8),rep("Control",9)),type="paired-end",row.names=colnames(fcx_als_gene_counts)[-1])
fcx_als_gene_counts.dds=DESeqDataSetFromMatrix(countData = fcx_als_gene_counts[,-1],colData = fcx_als_gene_counts.colData,design = ~condition)
fcx_als_gene_counts.dds=DESeq(fcx_als_gene_counts.dds)
#Remove intronic transcripts and uncharacterised gene symbols
rm_genes=grep(pattern = "_|\\.|^RP|-",fcx_als_gene_counts.normalised$fcx_als_gene_counts.symbol,invert = F)
#Normalise counts based on neg binomial model in DESeq2
fcx_als_gene_counts.normalised=data.frame(fcx_als_gene_counts$symbol,counts(fcx_als_gene_counts.dds,normalized=T))
fcx_als_gene_counts.agg=aggregate(x = fcx_als_gene_counts.normalised[-rm_genes,-1],by=list(GeneSymbol=fcx_als_gene_counts.normalised$fcx_als_gene_counts.symbol[-rm_genes]),mean)
fcx_als_genes_counts.final=fcx_als_gene_counts.agg[,-1]
rownames(fcx_als_genes_counts.final)=fcx_als_gene_counts.agg$GeneSymbol
saveRDS(object = fcx_als_genes_counts.final,"FCX_ALS_Gene_Counts_Final.RDS")
ProcessElement <- function(ic){
  A = ceiling((sqrt(8*(ic+1)-7)+1)/2)
  B = ic-choose(floor(1/2+sqrt(2*ic)),2)
  
  c_A = c_exprs_rank[A,]
  c_B = c_exprs_rank[B,]
  
  t_A = t_exprs_rank[A,]
  t_B = t_exprs_rank[B,]
  
  # get correlation between the summaries for the unique genes
  tmp = data.frame(Gene.A=gene.names[A],Gene.B=gene.names[B])
  c_cortest<-cor.test(unname(unlist(c_A)), unname(unlist(c_B)), method="spearman")
  t_cortest<-cor.test(unname(unlist(t_A)), unname(unlist(t_B)), method="spearman")
  rc<-c_cortest$estimate
  rt<-t_cortest$estimate
  diffcor<-cocor.indep.groups(rc, rt, n.c, n.t)
  tmp$r.c<-rc
  tmp$p.c<-c_cortest$p.value
  tmp$n.c<-n.c
  tmp$r.t<-rt
  tmp$p.t<-t_cortest$p.value
  tmp$n.t<-n.t
  tmp$p.cocor<-diffcor@fisher1925$p.value
  tmp$abs.corr.change<-abs(rt-rc)
  
  setTxtProgressBar(pb,ic)
  return(tmp)
}

fcx_als_genes_counts.deg=data.frame(results(fcx_als_gene_counts.dds),stringsAsFactors = F)
fcx_als_genes_counts.degGenes=select(x = org.Hs.eg.db,keys = rownames(fcx_als_genes_counts.deg[which(fcx_als_genes_counts.deg$padj<=0.05),]),columns = "SYMBOL",keytype = "ENSEMBL")[,2]

exprs_rank<-fcx_als_genes_counts.final[which(rownames(fcx_als_genes_counts.final)%in%na.omit(fcx_als_genes_counts.degGenes)[1:length(na.omit(fcx_als_genes_counts.degGenes))]),]
number_of_combinations<-choose(nrow(exprs_rank),2)
c_exprs_rank=exprs_rank[,which(colnames(exprs_rank)%in%fcx_als_gene_counts.Control)]
t_exprs_rank=exprs_rank[,which(colnames(exprs_rank)%in%fcx_als_gene_counts.C9)]
n.c<-ncol(c_exprs_rank)
n.t<-ncol(t_exprs_rank)
gene.names<-rownames(exprs_rank)

input = 1:number_of_combinations
pb = txtProgressBar(min=0,max=number_of_combinations,style=3,initial=0)
cat("\n")
res = mclapply(input,ProcessElement,mc.cores=nc)
close(pb)
