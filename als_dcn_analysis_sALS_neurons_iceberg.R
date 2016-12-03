library(DESeq2)
library(org.Hs.eg.db)
library(igraph)
library(cocor)
library(clusterProfiler)
library(data.table)
library(parallel)

setwd("/shared/hidelab2/shared/ALS_Claire_Datasets/Neurons_ALS_sALS")
neurons_sals_gene_counts_rawData=read.table("annotated_combined.counts",header = T,sep = "\t",row.names = 1,as.is = F)
neurons_sals_gene_counts=data.frame(neurons_sals_gene_counts_rawData[,c(22,grep(pattern = "sALS",x = colnames(neurons_sals_gene_counts_rawData)),grep(pattern = "control",x = colnames(neurons_sals_gene_counts_rawData)))])
neurons_sals_gene_counts.Control=grep(pattern = "control",x = colnames(neurons_sals_gene_counts_rawData),value = T)
neurons_sals_gene_counts.sALS=grep(pattern = "sALS",x = colnames(neurons_sals_gene_counts_rawData),value = T)
neurons_sals_gene_counts.colData=data.frame(condition=c(rep("sALS",length(neurons_sals_gene_counts.sALS)),rep("Control",length(neurons_sals_gene_counts.Control))),
                                       type="paired-end",row.names=colnames(neurons_sals_gene_counts)[-1])
neurons_sals_gene_counts.dds=DESeqDataSetFromMatrix(countData = neurons_sals_gene_counts[,-1],colData = neurons_sals_gene_counts.colData,design = ~condition)
neurons_sals_gene_counts.dds=DESeq(neurons_sals_gene_counts.dds)
neurons_sals_gene_counts.normalised=data.frame(neurons_sals_gene_counts_rawData$symbol,counts(neurons_sals_gene_counts.dds,normalized=T),stringsAsFactors = F)
#Remove intronic transcripts and uncharacterised gene symbols
rm_genes=grep(pattern = "_|\\.|^RP|-",neurons_sals_gene_counts.normalised$neurons_sals_gene_counts_rawData.symbol,invert = F)
#Normalise counts based on neg binomial model in DESeq2
neurons_sals_gene_counts.normalised=data.frame(neurons_sals_gene_counts$symbol,counts(neurons_sals_gene_counts.dds,normalized=T))
neurons_sals_gene_counts.agg=aggregate(x = neurons_sals_gene_counts.normalised[-rm_genes,-1],by=list(GeneSymbol=neurons_sals_gene_counts.normalised$neurons_sals_gene_counts.symbol[-rm_genes]),mean)
neurons_sals_genes_counts.final=neurons_sals_gene_counts.agg[,-1]
rownames(neurons_sals_genes_counts.final)=neurons_sals_gene_counts.agg$GeneSymbol
save.image("Neurons_sALS_genecounts_final.RData")

ncore = detectCores()
blocksize=100000
ProcessElement <- function(ic){
  A = ceiling((sqrt(8*(ic+1)-7)+1)/2)
  B = ic-choose(floor(1/2+sqrt(2*ic)),2)
  
  c_A = c_exprs_rank[A,]
  c_B = c_exprs_rank[B,]
  
  t_A = t_exprs_rank[A,]
  t_B = t_exprs_rank[B,]
  
  # get correlation between the summaries for the unique genes
  tmp = data.frame(Gene.A=gene.names[A],Gene.B=gene.names[B])
  c_cortest<-cor.test(unname(unlist(c_A)), unname(unlist(c_B)), method="pearson")
  t_cortest<-cor.test(unname(unlist(t_A)), unname(unlist(t_B)), method="pearson")
  rc<-unname(cor.test(unlist(c_A),unlist(c_B))$estimate)
  rt<-unname(cor.test(unlist(t_A),unlist(t_B))$estimate)
  # diffcor<-cocor.indep.groups(rc, rt, n.c, n.t)
  tmp$r.c<-rc
  tmp$p.c<-c_cortest$p.value
  tmp$n.c<-n.c
  tmp$r.t<-rt
  tmp$p.t<-t_cortest$p.value
  tmp$n.t<-n.t
  tmp$p.cocor<-NA
  tmp$abs.corr.change<-abs(rt-rc)
  if ( (!is.na(tmp$r.c)) && (!is.na(tmp$r.t)) )
  {
    diffcor<-cocor.indep.groups(tmp$r.c, tmp$r.t, tmp$n.c, tmp$n.t)
    tmp$p.cocor<-diffcor@fisher1925$p.value
  }
  setTxtProgressBar(pb,ic %% blocksize)
  return(tmp)
}
neurons_sals_genes_counts.final_keep=rowSums(neurons_sals_genes_counts.final >0) >=ncol(neurons_sals_genes_counts.final)/3
exprs_rank<-neurons_sals_genes_counts.final[neurons_sals_genes_counts.final_keep,]
number_of_combinations<-choose(nrow(exprs_rank),2)
sALS_samples=grep(pattern = "sALS",colnames(exprs_rank))
c_exprs_rank=exprs_rank[,-sALS_samples]
t_exprs_rank=exprs_rank[,sALS_samples]
n.c<-ncol(c_exprs_rank)
n.t<-ncol(t_exprs_rank)
gene.names<-rownames(exprs_rank)

i<-0
start<-i*blocksize+1
end<-min((i+1)*blocksize, number_of_combinations)
while( start < number_of_combinations)
{
  input<-start:end
  pb = txtProgressBar(min=0,max=length(input),style=3,initial=0)
  cat("\n")
  res = mclapply(input,ProcessElement,mc.cores=ncore)
  close(pb)	
  result <- rbindlist(res)
  result <- as.data.frame(result)
  result <- data.frame(result,stringsAsFactors = F)
  write.table(result, file=paste0("/shared/hidelab2/user/md4zsa/Work/Data/ALS/sALS_Neurons_DiffCorr/txt/neurons_sals_", i, ".t
xt"), sep="\t", row.names=FALSE, quote = FALSE)
  i<-i+1
  start<-i*blocksize+1
  end<-min((i+1)*blocksize, number_of_combinations)
}
cat(paste("Done!"))





