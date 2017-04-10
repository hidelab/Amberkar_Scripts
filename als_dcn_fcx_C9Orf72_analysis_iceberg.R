library(DESeq2)
library(org.Hs.eg.db)
library(igraph)
library(cocor)
library(clusterProfiler)
library(data.table)
library(parallel)

setwd("/shared/hidelab2/shared/ALS_Claire_Datasets/FCX_ALS_C9orf72_sALS/")
fcx_c9orf72_als_gene_counts_rawData=read.table("annotated_combined.counts",header = T,sep = "\t",row.names = 1,as.is = F)
fcx_c9orf72_als_gene_counts=data.frame(fcx_c9orf72_als_gene_counts_rawData[,c(28,grep(pattern = "C9",x = colnames(fcx_c9orf72_als_gene_counts_rawData)),grep(pattern = "control",x = colnames(fcx_c9orf72_als_gene_counts_rawData)))])
fcx_c9orf72_als_gene_counts.Control_Samples=grep(pattern = "control",x = colnames(fcx_c9orf72_als_gene_counts_rawData),value = T)
fcx_c9orf72_als_gene_counts.C9orf72_Samples=grep(pattern = "C9",x = colnames(fcx_c9orf72_als_gene_counts_rawData),value = T)
fcx_c9orf72_als_gene_counts.colData=data.frame(condition=c(rep("C9",length(fcx_c9orf72_als_gene_counts.C9orf72_Samples)),rep("Control",length(fcx_c9orf72_als_gene_counts.Control_Samples))),
                                           type="paired-end",row.names=colnames(fcx_c9orf72_als_gene_counts)[-1])
fcx_c9orf72_als_gene_counts.dds=DESeqDataSetFromMatrix(countData = fcx_c9orf72_als_gene_counts[,-1],colData = fcx_c9orf72_als_gene_counts.colData,design = ~condition)
fcx_c9orf72_als_gene_counts.dds=DESeq(fcx_c9orf72_als_gene_counts.dds)
fcx_c9orf72_als_gene_counts.normalised=data.frame(fcx_c9orf72_als_gene_counts_rawData$symbol,counts(fcx_c9orf72_als_gene_counts.dds,normalized=T),stringsAsFactors = F)
#Remove intronic transcripts and uncharacterised gene symbols
rm_genes=grep(pattern = "_|\\.|^RP|-",fcx_c9orf72_als_gene_counts.normalised$fcx_c9orf72_als_gene_counts_rawData.symbol,invert = F)
#Normalise counts based on neg binomial model in DESeq2
fcx_c9orf72_als_gene_counts.normalised=data.frame(fcx_c9orf72_als_gene_counts$symbol,counts(fcx_c9orf72_als_gene_counts.dds,normalized=T))
fcx_c9orf72_als_gene_counts.agg=aggregate(x = fcx_c9orf72_als_gene_counts.normalised[-rm_genes,-1],by=list(GeneSymbol=fcx_c9orf72_als_gene_counts.normalised$fcx_c9orf72_als_gene_counts.symbol[-rm_genes]),mean)
fcx_c9orf72_als_genes_counts.final=fcx_c9orf72_als_gene_counts.agg[,-1]
rownames(fcx_c9orf72_als_genes_counts.final)=fcx_c9orf72_als_gene_counts.agg$GeneSymbol
save.image("fcx_c9orf72_C9_genecounts_final.RData")

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
fcx_c9orf72_als_genes_counts.final_keep=rowSums(fcx_c9orf72_als_genes_counts.final >0) >=ncol(fcx_c9orf72_als_genes_counts.final)/3
exprs_rank<-fcx_c9orf72_als_genes_counts.final[fcx_c9orf72_als_genes_counts.final_keep,]
number_of_combinations<-choose(nrow(exprs_rank),2)
C9_samples=grep(pattern = "C9",colnames(exprs_rank))
c_exprs_rank=exprs_rank[,-C9_samples]
t_exprs_rank=exprs_rank[,C9_samples]
n.c<-ncol(c_exprs_rank)
n.t<-ncol(t_exprs_rank)
gene.names<-rownames(exprs_rank)

# input = 1:number_of_combinations
# pb = txtProgressBar(min=0,max=number_of_combinations,style=3,initial=0)
# cat("\n")
# res = mclapply(input,ProcessElement,mc.cores=ncore)
# close(pb)
# result <- rbindlist(res)
# result <- data.frame(result,stringsAsFactors = F)
# result$FDR<-p.adjust(result$p.cocor, method="fdr")
# result[,c(3:8,10:11)]=round(result[,c(3:8,10:11)],digits = 3)
# 
i<-478
start<-i*blocksize+1
end<-min((i+1)*blocksize, number_of_combinations)

while( start < number_of_combinations)
{
  input<-start:end
  pb = txtProgressBar(min=0,max=length(input),style=3,initial=0)
  cat("\n")
  res = mclapply(input,ProcessElement,mc.cores=ncore)
  close(pb)	
  # save results
  #saveRDS(res,paste0("/shared/hidelab2/user/md4zsa/Work/Data/ALS/fcx_c9orf72_als_", i, ".RDS"))
  #rbindlist is faster than rbind.fill
  result <- rbindlist(res)
  result <- as.data.frame(result)
  result <- data.frame(result,stringsAsFactors = F)
  #result$FDR<-p.adjust(result$p.cocor, method="fdr")
  #result[,c(3:8,10:11)]=round(result[,c(3:8,10:11)],digits = 3)
  write.table(result, file=paste0("/shared/hidelab2/user/md4zsa/Work/Data/ALS/C9orf72/txt/fcx_c9orf72_als_", i, ".txt"), sep="\t", row.names=FALSE, quote = FALSE)
  i<-i+1
  start<-i*blocksize+1
  end<-min((i+1)*blocksize, number_of_combinations)
}
cat(paste("Done!"))
