library(igraph)
library(org.Hs.eg.db)
library(cocor)
library(data.table)


setwd("/shared/hidelab2/user/md4zsa/Work/Data/MSMM_RNAseq/MSMM_RNAseq_FinalRelease/")
cat(paste("Reading data ...\n"))
load("msmm_rnaseq_FinalDataset_contPLQ.RData")
cat(paste("Applying significance threshold of 0.05 ...\n"))

msmm_rnaseq.cocor=vector(mode = "list",length = 3)
msmm_rnaseq.cocor[[1]]=msmm_rnaseq.cocor[[2]]=msmm_rnaseq.cocor[[3]]=vector(mode = "list",length = 2)
names(msmm_rnaseq.cocor)=names(msmm_rnaseq)
names(msmm_rnaseq.cocor[[1]])=names(msmm_rnaseq.cocor[[2]])=names(msmm_rnaseq.cocor[[3]])

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
nc=16
msmm_rnaseq_Plaque.Corr005.Genes=lapply(msmm_rnaseq_Plaque.Corr,function(x)x[which(x$Pval<=0.05),1])
msmm_rnaseq_Plaque.Corr001.Genes=lapply(msmm_rnaseq_Plaque.Corr,function(x)x[which(x$Pval<=0.01),1])
for (t in 1:length(names(msmm_rnaseq))){
  
  exprs_rank<-msmm_rnaseq[[t]][which(rownames(msmm_rnaseq[[t]])%in%msmm_rnaseq_Plaque.Corr005.Genes[[t]]),]
  number_of_combinations<-choose(nrow(exprs_rank),2)
  c_exprs_rank=exprs_rank[,which(colnames(exprs_rank)%in%low_Plaque_Samples[[t]])]
  t_exprs_rank=exprs_rank[,which(colnames(exprs_rank)%in%high_Plaque_Samples[[t]])]
  n.c<-ncol(c_exprs_rank)
  n.t<-ncol(t_exprs_rank)
  gene.names<-rownames(exprs_rank)
  
  input = 1:number_of_combinations
  pb = txtProgressBar(min=0,max=number_of_combinations,style=3,initial=0)
  cat("\n")
  res = mclapply(input,ProcessElement,mc.cores=nc)
  close(pb)
  
  result <- rbindlist(res)
  result <- data.frame(result,stringsAsFactors = F)
  result$FDR<-p.adjust(result$p.cocor, method="fdr")
  result[,c(3:8,10:11)]=round(result[,c(3:8,10:11)],digits = 3)
  msmm_rnaseq.cocor[[t]]=result
}

msmm_rnaseq.cocor_filtered=vector(mode = "list",length = 3)
names(msmm_rnaseq.cocor_filtered)=names(msmm_rnaseq.cocor)
for (t in 1:length(msmm_rnaseq.cocor_filtered)){
  msmm_rnaseq.cocor_filtered[[t]]=msmm_rnaseq.cocor[[t]][which((msmm_rnaseq.cocor[[t]]$p.cocor<=0.05)&(msmm_rnaseq.cocor[[t]]$FDR<=0.1)),]
  write.table(msmm_rnaseq.cocor_filtered[[t]],file = paste(names(msmm_rnaseq.cocor_filtered)[t],"p005_cocor_filtered_interactions.txt",sep="_"),sep = "\t",col.names = T,row.names = T,quote=F)
}
#[order(msmm_rnaseq.cocor[[t]]$abs.corr.change,decreasing = T),]
#which(msmm_rnaseq.cocor_filtered[[t]]$FDR<=0.1)
msmm_rnaseq.DCN=vector(mode = "list",length = 3)
names(msmm_rnaseq.DCN)=names(msmm_rnaseq.cocor)
msmm_rnaseq.DCN[[1]]=msmm_rnaseq.DCN[[2]]=msmm_rnaseq.DCN[[3]]=vector(mode = "list",length = 3)
names(msmm_rnaseq.DCN[[1]])=names(msmm_rnaseq.DCN[[2]])=names(msmm_rnaseq.DCN[[3]])=c("Low","High","Std")
for (t in 1:length(msmm_rnaseq.DCN)){
  msmm_rnaseq.DCN[[t]][[1]]=msmm_rnaseq.DCN[[t]][[2]]=msmm_rnaseq.DCN[[t]][[3]]=graph.data.frame(d = msmm_rnaseq.cocor_filtered[[t]][which(msmm_rnaseq.cocor_filtered[[t]]$FDR<=0.1),c(1:2)],directed = F)
  E(msmm_rnaseq.DCN[[t]][[1]])$weight=msmm_rnaseq.cocor_filtered[[t]]$r.c+1
  E(msmm_rnaseq.DCN[[t]][[2]])$weight=msmm_rnaseq.cocor_filtered[[t]]$r.t+1
  E(msmm_rnaseq.DCN[[t]][[3]])$weight=msmm_rnaseq.cocor_filtered[[t]]$abs.corr.change
  write(V(msmm_rnaseq.DCN[[t]][[1]])$name,paste(names(msmm_rnaseq.DCN)[t],"p005","DiffCoexp","Genes.txt",sep = "_"),sep = "\n")
  write.graph(msmm_rnaseq.DCN[[t]][[1]],paste(names(msmm_rnaseq.DCN)[t],names(msmm_rnaseq.DCN[[t]])[1],"p005","DiffCoexp","Subgraph.gml",sep = "_"),format = "gml")
  write.graph(msmm_rnaseq.DCN[[t]][[2]],paste(names(msmm_rnaseq.DCN)[t],names(msmm_rnaseq.DCN[[t]])[2],"p005","DiffCoexp","Subgraph.gml",sep = "_"),format = "gml")
  write.graph(msmm_rnaseq.DCN[[t]][[3]],paste(names(msmm_rnaseq.DCN)[t],names(msmm_rnaseq.DCN[[t]])[3],"p005","DiffCoexp","Subgraph.gml",sep = "_"),format = "gml")
}
par(mfrow=c(3,3))
for(t in 1:length(names(msmm_rnaseq.DCN))){
  hist(msmm_rnaseq.cocor_filtered[[t]]$r.c,col="darkolivegreen3",breaks = 100,xlab = "Spearmann transformed Rho",main = paste(names(msmm_rnaseq.DCN)[t],"_Low_PLQ Spearmann transformed Rho",sep = ""))
  hist(msmm_rnaseq.cocor_filtered[[t]]$r.t,col="firebrick3",breaks = 100,xlab = "Spearmann transformed Rho",main = paste(names(msmm_rnaseq.DCN)[t],"_High_PLQ Spearmann transformed Rho",sep = ""))
  hist(msmm_rnaseq.cocor_filtered[[t]]$abs.corr.change,col="deepskyblue3",breaks = 100,xlab = "abs(rho.Low_PLQ - rho.High_PLQ)",main = paste(names(msmm_rnaseq.DCN)[t]," Abs. change in correlation",sep = ""))
  
}


save.image("msmm_rnaseq_contScore_cocor_p001.RData")
save.image("msmm_rnaseq_contScore_cocor_p001.RData")
