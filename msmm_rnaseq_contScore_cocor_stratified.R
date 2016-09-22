library(igraph)
library(org.Hs.eg.db)
library(cocor)
library(data.table)

setwd("/shared/hidelab2/user/md4zsa/Work/Data/MSMM_RNAseq/MSMM_RNAseq_FinalRelease/Strat_DCN")
cat(paste("Reading data ...\n"))
load("../msmm_rnaseq_FinalDataset_contPLQ.RData")
cat(paste("Applying significance threshold of 0.05 ...\n"))
msmm_rnaseq.Sample_Strata=vector(mode = "list",length=3)
names(msmm_rnaseq.Sample_Strata)=names(msmm_rnaseq)
msmm_rnaseq.Sample_Strata$BM_10=msmm_rnaseq.Sample_Strata$BM_22=msmm_rnaseq.Sample_Strata$BM_36=vector(mode = "list",length = 2)
names(msmm_rnaseq.Sample_Strata$BM_10)=names(msmm_rnaseq.Sample_Strata$BM_22)=names(msmm_rnaseq.Sample_Strata$BM_36)=c("Stratum1","Stratum2")
msmm_rnaseq.Sample_Strata$BM_10$Stratum1=msmm_rnaseq.Sample_Strata$BM_10$Stratum2=msmm_rnaseq.Sample_Strata$BM_22$Stratum1=msmm_rnaseq.Sample_Strata$BM_22$Stratum2=msmm_rnaseq.Sample_Strata$BM_36$Stratum1=msmm_rnaseq.Sample_Strata$BM_36$Stratum2=vector(mode = "list",length = 2)
names(msmm_rnaseq.Sample_Strata$BM_10$Stratum1)=names(msmm_rnaseq.Sample_Strata$BM_10$Stratum2)=names(msmm_rnaseq.Sample_Strata$BM_22$Stratum1)=names(msmm_rnaseq.Sample_Strata$BM_22$Stratum2)=names(msmm_rnaseq.Sample_Strata$BM_36$Stratum1)=names(msmm_rnaseq.Sample_Strata$BM_36$Stratum2)=c("Low","High")

stratum1.Low=msmm_rnaseq_covariates$Sample.ID[which(msmm_rnaseq_covariates$PlaqueMean<1)]
stratum1.High=msmm_rnaseq_covariates$Sample.ID[which(msmm_rnaseq_covariates$PlaqueMean>=1&msmm_rnaseq_covariates$PlaqueMean<=10)]
stratum2.Low=msmm_rnaseq_covariates$Sample.ID[which(msmm_rnaseq_covariates$PlaqueMean<11)]
stratum2.High=msmm_rnaseq_covariates$Sample.ID[which(msmm_rnaseq_covariates$PlaqueMean>=11)]
# stratum3.Low=msmm_rnaseq_covariates$Sample.ID[which(msmm_rnaseq_covariates$PlaqueMean<21)]
# stratum3.High=msmm_rnaseq_covariates$Sample.ID[which(msmm_rnaseq_covariates$PlaqueMean>=21)]
# stratum4.Low=msmm_rnaseq_covariates$Sample.ID[which(msmm_rnaseq_covariates$PlaqueMean<31)]
# stratum4.High=msmm_rnaseq_covariates$Sample.ID[which(msmm_rnaseq_covariates$PlaqueMean>=31)]

msmm_rnaseq.Sample_Strata$BM_10$Stratum1$Low=grep(pattern = "BM_10",stratum1.Low,value = T)
msmm_rnaseq.Sample_Strata$BM_10$Stratum1$High=grep(pattern = "BM_10",stratum1.High,value = T)
msmm_rnaseq.Sample_Strata$BM_10$Stratum2$Low=grep(pattern = "BM_10",stratum2.Low,value = T)
msmm_rnaseq.Sample_Strata$BM_10$Stratum2$High=grep(pattern = "BM_10",stratum2.High,value = T)
# msmm_rnaseq.Sample_Strata$BM_10$Stratum3$Low=stratum3.Low=grep(pattern = "BM_10",stratum3.Low,value = T)
# msmm_rnaseq.Sample_Strata$BM_10$Stratum3$High=stratum3.High=grep(pattern = "BM_10",stratum3.High,value = T)

msmm_rnaseq.Sample_Strata$BM_22$Stratum1$Low=grep(pattern = "BM_22",stratum1.Low,value = T)
msmm_rnaseq.Sample_Strata$BM_22$Stratum1$High=grep(pattern = "BM_22",stratum1.High,value = T)
msmm_rnaseq.Sample_Strata$BM_22$Stratum2$Low=grep(pattern = "BM_22",stratum2.Low,value = T)
msmm_rnaseq.Sample_Strata$BM_22$Stratum2$High=grep(pattern = "BM_22",stratum2.High,value = T)
# msmm_rnaseq.Sample_Strata$BM_22$Stratum3$Low=stratum3.Low=grep(pattern = "BM_22",stratum3.Low,value = T)
# msmm_rnaseq.Sample_Strata$BM_22$Stratum3$High=stratum3.High=grep(pattern = "BM_22",stratum3.High,value = T)

msmm_rnaseq.Sample_Strata$BM_36$Stratum1$Low=grep(pattern = "BM_36",stratum1.Low,value = T)
msmm_rnaseq.Sample_Strata$BM_36$Stratum1$High=grep(pattern = "BM_36",stratum1.High,value = T)
msmm_rnaseq.Sample_Strata$BM_36$Stratum2$Low=grep(pattern = "BM_36",stratum2.Low,value = T)
msmm_rnaseq.Sample_Strata$BM_36$Stratum2$High=grep(pattern = "BM_36",stratum2.High,value = T)
# msmm_rnaseq.Sample_Strata$BM_36$Stratum3$Low=stratum3.Low=grep(pattern = "BM_36",stratum3.Low,value = T)
# msmm_rnaseq.Sample_Strata$BM_36$Stratum3$High=stratum3.High=grep(pattern = "BM_36",stratum3.High,value = T)


msmm_rnaseq.cocor=vector(mode = "list",length = 3)
msmm_rnaseq.cocor[[1]]=msmm_rnaseq.cocor[[2]]=msmm_rnaseq.cocor[[3]]=vector(mode = "list",length = 2)
names(msmm_rnaseq.cocor)=names(msmm_rnaseq)
names(msmm_rnaseq.cocor[[1]])=names(msmm_rnaseq.cocor[[2]])=names(msmm_rnaseq.cocor[[3]])=c("Stratum1","Stratum2")

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


  for (t in 1:length(names(msmm_rnaseq))){
    for (s in 1:2){
    exprs_rank<-msmm_rnaseq[[t]][which(rownames(msmm_rnaseq[[t]])%in%msmm_rnaseq_Plaque.Corr005.Genes[[t]]),]
    number_of_combinations<-choose(nrow(exprs_rank),2)
    c_exprs_rank=exprs_rank[,which(colnames(exprs_rank)%in%msmm_rnaseq.Sample_Strata[[t]][[s]]$Low)]
    t_exprs_rank=exprs_rank[,which(colnames(exprs_rank)%in%msmm_rnaseq.Sample_Strata[[t]][[s]]$High)]
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
    msmm_rnaseq.cocor[[t]][[s]]=result
  }
}  

msmm_rnaseq.cocor_filtered=vector(mode = "list",length = 3)
names(msmm_rnaseq.cocor_filtered)=names(msmm_rnaseq.cocor)
msmm_rnaseq.cocor_filtered$BM_10=msmm_rnaseq.cocor_filtered$BM_22=msmm_rnaseq.cocor_filtered$BM_36=vector(mode = "list",length = 2)
names(msmm_rnaseq.cocor_filtered$BM_10)=names(msmm_rnaseq.cocor_filtered$BM_22)=names(msmm_rnaseq.cocor_filtered$BM_36)=c("Stratum1","Stratum2")

cat(paste("Filtering differentially correlated genes ...\n"))
for (t in 1:length(msmm_rnaseq.cocor_filtered)){
  for (s in 1:2){
    msmm_rnaseq.cocor_filtered[[t]][[s]]=msmm_rnaseq.cocor[[t]][[s]][(msmm_rnaseq.cocor[[t]][[s]]$p.cocor<=0.05&msmm_rnaseq.cocor[[t]][[s]]$FDR<=0.1),]
    write.table(msmm_rnaseq.cocor_filtered[[t]][[s]],file = paste(names(msmm_rnaseq.cocor_filtered)[t],names(msmm_rnaseq.cocor_filtered[[t]]),"p005_cocor_filtered_interactions.txt",sep="_"),sep = "\t",col.names = T,row.names = T,quote=F)
  }
}
msmm_rnaseq.DCN=vector(mode = "list",length = 3)
names(msmm_rnaseq.DCN)=names(msmm_rnaseq.cocor)
msmm_rnaseq.DCN[[1]]=msmm_rnaseq.DCN[[2]]=msmm_rnaseq.DCN[[3]]=vector(mode = "list",length = 3)
names(msmm_rnaseq.DCN[[1]])=names(msmm_rnaseq.DCN[[2]])=names(msmm_rnaseq.DCN[[3]])=c("Low","High","Std")
msmm_rnaseq.DCN$BM_10$Low=msmm_rnaseq.DCN$BM_10$High=msmm_rnaseq.DCN$BM_10$Std=vector(mode = "list",length = 2)
msmm_rnaseq.DCN$BM_22$Low=msmm_rnaseq.DCN$BM_22$High=msmm_rnaseq.DCN$BM_22$Std=vector(mode = "list",length = 2)
msmm_rnaseq.DCN$BM_36$Low=msmm_rnaseq.DCN$BM_36$High=msmm_rnaseq.DCN$BM_36$Std=vector(mode = "list",length = 2)
names(msmm_rnaseq.DCN$BM_10$Low)=names(msmm_rnaseq.DCN$BM_10$High)=names(msmm_rnaseq.DCN$BM_10$Std)=c("Stratum1","Stratum2")
names(msmm_rnaseq.DCN$BM_22$Low)=names(msmm_rnaseq.DCN$BM_22$High)=names(msmm_rnaseq.DCN$BM_22$Std)=c("Stratum1","Stratum2")
names(msmm_rnaseq.DCN$BM_36$Low)=names(msmm_rnaseq.DCN$BM_36$High)=names(msmm_rnaseq.DCN$BM_36$Std)=c("Stratum1","Stratum2")

cat(paste("Building DCNs ..."))
for (t in 1:length(names(msmm_rnaseq))){
  for (y in 1:3){
    for (s in 1:2){
      msmm_rnaseq.DCN[[t]][[y]][[s]]=msmm_rnaseq.DCN[[t]][[y]][[s]]=msmm_rnaseq.DCN[[t]][[y]][[s]]=graph.data.frame(d = msmm_rnaseq.cocor_filtered[[t]][[s]][which(msmm_rnaseq.cocor_filtered[[t]][[s]]$FDR<=0.1),c(1:2)],directed = F)
      E(msmm_rnaseq.DCN[[t]][[y]][[s]])$weight=msmm_rnaseq.cocor_filtered[[t]][[s]]$r.c+1
      E(msmm_rnaseq.DCN[[t]][[y]][[s]])$weight=msmm_rnaseq.cocor_filtered[[t]][[s]]$r.t+1
      E(msmm_rnaseq.DCN[[t]][[y]][[s]])$weight=msmm_rnaseq.cocor_filtered[[t]][[s]]$abs.corr.change
      cat(paste("Writing DCNs for",names(msmm_rnaseq.DCN)[t],names(msmm_rnaseq.DCN[[t]])[[y]],names(msmm_rnaseq.DCN[[t]][[y]])[[s]],"...\n",sep = " "))
      write(V(msmm_rnaseq.DCN[[t]][[y]][[s]])$name,paste(names(msmm_rnaseq.DCN)[t],names(msmm_rnaseq.DCN[[t]])[[y]],names(msmm_rnaseq.DCN[[t]][[y]])[[s]],"p005","DiffCoexp","Genes.txt",sep = "_"),sep = "\n")
      write.graph(msmm_rnaseq.DCN[[t]][[y]][[s]],paste(names(msmm_rnaseq.DCN)[t],names(msmm_rnaseq.DCN[[t]])[[y]],names(msmm_rnaseq.DCN[[t]][[y]])[[s]],"p005","DiffCoexp","Subgraph.gml",sep = "_"),format = "gml")
      write.graph(msmm_rnaseq.DCN[[t]][[y]][[s]],paste(names(msmm_rnaseq.DCN)[t],names(msmm_rnaseq.DCN[[t]])[[y]],names(msmm_rnaseq.DCN[[t]][[y]])[[s]],"p005","DiffCoexp","Subgraph.gml",sep = "_"),format = "gml")
      write.graph(msmm_rnaseq.DCN[[t]][[y]][[s]],paste(names(msmm_rnaseq.DCN)[t],names(msmm_rnaseq.DCN[[t]])[[y]],names(msmm_rnaseq.DCN[[t]][[y]])[[s]],"p005","DiffCoexp","Subgraph.gml",sep = "_"),format = "gml")
    }
  }
}
cat(paste("Saving results in RData object ...\n"))
save.image("msmm_rnaseq_contScore_cocor_p005_Stratum.RData")
cat(paste("Done!"))