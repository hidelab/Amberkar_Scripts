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
  msmm_rnaseq.cocor_filtered[[t]]=msmm_rnaseq.cocor[[t]][(msmm_rnaseq.cocor[[t]]$p.cocor<=0.05),]
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

ffpe_wiesniewski=read.xls("../../../BrainExpression_Datasets/FFPE-Proteomics-Wiesniewski/srep15456-s6.xls",sheet = 2,header=T,as.is=T)
ffpe_wiesniewski2=ffpe_wiesniewski$Accession[which((ffpe_wiesniewski$Neuronal..more.stringent.database.=="Yes")&(ffpe_wiesniewski$Alzheimer.s.associated.protein=="Yes"))]
ffpe_wiesniewski2=gsub(pattern = "\\-[0-9]",replacement = "",x = ffpe_wiesniewski$Accession[which((ffpe_wiesniewski$Neuronal..more.stringent.database.=="Yes")&(ffpe_wiesniewski$Alzheimer.s.associated.protein=="Yes"))])
#Resolve dual IDs manually
ffpe_wiesniewski2[c(6,17,56,72,80,121,191)]=c("P68366","P09104","P38606","P30048","P18669","P48735","P27338")
seeds=intersect(V(msmm_rnaseq.DCN$BM_36$Std)$name,select(x = org.Hs.eg.db,keys = ffpe_wiesniewski2,columns = "SYMBOL",keytype = "UNIPROT")[,2])

ebc=cluster_edge_betweenness(graph = msmm_rnaseq.DCN$BM_36$Std,directed = F,modularity = T,weights = E(msmm_rnaseq.DCN$BM_36$Std)$weight)
comms=communities(ebc)
indices=names(unlist(lapply(communities(x = ebc),function(x)which(length(x)>=10))))
# msmm_rnaseq.PLQ_Corr005=vector(mode = "list",length = 3)
# names(msmm_rnaseq.PLQ_Corr005)=names(msmm_rnaseq.DCN)
# msmm_rnaseq.PLQ_Corr005[[1]]=scan("msmm_rnaseq_BM10_PLQ_Corr005_Genes.txt",what = "char",sep = "\n")
# msmm_rnaseq.PLQ_Corr005[[2]]=scan("msmm_rnaseq_BM22_PLQ_Corr005_Genes.txt",what = "char",sep = "\n")
# msmm_rnaseq.PLQ_Corr005[[3]]=scan("msmm_rnaseq_BM36_PLQ_Corr005_Genes.txt",what = "char",sep = "\n")

# 
# msmm_rnaseq.megena=vector(mode = "list",length = 3)
# names(msmm_rnaseq.megena)=names(msmm_rnaseq.DCN)
# megena_files=list.files("MEGENA_Networks/",pattern = "*.txt",full.names = T)
# for (t in 1:length(megena_files)){
#   dt=read.table(file = megena_files[t],sep="\t",header = T,as.is = T)
#   msmm_rnaseq.megena[[t]]=graph.data.frame(dt[,c(1:2)],directed = F)
#   E(msmm_rnaseq.megena[[t]])$weight=dt[,3]
# }
# 
# msmm_rnaseq.DCN$BM_22$Std=delete_vertices(graph = msmm_rnaseq.DCN$BM_22$Std,V(msmm_rnaseq.DCN$BM_22$Std)$name[which(V(msmm_rnaseq.DCN$BM_22$Std)$name%in%mapped_Ids[msmm_rnaseq.DCN_BM22_NA_Map$Index[which(msmm_rnaseq.DCN_BM22_NA_Map$ENTREZID=="<NA>")],1])])
save.image("msmm_rnaseq_contScore_cocor_p001.RData")
save.image("msmm_rnaseq_contScore_cocor_p001.RData")
