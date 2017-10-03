library(igraph)
library(foreach)
library(cocor)
library(synapseClient)
library(data.table)
library(org.Hs.eg.db)
library(doParallel)
library(DCGL)

synapseLogin(username = "s.amberkar@sheffield.ac.uk",apiKey = "evb/5m+/10KmKAOP2vS1G6+a20iWAQlDosD9UfoQhvvFUdip/R/kZCzuk3jYecQ7zti5F4ZePz8djJQ8PoRC6Q==",rememberMe = T)
cl=makeCluster(8)
registerDoParallel(cl)
#Define functions to be used
mapIds2<-function(IDs,IDFrom,IDTo){
  idmap=mapIds(x = org.Hs.eg.db,keys = IDs,column = IDTo,keytype = IDFrom,multiVals = "first")
  na_vec=names(idmap[is.na(idmap)==T])
  idmap=idmap[is.na(idmap)==F]
  idmap_df=data.frame("From"=names(idmap),"To"=unlist(unname(idmap)),stringsAsFactors = F)
  return(list(map=idmap_df,noMap=na_vec))
}
# Random_dCP=function(n_genes,all_exprs,c_exprs,t_exprs,blocksize){
#   #cat(paste("Iteration Nr ",iterations,"...\n",sep=""))
#   all_exprs=all_counts
#   #Randomly select between Control and Case samples
#   random_c_counts=all_counts[1:n_genes,colnames(all_exprs)%in%sample(colnames(all_exprs),size = length(colnames(c_exprs)))]
#   random_t_counts=all_counts[1:n_genes,colnames(all_exprs)%in%sample(colnames(all_exprs),size = length(colnames(t_exprs)))]
#   i=0
#   blocksize=blocksize
#   start<-i*blocksize+1
#   end<-min((i+1)*blocksize, length(rownames(random_c_counts)))
#   random_dCG.expr=c()
#   while(start<length(rownames(random_c_counts))){
#     #Initialise foreach loop to parallelise computations
#     random.res=foreach(j=start:end)%dopar%{
#       #Compute Spearman correlation for all genes
#       random_tmp_c.scc=apply(random_c_counts,1,function(x)cor.test(x = x,y = unname(unlist(random_c_counts[j,])),method = "spearman"))
#       random_tmp_t.scc=apply(random_t_counts,1,function(x)cor.test(x = x,y = unname(unlist(random_t_counts[j,])),method = "spearman"))
#       #Apply BH correction to p values
#       random_adj.p_c.scc=p.adjust(p = unname(unlist(lapply(apply(random_c_counts,1,function(x)cor.test(x = x,y = unname(unlist(random_c_counts[j,])),method = "spearman")),`[`,'p.value'))),method = "BH")
#       random_adj.p_t.scc=p.adjust(p = unname(unlist(lapply(apply(random_t_counts,1,function(x)cor.test(x = x,y = unname(unlist(random_t_counts[j,])),method = "spearman")),`[`,'p.value'))),method = "BH")
#       #Apply genenames to vector for convenient indexing
#       names(random_adj.p_c.scc)=names(random_tmp_c.scc)
#       names(random_adj.p_t.scc)=names(random_tmp_t.scc)
#       #Filter correlations at adj.p<0.05
#       random_filtered_c.scc=random_adj.p_c.scc<0.05
#       random_filtered_t.scc=random_adj.p_t.scc<0.05
#       random_keep.scc=mapply(FUN = function(a,b)a|b,random_filtered_c.scc,random_filtered_t.scc)
#       #Write the filtered correlations for gene j in a dataframe, separately for Control and Case
#       random_tmp_c.df=data.frame(cbind(rownames(random_c_counts[j,]),names(random_tmp_c.scc[random_keep.scc==T]),unlist(lapply(random_tmp_c.scc[random_keep.scc],function(a)round(a$estimate,digits = 3)))),stringsAsFactors = F)
#       random_tmp_t.df=data.frame(cbind(rownames(random_t_counts[j,]),names(random_tmp_t.scc[random_keep.scc==T]),unlist(lapply(random_tmp_t.scc[random_keep.scc],function(a)round(a$estimate,digits = 3)))),stringsAsFactors = F)
#       colnames(random_tmp_c.df)=colnames(random_tmp_t.df)=c("Gene.A","Gene.B","Spearman.Rho")
#       random_tmp_c.df$Spearman.Rho=as.numeric(random_tmp_c.df$Spearman.Rho)
#       random_tmp_t.df$Spearman.Rho=as.numeric(random_tmp_t.df$Spearman.Rho)
#       #Return the filtered correlation dataframes in a list
#       random_result=list(Corr_C=random_tmp_c.df,Corr_T=random_tmp_t.df)
#       
#     }
#     
#     random_Corr_C.df=data.frame(rbindlist(lapply(random.res,`[[`,1)),stringsAsFactors = F)
#     random_Corr_T.df=data.frame(rbindlist(lapply(random.res,`[[`,2)),stringsAsFactors = F)
#     # colnames(random_Corr_C.df)=colnames(random_Corr_T.df)=c("Gene.A","Gene.B","Spearman.Rho")
#     random_dCG.expr[start:end]=unlist(lapply(random.res,function(a)sqrt(sum((a$Corr_C$Spearman.Rho-a$Corr_T$Spearman.Rho)^2)/(length(unique(a$Corr_C$Gene.B))))))    
#     cat(paste("Processing block ...",i,"\n"))
#     fwrite(random_Corr_C.df,paste("random_tmp_scc_CorrC_",i,".txt",sep = ""),sep="\t",col.names = T,row.names = F,showProgress = T,buffMB = 10,nThread = 2)
#     fwrite(random_Corr_T.df,paste("random_tmp_scc_CorrT_",i,".txt",sep = ""),sep="\t",col.names = T,row.names = F,showProgress = T,buffMB = 10,nThread = 2)
#     i<-i+1
#     start<-i*blocksize+1
#     end<-min((i+1)*blocksize, length(rownames(random_c_counts)))
#     #Add conditional statement if number of genes are less than end block value
#     if(end>length(rownames(random_c_counts))){
#       end=(i*blocksize)+end-length(rownames(random_c_counts))
#     }
#     
#   }
#   return(random_dCG.expr)
# }
# dCP=function(n_genes,all_exprs,c_exprs,t_exprs,blocksize){
#   #cat(paste("Iteration Nr ",iterations,"...\n",sep=""))
#   all_exprs=all_counts
#   c_counts=c_exprs[1:n_genes,]
#   t_counts=t_exprs[1:n_genes,]
#   i=0
#   blocksize=blocksize
#   start<-i*blocksize+1
#   end<-min((i+1)*blocksize, length(rownames(c_counts)))
#   dCG.expr=c()
#   
#   while(start<length(rownames(c_counts))){
#     #Initialise foreach loop to parallelise computations
#     res=foreach(j=start:end)%dopar%{
#       #Compute Spearman correlation for all genes
#       
#       tmp_c.scc=apply(c_counts,1,function(x)cor.test(x = x,y = unname(unlist(c_counts[j,])),method = "spearman"))
#       tmp_t.scc=apply(t_counts,1,function(x)cor.test(x = x,y = unname(unlist(t_counts[j,])),method = "spearman"))
#       #Apply BH correction to p values
#       adj.p_c.scc=p.adjust(p = unname(unlist(lapply(apply(c_counts,1,function(x)cor.test(x = x,y = unname(unlist(c_counts[j,])),method = "spearman")),`[`,'p.value'))),method = "BH")
#       adj.p_t.scc=p.adjust(p = unname(unlist(lapply(apply(t_counts,1,function(x)cor.test(x = x,y = unname(unlist(t_counts[j,])),method = "spearman")),`[`,'p.value'))),method = "BH")
#       #Apply genenames to vector for convenient indexing
#       names(adj.p_c.scc)=names(tmp_c.scc)
#       names(adj.p_t.scc)=names(tmp_t.scc)
#       #Filter correlations at adj.p<0.05
#       filtered_c.scc=adj.p_c.scc<0.25
#       filtered_t.scc=adj.p_t.scc<0.25
#       keep.scc=mapply(FUN = function(a,b)a|b,filtered_c.scc,filtered_t.scc)
#       #Write the filtered correlations for gene j in a dataframe, separately for Control and Case
#       tmp_c.df=data.frame(cbind(rownames(c_counts[j,]),names(tmp_c.scc[keep.scc==T]),unlist(lapply(tmp_c.scc[keep.scc],function(a)round(a$estimate,digits = 3)))),stringsAsFactors = F)
#       tmp_t.df=data.frame(cbind(rownames(t_counts[j,]),names(tmp_t.scc[keep.scc==T]),unlist(lapply(tmp_t.scc[keep.scc],function(a)round(a$estimate,digits = 3)))),stringsAsFactors = F)
#       colnames(tmp_c.df)=colnames(tmp_t.df)=c("Gene.A","Gene.B","Spearman.Rho")
#       tmp_c.df$Spearman.Rho=as.numeric(tmp_c.df$Spearman.Rho)
#       tmp_t.df$Spearman.Rho=as.numeric(tmp_t.df$Spearman.Rho)
#       #Return the filtered correlation dataframes in a list
#       result=list(Corr_C=tmp_c.df,Corr_T=tmp_t.df)
#       
#     }
#     dCG.expr[start:end]=unlist(lapply(res,function(a)sqrt(sum((a$Corr_C$Spearman.Rho-a$Corr_T$Spearman.Rho)^2)/(length(unique(a$Corr_C$Gene.B))))))    
#     Corr_C.df=data.frame(rbindlist(lapply(res,`[[`,1)),stringsAsFactors = F)
#     Corr_T.df=data.frame(rbindlist(lapply(res,`[[`,2)),stringsAsFactors = F)
#     # colnames(Corr_C.df)=colnames(Corr_T.df)=c("Gene.A","Gene.B","Spearman.Rho")
#     cat(paste("Processing block ...",i,"\n"))
#     fwrite(Corr_C.df,paste("tmp_CorrC_block",i,".txt",sep = ""),sep="\t",col.names = T,row.names = F,showProgress = T,buffMB = 10,nThread = 2)
#     fwrite(Corr_T.df,paste("tmp_CorrT_block",i,".txt",sep = ""),sep="\t",col.names = T,row.names = F,showProgress = T,buffMB = 10,nThread = 2)
#     i<-i+1
#     start<-i*blocksize+1
#     end<-min((i+1)*blocksize, length(rownames(c_counts)))
#     #Add conditional statement if number of genes are less than end block value
#     if(end>length(rownames(c_counts))){
#       end=(i*blocksize)+end-length(rownames(c_counts))
#     }
#     
#   }
#   return(dCG.expr)
# }

#Set working directory
setwd("/shared/hidelab2/user/md4zsa/Work/Data/AMP-AD_RNAseq_ReSeq/Normalised_covariate_corrected_NoResiduals")
#Download data from Synapse
rosmap_reseq_data_pointer<-synGet(id='syn8456719')
rosmap_reseq_data=fread(rosmap_reseq_data_pointer@filePath,sep = "\t",header = T,stringsAsFactors = F,showProgress = T)

#Read ROSMAP covariates
rosmap_covariates=read.table("ROSMAP/ROSMAP_DLPFC_Covariates.tsv",header = T,sep = "\t",stringsAsFactors = F)

#Collapse Ensembl IDs to gene symbols. For multiple Ensembl IDs for same gene, compute average expression
ensembl_geneSymbol_map=mapIds2(IDs = rosmap_reseq_data$ensembl_gene_id,IDFrom = "ENSEMBL",IDTo = "SYMBOL")
rosmap_reseq_data2=rosmap_reseq_data[-which(rosmap_reseq_data$ensembl_gene_id%in%mapIds2(IDs = rosmap_reseq_data$ensembl_gene_id,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[2]]),]
rosmap_reseq_data2$gene_symbol=mapIds2(IDs = rosmap_reseq_data2$ensembl_gene_id,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[1]][,2]
rosmap_reseq_data2.agg=aggregate(x=rosmap_reseq_data2[,-c(1,634)],by=list(Symbol=rosmap_reseq_data2$gene_symbol),mean)
rownames(rosmap_reseq_data2.agg)=rosmap_reseq_data2.agg$Symbol
rosmap_reseq_data2.agg=rosmap_reseq_data2.agg[,-1]

#Read TF-Target interactions, preprocess data
regnet_tf2target=fread("/shared/hidelab2/user/md4zsa/Work/Data/TF_Databases/RegNetwork_human_regulators2.txt",header = T,sep = "\t",showProgress = T,data.table = F)%>%filter(evidence=="Experimental")%>%select(c(regulator_symbol,target_symbol))

#Segragate Control and AD samples
rosmap_c_counts=rosmap_reseq_data2.agg[,rosmap_covariates$SampleID[rosmap_covariates$Diagnosis=="CONTROL"]]
rosmap_t_counts=rosmap_reseq_data2.agg[,rosmap_covariates$SampleID[rosmap_covariates$Diagnosis=="AD"]]
rosmap_DCp=DCp(exprs.1 = rosmap_c_counts,exprs.2 = rosmap_t_counts,r.method = "spearman",link.method = "qth",cutoff = 0.05,N = 1000)
rosmap_DCe=DCe(exprs.1 = rosmap_c_counts,exprs.2 = rosmap_t_counts,r.method = "spearman",p = 0.05,link.method = "qth",cutoff = 0.05)
rosmap_DCsum.res=DCsum(rosmap_DCp,rosmap_DCe,DCpcutoff=0.05,DCecutoff=0.05)
rosmap_DRsort.res=DRsort(DCGs = rosmap_DCsum.res$DCGs,DCLs = rosmap_DCsum.res$DCLs,tf2target = regnet_tf2target,expGenes = rosmap_reseq_data2.agg)
saveRDS(rosmap_DCp,"ROSMAP_DCp_AnalysisResults.RDS")
saveRDS(rosmap_DCe,"ROSMAP_DCe_AnalysisResults.RDS")
saveRDS(rosmap_DCsum.res,"ROSMAP_DCsum_AnalysisResults.RDS")
saveRDS(rosmap_DRsort.res,"ROSMAP_DRsort_AnalysisResults.RDS")

proc.time()
cat(paste("Completed!"))
stopCluster(cl)