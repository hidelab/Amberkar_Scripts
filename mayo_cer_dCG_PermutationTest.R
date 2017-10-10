library(igraph)
library(foreach)
library(cocor)
library(synapseClient)
library(data.table)
library(org.Hs.eg.db)
library(doParallel)

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
Random_dCP=function(n_genes,all_exprs,c_exprs,t_exprs,blocksize){
  #cat(paste("Iteration Nr ",iterations,"...\n",sep=""))
  all_exprs=all_exprs
  #Randomly select between Control and Case samples
  random_c_counts=all_exprs[1:n_genes,colnames(all_exprs)%in%sample(colnames(all_exprs),size = length(colnames(c_exprs)))]
  random_t_counts=all_exprs[1:n_genes,colnames(all_exprs)%in%sample(colnames(all_exprs),size = length(colnames(t_exprs)))]
  i=0
  blocksize=blocksize
  start<-i*blocksize+1
  end<-min((i+1)*blocksize, length(rownames(random_c_counts)))
  random_dCG.expr=c()
  while(start<length(rownames(random_c_counts))){
    #Initialise foreach loop to parallelise computations
    random.res=foreach(j=start:end)%dopar%{
      #Compute Spearman correlation for all genes
      random_tmp_c.scc=apply(random_c_counts,1,function(x)cor.test(x = x,y = unname(unlist(random_c_counts[j,])),method = "spearman"))
      random_tmp_t.scc=apply(random_t_counts,1,function(x)cor.test(x = x,y = unname(unlist(random_t_counts[j,])),method = "spearman"))
      #Apply BH correction to p values
      random_adj.p_c.scc=p.adjust(p = unname(unlist(lapply(apply(random_c_counts,1,function(x)cor.test(x = x,y = unname(unlist(random_c_counts[j,])),method = "spearman")),`[`,'p.value'))),method = "BH")
      random_adj.p_t.scc=p.adjust(p = unname(unlist(lapply(apply(random_t_counts,1,function(x)cor.test(x = x,y = unname(unlist(random_t_counts[j,])),method = "spearman")),`[`,'p.value'))),method = "BH")
      #Apply genenames to vector for convenient indexing
      names(random_adj.p_c.scc)=names(random_tmp_c.scc)
      names(random_adj.p_t.scc)=names(random_tmp_t.scc)
      #Filter correlations at adj.p<0.05
      random_filtered_c.scc=random_adj.p_c.scc<0.05
      random_filtered_t.scc=random_adj.p_t.scc<0.05
      random_keep.scc=mapply(FUN = function(a,b)a|b,random_filtered_c.scc,random_filtered_t.scc)
      #Write the filtered correlations for gene j in a dataframe, separately for Control and Case
      random_tmp_c.df=data.frame(cbind(rownames(random_c_counts[j,]),names(random_tmp_c.scc[random_keep.scc==T]),unlist(lapply(random_tmp_c.scc[random_keep.scc],function(a)round(a$estimate,digits = 3))),random_adj.p_c.scc),stringsAsFactors = F)
      random_tmp_t.df=data.frame(cbind(rownames(random_t_counts[j,]),names(random_tmp_t.scc[random_keep.scc==T]),unlist(lapply(random_tmp_t.scc[random_keep.scc],function(a)round(a$estimate,digits = 3))),random_adj.p_t.scc),stringsAsFactors = F)
      colnames(random_tmp_c.df)=colnames(random_tmp_t.df)=c("Gene.A","Gene.B","Spearman.Rho","adj.p")
      random_tmp_c.df$Spearman.Rho=as.numeric(random_tmp_c.df$Spearman.Rho)
      random_tmp_t.df$Spearman.Rho=as.numeric(random_tmp_t.df$Spearman.Rho)
      random_tmp_c.df$adj.p=as.numeric(random_tmp_c.df$adj.p)
      random_tmp_t.df$adj.p=as.numeric(random_tmp_t.df$adj.p)
      #Return the filtered correlation dataframes in a list
      random_result=list(Corr_C=random_tmp_c.df,Corr_T=random_tmp_t.df)
      
    }
    
    random_Corr_C.df=data.frame(rbindlist(lapply(random.res,`[[`,1)),stringsAsFactors = F)
    random_Corr_T.df=data.frame(rbindlist(lapply(random.res,`[[`,2)),stringsAsFactors = F)
    # colnames(random_Corr_C.df)=colnames(random_Corr_T.df)=c("Gene.A","Gene.B","Spearman.Rho")
    random_dCG.expr[start:end]=unlist(lapply(random.res,function(a)sqrt(sum((a$Corr_C$Spearman.Rho-a$Corr_T$Spearman.Rho)^2)/(length(unique(a$Corr_C$Gene.B))))))    
    cat(paste("Processing block ...",i,"\n"))
    fwrite(random_Corr_C.df,paste("random_tmp_scc_CorrC_",i,".txt",sep = ""),sep="\t",col.names = T,row.names = F,showProgress = T,buffMB = 10,nThread = 2)
    fwrite(random_Corr_T.df,paste("random_tmp_scc_CorrT_",i,".txt",sep = ""),sep="\t",col.names = T,row.names = F,showProgress = T,buffMB = 10,nThread = 2)
    i<-i+1
    start<-i*blocksize+1
    end<-min((i+1)*blocksize, length(rownames(random_c_counts)))
    #Add conditional statement if number of genes are less than end block value
    if(end>length(rownames(random_c_counts))){
      end=(i*blocksize)+end-length(rownames(random_c_counts))
    }
    
  }
  return(random_dCG.expr)
}
dCP=function(n_genes,all_exprs,c_exprs,t_exprs,blocksize){
  #cat(paste("Iteration Nr ",iterations,"...\n",sep=""))
  all_exprs=all_exprs
  c_counts=c_exprs[1:n_genes,]
  t_counts=t_exprs[1:n_genes,]
  i=0
  blocksize=blocksize
  start<-i*blocksize+1
  end<-min((i+1)*blocksize, length(rownames(c_counts)))
  dCG.expr=c()
  
  while(start<length(rownames(c_counts))){
    #Initialise foreach loop to parallelise computations
    res=foreach(j=start:end)%dopar%{
      #Compute Spearman correlation for all genes
      
      tmp_c.scc=apply(c_counts,1,function(x)cor.test(x = x,y = unname(unlist(c_counts[j,])),method = "spearman"))
      tmp_t.scc=apply(t_counts,1,function(x)cor.test(x = x,y = unname(unlist(t_counts[j,])),method = "spearman"))
      #Apply BH correction to p values
      adj.p_c.scc=p.adjust(p = unname(unlist(lapply(apply(c_counts,1,function(x)cor.test(x = x,y = unname(unlist(c_counts[j,])),method = "spearman")),`[`,'p.value'))),method = "BH")
      adj.p_t.scc=p.adjust(p = unname(unlist(lapply(apply(t_counts,1,function(x)cor.test(x = x,y = unname(unlist(t_counts[j,])),method = "spearman")),`[`,'p.value'))),method = "BH")
      #Apply genenames to vector for convenient indexing
      names(adj.p_c.scc)=names(tmp_c.scc)
      names(adj.p_t.scc)=names(tmp_t.scc)
      #Filter correlations at adj.p<0.05
      filtered_c.scc=adj.p_c.scc<0.10
      filtered_t.scc=adj.p_t.scc<0.10
      keep.scc=mapply(FUN = function(a,b)a|b,filtered_c.scc,filtered_t.scc)
      #Write the filtered correlations for gene j in a dataframe, separately for Control and Case
      tmp_c.df=data.frame(cbind(rownames(c_counts[j,]),names(tmp_c.scc[keep.scc==T]),unlist(lapply(tmp_c.scc[keep.scc],function(a)round(a$estimate,digits = 3))),adj.p_c.scc),stringsAsFactors = F)
      tmp_t.df=data.frame(cbind(rownames(t_counts[j,]),names(tmp_t.scc[keep.scc==T]),unlist(lapply(tmp_t.scc[keep.scc],function(a)round(a$estimate,digits = 3))),adj.p_t.scc),stringsAsFactors = F)
      colnames(tmp_c.df)=colnames(tmp_t.df)=c("Gene.A","Gene.B","Spearman.Rho","adj.p")
      tmp_c.df$Spearman.Rho=as.numeric(tmp_c.df$Spearman.Rho)
      tmp_t.df$Spearman.Rho=as.numeric(tmp_t.df$Spearman.Rho)
      tmp_c.df$adj.p=as.numeric(tmp_c.df$adj.p)
      tmp_t.df$adj.p=as.numeric(tmp_t.df$adj.p)
      #Return the filtered correlation dataframes in a list
      result=list(Corr_C=tmp_c.df,Corr_T=tmp_t.df)
      
    }
    dCG.expr[start:end]=unlist(lapply(res,function(a)sqrt(sum((a$Corr_C$Spearman.Rho-a$Corr_T$Spearman.Rho)^2)/(length(unique(a$Corr_C$Gene.B))))))    
    Corr_C.df=data.frame(rbindlist(lapply(res,`[[`,1)),stringsAsFactors = F)
    Corr_T.df=data.frame(rbindlist(lapply(res,`[[`,2)),stringsAsFactors = F)
    # colnames(Corr_C.df)=colnames(Corr_T.df)=c("Gene.A","Gene.B","Spearman.Rho")
    cat(paste("Processing block ...",i,"\n"))
    fwrite(Corr_C.df,paste("tmp_CorrC_block",i,".txt",sep = ""),sep="\t",col.names = T,row.names = F,showProgress = T,quote = F,buffMB = 10,nThread = 2)
    fwrite(Corr_T.df,paste("tmp_CorrT_block",i,".txt",sep = ""),sep="\t",col.names = T,row.names = F,showProgress = T,quote = F,buffMB = 10,nThread = 2)
    i<-i+1
    start<-i*blocksize+1
    end<-min((i+1)*blocksize, length(rownames(c_counts)))
    #Add conditional statement if number of genes are less than end block value
    if(end>length(rownames(c_counts))){
      end=(i*blocksize)+end-length(rownames(c_counts))
    }
    
  }
  return(dCG.expr)
}

#Set working directory
setwd("/shared/hidelab2/user/md4zsa/Work/Data/AMP-AD_RNAseq_ReSeq/Normalised_covariate_corrected_NoResiduals")
#Download data from Synapse
mayo_reseq_data_pointer<-synGet(id='syn8466826')
mayo_reseq_data=fread(mayo_reseq_data_pointer@filePath,sep = "\t",header = T,stringsAsFactors = F,showProgress = T,data.table = F)
#Read ROSMAP covariates
mayo_covariates_pointer <- synGet(id='syn8466814')
mayo_covariates=read.table(mayo_covariates_pointer@filePath,header = T,sep = "\t",stringsAsFactors = F)

#Collapse Ensembl IDs to gene symbols. For multiple Ensembl IDs for same gene, compute average expression
ensembl_geneSymbol_map=mapIds2(IDs = mayo_reseq_data$ensembl_gene_id,IDFrom = "ENSEMBL",IDTo = "SYMBOL")
mayo_reseq_data2=mayo_reseq_data[-which(mayo_reseq_data$ensembl_gene_id%in%mapIds2(IDs = mayo_reseq_data$ensembl_gene_id,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[2]]),]
mayo_reseq_data2$gene_symbol=mapIds2(IDs = mayo_reseq_data2$ensembl_gene_id,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[1]][,2]
mayo_reseq_data2.agg=aggregate(x=mayo_reseq_data2[,-c(1,(length(colnames(mayo_reseq_data2))))],by=list(Symbol=mayo_reseq_data2$gene_symbol),mean)
rownames(mayo_reseq_data2.agg)=mayo_reseq_data2.agg$Symbol
mayo_reseq_data2.agg=mayo_reseq_data2.agg[,-1]

#Segregate data into TCX and CER brain regions
mayo_reseq_tcx_data=mayo_reseq_data2.agg[,which(colnames(mayo_reseq_data2.agg)%in%mayo_covariates$SampleID[grep(pattern = "TCX.AD|TCX.Control",mayo_covariates$BrainRegion.Diagnosis)])]
mayo_reseq_cer_data=mayo_reseq_data2.agg[,which(colnames(mayo_reseq_data2.agg)%in%mayo_covariates$SampleID[grep(pattern = "CER.AD|CER.Control",mayo_covariates$BrainRegion.Diagnosis)])]


#Segragate Control and AD samples
cer_c_counts=mayo_reseq_cer_data[,mayo_covariates$SampleID[grep(pattern = "CER.Control",mayo_covariates$BrainRegion.Diagnosis)]]
cer_t_counts=mayo_reseq_cer_data[,mayo_covariates$SampleID[grep(pattern = "CER.AD",mayo_covariates$BrainRegion.Diagnosis)]]
setwd("./MAYO/CER/results_dC/")
cer_all_counts=data.frame(cbind(cer_c_counts,cer_t_counts),stringsAsFactors = F)
AD_dCP=dCP(n_genes =55,all_exprs = cer_all_counts,c_exprs = cer_c_counts,t_exprs = cer_t_counts,blocksize = 10)
saveRDS(AD_dCP,"mayo_CER_AD_dCP.RDS")
system("awk 'FNR==1 && NR!=1 { while (/^<header>/) getline; }    1 {print}' random_tmp_scc_CorrC*.txt >mayo_CER_Control_random_dCL_allResults.txt")
system("awk 'FNR==1 && NR!=1 { while (/^<header>/) getline; }    1 {print}' tmp_CorrC*.txt >mayo_CER_Control_dCL_allResults.txt")
system("awk 'FNR==1 && NR!=1 { while (/^<header>/) getline; }    1 {print}' tmp_CorrT*.txt >mayo_CER_AD_dCL_allResults.txt")
mayo_random_cer_control_dCL=fread("mayo_CER_Control_random_dCL_allResults.txt",sep = "\t",header = T,showProgress = T,data.table = F,quote = "")
mayo_cer_control_dCL=fread("mayo_CER_Control_dCL_allResults.txt",sep = "\t",header = T,showProgress = T,data.table = F,quote = "")
mayo_cer_control_dCL=fread("mayo_CER_Control_dCL_allResults.txt",sep = "\t",header = T,showProgress = T,data.table = F,quote = "")
mayo_cer_AD_dCL=fread("mayo_CER_Control_dCL_allResults.txt",sep = "\t",header = T,showProgress = T,data.table = F)
mayo_cer_control_dCL.graph=graph.data.frame(d=mayo_cer_control_dCL[,c(1:2)],directed = F)
mayo_cer_AD_dCL.graph=graph.data.frame(d=mayo_cer_AD_dCL[,c(1:2)],directed = F)

#Initialise permutation matrix as list
ngenes=length(rownames(all_counts))
permut=10
random_dCG_exprList=lapply(AD_dCP,function(x){x<-c(rep(0,permut));x})
setwd("../random_results_dC/")

random_dCP_exprList=replicate(n = permut,Random_dCP(n_genes = 55,all_exprs = cer_all_counts,c_exprs = cer_c_counts,t_exprs = cer_t_counts,blocksize = 10))
colnames(random_dCP_exprList)=paste("Perm",1:permut,sep = "")
rownames(random_dCP_exprList)=rownames(cer_all_counts)[1:55]
saveRDS(random_dCG_exprList,"Random_dCP_PermuteList.RDS")
proc.time()
cat(paste("Completed!"))
stopCluster(cl)