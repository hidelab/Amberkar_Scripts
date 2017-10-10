library(igraph)
library(foreach)
library(cocor)
library(synapseClient)
synapseLogin(username = "s.amberkar@sheffield.ac.uk",apiKey = "evb/5m+/10KmKAOP2vS1G6+a20iWAQlDosD9UfoQhvvFUdip/R/kZCzuk3jYecQ7zti5F4ZePz8djJQ8PoRC6Q==",rememberMe = T)
mapIds2<-function(IDs,IDFrom,IDTo){
  idmap=mapIds(x = org.Hs.eg.db,keys = IDs,column = IDTo,keytype = IDFrom,multiVals = "first")
  na_vec=names(idmap[is.na(idmap)==T])
  idmap=idmap[is.na(idmap)==F]
  idmap_df=data.frame("From"=names(idmap),"To"=unlist(unname(idmap)),stringsAsFactors = F)
  return(list(map=idmap_df,noMap=na_vec))
}
setwd("/shared/hidelab2/user/md4zsa/Work/Data/AMP-AD_RNAseq_ReSeq/Normalised_covariate_corrected_NoResiduals")
rosmap_reseq_data_pointer<-synGet(id='syn8456719')
rosmap_reseq_data=fread(rosmap_reseq_data_pointer@filePath,sep = "\t",header = T,stringsAsFactors = F,showProgress = T)
rosmap_covariates=read.table("ROSMAP/ROSMAP_DLPFC_Covariates.tsv",header = T,sep = "\t",stringsAsFactors = F)
ensembl_geneSymbol_map=mapIds2(IDs = rosmap_reseq_data$ensembl_gene_id,IDFrom = "ENSEMBL",IDTo = "SYMBOL")
rosmap_reseq_data2=rosmap_reseq_data[-which(rosmap_reseq_data$ensembl_gene_id%in%mapIds2(IDs = rosmap_reseq_data$ensembl_gene_id,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[2]]),]
rosmap_reseq_data2$gene_symbol=mapIds2(IDs = rosmap_reseq_data2$ensembl_gene_id,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[1]][,2]
rosmap_reseq_data2.agg=aggregate(x=rosmap_reseq_data2[,-c(1,634)],by=list(Symbol=rosmap_reseq_data2$gene_symbol),mean)
rownames(rosmap_reseq_data2.agg)=rosmap_reseq_data2.agg$Symbol
rosmap_reseq_data2.agg=rosmap_reseq_data2.agg[,-1]
c_counts=rosmap_reseq_data2.agg[,rosmap_covariates$SampleID[rosmap_covariates$Diagnosis=="CONTROL"]]
t_counts=rosmap_reseq_data2.agg[,rosmap_covariates$SampleID[rosmap_covariates$Diagnosis=="AD"]]
#################################################################################################################################################
# mayo_reseq_data=fread("MAYO/MAYO_CBE_TCX_netResidualExpression.tsv",sep="\t",header=T,data.table=F)
# mayo_reseq_data=fread("MAYO/MAYO_CBE_TCX_netResidualExpression.tsv",sep="\t",header=T,data.table=F)
# rownames(mayo_reseq_data)=mayo_reseq_data$ensembl_gene_id
# ensembl_geneSymbol_map=mapIds2(IDs = mayo_reseq_data$ensembl_gene_id,IDFrom = "ENSEMBL",IDTo = "SYMBOL")
# mayo_reseq_data2=mayo_reseq_data[-which(mayo_reseq_data$ensembl_gene_id%in%mapIds2(IDs = mayo_reseq_data$ensembl_gene_id,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[2]]),]
# mayo_reseq_data2$gene_symbol=mapIds2(IDs = mayo_reseq_data$ensembl_gene_id,IDFrom = "ENSEMBL",IDTo = "SYMBOL")[[1]][,2]
# mayo_reseq_data2.agg=aggregate(x=mayo_reseq_data2[,-c(1,529)],by=list(Symbol=mayo_reseq_data2$gene_symbol),mean)
# rownames(mayo_reseq_data2.agg)=mayo_reseq_data2.agg$Symbol
# mayo_reseq_data2.agg=mayo_reseq_data2.agg[,-1]
# mayo_covariates_CER=fread("MAYO/MayoRNAseq_RNAseq_CBE_covariates.csv",sep="\t",header=T,data.table=F)

# exprs_rank=mayo_reseq_data2.agg
# c_counts=mayo_reseq_data2.agg[,grep(pattern = paste(mayo_covariates_CER$SampleID[grep(pattern = "Control",x = mayo_covariates_CER$Diagnosis)],collapse = "|"),x = colnames(mayo_reseq_data2.agg))]
# t_counts=mayo_reseq_data2.agg[,grep(pattern = paste(mayo_covariates_CER$SampleID[grep(pattern = "AD",x = mayo_covariates_CER$Diagnosis)],collapse = "|"),x = colnames(mayo_reseq_data2.agg))]



dir.create("./ROSMAP/results_dC/",mode = 0777)
setwd("./ROSMAP/results_dC/")
i=0
blocksize=10
start<-i*blocksize+1
end<-min((i+1)*blocksize, length(rownames(c_counts)))
while(start<length(rownames(test_c_counts))){
  #Initialise foreach loop to parallelise computations
  res=foreach(j=start:end)%dopar%{
    #Compute Spearman correlation for all genes
    tmp_c.scc=apply(test_c_counts,1,function(x)cor.test(x = x,y = unname(unlist(test_c_counts[j,])),method = "spearman"))
    tmp_t.scc=apply(test_t_counts,1,function(x)cor.test(x = x,y = unname(unlist(test_t_counts[j,])),method = "spearman"))
    #Apply BH correction to p values
    adj.p_c.scc=p.adjust(p = unname(unlist(lapply(apply(test_c_counts,1,function(x)cor.test(x = x,y = unname(unlist(test_c_counts[j,])),method = "spearman")),`[`,'p.value'))),method = "BH")
    adj.p_t.scc=p.adjust(p = unname(unlist(lapply(apply(test_c_counts,1,function(x)cor.test(x = x,y = unname(unlist(test_c_counts[j,])),method = "spearman")),`[`,'p.value'))),method = "BH")
    #Apply genenames to vector for convenient indexing
    names(adj.p_c.scc)=names(tmp_c.scc)
    names(adj.p_t.scc)=names(tmp_t.scc)
    #Filter correlations at adj.p<0.05
    filtered_c.scc=adj.p_c.scc<0.05
    filtered_t.scc=adj.p_t.scc<0.05
    keep.scc=mapply(FUN = function(a,b)a|b,filtered_c.scc,filtered_t.scc)
    #Write the filtered correlations for gene j in a dataframe, separately for Control and Case
    tmp_c.df=data.frame(cbind(rownames(test_c_counts[j,]),names(tmp_c.scc[keep.scc==T]),unlist(lapply(tmp_c.scc[keep.scc],function(a)round(a$estimate,digits = 3)))),stringsAsFactors = F)
    tmp_t.df=data.frame(cbind(rownames(test_t_counts[j,]),names(tmp_t.scc[keep.scc==T]),unlist(lapply(tmp_t.scc[keep.scc],function(a)round(a$estimate,digits = 3)))),stringsAsFactors = F)
    colnames(tmp_c.df)=colnames(tmp_t.df)=c("Gene.A","Gene.B","Spearman.Rho")
    tmp_c.df$Spearman.Rho=as.numeric(tmp_c.df$Spearman.Rho)
    tmp_t.df$Spearman.Rho=as.numeric(tmp_t.df$Spearman.Rho)
    #Return the filtered correlation dataframes in a list
    result=list(Corr_C=tmp_c.df,Corr_T=tmp_t.df)
    
  }
  Corr_C.df=data.frame(rbindlist(lapply(lapply(lapply(res,`[`,'Corr_C'),data.frame,stringsAsFactors=F),function(a)a%>%transform(Corr_C.3=as.numeric(Corr_C.3)))))
  Corr_T.df=data.frame(rbindlist(lapply(lapply(lapply(res,`[`,'Corr_T'),data.frame,stringsAsFactors=F),function(a)a%>%transform(Corr_T.3=as.numeric(Corr_T.3)))))
  colnames(Corr_C.df)=colnames(Corr_T.df)=c("Gene.A","Gene.B","Spearman.Rho")
  dCG.expr[start:end]=lapply(res,function(a)sqrt(sum((a$Corr_C$Spearman.Rho-a$Corr_T$Spearman.Rho)^2)/(length(unique(a$Corr_C$Gene.B)))))
  cat(paste("Processing block ...",i,"\n"))
  fwrite(Corr_C.df,paste("tmp_scc_CorrC_",i,".txt",sep = ""),sep="\t",col.names = T,row.names = F,showProgress = T,buffMB = 10,nThread = 2)
  fwrite(Corr_T.df,paste("tmp_scc_CorrT_",i,".txt",sep = ""),sep="\t",col.names = T,row.names = F,showProgress = T,buffMB = 10,nThread = 2)
  i<-i+1
  start<-i*blocksize+1
  
  end<-min((i+1)*blocksize, length(rownames(c_counts)))
  #Add conditional statement if number of genes are less than end block value
  if(end>length(rownames(test_c_counts))){
    end=(i*blocksize)+end-length(rownames(test_c_counts))
  }
}

dCG.expr=lapply(res,function(a)sqrt(sum((a$Corr_C$Spearman.Rho-a$Corr_T$Spearman.Rho)^2)/(length(unique(a$Corr_C$Gene.B)))))
names(dCG.expr)=rownames(test_c_counts)


# res=foreach(i=start:end)%dopar%{
#   #Compute Spearman correlation for all genes
#   tmp_c.scc=apply(test_c_counts,1,function(x)cor.test(x = x,y = unname(unlist(test_c_counts[i,])),method = "spearman"))
#   tmp_t.scc=apply(test_t_counts,1,function(x)cor.test(x = x,y = unname(unlist(test_t_counts[i,])),method = "spearman"))
#   #Add filter q-value step
#   adj.p_c.scc=p.adjust(p = unname(unlist(lapply(apply(test_c_counts,1,function(x)cor.test(x = x,y = unname(unlist(test_c_counts[i,])),method = "spearman")),`[`,'p.value'))),method = "BH")
#   adj.p_t.scc=p.adjust(p = unname(unlist(lapply(apply(test_c_counts,1,function(x)cor.test(x = x,y = unname(unlist(test_c_counts[i,])),method = "spearman")),`[`,'p.value'))),method = "BH")
#   
#   names(adj.p_c.scc)=names(tmp_c.scc)
#   names(adj.p_t.scc)=names(tmp_t.scc)
#   
#   filtered_c.scc=adj.p_c.scc<0.05
#   filtered_t.scc=adj.p_t.scc<0.05
#   
#   keep.scc=mapply(FUN = function(a,b)a|b,filtered_c.scc,filtered_t.scc)
#   tmp_c.df=cbind(rownames(test_c_counts[i,]),names(tmp_c.scc[keep.scc==T]),unlist(lapply(tmp_c.scc[keep.scc],function(a)round(a$estimate,digits = 3))))
#   tmp_t.df=cbind(rownames(test_t_counts[i,]),names(tmp_t.scc[keep.scc==T]),unlist(lapply(tmp_t.scc[keep.scc],function(a)round(a$estimate,digits = 3))))
#   result=list(Corr_C=tmp_c.df,Corr_T=tmp_t.df)
#   
# }



all_counts=cbind(c_counts,t_counts)

random_dCG_exprList=lapply(out,function(x){x<-c(rep(0,10));x})
c_exprs=all_counts[1:75,colnames(all_counts)%in%sample(colnames(all_counts),size = length(colnames(test_c_counts)))]
t_exprs=all_counts[1:75,colnames(all_counts)%in%sample(colnames(all_counts),size = length(colnames(test_t_counts)))]
Random_dCP=function(n_genes,all_exprs,c_exprs,t_exprs){
  #cat(paste("Iteration Nr ",iterations,"...\n",sep=""))
  all_exprs=all_counts
  random_c_counts=all_counts[1:n_genes,colnames(all_exprs)%in%sample(colnames(all_exprs),size = length(colnames(c_exprs)))]
  random_t_counts=all_counts[1:n_genes,colnames(all_exprs)%in%sample(colnames(all_exprs),size = length(colnames(t_exprs)))]
  i=0
  blocksize=10
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
      random_tmp_c.df=data.frame(cbind(rownames(random_c_counts[j,]),names(random_tmp_c.scc[random_keep.scc==T]),unlist(lapply(random_tmp_c.scc[random_keep.scc],function(a)round(a$estimate,digits = 3)))),stringsAsFactors = F)
      random_tmp_t.df=data.frame(cbind(rownames(random_t_counts[j,]),names(random_tmp_t.scc[random_keep.scc==T]),unlist(lapply(random_tmp_t.scc[random_keep.scc],function(a)round(a$estimate,digits = 3)))),stringsAsFactors = F)
      colnames(random_tmp_c.df)=colnames(random_tmp_t.df)=c("Gene.A","Gene.B","Spearman.Rho")
      random_tmp_c.df$Spearman.Rho=as.numeric(random_tmp_c.df$Spearman.Rho)
      random_tmp_t.df$Spearman.Rho=as.numeric(random_tmp_t.df$Spearman.Rho)
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

names(random_all_genes.res)=paste("Iteration",1:5,sep = "")


AD_wgs.lists[[1]]=scan("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Known_AD_genes_V2.txt",what = "char",sep = "\n")
AD_wgs.lists[[2]]=scan("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Exonic_control_genes.txt",what = "char",sep = "\n")
AD_wgs.lists[[3]]=scan("/Users/sandeepamberkar/Work/Collaborations/Tanzi_WGS/Exonic_case_genes.txt",what = "char",sep = "\n")
data.frame(enrichKEGG(gene = mapIds2(IDs = AD_wgs.lists[[1]],IDFrom = "SYMBOL",IDTo = "ENTREZID")[[1]][,2],organism = "hsa",pvalueCutoff = 0.05,pAdjustMethod = "BH"))[,2]
