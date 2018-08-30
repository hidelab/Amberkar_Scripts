library(dplyr)
library(magrittr)
library(data.table)
library(DCGL)
library(parallel)

setwd("/shared/hidelab2/user/md4zsa/Work/Data/MSBB_Array19/GSE84422/LateAD_diffcoexp/")
lateAD_diffcoexp_files=list.files(path = ".",pattern = "lateAD_diffcoexp",full.names = T)
lateAD_samples=readRDS("./msbb_gse84422_GPL96_97_lateAD_samplesToAnalyse.RDS")
lateAD_samples.exprs=readRDS("./msbb_gse84422_GPL96_97_lateAD_samplesToAnalyse_exprs.RDS")

humanRegnetwork=readRDS("humanRegnetwork.RDS")
regnet_tf2target_data=fread("/shared/hidelab2/user/md4zsa/Work/Data/TF_Databases/RegNetwork_human_regulators2.txt",header = T,sep = "\t",showProgress = T,data.table = F)
regnet_tf2target.HGNC=regnet_tf2target_data%>%dplyr::select(c(regulator_symbol,target_symbol))
humanRegnetwork=humanRegnetwork[unique(regnet_tf2target.HGNC$regulator_symbol)]
humanRegnetwork=humanRegnetwork[unlist(lapply(humanRegnetwork,function(x)length(x)>2))]

lateAD_diffcoexp_results=vector(mode = "list",length = length(lateAD_diffcoexp_files))
for(i in 1:length(lateAD_diffcoexp_files)){
  lateAD_diffcoexp_results[[i]]=readRDS(lateAD_diffcoexp_files[[i]])
}
names(lateAD_diffcoexp_results)=unlist(lapply(strsplit(split = "_lateAD",basename(lateAD_diffcoexp_files)),`[[`,1))

lateAD_DCGs=lapply(lateAD_diffcoexp_results,function(x)x$DCGs)
lateAD_DCGs.filtered=lapply(lateAD_diffcoexp_results,function(x)x$DCGs%>%dplyr::filter(q<=0.1))
lateAD_DCGs.filtered=Filter(f = function(x)dim(x)[1]>0,x = lateAD_DCGs.filtered)
lateAD_DCLs=lapply(lateAD_diffcoexp_results,function(x)x$DCLs)
lateAD_DCLs.filtered=lapply(lateAD_diffcoexp_results,function(x)x$DCLs%>%dplyr::filter(q.diffcor<=0.1))
lateAD_DCLs.filtered=Filter(f = function(x)dim(x)[1]>0,x = lateAD_DCLs.filtered)

lateAD_DRrank.TED=lateAD_DRrank.TED=vector(mode = "list",length = length(lateAD_DCGs.filtered))
names(lateAD_DRrank.TED)=names(lateAD_DRrank.TED)=names(lateAD_DCGs.filtered)
compute_time=list()
lateAD_samples.exprs=lateAD_samples.exprs[names(lateAD_DCGs.filtered)]

for(t in 1:length(lateAD_DCGs.filtered)){
  genes=lateAD_DCGs.filtered[[t]]$Gene
  dcls=lateAD_DCLs.filtered[[t]]
  prob<-nrow(dcls)/choose(nrow(lateAD_samples.exprs[[t]]), 2)
  
  
  ProcessElement<-function(ic) {
    regulator_gene_symbol<-names(humanRegnetwork)[ic]
    targets<-humanRegnetwork[[ic]]
    common.genes<-intersect(targets, genes)
    n<-length(common.genes)
    N<-n * (n-1)/2
    rowIntersect<-function(x) {
      length(intersect(x, common.genes))
    }
    y<-apply(dcls, 1, rowIntersect)
    e<-sum(y==2)
    #binom.test(e, N, prob, alternative="greater")
    p.value <- pbinom(e-1, N, prob, lower.tail = FALSE, log.p = FALSE)
    tmp <- data.frame(regulator_gene_symbol = regulator_gene_symbol,
                      n.target.gene.in.data = n,
                      n.target.pairs = e,
                      n.expected.target.pairs = N * prob,
                      p.value = p.value
    )
    setTxtProgressBar(pb,ic)
    return(tmp)
  }
  nc = 8
  input = 1:length(humanRegnetwork)
  pb = txtProgressBar(min=0,max=length(humanRegnetwork), style=3, initial=0)
  cat("\n")
  result<- mclapply(input, ProcessElement, mc.cores=nc)
  close(pb)
  library(data.table)
  result<-rbindlist(result)
  result<-as.data.frame(result)
  result$p.adjusted<-p.adjust(result$p.value)
  saveRDS(result,paste("earlyAD",names(lateAD_DCGs.filtered)[t],"allTFs_TED_new_pbinom_out.RDS",sep = "_"))
  compute_time[[t]]=proc.time()
  saveRDS(compute_time,paste("lateAD",names(lateAD_diffcoexp_results)[t],"TED_compute_time.RDS",sep = "_"))
}
