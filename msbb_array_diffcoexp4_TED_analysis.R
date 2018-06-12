library(dplyr)
library(magrittr)
library(data.table)
library(DCGL)
library(parallel)

setwd("/shared/hidelab2/user/md4zsa/Work/Data/MSBB_Array19/GSE84422/")
earlyAD_diffcoexp_files=list.files(path = ".",pattern = "earlyAD_diffcoexp",full.names = T)
earlyAD_samples=readRDS("./msbb_gse84422_GPL96_97_earlyAD_samplesToAnalyse.RDS")
earlyAD_samples.exprs=readRDS("./msbb_gse84422_GPL96_97_earlyAD_samplesToAnalyse_exprs.RDS")

humanRegnetwork=readRDS("humanRegnetwork.RDS")

earlyAD_diffcoexp_results=vector(mode = "list",length = length(earlyAD_diffcoexp_files))
for(i in 1:length(earlyAD_diffcoexp_files)){
  earlyAD_diffcoexp_results[[i]]=readRDS(earlyAD_diffcoexp_files[[i]])
}
names(earlyAD_diffcoexp_results)=names(earlyAD_samples.exprs)

earlyAD_DCGs=lapply(earlyAD_diffcoexp_results,function(x)x$DCGs)
earlyAD_DCLs=lapply(earlyAD_diffcoexp_results,function(x)x$DCLs)

earlyAD_DRrank.TED=earlyAD_DRrank.TED=vector(mode = "list",length = length(earlyAD_samples.exprs))
names(earlyAD_DRrank.TED)=names(earlyAD_DRrank.TED)=names(earlyAD_diffcoexp_results)
compute_time=list()
for(t in 9:length(earlyAD_diffcoexp_results)){
  genes=earlyAD_DCGs[[t]]$Gene
  dcls=earlyAD_DCLs[[t]]
  prob<-nrow(dcls)/choose(nrow(earlyAD_samples.exprs[[t]]), 2)
  
  
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
  nc = 4
  input = 1:length(humanRegnetwork)
  pb = txtProgressBar(min=0,max=length(humanRegnetwork), style=3, initial=0)
  cat("\n")
  result<- mclapply(input, ProcessElement, mc.cores=nc)
  close(pb)
  library(data.table)
  result<-rbindlist(result)
  result<-as.data.frame(result)
  result$p.adjusted<-p.adjust(result$p.value)
  saveRDS(result,paste("earlyAD",names(earlyAD_diffcoexp_results)[t],"TED_new_pbinom_out.RDS",sep = "_"))
  compute_time[[t]]=proc.time()
}
