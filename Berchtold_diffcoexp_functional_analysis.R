library(org.Hs.eg.db)
library(DCGL)
library(dplyr)
library(magrittr)
library(enrichR)


#Read diffcoexp resultslist
setwd("/Users/sandeepamberkar/Work/Data/AD_GSE48350/")
gse48350_diffcoexp_files=list.files(path = ".",pattern = "*diffcoexp.RDS",full.names = T)
regnet_tf2target.HGNC=fread("/Users/sandeepamberkar/Work/Data/TF_Databases/RegNetwork_human_regulators2.txt",header = T,sep = "\t",showProgress = T,data.table = F)%>%filter(evidence=="Experimental")%>%dplyr::select(c(regulator_symbol,target_symbol))
gse48350_brain_regions.diffcoexp=vector(mode = "list",length = 4)
names(gse48350_brain_regions.diffcoexp)=gsub(pattern = "./|diffcoexp.RDS",replacement = "",x = gse48350_diffcoexp_files)
for(f in 1:4){
  gse48350_brain_regions.diffcoexp[[f]]=readRDS(gse48350_diffcoexp_files[f])  
}
gse48350.DCGs=lapply(gse48350_brain_regions.diffcoexp,function(x)x$DCGs)
gse48350.DCLs=lapply(gse48350_brain_regions.diffcoexp,function(x)x$DCLs)
gse48350.DRsort=vector(mode = "list",length = 4)
names(gse48350.DRsort)=names(gse48350_brain_regions.diffcoexp)
gse48350.DRsort=mapply(FUN = function(a,b)DRsort(DCGs = a,DCLs = b,tf2target = regnet_tf2target.HGNC,expGenes = rownames(gse48350_exprs.agg)),gse48350.DCGs,gse48350.DCLs)

