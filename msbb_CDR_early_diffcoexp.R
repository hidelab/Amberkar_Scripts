library(org.Hs.eg.db)
library(dplyr)
library(magrittr)
library(data.table)
library(diffcoexp)

setwd("/shared/hidelab2/user/md4zsa/Work/Data/MSBB_Array19/GSE84422/")
load("msbb_gse84422_Lobe_DEG_DEP.RData")
for(i in 1:length(msbb_gse84422_GPL96_97_byRegion.exprs)){
  c_exprs=msbb_gse84422_GPL96_97_byRegion.exprs[[i]][,msbb_gse84422.pData$GPL96$pseudoSampleID[which(msbb_gse84422.pData$GPL96$`brain region:ch1`==names(msbb_gse84422_GPL96_97_byRegion.exprs)[i]&msbb_gse84422.pData$GPL96$SampleTypeCDR=="CDR0")]]
  d_exprs=msbb_gse84422_GPL96_97_byRegion.exprs[[i]][,msbb_gse84422.pData$GPL96$pseudoSampleID[which(msbb_gse84422.pData$GPL96$`brain region:ch1`==names(msbb_gse84422_GPL96_97_byRegion.exprs)[i]&msbb_gse84422.pData$GPL96$SampleTypeCDR=="CDR1")]]
  diffcoexp_out=diffcoexp(exprs.1 = c_exprs[,],exprs.2 = d_exprs[,],r.method="spearman",rth=0.6,q.diffth=0.1,q.dcgth=0.1)
  fn_prefix=gsub(pattern = " ",replacement = "_",names(msbb_gse84422_GPL96_97_byRegion.exprs)[i])
  saveRDS(diffcoexp_out,paste(fn_prefix,"diffcoexp_CDR0_CDR1.RDS",sep="_"))
  proc.time()
}
