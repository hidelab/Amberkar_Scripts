library(data.table)
library(dplyr)
library(magrittr)
library(devtools)
library(DESeq2)
library(EnsDb.Hsapiens.v79)
library(tibble)

rosmap_reseq_counts.res=readRDS("rosmap_reseq_counts_res.RDS")
rosmap_reseq_DEGs.list=lapply(rosmap_reseq_counts.res,function(x)data.frame(x))
rosmap_reseq_DEGs.list=lapply(rosmap_reseq_counts.res,function(x)data.frame(x)%>%rownames_to_column("Gene"))

rosmap_reseq_filtered_DEGs.list=lapply(rosmap_reseq_DEGs.list,function(x)x%>%dplyr::filter(padj<=0.1))
for(i in 1:3){
  fwrite(rosmap_reseq_filtered_DEGs.list[[i]],paste("rosmap",names(rosmap_reseq_filtered_DEGs.list)[i],"DEGs.txt",sep = "_"),sep = "\t",col.names = T)
}
