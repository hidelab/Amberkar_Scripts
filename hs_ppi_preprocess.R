library(org.Hs.eg.db)
library(data.table)
library(igraph)

#Read PPI data, Irefindex v14.0 (http://irefindex.org/wiki/index.php?title=README_MITAB2.6_for_iRefIndex)
#PPI data has been filtered for duplicates, non-human PPI and interactions without a curated gene ID
setwd("/shared/hidelab2/shared/")
hs_ppi_table=fread('iref14_Human_UP_noDup_table.txt',sep = '\t',header = T,data.table = T,showProgress = T)
#Remove isoforms and collapse to gene level
hs_ppi.Isoform_A=grep(pattern = '-',x = hs_ppi_table$V1,value = F)
names(hs_ppi.Isoform_A)=grep(pattern = '-',x = hs_ppi_table$V1,value = T)
hs_ppi.Isoform_B=grep(pattern = '-',x = hs_ppi_table$V2,value = F)
names(hs_ppi.Isoform_B)=grep(pattern = '-',x = hs_ppi_table$V2,value = T)
hs_ppi_table$V1[hs_ppi.Isoform_A]=gsub(pattern = '-\\d+',replacement='',x = grep(pattern = '-',x = names(hs_ppi.Isoform_A),value = T))
hs_ppi_table$V2[hs_ppi.Isoform_B]=gsub(pattern = '-\\d+',replacement='',x = grep(pattern = '-',x = names(hs_ppi.Isoform_B),value = T))
#Map Uniprot IDs to Gene symbols
hs_ppi_table$Gene.A=unname(mapIds(x = org.Hs.eg.db,keys = hs_ppi_table$V1,column = 'SYMBOL',keytype = 'UNIPROT'))
hs_ppi_table$Gene.B=unname(mapIds(x = org.Hs.eg.db,keys = hs_ppi_table$V2,column = 'SYMBOL',keytype = 'UNIPROT'))
hs_ppi=graph.data.frame(d = hs_ppi_table[which((is.na(hs_ppi_table$Gene.A)==F)&(is.na(hs_ppi_table$Gene.B)==F)),c(4:5)],directed = F)