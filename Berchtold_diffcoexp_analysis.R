library(GEOquery)
library(diffcoexp)
library(org.Hs.eg.db)
library(DCGL)



mapIds2<-function(IDs,IDFrom,IDTo){
  require(org.Hs.eg.db)
  idmap=mapIds(x = org.Hs.eg.db,keys = IDs,column = IDTo,keytype = IDFrom,multiVals = "first")
  na_vec=names(idmap[is.na(idmap)==T])
  idmap=idmap[is.na(idmap)==F]
  idmap_df=data.frame("From"=names(idmap),"To"=unlist(unname(idmap)),stringsAsFactors = F)
  return(list(map=idmap_df,noMap=na_vec))
}
jaccard=function(A,B){
  jc=set_cardinality(intersect(A,B))/set_cardinality(union(A,B))
  return(jc)
}

setwd("/shared/hidelab2/user/md4zsa/Work/Data/AD_GSE48350")
# gse48350_series_matrix=getGEO(GEO = "GSE48350",GSEMatrix = T)
# gse48350_exprs=cbind.data.frame(pData(gse48350_series_matrix@featureData[,c(1,11:12)]),exprs(gse48350_series_matrix))
# gse48350_exprs.agg=aggregate(x = gse48350_exprs[,-c(1:3)],by=list(symbol=gse48350_exprs$`Gene Symbol`),mean)
# rownames(gse48350_exprs.agg)=gse48350_exprs.agg$symbol
# gse48350_exprs.agg=gse48350_exprs.agg[-1,]
# gse48350_exprs.agg=gse48350_exprs.agg[,-1]
gse48350_exprs.agg=readRDS("GSE48350_exprs.RDS")
gse48350_phenodata=readRDS("GSE48350_pData.RDS")

gse48350_phenodata$`brain region:ch1`=gsub(pattern = "^ ",replacement = "",x = gse48350_phenodata$`brain region:ch1`)
gse48350_phenodata$`brain region:ch1`=gsub(pattern = "-",replacement = "",x = gse48350_phenodata$`brain region:ch1`)
gse48350_AD_samples=gse48350_phenodata$geo_accession[grep(pattern = "AD",gse48350_phenodata$title,value = F,invert = F)]
gse48350_Control_samples=gse48350_phenodata$geo_accession[grep(pattern = "AD",gse48350_phenodata$title,value = F,invert = T)]
gse48350_total_samples=c(gse48350_Control_samples,gse48350_AD_samples)
names(gse48350_total_samples)=gse48350_phenodata$`brain region:ch1`
gse48350_brain_regions=unique(gse48350_phenodata$`brain region:ch1`)

gse48350_brain_regions.diffcoexp=gse48350_brain_regions.exprs=vector(mode = "list",length = 4)
names(gse48350_brain_regions.diffcoexp)=names(gse48350_brain_regions.exprs)=gse48350_brain_regions
names(gse48350_brain_regions.diffcoexp)=gsub(pattern = " ",replacement = "_",x = names(gse48350_brain_regions.diffcoexp))
for(i in 1:4){
  gse48350_brain_regions.exprs[[i]]=gse48350_exprs.agg[,(gse48350_total_samples[names(gse48350_total_samples)==gse48350_brain_regions[i]])]
  c_exprs=gse48350_brain_regions.exprs[[i]][,colnames(gse48350_brain_regions.exprs[[i]])%in%gse48350_Control_samples]
  t_exprs=gse48350_brain_regions.exprs[[1]][,colnames(gse48350_brain_regions.exprs[[1]])%in%gse48350_AD_samples]
  gse48350_brain_regions.diffcoexp[[i]]==diffcoexp(exprs.1=c_exprs,exprs.2=t_exprs,r.method="spearman",rth=0.6,q.diffth=0.1,q.dcgth=0.1)
  cat(paste("Processing region ",names(gse48350_brain_regions.diffcoexp)[i],"...\n",sep = ""))
  saveRDS(gse48350_brain_regions.diffcoexp[[i]],paste(names(gse48350_brain_regions.diffcoexp)[i],"diffcoexp.RDS",sep = ""))
  cat(paste("Done!\n"))
}
proc.time()




