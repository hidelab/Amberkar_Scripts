

emory_data=read.delim("~/Work/Data/EMORY/AMP-AD_Emory_Emory_Protein.tsv",sep="\t",header=T,as.is = T)
refseq_map=read.table("~/Work/Data/EMORY/RefseqP_mapping.txt",sep = "\t",header = T,as.is = T)
emory_data=emory_data[grep("^NP",emory_data$Reference),]
emory_data$Reference=gsub(pattern = "\\.[[:digit:]]{1,}",replacement = "",emory_data$Reference)
index=match(emory_data$Reference,refseq_map$RefSeq.Protein.ID..e.g..NP_001005353.)[is.na(match(emory_data$Reference,refseq_map$RefSeq.Protein.ID..e.g..NP_001005353.))==F]
emory_data2=data.frame(emory_data$Group[1:length(index)],refseq_map[index,1],emory_data[1:length(index),2:21],stringsAsFactors = F)
rsub=apply(pd,2,function(column)all(column>1))

# Analysis with modified dataset
emory_data=read.delim("AMP-AD_Emory_Emory_Protein_modified.tsv",sep = "\t",header = T,as.is = T)
#Subset AD samples
emory_ad=emory_data[,c(3,15,16,grep("AD",colnames(emory_data)))]
emory_ad2=emory_ad[,-1]
rownames(emory_ad2)=emory_ad$Reference
colnames(emory_ad2)=c("Control.Pool1","Control.Pool2","AD.PD.Pool1","AD.PD.Pool2","AD.Pool1","AD.Pool2")
emory_design=data.frame(Sample=colnames(emory_ad2),
                        NeuroPath=factor(c(rep("Control",2),rep("AD.PD",2),rep("AD",2))),
                        Pool=factor(c(rep(c("Pool1","Pool2"),3))))
emory_y=DGEList(emory_ad2,genes = rownames(emory_ad2))
NeuroPath=emory_design$NeuroPath
Pool=emory_design$Pool
emory_design2=model.matrix(~NeuroPath+Pool)
rownames(emory_design2)=colnames(emory_ad2)
#Compute normalised counts using edgeR, as spectral peptide counts are 
#analogous to raw counts from RNAseq experiment
emory_y=DGEList(emory_ad2,genes = rownames(emory_ad2))
emory_y=calcNormFactors(emory_y,method = "TMM")
emory_ad2.counts=cpm(emory_y,normalized.lib.sizes = T,log = T)
emory_y=estimateDisp(emory_y,design = emory_design2)
plotMDS(emory_y,col=c("Green","Green","Red","Red","Blue","Blue"))
emory_ad2.fit=glmFit(emory_y,design = emory_design2)
emory_ad2.lrt=glmLRT(emory_ad2.fit)

#BEgin WGCNA analysis
emory_multiExpr=vector(mode="list",length=3)
emory_multiExpr=vector(mode="list",length=3)
emory_multiExpr[[1]]=list(data=t(emory_ad2.counts[,c(1:2)]))
emory_multiExpr[[2]]=list(data=t(emory_ad2.counts[,c(3:4)]))
emory_multiExpr[[3]]=list(data=t(emory_ad2.counts[,c(5:6)]))
names(emory_multiExpr)=c("Control","AD.PD","AD")
emory_exprSize=checkSets(emory_multiExpr)
ylim = matrix(NA, nrow = 3, ncol = 4);
for (set in 1:nSets){
  for (col in 1:length(plotCols)){
    ylim[1, col] = min(ylim[1, col], emory_powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE)
    ylim[2, col] = max(ylim[2, col], emory_powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
  }
}
for (col in 1:length(plotCols)) for (set in 1:3){
  if (set==1){
    plot(emory_powerTables[[set]]$data[,1], -sign(emory_powerTables[[set]]$data[,3])*emory_powerTables[[set]]$data[,2],xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],main = colNames[col]);
    addGrid();
    }
  if (col==1){
    text(emory_powerTables[[set]]$data[,1], -sign(emory_powerTables[[set]]$data[,3])*emory_powerTables[[set]]$data[,2],labels=powers,cex=cex1,col=colors[set]);
  }
  else
    text(emory_powerTables[[set]]$data[,1], emory_powerTables[[set]]$data[,plotCols[col]],labels=powers,cex=cex1,col=colors[set]);
  if (col==1)
  {
    legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
  } else
    legend("topright", legend = setLabels, col = colors, pch = 20) ;
}





















































pd2_GS1=as.numeric(cor(pd2_trait,pd2,use="p"))
pd2_p.std=corPvalueFisher(pd2_GS1,nSamples = length(pd2_trait))
pd2_p.std2=pd2_p.std
pd2_p.std2[is.na(pd2_p.std)]=1
pd2_q.std=qvalue(pd2_p.std2)$qvalues
pd2_StdGeneScreeningResults=data.frame(colnames(pd2),PearsonCorrelation=pd2_GS1,pd2_p.std,pd2_q.std)
pd3=as.data.frame(lapply(pd2[,1:dim(pd2)[2]],as.numeric))
pd3_ADJ=abs(cor(pd3,use = "p"))^6



