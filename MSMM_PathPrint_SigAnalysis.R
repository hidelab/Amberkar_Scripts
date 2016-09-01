getEntropy <- function(mat, index){
  if (index > 2 | index < 1)
    stop("Indicate 1 for rows or 2 for columns")
  d <- apply(as.matrix(mat), index, function(x){discretize(x, numBins = 3, r=c(-1,1))})
  entropy.vec <- apply(d, 2, entropy)
  return(entropy.vec)
}
pvals=matrix(NA,nrow=nrow(POE.msmm.rnaSeq.matrix.STG),ncol=500)
sample_reps=replicate(500,sample(1:60,60,replace=F))
for (i in 1:dim(pvals)[1]){
  for (j in 1:dim(sample_reps)[2]){
    pvals[i,j]=wilcox.test(POE.msmm.expr.stg[i,],POE.msmm.rnaSeq.matrix.STG[i,sample_reps[j]],paired = F,correct = T,conf.level = 0.05)$p.value
  }
}

entropy.vector.PHG.col=getEntropy(POE.msmm.rnaSeq.matrix.PHG,2)
entropy.vector.PHG.row=getEntropy(POE.msmm.rnaSeq.matrix.PHG,1)
entropy.vector.PHG.col=getEntropy(POE.msmm.rnaSeq.matrix.PHG,2)
entropy.vector.PHG.row=getEntropy(POE.msmm.rnaSeq.matrix.PHG,1)
entropy.vector.PHG.col=getEntropy(POE.msmm.rnaSeq.matrix.PHG,2)
entropy.vector.PHG.row=getEntropy(POE.msmm.rnaSeq.matrix.PHG,1)
high.entropy.PHG.col=names(entropy.vector.PHG.col)[which(entropy.vector.PHG.col >= 0.25)]
high.entropy.PHG.row=names(entropy.vector.PHG.row)[which(entropy.vector.PHG.row >= 0.4)]
phm_msmm_rnaseq_PHG=pheatmap(POE.msmm.rnaSeq.matrix.PHG[high.entropy.PHG.row,high.entropy.PHG.col],
                                                          color=c("Blue","White","Red"),
                                                          cluster_rows=T,cluster_cols=T,
                                                          clustering_distance_rows = "manhattan",
                                                          clustering_distance_cols = "manhattan",
                                                          legend_breaks = c(-1, 0, 1),
                                                          fontsize = 9,fontsize_row = 9,fontsize_col = 9,
                                                          cellwidth = 9,cellheight = 9,
                                                          main="PathPrint - MSMM.PHG samples",
                                                          annotation=msmm_annotation_PHG)
entropy.STG.row=getEntropy(POE.msmm.rnaSeq.matrix.STG,1)
entropy.STG.col=getEntropy(POE.msmm.rnaSeq.matrix.STG,2)
high.entropy.STG.row=names(entropy.STG.row)[which(entropy.STG.row >= 0.5)]
high.entropy.STG.col=names(entropy.STG.col)[which(entropy.STG.col >= 0.2)]
phm_msmm_rnaseq_STG=pheatmap(POE.msmm.rnaSeq.matrix.STG[high.entropy.STG.row,],
                             color=c("Blue","White","Red"),
                             cluster_rows=T,cluster_cols=T,
                             clustering_distance_rows = "manhattan",
                             clustering_distance_cols = "manhattan",
                             legend_breaks = c(-1, 0, 1),
                             fontsize = 8.5,fontsize_row = 8.5,fontsize_col = 8.5,
                             cellwidth = 8.5,cellheight = 8.5,
                             main="PathPrint - MSMM.STG samples",
                             annotation=msmm_annotation_STG)
entropy.FP.row=getEntropy(POE.msmm.rnaSeq.matrix.FP,1)
entropy.FP.col=getEntropy(POE.msmm.rnaSeq.matrix.FP,2)
high.entropy.FP.row=names(entropy.FP.row)[which(entropy.FP.row >= 0.3)]
high.entropy.FP.col=names(entropy.FP.col)[which(entropy.FP.col >= 0.3)]
phm_msmm_rnaseq_FP=pheatmap(POE.msmm.rnaSeq.matrix.FP[high.entropy.FP.row,high.entropy.FP.col],
                             color=c("Blue","White","Red"),
                             cluster_rows=T,cluster_cols=T,
                             clustering_distance_rows = "manhattan",
                             clustering_distance_cols = "manhattan",
                             legend_breaks = c(-1, 0, 1),
                             fontsize = 9,fontsize_row = 9,fontsize_col = 9,
                             cellwidth = 9,cellheight = 9,
                             main="PathPrint - MSMM.FP samples",
                             annotation=msmm_annotation_FP)
