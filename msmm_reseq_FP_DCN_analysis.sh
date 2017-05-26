#!/bin/bash 
#Job time 
#-l h_rt=120:00:00 
#Request resources 
#-l rmem=12G 
#Cores 
#$ -pe openmp 8 
#Queue
# -P hidelab
#Job Name 
#-N "MSMM_ReSeq_FP_DCN"
module load apps/R/3.3.1
module load compilers/gcc/6.2
R CMD BATCH --no-save --no-restore /shared/hidelab2/user/md4zsa/Work/Amberkar_Scripts/msmm_reseq_FP_DCN_analysis.R