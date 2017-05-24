#!/bin/bash 
#Job time 
#-l h_rt=120:00:00 
#Request resources 
#-l rmem=16G 
#Cores 
#$ -pe openmp 8 
#Job Name 
#-N MSMM_ReSeq_PHG_DCN
module load module load apps/R/3.3.1
module load module load compilers/gcc/6.2
R CMD BATCH --no-save --no-restore /shared/hidelab2/user/md4zsa/Work/Amberkar_Scripts/mssm_reseq_PHG_DCN_analysis.R