#!/bin/bash 
#Job time 
#-l h_rt=120:00:00 
#Request resources 
#-l rmem=12G 
#Cores 
#$ -pe openmp 8 
#Job Name 
#-N MSMM_ReSeq_PHG_DCN
module load apps/R/3.4.0/intel-17.0-parallel
module load dev/intel-compilers/17.0.0
R CMD BATCH --no-save --no-restore /shared/hidelab2/user/md4zsa/Work/Amberkar_Scripts/mssm_reseq_PHG_DCN_analysis.R