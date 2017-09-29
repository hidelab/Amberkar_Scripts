#!/bin/bash 
#Job time 
#-l h_rt=96:00:00 
#Request resources 
#-l rmem=12G 
#Cores 
#$-pe openmp 12 
#Job Name 
#$ -N MAYO_TCX_DCe

module load apps/R/3.3.1
module load compilers/gcc/6.2
R CMD BATCH --no-save --no-restore /shared/hidelab2/user/md4zsa/Work/Amberkar_Scripts/mayo_tcx_DCL_analysis.R