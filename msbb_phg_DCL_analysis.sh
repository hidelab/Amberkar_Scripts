#!/bin/bash 
#Job time 
#$-l h_rt=08:00:00 
#Request resources 
#$-l rmem=12G 
#$-l mem=12G
#Cores 
#$-pe openmp 8 
#Job Name 
#$-N MSBB_PHG_DCe


module load apps/R/3.3.1
module load compilers/gcc/6.2
R CMD BATCH --no-save --no-restore /shared/hidelab2/user/md4zsa/Work/Amberkar_Scripts/msbb_phg_DCL_analysis.R