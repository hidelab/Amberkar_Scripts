#!/bin/bash 
#Job time 
#-l h_rt=96:00:00 
#Request resources 
#-l rmem=8G 
#Cores 
#$ -pe openmp 16 
#Job Name 
#-N MSBB_SpTp
module load dev/intel-compilers/17.0.0
module load apps/R/3.4.0/intel-17.0-parallel
R CMD BATCH --no-save --no-restore /shared/hidelab2/user/md4zsa/Work/Amberkar_Scripts/msmm_array19_pathprint_fingerprint_SpTp.R