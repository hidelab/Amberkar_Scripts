#!/bin/bash 
#Job time 
#$-l h_rt=96:00:00 
#Request resources 
#$-l mem=12G -l rmem=8G 
#Cores 
#$ -pe openmp 12
#Job Name 
#$-N MSBB2016_DiffCorr
module load apps/R/3.3.1
module load compilers/gcc/5.3 
R CMD BATCH --no-save --no-restore /shared/hidelab2/user/md4zsa/Work/Amberkar_Scripts/msbb_rnaseq2016_diffcorr_analysis.R