#!/bin/bash 
#Job time 
#$-l h_rt=120:00:00 
#Request resources 
#$-l mem=12G -l rmem=12G 
#Cores 
#$ -pe openmp   10
#Queue
#$ -P hidelab
#Job Name 
#$-N MSBB_PLQ_DiffCorr
#Email
#$-M s.amberkar@sheffield.ac.uk

module load apps/R/3.3.1
module load compilers/gcc/5.3 
R CMD BATCH --no-save --no-restore /shared/hidelab2/user/md4zsa/Work/Amberkar_Scripts/msbb_rnaseq2016_PLQGenes_diffcorr_analysis.R