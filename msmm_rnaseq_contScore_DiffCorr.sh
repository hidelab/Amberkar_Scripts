#!/bin/bash
#Job time
#-l h_rt=48:00:00
#Request resources
#-l mem=64G -l rmem=32G
#Queue
#$ -P hidelab
#Cores
#$ -pe openmp 4
#Job Name
#-N msmm_contScore_cocor

module load apps/R/3.3.1
module load compilers/gcc/5.3
R CMD BATCH --no-save --no-restore /shared/hidelab2/user/md4zsa/Work/Amberkar_Scripts/msmm_rnaseq_contScore_DiffCorr.R
