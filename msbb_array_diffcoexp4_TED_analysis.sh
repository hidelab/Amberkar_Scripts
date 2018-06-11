#!/bin/bash 
#Job time 
#$ -l h_rt=00:01:00 
#Request resources 
#$ -l rmem=64G
#$ -p smp 4
#Job Name 
#$ -N EarlyCDR_TED
#$ -o msbb_array_diffcoexp4_TED_analysis.out
#Queue
#$ -P rse



module load apps/R/3.4.0/gcc-4.8.5
module load mpi/openmpi/2.1.1/gcc-6.2
R CMD BATCH --no-save --no-restore /shared/hidelab2/user/md4zsa/Work/Amberkar_Scripts/msbb_array_diffcoexp4_TED_analysis.R
