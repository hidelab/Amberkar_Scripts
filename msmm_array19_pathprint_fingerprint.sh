#!/bin/bash 
#Job time 
#-l h_rt=24:00:00 
#Request resources 
#-l mem=8G -l rmem=8G 
#Cores 
#$ -pe openmp 8
#Job Name 
#-N MSBB_CDR_Pathprint
module load apps/R/3.3.2/intel-17.0-parallel
module load dev/intel-compilers/17.0.0
module load mpi/openmpi/2.0.1/intel-17.0.0
R CMD BATCH --no-save --no-restore /shared/hidelab2/user/md4zsa/Work/Amberkar_Scripts/msmm_array19_pathprint_fingerprint.R