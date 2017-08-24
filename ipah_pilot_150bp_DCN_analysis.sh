#!/bin/bash 
#Job time 
#-l h_rt=120:00:00 
#Request resources 
#-l mem=8G 
#Cores 
#$-pe openmp 12 
#Job Name 
#-N IPAH_DCN2

module load mpi/openmpi/2.1.1/gcc-4.8.5
module load apps/R/3.4.0/gcc-4.8.5
R CMD BATCH --no-save --no-restore /shared/hidelab2/user/md4zsa/Work/Amberkar_Scripts/ipah_pilot_150bp_DCN_analysis.R