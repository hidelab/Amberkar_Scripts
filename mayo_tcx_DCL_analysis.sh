#!/bin/bash 
#Job time 
#$-l h_rt=120:00:00 
#Request resources 
#$-l rmem=12G 
#Cores 
#$-pe openmp 12 
#Job Name 
#$-N MAYO_TCX_DCe
#Queue
#$-P rse

module load apps/R/3.4.0/gcc-4.8.5
module load mpi/openmpi/2.1.1/gcc-4.8.5
R CMD BATCH --no-save --no-restore /shared/hidelab2/user/md4zsa/Work/Amberkar_Scripts/mayo_tcx_DCL_analysis.R