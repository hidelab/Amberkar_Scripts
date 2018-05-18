#!/bin/bash 
#Job time 
#$-l h_rt=24:00:00 
#Request resources 
#$-l mem=128G
#$-l rmem=128G
#Job Name 
#$-N EarlyCDR_TED
#Queue
#$-P rse


module load apps/R/3.4.0/gcc-4.8.5
module load mpi/openmpi/2.1.1/gcc-6.2
R CMD BATCH --no-save --no-restore /shared/hidelab2/user/md4zsa/Work/Amberkar_Scripts/msbb_array_diffcoexp_earlyCDR_TDD.sh
