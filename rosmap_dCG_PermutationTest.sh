#!/bin/bash 
#Job time 
#$-l h_rt=240:00:00 
#Request resources 
#$-l mem=12G 
#Cores 
#$-pe openmp 12
#Job Name 
#$-N ROSMAP_dCP
#Job queue
#-P hidelab

module load compilers/gcc/6.2
module load apps/R/3.3.1
R CMD BATCH --no-save --no-restore /shared/hidelab2/user/md4zsa/Work/Amberkar_Scripts/rosmap_dCG_PermutationTest.R