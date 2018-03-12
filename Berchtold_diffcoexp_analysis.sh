#!/bin/bash 
#Job time 
#$-l h_rt=12:00:00 
#Request resources 
#$-l mem=64G
#$-l rmem=32G
#$-pe openmp 4
#Job Name 
#$-N Berchtold_diffcoep
#Queue
#-P hidelab

module load apps/R/3.3.1
module load compilers/gcc/6.2
R CMD BATCH --no-save --no-restore /shared/hidelab2/user/md4zsa/Work/Amberkar_Scripts/Berchtold_diffcoexp_analysis.R
