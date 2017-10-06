#!/bin/bash 
#Job time 
#$-l h_rt=08:00:00 
#Request resources 
#$-l rmem=64G 
#$-l mem=64G
#Job Name 
#$-N MSBB_IFG_DCe
#Queue
#$-P hidelab

module load apps/R/3.3.1
module load compilers/gcc/6.2
R CMD BATCH --no-save --no-restore /shared/hidelab2/user/md4zsa/Work/Amberkar_Scripts/msbb_ifg_DCL_analysis.R