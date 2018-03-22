#!/bin/bash 
#Job time 
#$-l h_rt=24:00:00 
#Request resources 
#$-l mem=256G
#$-l rmem=128G
#Job Name 
#$-N MSBB_TED
#Queue
#$-P hidelab


module load apps/R/3.3.1
module load compilers/gcc/6.2
R CMD BATCH --no-save --no-restore /shared/hidelab2/user/md4zsa/Work/Amberkar_Scripts/msbb_array_diffcoexp2_TED_analysis.R
