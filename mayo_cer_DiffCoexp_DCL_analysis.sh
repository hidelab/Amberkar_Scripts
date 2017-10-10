#!/bin/bash 
#Job time 
#$-l h_rt=12:00:00 
#Request resources 
#$-l rmem=32G
#$-l mem=32G
#Job Name 
#$-N Mayo_CER_DiffCoexp


module load apps/R/3.3.1
module load compilers/gcc/6.2
R CMD BATCH --no-save --no-restore /shared/hidelab2/user/md4zsa/Work/Amberkar_Scripts/mayo_cer_DiffCoexp_DCL_analysis.R