#!/bin/bash 
#Job time 
#$-l h_rt=240:00:00 
#Request resources 
#$-l mem=256G
#$-l rmem=256G
#Job Name 
#$-N IPAH_150bp_DC

module load apps/R/3.3.1
module load compilers/gcc/6.2
R CMD BATCH --no-save --no-restore /shared/hidelab2/user/md4zsa/Work/Amberkar_Scripts/ipah_pilot_150bp_DCN_analysis.R