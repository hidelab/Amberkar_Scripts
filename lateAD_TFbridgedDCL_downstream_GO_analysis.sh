#!/bin/bash 
#Job time 
#$ -l h_rt=08:00:00 
#Request resources 
#$ -l mem=32G
#$ -l rmem=32G
#Job Name 
#$ -N lateAD_TFbridgedDCL_GO
#$ -o lateAD_TFbridgedDCL_downstream_GO_analysis.out
#Queue
#$ -m bea # send mails at beginning, end and if aborted unexpectedly
#$ -M s.amberkar@Sheffield.ac.uk # mail sent to this address
#$ -o lateAD_TFbridgedDCL_GO.log # output file
#$ -j y # send output and error to same file



module load apps/R/3.3.1
module load compilers/gcc/6.2
R CMD BATCH --no-save --no-restore /shared/hidelab2/user/md4zsa/Work/Amberkar_Scripts/lateAD_TFbridgedDCL_downstream_GO_analysis.R