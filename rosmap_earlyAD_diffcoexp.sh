#!/bin/bash 
#Job time 
#$-l h_rt=12:00:00 
#Request resources 
#$-l mem=512G
#$-l rmem=512G
#Job Name 
#$-N ROSMAP_EarlyAD_DCN
#Queue
#$-P rse
#$ -m bea # send mails at beginning, end and if aborted unexpectedly
#$ -M s.amberkar@Sheffield.ac.uk # mail sent to this address
#$ -o ROSMAP_EarlyAD_DCN.log # output file
#$ -j y # send output and error to same file

module load apps/R/3.4.0/gcc-4.8.5
module load mpi/openmpi/2.1.1/gcc-6.2
R CMD BATCH --no-save --no-restore /shared/hidelab2/user/md4zsa/Work/Amberkar_Scripts/rosmap_earlyAD_diffcoexp_analysis.R