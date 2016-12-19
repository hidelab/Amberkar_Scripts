#!/bin/bash
#$ -q long   # <- the name of the Q you want to submit to
#$ -rmem 32G # <- Real memory request
#$ -mem 64G # <- Virtual memory request
#$ -binding linear:16 #  <- seek 16 cores
#$ -S /bin/bash   # <- run the job under bash
#$ -N TDP43_RNAseq # <- name of the job in the qstat output
#$ -o TDP43_RNAseq.out # <- name of the output file.
#$ -e TDP43_RNAseq.stderr # <- name of the stderr file.
module load apps/R/3.2.0
module load compilers/intel/15.0.3
R CMD BATCH TDP_Align.R
