#!/bin/bash -l

#SBATCH -A g2020012
#SBATCH -t 1:00:00

module load gcc openmpi
mpirun ./matmult /proj/g2020012/nobackup/matmul_indata/input5716.txt output.txt
