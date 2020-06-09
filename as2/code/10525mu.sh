#!/bin/bash -l

#SBATCH -A g2020012
#SBATCH -t 20:00:00

module load gcc openmpi
mpirun ./matmult /proj/g2020012/nobackup/matmul_indata/input10525.txt output.txt

