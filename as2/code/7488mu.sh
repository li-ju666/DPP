#!/bin/bash -l

#SBATCH -A g2020012
#SBATCH -t 00:15:00 --qos=short

module load gcc openmpi
mpirun ./matmult /proj/g2020012/nobackup/matmul_indata/input7488.txt output.txt

