#!/bin/bash -l

#SBATCH -A g2020012
#SBATCH -t 05:00 --qos=short

module load gcc openmpi
mpirun ./test /proj/g2020012/nobackup/matmul_indata/input3600.txt



