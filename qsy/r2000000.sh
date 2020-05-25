#!/bin/bash -l

#SBATCH -A g2020012
#SBATCH -t 15:00 --qos=short

module load gcc openmpi
mpirun ./mc 2000000 output.txt

