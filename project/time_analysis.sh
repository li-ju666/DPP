#!/bin/bash -l

#SBATCH -A g2020012
#SBATCH -t 17:40:00

module load gcc openmpi
export OMPI_MCA_btl_openib_allow_ib=1

mpirun ./test 10000 result100000
# mpirun ./test 120000 result1200000
# mpirun ./test 150000 result1500000
