#!/bin/bash -l

#SBATCH -A g2020012
#SBATCH -t 00:15:00
#SBATCH -p core
#SBATCH -n 2

module load gcc openmpi
export OMPI_MCA_btl_openib_allow_ib=1

echo strong2
mpirun ./simulation 50000 result
mpirun ./simulation 50000 result
mpirun ./simulation 50000 result
mpirun ./simulation 50000 result
mpirun ./simulation 50000 result

# mpirun ./test 120000 result1200000
# mpirun ./test 150000 result1500000
