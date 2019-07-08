#!/bin/bash
#PBS -l nodes=1:ppn=24,walltime=00:10:00
#PBS -q molssi
#PBS -N example1
#PBS -j oe

source /gpfs/projects/molssi/modules-intel
cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=4

mpirun -n 1 ./md > output
