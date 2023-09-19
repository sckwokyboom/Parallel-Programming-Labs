#!/bin/bash

#PBS -l select=2:ncpus=8:mpiprocs=4:mem=16000m
#PBS -l walltime=00:05:00
#PBS -l place=scatter
#PBS -m n
#PBS -o log.txt
#PBS -e err.txt

MPI_NP=$(wc -l $PBS_NODEFILE | awk '{ print $1 }')

cd $PBS_O_WORKDIR

echo "Node file path: $PBS_NODEFILE"
echo "Node file contents:"
cat $PBS_NODEFILE

echo "Using mpirun at `which mpirun`"
echo "Running $MPI_NP MPI processes"

mpirun -machinefile $PBS_NODEFILE -np $MPI_NP ./linpack
