#!/bin/bash
#SBATCH -A m1709
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1

export NUMEXPR_MAX_THREADS=64
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4

srun python3 process_output.py
