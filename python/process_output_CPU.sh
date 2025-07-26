#!/bin/bash
#SBATCH -A m1709
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -t 10:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1

srun python3 process_output.py
