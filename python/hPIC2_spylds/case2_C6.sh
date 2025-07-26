#!/bin/bash
#SBATCH -A m1709
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -t 03:30:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1


srun python3 case2_C6.py
