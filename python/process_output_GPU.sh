#!/bin/bash
#SBATCH -A m1709
#SBATCH -C gpu
#SBATCH -q regular
#SBATCH -t 09:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -G 1
#SBATCH --gpus-per-task=1
#SBATCH --gpu-bind=none

srun python3 process_output.py
