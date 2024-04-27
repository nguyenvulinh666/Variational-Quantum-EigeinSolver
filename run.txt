#!/bin/sh

#SBATCH --job-name=VQE_Ising
#SBATCH --partition=small
#SBATCH --ntasks=8
#SBATCH --nodes=4
#SBATCH --tasks-per-node=2
#SBATCH --cpus-per-task 16
#SBATCH --output=/home/ctv.linnv/out/%x.%j.out
#SBATCH --error=/home/ctv.linnv/err/%x.%j.err

module load python/anaconda3

python run_VQE_modified_on_HPC_sample.py

