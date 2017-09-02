#!/bin/bash

#SBATCH -p t2standard
#SBATCH --ntasks=72
#SBATCH --ntasks-per-node=24
#SBATCH -t 20

#SBATCH --output=OUTPUT_FILES/%j.o
#SBATCH --job-name=go_wfdiff

cd $SLURM_SUBMIT_DIR

source activate sln2
mpirun -np $SLURM_NTASKS python run_wfdiff_alaska.py