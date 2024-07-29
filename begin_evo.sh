#!/bin/bash
# Specify a partition 
#SBATCH --partition=bluemoon
# Request nodes 
#SBATCH --nodes=1
# Request some processor cores 
#SBATCH --ntasks=20
# Maximum runtime of 2 hours
#SBATCH --time=30:00:00
# Memory per CPU core (optional, if needed)
#SBATCH --mem-per-cpu=4G  # Adjust based on your memory requirements
python begin_evolving.py