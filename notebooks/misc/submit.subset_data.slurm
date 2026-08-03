#!/bin/bash
#SBATCH --job-name="subset_WeddellGyre"
#SBATCH --time=00:30:00
#SBATCH --mem=128G
#SBATCH --account=co2clim
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH -o %j.out
#SBATCH -e %j.err

script_args=("$@")
set --

module load jaspy
python subset_data.py