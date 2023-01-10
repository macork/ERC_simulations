#!/bin/bash

#SBATCH --partition=fasse
#SBATCH --qos=normal
#SBATCH --account=dominici_lab
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name erc_data_app
#SBATCH --output data_post_test.out
#SBATCH --mem=20GB
#SBATCH --time=01-00:00
#SBATCH --mail-user=mcork@g.harvard.edu
#SBATCH --mail-type=ALL

module load gcc/9.3.0-fasrc01 R/4.0.5-fasrc02

R CMD BATCH --no-save "--args $1 $2" post_processing.R Logs/data_post_test.Rout
