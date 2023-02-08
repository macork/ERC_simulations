#!/bin/bash

#SBATCH --partition=fasse
#SBATCH --qos=normal
#SBATCH --account=dominici_lab
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --ntasks-per-node=1
#SBATCH --job-name erc_data_app_fasse
#SBATCH --output data_app_fasse.out
#SBATCH --mem=184GB
#SBATCH --time=07-00:00
#SBATCH --mail-user=mcork@g.harvard.edu
#SBATCH --mail-type=ALL

module load gcc/9.3.0-fasrc01 R/4.0.5-fasrc02

R CMD BATCH --no-save "--args $1 $2 $3 $4" 1_model_fit.R Logs/data_app_fasse_$2.Rout
