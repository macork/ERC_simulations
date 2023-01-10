#!/bin/bash

#SBATCH --partition=fasse_bigmem
#SBATCH --qos=normal
#SBATCH --account=dominici_lab
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name erc_data_app
#SBATCH --output data_app.out
#SBATCH --mem=499GB
#SBATCH --time=07-00:00
#SBATCH --mail-user=mcork@g.harvard.edu
#SBATCH --mail-type=ALL

module load gcc/9.3.0-fasrc01 R/4.0.5-fasrc02

R CMD BATCH --no-save "--args $1 $2" model_fit.R Logs/data_app.Rout
