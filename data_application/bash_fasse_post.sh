#!/bin/bash

#SBATCH --partition=fasse
#SBATCH --qos=normal
#SBATCH --account=dominici_lab
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=1
#SBATCH --job-name erc_data_app_fasse_post
#SBATCH --output data_app_fasse_post.out
#SBATCH --mem=184GB
#SBATCH --time=07-00:00
#SBATCH --mail-user=mcork@g.harvard.edu
#SBATCH --mail-type=ALL

module load gcc/9.5.0-fasrc01 R/4.2.2-fasrc01

R CMD BATCH --no-save "--args $1 $2" 2_post_processing.R Logs/data_app_fasse_post.Rout
