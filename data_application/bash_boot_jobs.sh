#!/bin/bash

#SBATCH --partition=fasse
#SBATCH --qos=normal
#SBATCH --account=dominici_lab
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name erc_data_boot
#SBATCH --output boot.out
#SBATCH --mem-per-cpu=25GB
#SBATCH --time=09:00:00
#SBATCH --array=101-1000%100
#SBATCH --mail-user=mcork@g.harvard.edu
#SBATCH --mail-type=ALL

module load gcc/9.5.0-fasrc01 R/4.2.2-fasrc01

R CMD BATCH --no-save "--args $1 $2" 1_model_fit_boot.R Logs/${1}_${2}_run${SLURM_ARRAY_TASK_ID}.Rout