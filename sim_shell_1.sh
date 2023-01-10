#!/bin/bash

#SBATCH --partition=fasse
#SBATCH --qos=normal
#SBATCH --account=dominici_lab
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --job-name erc_sim
#SBATCH --output sim.out
#SBATCH --mem-per-cpu=4000MB
#SBATCH --time=09:00:00
#SBATCH --array=1-100
#SBATCH --mail-user=mcork@g.harvard.edu
#SBATCH --mail-type=ALL

module load gcc/9.3.0-fasrc01 R/4.0.5-fasrc02

R CMD BATCH --no-save "--args $1 $2" sim_fasse_1.R Logs/${1}_${2}_run${SLURM_ARRAY_TASK_ID}.Rout