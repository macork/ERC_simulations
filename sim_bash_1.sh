#!/bin/bash

#SBATCH --partition=fasse
#SBATCH --qos=normal
#SBATCH --account=dominici_lab
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --job-name erc_sim
#SBATCH --output sim.out
#SBATCH --mem-per-cpu=3000MB
#SBATCH --time=02:00:00
#SBATCH --array=1-50
#SBATCH --mail-user=mcork@g.harvard.edu
#SBATCH --mail-type=ALL

ml_r4

R CMD BATCH --no-save "--args $1 $2" sim_fasse1.R Logs/exp_${1}_adjust${2}_run${SLURM_ARRAY_TASK_ID}.Rout
