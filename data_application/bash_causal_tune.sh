#!/bin/bash

#SBATCH --partition=fasse
#SBATCH --qos=normal
#SBATCH --account=dominici_lab
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --ntasks-per-node=1
#SBATCH --job-name causal_tune_fasse
#SBATCH --output causal_tune.out
#SBATCH --mem=184GB
#SBATCH --time=02-00:00
#SBATCH --mail-user=mcork@g.harvard.edu
#SBATCH --mail-type=ALL

module load gcc/9.3.0-fasrc01 R/4.0.5-fasrc02

R CMD BATCH --no-save tune_CausalGPS.R Logs/causal_tune.Rout
