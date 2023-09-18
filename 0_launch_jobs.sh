#!/bin/bash

#SBATCH --partition=fasse
#SBATCH --qos=normal
#SBATCH --account=dominici_lab
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name erc_sim
#SBATCH --output sim.out
#SBATCH --mem=9GB
#SBATCH --time=4:00:00
#SBATCH --array=1-1000%50
#SBATCH --mail-user=mcork@g.harvard.edu
#SBATCH --mail-type=ALL
#SBATCH --error=sim.err

echo "Job script started"

module load gcc/9.5.0-fasrc01 R/4.2.2-fasrc01
if [ $? -ne 0 ]; then
  echo "Module load failed"
  exit 1
fi

echo "Modules loaded successfully"

R CMD BATCH --no-save "--args $1 $2" 1_run_simulation.R Logs/${1}_${2}_run${SLURM_ARRAY_TASK_ID}.Rout