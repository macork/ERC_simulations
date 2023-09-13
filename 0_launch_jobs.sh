#!/bin/bash

#SBATCH --partition=test
#SBATCH --qos=normal
#SBATCH --account=dominici_lab
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name erc_sim
#SBATCH --output sim.out
#SBATCH --mem-per-cpu=4000MB
#SBATCH --time=09:00:00
#SBATCH --array=1-2
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