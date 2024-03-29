#!/bin/bash
#SBATCH -J BV1
#SBATCH --account=def-ubcxzh
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH -t 0-3:00

module load r/4.1.2
Rscript --max-ppsize=500000  $1.r $2 $3