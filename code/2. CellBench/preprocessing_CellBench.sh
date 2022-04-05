#!/bin/bash
#SBATCH -J BV1
#SBATCH --account=rrg-ubcxzh
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH -t 0-1:00

module load r/4.1.2
Rscript --max-ppsize=500000  preprocessing_CellBench.R