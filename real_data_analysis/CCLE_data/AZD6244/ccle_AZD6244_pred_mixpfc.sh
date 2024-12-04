#!/bin/bash

#SBATCH --job-name="ccle_AZD6244_pred_mixpfc"
#SBATCH -N 1
#SBATCH -n 25
#SBATCH --cores-per-socket=8
#SBATCH -p backfill2
#SBATCH --mail-type="END, FAIL"

#SBATCH -o ccle_AZD6244_pred_mixpfc.out
#SBATCH -e ccle_AZD6244_pred_mixpfc.err
module load R
Rscript ccle_AZD6244_pred_mixpfc.R