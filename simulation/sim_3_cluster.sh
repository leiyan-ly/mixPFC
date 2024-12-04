#!/bin/bash

#SBATCH --job-name="mixpfc_sim3_cluster"
#SBATCH -N 1
#SBATCH -n 25
#SBATCH --cores-per-socket=8
#SBATCH -p statistics_q
#SBATCH --mail-type="END, FAIL"

#SBATCH -o mixpfc_sim3_cluster.out
#SBATCH -e mixpfc_sim3_cluster.err
module load R
Rscript sim_3_cluster.R