#!/bin/bash

#SBATCH --time=36:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joycehsiao@uchicago.edu
#SBATCH --partition=broadwl
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=8GB
#SBATCH --nodes=5

set -e

dir_code=$1
#ncores=$2

Rscript ${dir_code}test-bpsc.R
