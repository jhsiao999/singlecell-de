#!/bin/bash

#SBATCH --time=36:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joycehsiao@uchicago.edu
#SBATCH --partition=broadwl
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=8GB

set -e

dir_code=$1
data=$2
dir_output=$3

Rscript ${dir_code}/makeData.R ${data} ${dir_output}
