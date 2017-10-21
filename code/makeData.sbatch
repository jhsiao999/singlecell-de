#!/bin/bash

#SBATCH --time=36:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joycehsiao@uchicago.edu
#SBATCH --mem=15G
#SBATCH --partition=broadwl

set -e

dir_code=$1
data=$2
dir_output=$3

Rscript ${dir_code}/makeData.R ${data} ${dir_output}
