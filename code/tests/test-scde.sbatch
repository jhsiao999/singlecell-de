#!/bin/bash

#SBATCH --time=36:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joycehsiao@uchicago.edu
#SBATCH --partition=broadwl

set -e

dir_code=$1
#ncores=$2

Rscript ${dir_code}test-scde.R
