#!/bin/bash

#SBATCH --time=36:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joycehsiao@uchicago.edu
#SBATCH --partition=broadwl
#SBATCH --mem=20GB

set -e

dir_code=$1
dir_output_dt=$2
dir_rocavg_dt=$3
data=$4

Rscript ${dir_code}/combineResults.R ${dir_output_dt} ${dir_rocavg_dt} ${data}
