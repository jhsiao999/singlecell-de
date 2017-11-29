#!/bin/bash

#SBATCH --time=36:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joycehsiao@uchicago.edu
#SBATCH --partition=gilad
#SBATCH --account=pi-gilad
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=1G

set -e

dir_code=$1
dir_output_dt=$2
dir_rocavg_dt=$3
data=$4

Rscript ${dir_code}/combineResults.R ${dir_output_dt} ${dir_rocavg_dt} ${data}
