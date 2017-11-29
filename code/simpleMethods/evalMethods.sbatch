#!/bin/bash

#SBATCH --time=36:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joycehsiao@uchicago.edu
#SBATCH --partition=gilad
#SBATCH --account=pi-gilad
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=3G

set -e

dir_code=$1
dir_simdata_dt=$2
dir_output_dt=$3
data=$4
#methodsNormalize=$5
#methodsMeanExpression=$6

Rscript ${dir_code}/evalMethods.R ${dir_simdata_dt} ${dir_output_dt} ${data}
