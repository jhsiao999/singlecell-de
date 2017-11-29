#!/bin/bash
set -e
#
# Description:
#   Submit a batch job for one experimental data
#
# Usage:
#
# ./submit-evalMethods
#
#
#datasets=("mousezeiselbrain" "mousekleinesc" "humantungipsc" \
#          "gtex001")

datasets=("mousekleinesc" "mouseengeltcell" "gtex01")
#datasets=("gtex01")
#datasets=("mousekleinesc")
#datasets=("mouseengeltcell")

#methodsNormalize=$1
#methodsMeanExpression=$2

dir_code="/project2/gilad/joycehsiao/singlecell-de/code/simpleMethods/"
dir_simdata="/project2/gilad/joycehsiao/singlecell-de/data_simulated/"
dir_output="/project2/gilad/joycehsiao/singlecell-de/output_eval/"

for data in ${datasets[@]}; do

  dir_simdata_dt=${dir_simdata}${data}/
  dir_output_dt=${dir_output}${data}/
  sbatch evalMethods.sbatch ${dir_code} ${dir_simdata_dt} ${dir_output_dt} ${data}

done
