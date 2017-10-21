#!/bin/bash
set -e
#
# Description:
#   Submit a batch job for one experimental data
#
# Usage:
#
# ./script.sh ex1 ex2 ex3
#
#
#datasets=("mousezeiselbrain" "mousekleinesc" "humantungipsc" \
#          "gtex001")
#datasets=("mousekleinesc" "humantungipsc")
datasets=("mousezeiselbrain")

dir_code="/project2/gilad/joycehsiao/singlecell-de/code/"
dir_simdata="/scratch/midway2/joycehsiao/singlecell-de/simulated_data/"
dir_output="/project2/gilad/joycehsiao/singlecell-de/output_eval/"

for data in ${datasets[@]}; do

  dir_simdata_dt=${dir_simdata}${data}/
  dir_output_dt=${dir_output}${data}/
  sbatch evalMethods.sbatch ${dir_code} ${dir_simdata_dt} ${dir_output_dt} ${data}

done
