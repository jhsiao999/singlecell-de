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
datasets=("gtex01" "mousekleinesc" "mouseengeltcell")
#datasets=("gtex01")
#datasets=("mousekleinesc")

dir_code="/project2/gilad/joycehsiao/singlecell-de/code/simpleMethods/"
dir_output="/project2/gilad/joycehsiao/singlecell-de/output_eval/"
dir_rocavg="/project2/gilad/joycehsiao/singlecell-de/output_rocavg/"

for data in ${datasets[@]}; do

  dir_rocavg_dt=${dir_rocavg}${data}/
  dir_output_dt=${dir_output}${data}/
  sbatch combineResults.sbatch ${dir_code} ${dir_output_dt} ${dir_rocavg_dt} ${data}

done
