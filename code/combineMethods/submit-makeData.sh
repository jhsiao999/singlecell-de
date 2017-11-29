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
#datasets=("mousezeiselbrain.rds" "mousekleinesc.rds" "humantungipsc.rds" \
#            "gtex01.rds")

datasets=("mousekleinesc.rds" "mouseengeltcell" "gtex01.rds")


dir_datasets="/project2/gilad/joycehsiao/singlecell-de/data/"
dir_code="/project2/gilad/joycehsiao/singlecell-de/code/"
dir_output="/scratch/midway2/joycehsiao/singlecell-de/simulated_data/"

for data in ${datasets[@]}; do

  data=${dir_datasets}${data}
  sbatch makeData.sbatch ${dir_code} ${data} ${dir_output}

done
