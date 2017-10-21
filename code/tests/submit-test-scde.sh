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

dir_code="/project2/gilad/joycehsiao/singlecell-de/code/tests/"

sbatch ${dir_code}test-code.sbatch --ntasks=20 --mem-per-cpu=8GB --nodes=5 ${dir_code}

sbatch ${dir_code}test-code.sbatch --ntasks=10 --mem-per-cpu=8GB --nodes=5 ${dir_code}
