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

sbatch ${dir_code}test-bpsc-1.sbatch ${dir_code}

sbatch ${dir_code}test-bpsc-2.sbatch ${dir_code}
