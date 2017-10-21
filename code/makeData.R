#!/usr/bin/env Rscript

# Command line arguments can be read with commandArgs.
#
# Usage:
#
# dir_output="/scratch/midway2/joycehsiao/singlecell-de/simulated_data/"
# data="/project2/gilad/joycehsiao/singlecell-de/data/mousezeiselbrain.rds"
#
# Rscript ./makeData.R ${data} ${dir_output}
#

args <- commandArgs(trailingOnly = TRUE)

for (i in seq_along(args)) {
  cat(sprintf("Argument %d: %s", i, args[i]), sep = "\n")
}

cat(5)

#dataset <- args[1]
#dir_output <- args[2]

#data <- readRDS(as.character(args[1]))

#print(args[1])
