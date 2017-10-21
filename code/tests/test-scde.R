#!/usr/bin/env Rscript

# Command line arguments can be read with commandArgs.
#
# Usage:
# dir_code="/project2/gilad/joycehsiao/singlecell-de/code/tests/"
# Rscript ${dir_code}test-scde.R
#

args <- commandArgs(trailingOnly = TRUE)

for (i in seq_along(args)) {
  cat(sprintf("Argument %d: %s", i, args[i]), sep = "\n")
}

#args[1]=ncores

ncores <- 10

library(ashbun)
library(scde)

# read in one simulated datasets
dir_simdata <- "/scratch/midway2/joycehsiao/singlecell-de/simulated_data/mousezeiselbrain/"
data <- readRDS(paste0(dir_simdata,"mousezeiselbrain.allgenes.pi05.n50.bignormal.data.rds"))

counts <- data[[1]]$counts
condition <- data[[1]]$condition

#sessionInfo()
#.libPaths()

res <- methodWrapper.scde(counts, condition,
                   control = list(save_modelFit = FALSE,
                                   n.randomizations = 100,
                                   n.cores = ncores,
                                   min.size.entries = ncol(counts)))

print(head(res))
