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


#dataset <- args[1]
#dir_output <- args[2]

# manually set 20 simulated datsets
Nsim <- 20


# call libraries
library(ashbun)

data <- readRDS(args[1])

for (sam_method in c("all_genes", "per_gene")) {
  for (pi0 in c(.5, .9)) {
    for (Nsam in c(50)) {

    # make data labels
    data_name <- strsplit( basename(args[1]), split = ".", fix = TRUE)[[1]][1]
    pi0_lab <- paste0("pi0", pi0*10)
    Nsamples_lab <- paste0("n", Nsam)

    tmp <- strsplit(sam_method, split = "_")[[1]]
    sam_label <- paste0(tmp[1],tmp[2])

    message(cat( data_name, ",", sam_label, ",",
               "pi0=", pi0, ",", Nsam, "samples/condition", "\n"))

    # simulate data version 1
    # this version selects samples and genes, and then do filtering
#    simdata <- simulationWrapper(data, Nsim = Nsim, Ngenes = 1000,
#                                Nsamples = Nsam, sample_method = sam_method,
#                                pi0 = pi0,
#                                beta_args = args.big_normal(betapi = 1,
#                                                            betamu = 0,
#                                                            betasd = .8))

    # simulate data version 2
    # this version do filtering of samples, then select samples, then filter genes,
    # then select genes
    simdata <- simulationWrapper.filter(data, Nsim = Nsim, Ngenes = 1000,
                                Nsamples = Nsam, sample_method = sam_method,
                                pi0 = pi0,
                                samplesFractionExpressed=.25,
                                featuresFractionExpressed=.25,
                                beta_args = args.big_normal(betapi = 1,
                                                            betamu = 0,
                                                            betasd = .8))


    saveRDS(simdata,
      file = paste0(args[2], data_name, "/", data_name,".", sam_label, ".",
                    pi0_lab, ".", Nsamples_lab, ".bignormal.data.rds"))
    }
   }
 }
