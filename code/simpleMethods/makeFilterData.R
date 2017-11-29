#!/usr/bin/env Rscript

# Description:
#  filter samples and genes and save the filtered data to rds
# Usage:
#
# dir_output="/scratch/midway2/joycehsiao/singlecell-de/data_n50/"
# dir_input="/project2/gilad/joycehsiao/singlecell-de/data/"

dataset <- c("mousezeiselbrain.rds", "mousekleinesc.rds",
             "humantungipsc.rds", "gtex001.rds")

library(ashbun)

#----- filtering

for (index in 1:length(dataset)) {
 tmp <- dataset[index]
 data <- readRDS(paste0(dir_input, tmp))

 
 fl_filtered <- filter.Wrapper(counts = fl,
                                 condition = condition,
                                 thresholdDetection = thresholdDetection,
                                 fractionExpressed = fractionExpressed,
                                 is_nullgene = is_nullgene)
}
data_filtered <- filter.Wrapper(counts = counts,
                                condition = condition,
                                thresholdDetection = thresholdDetection,
                                fractionExpressed = fractionExpressed,
                                is_nullgene = is_nullgene)





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

    # simulate data
    simdata <- simulationWrapper(data, Nsim = Nsim, Ngenes = 1000,
                                Nsamples = Nsam, sample_method = sam_method,
                                pi0 = pi0,
                                beta_args = args.big_normal(betapi = 1,
                                                            betamu = 0,
                                                            betasd = .8))


    saveRDS(simdata,
      file = paste0(args[2], data_name, "/", data_name,".", sam_label, ".",
                    pi0_lab, ".", Nsamples_lab, ".bignormal.data.rds"))
    }
   }
 }
