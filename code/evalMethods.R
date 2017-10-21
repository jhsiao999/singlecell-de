#!/usr/bin/env Rscript

# Command line arguments can be read with commandArgs.
#
# Usage:
#
# dir_simdata_dt="/scratch/midway2/joycehsiao/singlecell-de/simulated_data/mousezeiselbrain/"
# dir_output_dt="/project2/gilad/joycehsiao/singlecell-de/output_eval/mousezeiselbrain/"
# data="mousezeiselbrain"
#
# args=c(dir_simdata_dt, dir_output_dt, data)
#
# Rscript ./evalMethods.R ${dir_simdata_dt} ${dir_output_dt} ${data}
#


args <- commandArgs(trailingOnly = TRUE)

for (i in seq_along(args)) {
  cat(sprintf("Argument %d: %s", i, args[i]), sep = "\n")
}

library(ashbun)

# --- register parallel computing clusters
library(doParallel)
cl <- makeCluster(8)
registerDoParallel(cl)

for (sam_method in c("all_genes", "per_gene")) {
  for (pi0 in c(.5, .9)) {
    for (Nsam in c(50)) {

     # make data labels
     data_name <- args[3]
     pi0_lab <- paste0("pi0", pi0*10)
     Nsamples_lab <- paste0("n", Nsam)

     tmp <- strsplit(sam_method, split = "_")[[1]]
     sam_label <- paste0(tmp[1],tmp[2])

     message(cat( data_name, ",", sam_label, ",",
                 "pi0=", pi0, ",", Nsam, "samples/condition", "\n"))

     simdata <- readRDS(paste0(args[1], args[3],".", sam_label, ".",
                        pi0_lab, ".", Nsamples_lab, ".bignormal.data.rds"))

     # apply methods
     eval_output <- vector("list", length(simdata))
     for (index in 1:length(simdata)) {
     
      message("------- simulated data ",index,"\n")
      eval_output[[index]]  <- query.evaluation(counts = simdata[[index]]$counts,
               condition = simdata[[index]]$counts,
               is_nullgene = simdata[[index]]$counts,
               thresholdDetection = 1, fractionExpressed = .01,
               methodsNormalize = c("LIB", "TMM", "census"),
               methodsMeanExpression = c("MAST", "BPSC",
                                         "DESeq2", "edgeR", "limmaVoom"),
#               methodsMeanExpression = c("BPSC", "MAST", "ROTS", "SCDE",
#                                         "DESeq2", "edgeR", "limmaVoom"),
               nsim = index)
     }
     saveRDS(eval_output,
        file = paste0(args[2], args[3],".", sam_label, ".",
                      pi0_lab, ".", Nsamples_lab, ".bignormal.output.rds"))
    }
   }
 }
