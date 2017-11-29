#!/usr/bin/env Rscript

# Description:
#  Compute ROC average
#
# Usage:
#
# dir_output_dt="/project2/gilad/joycehsiao/singlecell-de/output_eval/mousekleinesc/"
# dir_rocavg_dt="/project2/gilad/joycehsiao/singlecell-de/output_rocavg/mousekleinesc/"
# data="mousekleinesc"
#
# args=c(dir_output_dt, dir_rocavg_dt, data)
#
# Rscript ./combineResults.R ${dir_output_dt} ${dir_rocavg_dt} ${data}
#


args <- commandArgs(trailingOnly = TRUE)

for (i in seq_along(args)) {
  cat(sprintf("Argument %d: %s", i, args[i]), sep = "\n")
}

library(ashbun)

for (sam_method in c("all_genes", "per_gene")) {
#  for (gene_method in c("high", "medium", "low")) {
  for (gene_method in c("high")) {

  for (pi0 in c(.5, .9)) {
    for (Nsam in c(40)) {

     # make data labels
     data_name <- args[3]
     pi0_lab <- paste0("pi0", pi0*10)
     Nsamples_lab <- paste0("n", Nsam)

     tmp <- strsplit(sam_method, split = "_")[[1]]
     sam_label <- paste0(tmp[1],tmp[2])

     message(cat( data_name, ",", sam_label, ",",
                 "pi0=", pi0, ",", Nsam, "samples/condition", "\n"))

     output <- readRDS(paste0(args[1], args[3],".", sam_label, ".", gene_method, ".",
                      pi0_lab, ".", Nsamples_lab, ".bignormal.output.rds"))

     # combine results
     roc_combined <- do.call(rbind, lapply(output, "[[", 2))
     # compute ROC average
     roc_average <- do.call(rbind, lapply(1:nlevels(roc_combined$methodsMeanExpression),
                                 function(index_methodsMeanExpression) {
             meanMethod <- levels(roc_combined$methodsMeanExpression)[index_methodsMeanExpression]
             roc_combined_subset <- subset(roc_combined,
                 subset = methodsMeanExpression %in%  meanMethod )
             foo <- getROC.average(roc_combined_subset)
             foo$methodsMeanExpression <- meanMethod
             return(foo)
       }) )

    saveRDS(roc_average,
       file = paste0(args[2], args[3],".", sam_label, ".", gene_method, ".",
                     pi0_lab, ".", Nsamples_lab, ".bignormal.rocavg.rds"))
    }
  }
 }
}
