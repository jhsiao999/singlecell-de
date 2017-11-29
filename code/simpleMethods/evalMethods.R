#!/usr/bin/env Rscript

# Command line arguments can be read with commandArgs.
#
# Usage:
#
# dir_simdata_dt="/project2/gilad/joycehsiao/singlecell-de/data_simulated/mousekleinesc/"
# dir_output_dt="/project2/gilad/joycehsiao/singlecell-de/output_eval/mousekleinesc/"
# data="mousekleinesc"
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
#library(methods)

# --- register parallel computing clusters
#library(doParallel)
#cl <- makeCluster(5)
#registerDoParallel(cl)

for (sam_method in c("all_genes", "per_gene")) {
#for (sam_method in c("per_gene")) {
  tmp <- strsplit(sam_method, split = "_")[[1]]
  sam_label <- paste0(tmp[1],tmp[2])

#  for (gene_method in c("high", "medium", "low")) {
 for (gene_method in c("high")) {

  for (pi0 in c(.5, .9)) {
#for (pi0 in c(.9)) {
    for (Nsam in c(40)) {

     # make data labels
     data_name <- args[3]
     pi0_lab <- paste0("pi0", pi0*10)
     Nsamples_lab <- paste0("n", Nsam)

     message(cat( data_name, ",", sam_label, ",", gene_method, ",",
                 "pi0=", pi0, ",", Nsam, "samples/condition", "\n"))

     simdata <- readRDS(paste0(args[1], args[3],".", sam_label, ".", gene_method, ".",
                        pi0_lab, ".", Nsamples_lab, ".bignormal.data.rds"))

     # apply methods
     eval_output <- vector("list", length(simdata))
     for (index in 1:length(simdata)) {

      message("------- simulated data ",index,"\n")

      eval_output[[index]]  <- query.evaluation.simple(counts = simdata[[index]]$counts,
               condition = simdata[[index]]$condition,
               is_nullgene = simdata[[index]]$is_nullgene,
               methodsMeanExpression = c("MAST", "limmaVoom",
                                        "DESeq2"),
               nsim = index)
     }

     message("saving: ", paste0(args[2], args[3],".", sam_label, ".", gene_method, ".",
                              pi0_lab, ".", Nsamples_lab, ".bignormal.output.rds"), "\n")
     saveRDS(eval_output,
        file = paste0(args[2], args[3],".", sam_label, ".", gene_method, ".",
                      pi0_lab, ".", Nsamples_lab, ".bignormal.output.rds"))
    }
   }
 }
}
