library(ashbun)

Nsim <- 10
data <- c("mouseengeltcell.rds", "mousekleinesc.rds")
dir_output <- "output_tmp/"


for (index in 1:length(data)) {

  df <- readRDS(paste0("data/",data[index]))

for (sam_method in c("all_genes", "per_gene")) {

    for (pi0 in c(.5, .9)) {
      for (Nsam in c(40)) {

        # make data labels
        data_name <- strsplit( basename( data[index] ), split = ".", fix = TRUE)[[1]][1]
        pi0_lab <- paste0("pi0", pi0*10)
        Nsamples_lab <- paste0("n", Nsam)

        tmp <- strsplit(sam_method, split = "_")[[1]]
        sam_label <- paste0(tmp[1],tmp[2])

        message(cat( data_name, ",", sam_label, ",",
                     "pi0=", pi0, ",", Nsam, "samples/condition", "\n"))

        simdata <- simulationWrapper(df, Nsim = Nsim,
                                          Ngenes = 1000,
                                          Nsamples = Nsam,
                                          sample_method =  sam_method,
                                          pi0 = pi0,
                                          beta_args = args.big_normal(betapi = 1,
                                                                      betamu = 0, betasd = .8))


        saveRDS(simdata,
                file = paste0(dir_output, data_name,".", sam_label, ".",
                              pi0_lab, ".", Nsamples_lab, ".bignormal.data.rds"))
      }
    }
  }
}



#--------- compute simulation results

library(ashbun)
library(methods)

# --- register parallel computing clusters
library(doParallel)
cl <- makeCluster(8)
registerDoParallel(cl)

for (index in 1:length(data)) {
#for (sam_method in c("all_genes", "per_gene")) {
for (sam_method in c("all_genes")) {
  for (gene_method in c("high", "medium")) {

  for (pi0 in c(.5, .9)) {
    for (Nsam in c(50)) {

      # make data labels
      data_name <- strsplit( basename( data[index] ), split = ".", fix = TRUE)[[1]][1]
      pi0_lab <- paste0("pi0", pi0*10)
      Nsamples_lab <- paste0("n", Nsam)

      tmp <- strsplit(sam_method, split = "_")[[1]]
      sam_label <- paste0(tmp[1],tmp[2])

      message(cat( data_name, ",", sam_label, ",", gene_method, ",",
                   "pi0=", pi0, ",", Nsam, "samples/condition", "\n"))

      simdata <- readRDS(paste0(dir_output, data_name,".", sam_label, ".", gene_method, ".",
                                pi0_lab, ".", Nsamples_lab, ".bignormal.data.rds"))

      # apply methods
      eval_output <- vector("list", length(simdata))
      for (index in 1:length(simdata)) {

        message("------- simulated data ",index,"\n")
        eval_output[[index]]  <- query.evaluation.simple(counts = simdata[[index]]$counts,
                                      condition = simdata[[index]]$condition,
                                      is_nullgene = simdata[[index]]$is_nullgene,
                                      methodsMeanExpression = c("limmaVoom", "MAST"),
                                      report.control = list(fdr_control_threshold = .05),
                                      nsim = index)
      }
      saveRDS(eval_output,
              file = paste0(dir_output, data_name,".", sam_label, ".", gene_method, ".",
                            pi0_lab, ".", Nsamples_lab, ".bignormal.output.rds"))
    }
  }
}
}
}




output <- vector("list", length(simdata_list))

for (index in 1:length(output)) {
  output[[index]]  <- query.evaluation.simple(counts = simdata_list[[index]]$counts,
                                       condition = simdata_list[[index]]$condition,
                                       is_nullgene = simdata_list[[index]]$is_nullgene,
                                       methodsMeanExpression = c("DESeq2", "limmaVoom",
                                                                 "MAST"),
                                       report.control = list(fdr_control_threshold = .05),
                                       nsim = index)
}
results.multiple.data <- saveRDS(output, file = "tests.results.multiple.data.rds")


