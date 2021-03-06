# --- Title: Testing ashbun package


# --- Load packages

library(Biobase)
library(ashbun)
library(singleCellRNASeqMouseZeiselBrain)


# --- Test for one simulated single cell data.

# extract count table
eset <- get(data("MouseZeiselBrain"))
counts <- exprs(eset)

# filter samples with zero depth
samples_any_gene <- which( colSums(counts) != 0)
counts <- counts[,samples_any_gene]

# filter genes with zero count
genes_to_include <- which(apply(counts,1,sum)>0)
counts <- counts[genes_to_include, ]


# --- generat one simulated data
simdata_list <- simulationWrapper(counts, Nsim = 1,
                                  Ngenes = 1000,
                                  Nsamples = 50,
                                  sample_method = "all_genes",
                                  pi0 = .9,
                                  beta_args = args.big_normal(betapi = 1,
                                                              betamu = 0, betasd = .8))

simdata <- simdata_list[[1]]

# --- register parallel computing clusters
library(doParallel)
cl <- makeCluster(8)
registerDoParallel(cl)

# ---- gather evaluation results
output  <- query.evaluation(counts = simdata$counts,
                            condition = simdata$condition,
                            is_nullgene = simdata$is_nullgene,
                            methodsNormalize = c("LIB", "TMM", "census"),
                            methodsMeanExpression = c("BPSC", "MAST", "ROTS", "SCDE",
                                                      "DESeq2", "edgeR", "limmaVoom"))

# library(ggplot2)
# qplot(x = FPR, y = TPR, data = output$roc,
#       colour = methodsMeanExpression,
#       facets = ~ methodsNormalize,
#       geom = "path")


Run on 5 simulated data.

```{r}
simdata_list <- simulationWrapper(counts, Nsim = 5,
                                  Ngenes = 1000,
                                  Nsamples = 50,
                                  sample_method = "all_genes",
                                  pi0 = .9,
                                  beta_args = args.big_normal(betapi = 1,
                                                              betamu = 0, betasd = .8))

output <- vector("list", 5)
for (index in 1:length(eval_output)) {
  output[[index]]  <- query.evaluation(counts = simdata_list[[index]]$counts,
                                       condition = simdata_list[[index]]$condition,
                                       is_nullgene = simdata_list[[index]]$is_nullgene,
                                       methodsNormalize = c("LIB", "TMM", "census","scran"),
                                       methodsMeanExpression = c("BPSC", "MAST", "ROTS", "SCDE",
                                                                 "DESeq2", "edgeR", "limmaVoom"),
                                       nsim = index)
}

# summarize ROC results
roc_combined <- do.call(rbind, lapply(output, "[[", 2))
roc_average <- do.call(rbind, lapply(1:nlevels(roc_combined$methodsNormalize),
                                     function(index_methodsNormalize) {
                                       foo2 <-
                                         do.call(rbind, lapply(1:nlevels(roc_combined$methodsMeanExpression),
                                                               function(index_methodsMeanExpression) {
                                                                 roc_combined_subset <- subset(roc_combined,
                                                                                               subset = methodsMeanExpression %in% methodsMeanExpression[index_methodsMeanExpression] & methodsNormalize %in% methodsNormalize[index_methodsNormalize] )
                                                                 foo <- getROC.average(roc_combined_subset)
                                                                 foo$methodsMeanExpression <- unique(methodsMeanExpression[index_methodsMeanExpression])
                                                                 return(foo)
                                                               }) )
                                       foo2$methodsNormalize <- unique(methodsNormalize[index_methodsNormalize])
                                       return(foo2)
                                     }) )


library(ggplot2)
qplot(data = roc_average,
      x = FPR, y = TPR,
      colour = methodsMeanExpression,
      facets = ~ methodsNormalize,
      geom = "path")
```
